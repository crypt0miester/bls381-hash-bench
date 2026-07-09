//! Syscall-assisted RFC 9380 hash_to_G1 for BLS12-381 (min-sig suite
//! BLS_SIG_BLS12381G1_XMD:SHA-256_SSWU_RO_POP_).
//!
//! Strategy: expand_message_xmd through the sha256 syscall, wide reduction /
//! inversion / Legendre / sqrt through big_mod_exp, SSWU + iso-11 polynomial
//! evaluation with an in-program variable-time Montgomery multiplier, the two
//! mapped points added on the isogenous curve so the isogeny runs once, and
//! cofactor clearing as a double-and-add chain over the g1 add syscall.

use solana_program_error::ProgramError;
use alloc::vec::Vec;

use crate::consts_g1::{
    ISO11A_XDEN, ISO11A_XNUM, ISO11A_YDEN, ISO11A_YNUM, ISO11_XDEN, ISO11_XNUM, ISO11_YDEN,
    ISO11_YNUM, R, R2, SSWU_ELLP_A, SSWU_ELLP_B,
};
use crate::consts_g2::C256_MONT;
use crate::fp::*;

const H_EFF: u64 = 0xd201000000010001;

fn expand_message_xmd(dst: &[u8], msg: &[u8]) -> [[u8; 32]; 4] {
    use solana_sha256_hasher::hashv;

    let z_pad = [0u8; 64];
    let l_i_b = [0u8, 128];
    let dst_len = [dst.len() as u8];

    let b0 = hashv(&[&z_pad, msg, &l_i_b, &[0u8], dst, &dst_len]).to_bytes();

    let mut blocks = [[0u8; 32]; 4];
    blocks[0] = hashv(&[&b0, &[1u8], dst, &dst_len]).to_bytes();
    for i in 1..4 {
        let mut x = [0u8; 32];
        for j in 0..32 {
            x[j] = b0[j] ^ blocks[i - 1][j];
        }
        blocks[i] = hashv(&[&x, &[i as u8 + 1], dst, &dst_len]).to_bytes();
    }
    blocks
}

/// hash_to_field for two Fp elements, canonical form.
fn hash_to_field(dst: &[u8], msg: &[u8]) -> Result<[Fp; 2], ProgramError> {
    let blocks = expand_message_xmd(dst, msg);
    let one = [1u8];

    let mut wide = [0u8; 64];
    wide[..32].copy_from_slice(&blocks[0]);
    wide[32..].copy_from_slice(&blocks[1]);
    let u0 = be_to_limbs(&modexp_bytes(&wide, &one)?);

    wide[..32].copy_from_slice(&blocks[2]);
    wide[32..].copy_from_slice(&blocks[3]);
    let u1 = be_to_limbs(&modexp_bytes(&wide, &one)?);

    Ok([u0, u1])
}

struct Exps {
    inv: [u8; 48],
    legendre: [u8; 48],
    sqrt: [u8; 48],
}

impl Exps {
    fn new() -> Self {
        Self {
            inv: exp_inverse(),
            legendre: exp_legendre(),
            sqrt: exp_sqrt(),
        }
    }
}

/// Affine point on the 11-isogenous curve E', Montgomery form coordinates.
struct PointPrime {
    x: Fp,
    y: Fp,
}

fn is_one_canonical(a: &Fp) -> bool {
    a[0] == 1 && a[1..].iter().all(|&l| l == 0)
}

/// xi = 11: multiply by the SSWU non-residue with an addition chain.
fn mul_by_xi(a: &Fp) -> Fp {
    let a2 = add_mod(a, a);
    let a4 = add_mod(&a2, &a2);
    let a8 = add_mod(&a4, &a4);
    add_mod(&add_mod(&a8, &a2), a)
}

/// Simplified SWU map onto E', per RFC 9380 section 6.6.2 (variable time).
fn map_to_curve_sswu(u: &Fp, c_neg_b_over_a: &Fp, exps: &Exps) -> Result<PointPrime, ProgramError> {
    let um = to_mont(u);

    let usq = mont_sqr(&um);
    let xi_usq = mul_by_xi(&usq);
    let tv2 = add_mod(&mont_sqr(&xi_usq), &xi_usq);

    // tv2 == 0 has probability ~2^-381; the exceptional-case branch is omitted
    if is_zero(&tv2) {
        return Err(ProgramError::InvalidInstructionData);
    }

    let inv_tv2 = inverse_mont(&tv2, &exps.inv)?;
    let x1 = mont_mul(c_neg_b_over_a, &add_mod(&R, &inv_tv2));

    let gx = |x: &Fp| -> Fp {
        let xsq = mont_sqr(x);
        let x3 = mont_mul(&xsq, x);
        add_mod(&add_mod(&x3, &mont_mul(&SSWU_ELLP_A, x)), &SSWU_ELLP_B)
    };

    let gx1 = gx(&x1);
    let legendre = modexp(&from_mont(&gx1), &exps.legendre)?;

    let (x, gx_val) = if is_one_canonical(&legendre) {
        (x1, gx1)
    } else {
        let x2 = mont_mul(&xi_usq, &x1);
        let gx2 = gx(&x2);
        (x2, gx2)
    };

    let mut y = modexp(&from_mont(&gx_val), &exps.sqrt)?;

    // sgn0 correction: parity of y must match parity of u
    if (y[0] & 1) != (u[0] & 1) {
        y = neg_mod(&y);
    }

    Ok(PointPrime { x, y: to_mont(&y) })
}

/// Affine addition on E' (variable time; errors on the infinity outcome).
fn add_prime(p: &PointPrime, q: &PointPrime, exps: &Exps) -> Result<PointPrime, ProgramError> {
    let lambda = if p.x == q.x {
        if p.y != q.y || is_zero(&p.y) {
            return Err(ProgramError::InvalidInstructionData);
        }
        // doubling: (3x^2 + A) / 2y
        let xsq = mont_sqr(&p.x);
        let num = add_mod(&add_mod(&add_mod(&xsq, &xsq), &xsq), &SSWU_ELLP_A);
        let den = add_mod(&p.y, &p.y);
        mont_mul(&num, &inverse_mont(&den, &exps.inv)?)
    } else {
        let num = sub_mod(&q.y, &p.y);
        let den = sub_mod(&q.x, &p.x);
        mont_mul(&num, &inverse_mont(&den, &exps.inv)?)
    };

    let x3 = sub_mod(&sub_mod(&mont_sqr(&lambda), &p.x), &q.x);
    let y3 = sub_mod(&mont_mul(&lambda, &sub_mod(&p.x, &x3)), &p.y);
    Ok(PointPrime { x: x3, y: y3 })
}

/// Evaluate the four iso-11 polynomials with Knuth-adapted constants:
/// 27 multiplications against 51 for shared-nothing Horner. 
pub(crate) fn iso11_adapted(x: &Fp) -> (Fp, Fp, Fp, Fp) {
    let w = mont_mul(x, x);
    let mut xnum = add_mod(
        &mont_mul(&add_mod(x, &ISO11A_XNUM[0]), &add_mod(&w, &ISO11A_XNUM[1])),
        &ISO11A_XNUM[2],
    );
    let t = add_mod(&w, &ISO11A_XNUM[3]);
    xnum = add_mod(&mont_mul(&xnum, &t), &ISO11A_XNUM[4]);
    let mut t = add_mod(&w, &ISO11A_XNUM[5]);
    t = add_mod(&t, x);
    xnum = add_mod(&mont_mul(&xnum, &t), &ISO11A_XNUM[6]);
    let mut t = add_mod(&w, &ISO11A_XNUM[7]);
    t = add_mod(&t, x);
    t = add_mod(&t, x);
    t = add_mod(&t, x);
    xnum = add_mod(&mont_mul(&xnum, &t), &ISO11A_XNUM[8]);
    let t = add_mod(&w, &ISO11A_XNUM[9]);
    xnum = add_mod(&mont_mul(&xnum, &t), &ISO11A_XNUM[10]);
    let xnum = mont_mul(&xnum, &ISO11A_XNUM[11]);

    let t = add_mod(&add_mod(&w, &mont_mul(&ISO11A_XDEN[0], x)), &ISO11A_XDEN[1]);
    let mut xden = add_mod(&mont_mul(&t, &add_mod(&w, &ISO11A_XDEN[2])), &ISO11A_XDEN[3]);
    let mut t = add_mod(&w, &ISO11A_XDEN[4]);
    t = add_mod(&t, x);
    t = add_mod(&t, x);
    xden = add_mod(&mont_mul(&xden, &t), &ISO11A_XDEN[5]);
    let mut t = add_mod(&w, &ISO11A_XDEN[6]);
    t = add_mod(&t, x);
    xden = add_mod(&mont_mul(&xden, &t), &ISO11A_XDEN[7]);
    let mut t = add_mod(&w, &ISO11A_XDEN[8]);
    t = add_mod(&t, x);
    t = add_mod(&t, x);
    xden = add_mod(&mont_mul(&xden, &t), &ISO11A_XDEN[9]);

    let mut ynum = add_mod(
        &mont_mul(&add_mod(x, &ISO11A_YNUM[0]), &add_mod(&w, &ISO11A_YNUM[1])),
        &ISO11A_YNUM[2],
    );
    let mut t = add_mod(&w, &ISO11A_YNUM[3]);
    t = add_mod(&t, x);
    ynum = add_mod(&mont_mul(&ynum, &t), &ISO11A_YNUM[4]);
    let t = add_mod(&w, &ISO11A_YNUM[5]);
    ynum = add_mod(&mont_mul(&ynum, &t), &ISO11A_YNUM[6]);
    let mut t = add_mod(&w, &ISO11A_YNUM[7]);
    t = add_mod(&t, x);
    ynum = add_mod(&mont_mul(&ynum, &t), &ISO11A_YNUM[8]);
    let t = add_mod(&w, &ISO11A_YNUM[9]);
    ynum = add_mod(&mont_mul(&ynum, &t), &ISO11A_YNUM[10]);
    let mut t = add_mod(&w, &ISO11A_YNUM[11]);
    t = add_mod(&t, x);
    t = add_mod(&t, x);
    ynum = add_mod(&mont_mul(&ynum, &t), &ISO11A_YNUM[12]);
    let t = add_mod(&w, &ISO11A_YNUM[13]);
    ynum = add_mod(&mont_mul(&ynum, &t), &ISO11A_YNUM[14]);
    let ynum = mont_mul(&ynum, &ISO11A_YNUM[15]);

    let mut yden = add_mod(
        &mont_mul(&add_mod(x, &ISO11A_YDEN[0]), &add_mod(&w, &ISO11A_YDEN[1])),
        &ISO11A_YDEN[2],
    );
    let mut t = add_mod(&w, &ISO11A_YDEN[3]);
    t = add_mod(&t, x);
    yden = add_mod(&mont_mul(&yden, &t), &ISO11A_YDEN[4]);
    let mut t = add_mod(&w, &ISO11A_YDEN[5]);
    t = add_mod(&t, x);
    t = add_mod(&t, x);
    yden = add_mod(&mont_mul(&yden, &t), &ISO11A_YDEN[6]);
    let t = add_mod(&w, &ISO11A_YDEN[7]);
    yden = add_mod(&mont_mul(&yden, &t), &ISO11A_YDEN[8]);
    let mut t = add_mod(&w, &ISO11A_YDEN[9]);
    t = add_mod(&t, x);
    yden = add_mod(&mont_mul(&yden, &t), &ISO11A_YDEN[10]);
    let t = add_mod(&w, &ISO11A_YDEN[11]);
    yden = add_mod(&mont_mul(&yden, &t), &ISO11A_YDEN[12]);
    let t = add_mod(&w, &ISO11A_YDEN[13]);
    yden = add_mod(&mont_mul(&yden, &t), &ISO11A_YDEN[14]);

    (xnum, xden, ynum, yden)
}

fn horner(coeffs: &[[u64; 6]], x: &Fp) -> Fp {
    let mut acc = coeffs[coeffs.len() - 1];
    for c in coeffs[..coeffs.len() - 1].iter().rev() {
        acc = add_mod(&mont_mul(&acc, x), c);
    }
    acc
}

/// 11-isogeny from E' to E, affine in and out (Montgomery form).
fn iso_map(p: &PointPrime, exps: &Exps) -> Result<([u8; 48], [u8; 48]), ProgramError> {
    let (x_num, x_den, y_num, y_den) = iso11_adapted(&p.x);

    // batch inversion of both denominators
    let t = mont_mul(&x_den, &y_den);
    let t_inv = inverse_mont(&t, &exps.inv)?;
    let x_den_inv = mont_mul(&t_inv, &y_den);
    let y_den_inv = mont_mul(&t_inv, &x_den);

    let x = mont_mul(&x_num, &x_den_inv);
    let y = mont_mul(&p.y, &mont_mul(&y_num, &y_den_inv));

    Ok((limbs_to_be(&from_mont(&x)), limbs_to_be(&from_mont(&y))))
}

const BLS12_381_G1_BE: u64 = 5 | 0x80;
const OP_ADD: u64 = 0;

fn g1_add(a: &[u8; 96], b: &[u8; 96]) -> Result<[u8; 96], ProgramError> {
    let mut out = [0u8; 96];
    let rc = unsafe {
        sys::sol_curve_group_op(
            BLS12_381_G1_BE,
            OP_ADD,
            a.as_ptr(),
            b.as_ptr(),
            out.as_mut_ptr(),
        )
    };
    if rc != 0 {
        return Err(ProgramError::InvalidInstructionData);
    }
    Ok(out)
}

/// Multiplies by the effective cofactor with double-and-add over the g1 add
/// syscall, which skips the subgroup check that blocks the mul syscall here.
pub(crate) fn clear_cofactor(p: &[u8; 96]) -> Result<[u8; 96], ProgramError> {
    let mut acc = *p;
    for bit in (0..63).rev() {
        acc = g1_add(&acc, &acc)?;
        if (H_EFF >> bit) & 1 == 1 {
            acc = g1_add(&acc, p)?;
        }
    }
    Ok(acc)
}

pub(crate) fn validate(p: &[u8; 96]) -> Result<(), ProgramError> {
    let mut out = 0u8;
    let rc = unsafe { sys::sol_curve_validate_point(BLS12_381_G1_BE, p.as_ptr(), &mut out) };
    if rc != 0 {
        return Err(ProgramError::InvalidInstructionData);
    }
    Ok(())
}

pub(crate) fn point_bytes(x: &[u8; 48], y: &[u8; 48]) -> [u8; 96] {
    let mut out = [0u8; 96];
    out[..48].copy_from_slice(x);
    out[48..].copy_from_slice(y);
    out
}

fn c_neg_b_over_a(exps: &Exps) -> Result<Fp, ProgramError> {
    let a_inv = inverse_mont(&SSWU_ELLP_A, &exps.inv)?;
    Ok(neg_mod(&mont_mul(&SSWU_ELLP_B, &a_inv)))
}

/// Stages, cumulative: 0 = hash_to_field, 1 = + both SSWU maps,
/// 2 = + E' add + isogeny, 3 = full with cofactor clearing and validation.
pub fn run(dst: &[u8], stage: u8, msg: &[u8]) -> Result<Vec<u8>, ProgramError> {
    let u = hash_to_field(dst, msg)?;
    if stage == 0 {
        let mut out = Vec::with_capacity(96);
        out.extend_from_slice(&limbs_to_be(&u[0]));
        out.extend_from_slice(&limbs_to_be(&u[1]));
        return Ok(out);
    }

    let exps = Exps::new();
    let c = SSWU_C1_NEG_B_OVER_A;
    let p0 = map_to_curve_sswu(&u[0], &c, &exps)?;
    let p1 = map_to_curve_sswu(&u[1], &c, &exps)?;
    if stage == 1 {
        return Ok(limbs_to_be(&from_mont(&p0.x)).to_vec());
    }

    let sum = add_prime(&p0, &p1, &exps)?;
    let (x, y) = iso_map(&sum, &exps)?;
    let uncleared = point_bytes(&x, &y);
    if stage == 2 {
        return Ok(uncleared.to_vec());
    }

    let cleared = clear_cofactor(&uncleared)?;
    validate(&cleared)?;
    Ok(cleared.to_vec())
}


/// CIOS variant: single interleaved pass, less memory traffic than SOS.
const W_MAP: usize = 1 + 48 + 48;
const W_TOTAL: usize = 2 * W_MAP + 3 * 48;

use crate::consts_g1::SSWU_C1_NEG_B_OVER_A;

struct FieldElem {
    canonical: Fp,
    mont: Fp,
}

/// hash_to_field without modexp: split the 64-byte value at bit 256 and fold
/// with a precomputed 2^256 mod p.
fn hash_to_field_folded(dst: &[u8], msg: &[u8]) -> [FieldElem; 2] {
    let blocks = expand_message_xmd(dst, msg);
    let mut out = [
        FieldElem { canonical: [0; 6], mont: [0; 6] },
        FieldElem { canonical: [0; 6], mont: [0; 6] },
    ];
    for (i, elem) in out.iter_mut().enumerate() {
        let mut hi = [0u8; 48];
        let mut lo = [0u8; 48];
        hi[16..].copy_from_slice(&blocks[i * 2]);
        lo[16..].copy_from_slice(&blocks[i * 2 + 1]);
        // canonical * Montgomery-form constant gives a canonical product
        let t = mont_mul(&be_to_limbs(&hi), &C256_MONT);
        let canonical = add_mod(&t, &be_to_limbs(&lo));
        elem.canonical = canonical;
        elem.mont = to_mont(&canonical);
    }
    out
}

fn gx_at(x: &Fp) -> Fp {
    let xsq = mont_sqr(x);
    let x3 = mont_mul(&xsq, x);
    add_mod(&add_mod(&x3, &mont_mul(&SSWU_ELLP_A, x)), &SSWU_ELLP_B)
}

/// The witness arrives in Montgomery form, so the check is one multiply.
pub(crate) fn check_inverse(v_m: &Fp, witness_m: &Fp) -> Result<Fp, ProgramError> {
    if mont_mul(v_m, witness_m) != R {
        return Err(ProgramError::InvalidInstructionData);
    }
    Ok(*witness_m)
}

fn map_to_curve_witnessed(u: &FieldElem, wit: &[u8]) -> Result<PointPrime, ProgramError> {
    let flag = wit[0];
    let w_inv = wit48(&wit[1..49])?;
    let y_w = wit48(&wit[49..97])?;

    let usq = mont_sqr(&u.mont);
    let xi_usq = mul_by_xi(&usq);
    let tv2 = add_mod(&mont_sqr(&xi_usq), &xi_usq);
    if is_zero(&tv2) {
        return Err(ProgramError::InvalidInstructionData);
    }

    let inv_m = check_inverse(&tv2, &w_inv)?;
    let x1 = mont_mul(&SSWU_C1_NEG_B_OVER_A, &add_mod(&R, &inv_m));

    // gx2 = (Z u^2)^3 gx1 with Z a non-residue, so gx1 and gx2 always have
    // opposite quadratic characters: one sqrt witness proves its own branch.
    // A flag above 1 is non-canonical. The lone ambiguity is gx1 == 0 (which
    // forces gx2 == 0): both branches would then accept y = 0 and let the flag
    // steer the output between (x1,0) and (x2,0). blst takes x1 (is_square(0)
    // is true), so reject the x2 branch when gx == 0.
    if flag > 1 {
        return Err(ProgramError::InvalidInstructionData);
    }
    let (x, gx) = if flag == 0 {
        (x1, gx_at(&x1))
    } else {
        let x2 = mont_mul(&xi_usq, &x1);
        (x2, gx_at(&x2))
    };
    if flag == 1 && is_zero(&gx) {
        return Err(ProgramError::InvalidInstructionData);
    }

    let yw_m = to_mont(&y_w);
    if mont_sqr(&yw_m) != gx {
        return Err(ProgramError::InvalidInstructionData);
    }

    // Negation commutes with the Montgomery map, so flip the square root we
    // already converted rather than converting the canonical value a second time.
    let y = if (y_w[0] & 1) != (u.canonical[0] & 1) {
        neg_mod(&yw_m)
    } else {
        yw_m
    };

    Ok(PointPrime { x, y })
}

fn add_prime_witnessed(
    p: &PointPrime,
    q: &PointPrime,
    w_dx: &Fp,
) -> Result<PointPrime, ProgramError> {
    if p.x == q.x {
        return Err(ProgramError::InvalidInstructionData);
    }
    let dx = sub_mod(&q.x, &p.x);
    let inv_m = check_inverse(&dx, w_dx)?;
    let lambda = mont_mul(&sub_mod(&q.y, &p.y), &inv_m);
    let x3 = sub_mod(&sub_mod(&mont_sqr(&lambda), &p.x), &q.x);
    let y3 = sub_mod(&mont_mul(&lambda, &sub_mod(&p.x, &x3)), &p.y);
    Ok(PointPrime { x: x3, y: y3 })
}

fn iso_map_witnessed(
    p: &PointPrime,
    w_xd: &Fp,
    w_yd: &Fp,
) -> Result<([u8; 48], [u8; 48]), ProgramError> {
    let (x_num, x_den, y_num, y_den) = iso11_adapted(&p.x);

    let xd_inv = check_inverse(&x_den, w_xd)?;
    let yd_inv = check_inverse(&y_den, w_yd)?;

    let x = mont_mul(&x_num, &xd_inv);
    let y = mont_mul(&p.y, &mont_mul(&y_num, &yd_inv));
    Ok((limbs_to_be(&from_mont(&x)), limbs_to_be(&from_mont(&y))))
}



/// Single-element hash_to_field for the NU (encode_to_curve) suite.
fn hash_to_field_nu(dst: &[u8], msg: &[u8]) -> FieldElem {
    use solana_sha256_hasher::hashv;
    let z_pad = [0u8; 64];
    let l_i_b = [0u8, 64];
    let dst_len = [dst.len() as u8];
    let b0 = hashv(&[&z_pad, msg, &l_i_b, &[0u8], dst, &dst_len]).to_bytes();
    let b1 = hashv(&[&b0, &[1u8], dst, &dst_len]).to_bytes();
    let mut x = [0u8; 32];
    for j in 0..32 {
        x[j] = b0[j] ^ b1[j];
    }
    let b2 = hashv(&[&x, &[2u8], dst, &dst_len]).to_bytes();

    let mut hi = [0u8; 48];
    let mut lo = [0u8; 48];
    hi[16..].copy_from_slice(&b1);
    lo[16..].copy_from_slice(&b2);
    let t = mont_mul(&be_to_limbs(&hi), &C256_MONT);
    let canonical = add_mod(&t, &be_to_limbs(&lo));
    FieldElem { canonical, mont: to_mont(&canonical) }
}

/// Witnessed encode_to_curve (RFC 9380 NU): one map, no addition.
/// Blob: flag, w_inv, y, w_xd, w_yd then the message.
pub fn encode_to_g1(dst: &[u8], payload: &[u8]) -> Result<Vec<u8>, ProgramError> {
    const NU_TOTAL: usize = W_MAP + 2 * 48;
    let (wits, msg) = split_witness(payload, NU_TOTAL)?;
    let w_xd = wit48(&wits[W_MAP..W_MAP + 48])?;
    let w_yd = wit48(&wits[W_MAP + 48..])?;

    let u = hash_to_field_nu(dst, msg);
    let p = map_to_curve_witnessed(&u, &wits[..W_MAP])?;
    let (x, y) = iso_map_witnessed(&p, &w_xd, &w_yd)?;
    let cleared = clear_cofactor(&point_bytes(&x, &y))?;
    validate(&cleared)?;
    Ok(cleared.to_vec())
}

pub fn hash_to_g1(dst: &[u8], payload: &[u8]) -> Result<Vec<u8>, ProgramError> {
    let (wits, msg) = split_witness(payload, W_TOTAL)?;

    let u = hash_to_field_folded(dst, msg);
    let p0 = map_to_curve_witnessed(&u[0], &wits[..W_MAP])?;
    let p1 = map_to_curve_witnessed(&u[1], &wits[W_MAP..2 * W_MAP])?;

    let base = 2 * W_MAP;
    let w_dx = wit48(&wits[base..base + 48])?;
    let w_xd = wit48(&wits[base + 48..base + 96])?;
    let w_yd = wit48(&wits[base + 96..base + 144])?;

    let sum = add_prime_witnessed(&p0, &p1, &w_dx)?;
    let (x, y) = iso_map_witnessed(&sum, &w_xd, &w_yd)?;

    let cleared = clear_cofactor(&point_bytes(&x, &y))?;
    validate(&cleared)?;
    Ok(cleared.to_vec())
}

/// Host-side witness generation, mirroring the on-chain pipeline with the
/// expensive results computed via square-and-multiply.
#[cfg(not(target_os = "solana"))]
pub mod witness {
    use super::*;

    /// The adapted chains must agree with Horner over the original
    /// coefficient tables at arbitrary points.
    pub fn iso11_adapted_selftest() {
        let mut x = R2;
        for i in 0..50u64 {
            let expect = (
                horner(&ISO11_XNUM, &x),
                horner(&ISO11_XDEN, &x),
                horner(&ISO11_YNUM, &x),
                horner(&ISO11_YDEN, &x),
            );
            assert_eq!(iso11_adapted(&x), expect, "adapted iso11 diverged");
            x = add_mod(&mont_mul(&x, &R2), &[i, 1, 0, 0, 0, 0]);
        }
    }

    /// Correctness guard: the bare Montgomery reduction must equal the general
    /// multiply by ONE at arbitrary points, edge cases included.
    pub fn redc_selftest() {
        use crate::consts_g1::MODULUS;
        let one: Fp = [1, 0, 0, 0, 0, 0];
        let mut samples: Vec<Fp> =
            alloc::vec![[0u64; 6], one, R, R2, sub_nocheck(&MODULUS, &one)];
        let mut x = R2;
        for i in 0..300u64 {
            x = add_mod(&mont_mul(&x, &R2), &[i, 1, 0, 0, 0, 0]);
            samples.push(x);
        }
        for s in &samples {
            assert_eq!(mont_redc_cios32(s), mont_mul(s, &one), "redc != from_mont");
        }
    }

    fn pow_mont(base: &Fp, exp_be: &[u8; 48]) -> Fp {
        let mut acc = R;
        for byte in exp_be {
            for bit in (0..8).rev() {
                acc = mont_mul(&acc, &acc);
                if (byte >> bit) & 1 == 1 {
                    acc = mont_mul(&acc, base);
                }
            }
        }
        acc
    }

    fn inverse(v_m: &Fp) -> [u8; 48] {
        limbs_to_be(&pow_mont(v_m, &exp_inverse()))
    }

    // The other square root of gx: an equally valid witness that the sign
    // correction must resolve to the same output point.
    pub fn flip_first_sqrt(blob: &[u8]) -> Vec<u8> {
        let y = wit48(&blob[49..97]).unwrap();
        let mut out = blob[..49].to_vec();
        out.extend_from_slice(&limbs_to_be(&neg_mod(&y)));
        out.extend_from_slice(&blob[97..]);
        out
    }

    pub fn generate_nu(msg: &[u8]) -> Vec<u8> {
        let elem = hash_to_field_nu(crate::dst::G1_NU, msg);
        let (blob_map, point) = map_blob(&elem);
        let mut blob = blob_map;
        let x_den = horner(&ISO11_XDEN, &point.x);
        let y_den = horner(&ISO11_YDEN, &point.x);
        blob.extend_from_slice(&inverse(&x_den));
        blob.extend_from_slice(&inverse(&y_den));
        blob.extend_from_slice(msg);
        blob
    }

    /// One map's witness blob plus the mapped E' point.
    fn map_blob(elem: &FieldElem) -> (Vec<u8>, PointPrime) {
        let usq = mont_sqr(&elem.mont);
        let xi_usq = mul_by_xi(&usq);
        let tv2 = add_mod(&mont_sqr(&xi_usq), &xi_usq);
        assert!(!is_zero(&tv2));
        let w_inv = inverse(&tv2);
        let inv_m = be_to_limbs(&w_inv);
        let x1 = mont_mul(&SSWU_C1_NEG_B_OVER_A, &add_mod(&R, &inv_m));
        let gx1 = gx_at(&x1);
        let legendre = pow_mont(&gx1, &exp_legendre());
        let (flag, x, gx) = if is_zero(&gx1) || legendre == R {
            (0u8, x1, gx1)
        } else {
            let x2 = mont_mul(&xi_usq, &x1);
            (1u8, x2, gx_at(&x2))
        };
        let y_m = pow_mont(&gx, &exp_sqrt());
        assert_eq!(mont_sqr(&y_m), gx);
        let y_c = from_mont(&y_m);
        let mut y_final = y_c;
        if (y_final[0] & 1) != (elem.canonical[0] & 1) {
            y_final = neg_mod(&y_final);
        }
        let mut blob = vec![flag];
        blob.extend_from_slice(&w_inv);
        blob.extend_from_slice(&limbs_to_be(&y_c));
        (blob, PointPrime { x, y: to_mont(&y_final) })
    }

    pub fn generate(msg: &[u8]) -> Vec<u8> {
        let u = hash_to_field_folded(crate::dst::G1_RO, msg);
        let mut blob = Vec::with_capacity(W_TOTAL);
        let mut points = Vec::new();

        for elem in &u {
            let usq = mont_sqr(&elem.mont);
            let xi_usq = mul_by_xi(&usq);
            let tv2 = add_mod(&mont_sqr(&xi_usq), &xi_usq);
            assert!(!is_zero(&tv2));

            let w_inv = inverse(&tv2);
            let inv_m = be_to_limbs(&w_inv);
            let x1 = mont_mul(&SSWU_C1_NEG_B_OVER_A, &add_mod(&R, &inv_m));
            let gx1 = gx_at(&x1);

            // is_square(0) is true (matches blst / the on-chain zero guard).
            let legendre = pow_mont(&gx1, &exp_legendre());
            let (flag, x, gx) = if is_zero(&gx1) || legendre == R {
                (0u8, x1, gx1)
            } else {
                let x2 = mont_mul(&xi_usq, &x1);
                let gx2 = gx_at(&x2);
                (1u8, x2, gx2)
            };
            let y_m = pow_mont(&gx, &exp_sqrt());
            assert_eq!(mont_mul(&y_m, &y_m), gx);
            let y_c = from_mont(&y_m);

            let mut y_final = y_c;
            if (y_final[0] & 1) != (elem.canonical[0] & 1) {
                y_final = neg_mod(&y_final);
            }

            blob.push(flag);
            blob.extend_from_slice(&w_inv);
            blob.extend_from_slice(&limbs_to_be(&y_c));

            points.push(PointPrime { x, y: to_mont(&y_final) });
        }

        let dx = sub_mod(&points[1].x, &points[0].x);
        blob.extend_from_slice(&inverse(&dx));

        let inv_m = be_to_limbs(&inverse(&dx));
        let lambda = mont_mul(&sub_mod(&points[1].y, &points[0].y), &inv_m);
        let x3 = sub_mod(
            &sub_mod(&mont_mul(&lambda, &lambda), &points[0].x),
            &points[1].x,
        );
        let sum_x = x3;

        let x_den = horner(&ISO11_XDEN, &sum_x);
        let y_den = horner(&ISO11_YDEN, &sum_x);
        blob.extend_from_slice(&inverse(&x_den));
        blob.extend_from_slice(&inverse(&y_den));

        assert_eq!(blob.len(), W_TOTAL);
        blob
    }
}
