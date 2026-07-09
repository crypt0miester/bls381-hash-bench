//! BLS12-381 base-field (Fp) arithmetic: 32-bit CIOS Montgomery plus the
//! big_mod_exp helpers used by the host and the modexp-assisted path.

use solana_program_error::ProgramError;

use crate::consts_g1::{INV, MODULUS, R2};

pub(crate) type Fp = [u64; 6];

// Split a tag payload into (witness, message), erroring if it is too short.
pub(crate) fn split_witness(payload: &[u8], total: usize) -> Result<(&[u8], &[u8]), ProgramError> {
    if payload.len() < total {
        return Err(ProgramError::InvalidInstructionData);
    }
    Ok(payload.split_at(total))
}

const ONE: Fp = [1, 0, 0, 0, 0, 0];

#[inline(always)]
fn mul_64x64(a: u64, b: u64) -> (u64, u64) {
    let a_lo = a & 0xffff_ffff;
    let a_hi = a >> 32;
    let b_lo = b & 0xffff_ffff;
    let b_hi = b >> 32;
    let p0 = a_lo * b_lo;
    let p1 = a_lo * b_hi;
    let p2 = a_hi * b_lo;
    let p3 = a_hi * b_hi;
    let mid = (p0 >> 32) + (p1 & 0xffff_ffff) + (p2 & 0xffff_ffff);
    let lo = (p0 & 0xffff_ffff) | (mid << 32);
    let hi = p3 + (p1 >> 32) + (p2 >> 32) + (mid >> 32);
    (lo, hi)
}

#[inline(always)]
fn mac(acc: u64, a: u64, b: u64, carry: u64) -> (u64, u64) {
    let (lo, hi) = mul_64x64(a, b);
    let (lo, c1) = lo.overflowing_add(acc);
    let (lo, c2) = lo.overflowing_add(carry);
    (lo, hi + c1 as u64 + c2 as u64)
}

#[inline(always)]
fn adc(a: u64, b: u64, carry: u64) -> (u64, u64) {
    let (s, c1) = a.overflowing_add(b);
    let (s, c2) = s.overflowing_add(carry);
    (s, c1 as u64 + c2 as u64)
}

#[inline(always)]
fn sbb(a: u64, b: u64, borrow: u64) -> (u64, u64) {
    let (d, b1) = a.overflowing_sub(b);
    let (d, b2) = d.overflowing_sub(borrow);
    (d, b1 as u64 + b2 as u64)
}

pub(crate) fn geq(a: &Fp, b: &Fp) -> bool {
    for i in (0..6).rev() {
        if a[i] > b[i] {
            return true;
        }
        if a[i] < b[i] {
            return false;
        }
    }
    true
}

pub(crate) fn sub_nocheck(a: &Fp, b: &Fp) -> Fp {
    let mut r = [0u64; 6];
    let mut borrow = 0u64;
    for i in 0..6 {
        let (d, br) = sbb(a[i], b[i], borrow);
        r[i] = d;
        borrow = br;
    }
    r
}

#[inline(always)]
pub(crate) fn add_mod(a: &Fp, b: &Fp) -> Fp {
    let mut r = [0u64; 6];
    let mut carry = 0u64;
    for i in 0..6 {
        let (s, c) = adc(a[i], b[i], carry);
        r[i] = s;
        carry = c;
    }
    if carry != 0 || geq(&r, &MODULUS) {
        r = sub_nocheck(&r, &MODULUS);
    }
    r
}

#[inline(always)]
pub(crate) fn sub_mod(a: &Fp, b: &Fp) -> Fp {
    if geq(a, b) {
        sub_nocheck(a, b)
    } else {
        add_carryless(&sub_nocheck(a, b))
    }
}

pub(crate) fn add_carryless(r: &Fp) -> Fp {
    // wrapped subtraction result plus p restores the field value
    let mut out = [0u64; 6];
    let mut carry = 0u64;
    for i in 0..6 {
        let (s, c) = adc(r[i], MODULUS[i], carry);
        out[i] = s;
        carry = c;
    }
    out
}

pub(crate) fn neg_mod(a: &Fp) -> Fp {
    if a.iter().all(|&l| l == 0) {
        return [0u64; 6];
    }
    sub_nocheck(&MODULUS, a)
}

pub(crate) fn is_zero(a: &Fp) -> bool {
    a.iter().all(|&l| l == 0)
}

pub(crate) fn mont_mul(a: &Fp, b: &Fp) -> Fp {
    mont_mul_cios32(a, b)
}

/// Square through the general multiply; a dedicated squaring does not
/// pay for itself on sbpf.
pub(crate) fn mont_sqr(a: &Fp) -> Fp {
    mont_mul_cios32(a, a)
}


pub(crate) fn to_mont(a: &Fp) -> Fp {
    mont_mul(a, &R2)
}

/// `a * R^-1 mod p`, via the bare Montgomery reduction.
pub(crate) fn from_mont(a: &Fp) -> Fp {
    mont_redc_cios32(a)
}

pub(crate) fn limbs_to_be(a: &Fp) -> [u8; 48] {
    let mut out = [0u8; 48];
    for i in 0..6 {
        out[i * 8..i * 8 + 8].copy_from_slice(&a[5 - i].to_be_bytes());
    }
    out
}

pub(crate) fn be_to_limbs(b: &[u8; 48]) -> Fp {
    let mut r = [0u64; 6];
    for i in 0..6 {
        r[5 - i] = u64::from_be_bytes(b[i * 8..i * 8 + 8].try_into().unwrap());
    }
    r
}

pub(crate) fn shr1(a: &Fp) -> Fp {
    let mut r = [0u64; 6];
    for i in 0..6 {
        r[i] = a[i] >> 1;
        if i < 5 {
            r[i] |= a[i + 1] << 63;
        }
    }
    r
}

pub(crate) fn exp_inverse() -> [u8; 48] {
    limbs_to_be(&sub_nocheck(&MODULUS, &[2, 0, 0, 0, 0, 0]))
}

pub(crate) fn exp_legendre() -> [u8; 48] {
    limbs_to_be(&shr1(&sub_nocheck(&MODULUS, &ONE)))
}

pub(crate) fn exp_sqrt() -> [u8; 48] {
    let mut p1 = MODULUS;
    p1[0] += 1;
    limbs_to_be(&shr1(&shr1(&p1)))
}

#[cfg(target_os = "solana")]
pub(crate) mod sys {
    use solana_define_syscall::define_syscall;

    define_syscall!(fn sol_curve_validate_point(curve_id: u64, point_addr: *const u8, result: *mut u8) -> u64);
    define_syscall!(fn sol_curve_group_op(curve_id: u64, group_op: u64, left_input_addr: *const u8, right_input_addr: *const u8, result_point_addr: *mut u8) -> u64);
    define_syscall!(fn sol_big_mod_exp(params: *const u8, result: *mut u8) -> u64);
}

#[cfg(not(target_os = "solana"))]
#[allow(clippy::missing_safety_doc)]
pub(crate) mod sys {
    pub unsafe fn sol_curve_validate_point(_: u64, _: *const u8, _: *mut u8) -> u64 {
        unimplemented!()
    }
    pub unsafe fn sol_curve_group_op(_: u64, _: u64, _: *const u8, _: *const u8, _: *mut u8) -> u64 {
        unimplemented!()
    }
    pub unsafe fn sol_big_mod_exp(_: *const u8, _: *mut u8) -> u64 {
        unimplemented!()
    }
}

#[repr(C)]
struct BigModExpParams {
    base: *const u8,
    base_len: u64,
    exponent: *const u8,
    exponent_len: u64,
    modulus: *const u8,
    modulus_len: u64,
}

pub(crate) fn modexp_bytes(base: &[u8], exp: &[u8]) -> Result<[u8; 48], ProgramError> {
    let modulus = limbs_to_be(&MODULUS);
    let params = BigModExpParams {
        base: base.as_ptr(),
        base_len: base.len() as u64,
        exponent: exp.as_ptr(),
        exponent_len: exp.len() as u64,
        modulus: modulus.as_ptr(),
        modulus_len: 48,
    };
    let mut out = [0u8; 48];
    let rc = unsafe {
        sys::sol_big_mod_exp(
            &params as *const BigModExpParams as *const u8,
            out.as_mut_ptr(),
        )
    };
    if rc != 0 {
        return Err(ProgramError::InvalidInstructionData);
    }
    Ok(out)
}

/// Modular exponentiation of a canonical-form element, returning canonical.
pub(crate) fn modexp(base: &Fp, exp: &[u8; 48]) -> Result<Fp, ProgramError> {
    Ok(be_to_limbs(&modexp_bytes(&limbs_to_be(base), exp)?))
}

/// Inverse of a Montgomery-form element, returned in Montgomery form.
pub(crate) fn inverse_mont(a: &Fp, exp_inv: &[u8; 48]) -> Result<Fp, ProgramError> {
    Ok(to_mont(&modexp(&from_mont(a), exp_inv)?))
}

const INV32: u64 = INV & 0xffff_ffff;

const fn split32(x: &Fp) -> [u64; 12] {
    let mut out = [0u64; 12];
    let mut i = 0;
    while i < 6 {
        out[i * 2] = x[i] & 0xffff_ffff;
        out[i * 2 + 1] = x[i] >> 32;
        i += 1;
    }
    out
}

// The modulus never changes, so its 32-bit lanes belong in a constant rather
// than being rebuilt on every multiply.
const P32: [u64; 12] = split32(&MODULUS);

#[inline(always)]
fn mac32(acc: u64, a: u64, b: u64, carry: u64) -> (u64, u64) {
    let t = acc + a * b + carry;
    (t & 0xffff_ffff, t >> 32)
}

/// CIOS with 32-bit limbs: the multiply-accumulate needs no wide arithmetic.
// Indexed loops are load-bearing here: iterator rewrites measurably regress CU
// on sBPF, so the range-loop lint is silenced rather than followed. For the same
// reason the reduction round is kept inline (not shared with mont_redc_cios32):
// factoring it into an #[inline(always)] helper measurably regressed CU.
#[inline(always)]
#[allow(clippy::needless_range_loop)]
fn mont_mul_cios32(a: &Fp, b: &Fp) -> Fp {
    let a32 = split32(a);
    let b32 = split32(b);

    let mut t = [0u64; 14];
    for i in 0..12 {
        let ai = a32[i];
        let mut carry = 0u64;
        for j in 0..12 {
            let (lo, hi) = mac32(t[j], ai, b32[j], carry);
            t[j] = lo;
            carry = hi;
        }
        let s = t[12] + carry;
        t[12] = s & 0xffff_ffff;
        t[13] = s >> 32;

        let m = (t[0].wrapping_mul(INV32)) & 0xffff_ffff;
        let (_, mut carry) = mac32(t[0], m, P32[0], 0);
        for j in 1..12 {
            let (lo, hi) = mac32(t[j], m, P32[j], carry);
            t[j - 1] = lo;
            carry = hi;
        }
        let s = t[12] + carry;
        t[11] = s & 0xffff_ffff;
        t[12] = t[13] + (s >> 32);
        t[13] = 0;
    }

    let mut r = [0u64; 6];
    for i in 0..6 {
        r[i] = t[i * 2] | (t[i * 2 + 1] << 32);
    }
    if t[12] != 0 || geq(&r, &MODULUS) {
        r = sub_nocheck(&r, &MODULUS);
    }
    r
}

/// `from_mont`: `a * R^-1 mod p`. The Montgomery reduction step of a multiply,
/// without the product loop.
#[allow(clippy::needless_range_loop)]
pub(crate) fn mont_redc_cios32(a: &Fp) -> Fp {
    let a32 = split32(a);
    let mut t = [0u64; 14];
    for i in 0..12 {
        t[i] = a32[i];
    }
    for _ in 0..12 {
        let m = (t[0].wrapping_mul(INV32)) & 0xffff_ffff;
        let (_, mut carry) = mac32(t[0], m, P32[0], 0);
        for j in 1..12 {
            let (lo, hi) = mac32(t[j], m, P32[j], carry);
            t[j - 1] = lo;
            carry = hi;
        }
        let s = t[12] + carry;
        t[11] = s & 0xffff_ffff;
        t[12] = t[13] + (s >> 32);
        t[13] = 0;
    }
    let mut r = [0u64; 6];
    for i in 0..6 {
        r[i] = t[i * 2] | (t[i * 2 + 1] << 32);
    }
    if t[12] != 0 || geq(&r, &MODULUS) {
        r = sub_nocheck(&r, &MODULUS);
    }
    r
}

pub(crate) fn wit48(bytes: &[u8]) -> Result<Fp, ProgramError> {
    let arr: &[u8; 48] = bytes
        .try_into()
        .map_err(|_| ProgramError::InvalidInstructionData)?;
    let limbs = be_to_limbs(arr);
    if geq(&limbs, &MODULUS) {
        return Err(ProgramError::InvalidInstructionData);
    }
    Ok(limbs)
}

