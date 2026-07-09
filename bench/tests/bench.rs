use bls12_381::{
    hash_to_curve::{ExpandMsgXmd, HashToCurve},
    G1Affine, G1Projective, G2Affine, G2Projective, Scalar,
};
use mollusk_svm::{program::loader_keys::LOADER_V3, result::InstructionResult, Mollusk};
use solana_instruction::Instruction;
use solana_pubkey::Pubkey;

const ID: Pubkey = Pubkey::new_from_array([7u8; 32]);

const DST_G2: &[u8] = b"BLS_SIG_BLS12381G2_XMD:SHA-256_SSWU_RO_POP_";
const DST_G1: &[u8] = b"BLS_SIG_BLS12381G1_XMD:SHA-256_SSWU_RO_POP_";

const MESSAGE: &[u8] = b"tapedrive vote payload: epoch 42, slot 1337, snapshot root cafebabe";

fn mollusk() -> Mollusk {
    let elf = std::fs::read(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../target/deploy/bls381_bench.so"
    ))
    .expect("build the program first: cd program && cargo build-sbf");
    let mut mollusk = Mollusk::default();
    mollusk.add_program_with_loader_and_elf(&ID, &LOADER_V3, &elf);
    mollusk.compute_budget.compute_unit_limit = 500_000_000;
    mollusk
}

fn run(mollusk: &Mollusk, tag: u8, payload: &[u8]) -> InstructionResult {
    let mut data = vec![tag];
    data.extend_from_slice(payload);
    let instruction = Instruction::new_with_bytes(ID, &data, vec![]);
    mollusk.process_instruction(&instruction, &[])
}

fn cu(mollusk: &Mollusk, tag: u8, payload: &[u8], label: &str) -> u64 {
    let result = run(mollusk, tag, payload);
    assert!(
        !result.program_result.is_err(),
        "{label} failed: {:?}",
        result.program_result
    );
    println!("{label}: {} CU", result.compute_units_consumed);
    result.compute_units_consumed
}

#[test]
fn bench_hash_to_curve_pipeline() {
    let mollusk = mollusk();

    let full_g2 = {
        let result = run(&mollusk, 0, MESSAGE);
        assert!(!result.program_result.is_err(), "hash_to_g2 failed: {:?}", result.program_result);
        let expected = G2Affine::from(
            <G2Projective as HashToCurve<ExpandMsgXmd<sha2::Sha256>>>::hash_to_curve(
                MESSAGE, DST_G2,
            ),
        )
        .to_compressed();
        assert_eq!(result.return_data, expected.to_vec(), "hash_to_g2 output mismatch");
        println!("hash_to_g2 full (with compress): {} CU", result.compute_units_consumed);
        result.compute_units_consumed
    };

    let full_g1 = {
        let result = run(&mollusk, 1, MESSAGE);
        assert!(!result.program_result.is_err(), "hash_to_g1 failed: {:?}", result.program_result);
        let expected = G1Affine::from(
            <G1Projective as HashToCurve<ExpandMsgXmd<sha2::Sha256>>>::hash_to_curve(
                MESSAGE, DST_G1,
            ),
        )
        .to_compressed();
        assert_eq!(result.return_data, expected.to_vec(), "hash_to_g1 output mismatch");
        println!("hash_to_g1 full (with compress): {} CU", result.compute_units_consumed);
        result.compute_units_consumed
    };

    let field = cu(&mollusk, 2, MESSAGE, "hash_to_field only");
    let mapped = cu(&mollusk, 3, MESSAGE, "field + 2x map_to_curve + add");
    let cleared = cu(&mollusk, 4, MESSAGE, "field + map + add + clear_h");

    println!();
    println!("phase breakdown (G2):");
    println!("  hash_to_field:      {field} CU");
    println!("  2x map_to_curve:    {} CU", mapped.saturating_sub(field));
    println!("  clear_cofactor:     {} CU", cleared.saturating_sub(mapped));
    println!("  affine + compress:  {} CU", full_g2.saturating_sub(cleared));
    println!("  total G2:           {full_g2} CU");
    println!("  total G1:           {full_g1} CU");
}

#[test]
fn bench_hash_to_g2_matches_blst() {
    let mollusk = mollusk();
    let result = run(&mollusk, 0, MESSAGE);
    assert!(!result.program_result.is_err());

    let mut point = blst::blst_p2::default();
    let mut compressed = [0u8; 96];
    unsafe {
        blst::blst_hash_to_g2(
            &mut point,
            MESSAGE.as_ptr(),
            MESSAGE.len(),
            DST_G2.as_ptr(),
            DST_G2.len(),
            std::ptr::null(),
            0,
        );
        blst::blst_p2_compress(compressed.as_mut_ptr(), &point);
    }
    assert_eq!(result.return_data, compressed.to_vec(), "SBF output differs from blst");
}

#[test]
fn bench_syscalls() {
    let mollusk = mollusk();

    let g2_gen = G2Affine::generator().to_uncompressed();
    let g1_gen = G1Affine::generator().to_uncompressed();

    // Validate.
    let validate_g2 = run(&mollusk, 10, &g2_gen);
    println!(
        "g2 validate: rc={:?} {} CU",
        validate_g2.return_data.first(),
        validate_g2.compute_units_consumed
    );
    let validate_g1 = run(&mollusk, 18, &g1_gen);
    println!(
        "g1 validate: rc={:?} {} CU",
        validate_g1.return_data.first(),
        validate_g1.compute_units_consumed
    );

    // Add: gen + gen == 2*gen.
    let mut payload = g2_gen.to_vec();
    payload.extend_from_slice(&g2_gen);
    let add_g2 = run(&mollusk, 11, &payload);
    let expected = G2Affine::from(G2Projective::generator().double()).to_uncompressed();
    println!(
        "g2 add: rc={:?} match={} {} CU",
        add_g2.return_data.first(),
        add_g2.return_data.get(1..) == Some(expected.as_ref()),
        add_g2.compute_units_consumed
    );

    let mut payload = g1_gen.to_vec();
    payload.extend_from_slice(&g1_gen);
    let add_g1 = run(&mollusk, 15, &payload);
    let expected = G1Affine::from(G1Projective::generator().double()).to_uncompressed();
    println!(
        "g1 add: rc={:?} match={} {} CU",
        add_g1.return_data.first(),
        add_g1.return_data.get(1..) == Some(expected.as_ref()),
        add_g1.compute_units_consumed
    );

    // Mul: 7 * gen, scalar big-endian.
    let scalar_be = {
        let mut b = [0u8; 32];
        b[31] = 7;
        b
    };
    let mut payload = scalar_be.to_vec();
    payload.extend_from_slice(&g2_gen);
    let mul_g2 = run(&mollusk, 12, &payload);
    let expected =
        G2Affine::from(G2Projective::generator() * Scalar::from(7u64)).to_uncompressed();
    println!(
        "g2 mul: rc={:?} match={} {} CU",
        mul_g2.return_data.first(),
        mul_g2.return_data.get(1..) == Some(expected.as_ref()),
        mul_g2.compute_units_consumed
    );

    let mut payload = scalar_be.to_vec();
    payload.extend_from_slice(&g1_gen);
    let mul_g1 = run(&mollusk, 16, &payload);
    let expected =
        G1Affine::from(G1Projective::generator() * Scalar::from(7u64)).to_uncompressed();
    println!(
        "g1 mul: rc={:?} match={} {} CU",
        mul_g1.return_data.first(),
        mul_g1.return_data.get(1..) == Some(expected.as_ref()),
        mul_g1.compute_units_consumed
    );

    // Decompress.
    let decompress_g2 = run(&mollusk, 13, &G2Affine::generator().to_compressed());
    println!(
        "g2 decompress: rc={:?} match={} {} CU",
        decompress_g2.return_data.first(),
        decompress_g2.return_data.get(1..) == Some(g2_gen.as_ref()),
        decompress_g2.compute_units_consumed
    );
    let decompress_g1 = run(&mollusk, 17, &G1Affine::generator().to_compressed());
    println!(
        "g1 decompress: rc={:?} match={} {} CU",
        decompress_g1.return_data.first(),
        decompress_g1.return_data.get(1..) == Some(g1_gen.as_ref()),
        decompress_g1.compute_units_consumed
    );

    // Pairing, 1 and 2 pairs.
    let mut payload = g1_gen.to_vec();
    payload.extend_from_slice(&g2_gen);
    let pair_one = run(&mollusk, 14, &payload);
    println!(
        "pairing 1 pair: rc={:?} {} CU",
        pair_one.return_data.first(),
        pair_one.compute_units_consumed
    );

    let mut payload = g1_gen.to_vec();
    payload.extend_from_slice(&g1_gen);
    payload.extend_from_slice(&g2_gen);
    payload.extend_from_slice(&g2_gen);
    let pair_two = run(&mollusk, 14, &payload);
    println!(
        "pairing 2 pairs: rc={:?} {} CU",
        pair_two.return_data.first(),
        pair_two.compute_units_consumed
    );

    // big_mod_exp with 48-byte operands: base^exp mod p.
    let p_be = hex::decode(
        "1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaab",
    )
    .unwrap();
    let mut payload = vec![0u8; 48];
    payload[47] = 5;
    let mut exp = vec![0u8; 48];
    exp[47] = 3;
    payload.extend_from_slice(&exp);
    payload.extend_from_slice(&p_be);
    let modexp = run(&mollusk, 20, &payload);
    println!(
        "big_mod_exp 48B: rc={:?} out_tail={:?} {} CU",
        modexp.return_data.first(),
        modexp.return_data.last(),
        modexp.compute_units_consumed
    );
}

#[test]
fn bench_u128_mac_loop() {
    let mollusk = mollusk();
    let base = run(&mollusk, 21, &0u64.to_le_bytes());
    let loop_100k = run(&mollusk, 21, &100_000u64.to_le_bytes());
    assert!(!base.program_result.is_err());
    assert!(!loop_100k.program_result.is_err());
    let per_op = (loop_100k.compute_units_consumed - base.compute_units_consumed) as f64 / 100_000.0;
    println!(
        "u128 mac: base={} loop={} per_op={:.2} CU",
        base.compute_units_consumed, loop_100k.compute_units_consumed, per_op
    );
}

#[test]
fn bench_syscall_assisted_hash_to_g1() {
    use bls12_381::hash_to_curve::{HashToField, MapToCurve};

    let mollusk = mollusk();

    let field = run(&mollusk, 30, MESSAGE);
    assert!(!field.program_result.is_err(), "stage field: {:?}", field.program_result);

    let mapped = run(&mollusk, 31, MESSAGE);
    assert!(!mapped.program_result.is_err(), "stage maps: {:?}", mapped.program_result);

    let iso = run(&mollusk, 32, MESSAGE);
    assert!(!iso.program_result.is_err(), "stage iso: {:?}", iso.program_result);

    let full = run(&mollusk, 33, MESSAGE);
    assert!(!full.program_result.is_err(), "stage full: {:?}", full.program_result);

    // Reference: sum of the two mapped points before cofactor clearing.
    type F = <G1Projective as MapToCurve>::Field;
    let mut u = [F::default(); 2];
    F::hash_to_field::<ExpandMsgXmd<sha2::Sha256>>(MESSAGE, DST_G1, &mut u);
    let sum = G1Projective::map_to_curve(&u[0]) + G1Projective::map_to_curve(&u[1]);
    let expected_uncleared = G1Affine::from(sum).to_uncompressed();
    assert_eq!(
        iso.return_data,
        expected_uncleared.to_vec(),
        "pre-clearing point differs from zkcrypto"
    );

    // Reference: full hash_to_curve, zkcrypto and blst.
    let expected_full = G1Affine::from(
        <G1Projective as HashToCurve<ExpandMsgXmd<sha2::Sha256>>>::hash_to_curve(MESSAGE, DST_G1),
    )
    .to_uncompressed();
    assert_eq!(full.return_data, expected_full.to_vec(), "final point differs from zkcrypto");

    let mut point = blst::blst_p1::default();
    let mut serialized = [0u8; 96];
    unsafe {
        blst::blst_hash_to_g1(
            &mut point,
            MESSAGE.as_ptr(),
            MESSAGE.len(),
            DST_G1.as_ptr(),
            DST_G1.len(),
            std::ptr::null(),
            0,
        );
        blst::blst_p1_serialize(serialized.as_mut_ptr(), &point);
    }
    assert_eq!(full.return_data, serialized.to_vec(), "final point differs from blst");

    let f = field.compute_units_consumed;
    let m = mapped.compute_units_consumed;
    let i = iso.compute_units_consumed;
    let t = full.compute_units_consumed;
    println!();
    println!("syscall-assisted hash_to_G1 (min-sig):");
    println!("  hash_to_field:          {f} CU");
    println!("  2x sswu map:            {} CU", m.saturating_sub(f));
    println!("  E' add + iso-11:        {} CU", i.saturating_sub(m));
    println!("  clear_cofactor + check: {} CU", t.saturating_sub(i));
    println!("  TOTAL:                  {t} CU");
}

#[test]
fn bench_witness_hash_to_g1() {
    let mollusk = mollusk();

    let witnesses = bls381_hash::witness::g1::generate(MESSAGE);
    let mut payload = witnesses.clone();
    payload.extend_from_slice(MESSAGE);

    let result = run(&mollusk, 40, &payload);
    assert!(
        !result.program_result.is_err(),
        "witnessed hash_to_g1 failed: {:?}",
        result.program_result
    );

    let mut point = blst::blst_p1::default();
    let mut serialized = [0u8; 96];
    unsafe {
        blst::blst_hash_to_g1(
            &mut point,
            MESSAGE.as_ptr(),
            MESSAGE.len(),
            DST_G1.as_ptr(),
            DST_G1.len(),
            std::ptr::null(),
            0,
        );
        blst::blst_p1_serialize(serialized.as_mut_ptr(), &point);
    }
    assert_eq!(result.return_data, serialized.to_vec(), "differs from blst");

    println!(
        "witness-assisted hash_to_G1: {} CU ({} witness bytes)",
        result.compute_units_consumed,
        witnesses.len()
    );

    // corrupted witness must abort, not produce a different point
    let mut bad = payload.clone();
    bad[60] ^= 1;
    let rejected = run(&mollusk, 40, &bad);
    assert!(rejected.program_result.is_err(), "corrupt witness was accepted");
}

#[test]
fn bench_witness_svdw_hash_to_g1() {
    let mollusk = mollusk();

    let witnesses = bls381_hash::witness::g1_svdw::generate(MESSAGE);
    let mut payload = witnesses.clone();
    payload.extend_from_slice(MESSAGE);

    let result = run(&mollusk, 42, &payload);
    assert!(
        !result.program_result.is_err(),
        "witnessed svdw hash_to_g1 failed: {:?}",
        result.program_result
    );

    // Reference: host-side pre-clearing sum, effective cofactor applied
    // through zkcrypto scalar multiplication.
    let pre = bls381_hash::witness::g1_svdw::reference_preclear(MESSAGE);
    let aff = Option::<G1Affine>::from(G1Affine::from_uncompressed_unchecked(&pre))
        .expect("reference point parses");
    let expected = G1Affine::from(G1Projective::from(aff) * Scalar::from(0xd201000000010001u64))
        .to_uncompressed();
    assert_eq!(result.return_data, expected.to_vec(), "differs from host reference");

    // Different map, different suite: must NOT match the SSWU hash.
    let mut point = blst::blst_p1::default();
    let mut sswu = [0u8; 96];
    unsafe {
        blst::blst_hash_to_g1(
            &mut point,
            MESSAGE.as_ptr(),
            MESSAGE.len(),
            DST_G1.as_ptr(),
            DST_G1.len(),
            std::ptr::null(),
            0,
        );
        blst::blst_p1_serialize(sswu.as_mut_ptr(), &point);
    }
    assert_ne!(result.return_data, sswu.to_vec(), "svdw output cannot equal the sswu suite");

    println!(
        "witness-assisted SvdW hash_to_G1: {} CU ({} witness bytes)",
        result.compute_units_consumed,
        witnesses.len()
    );

    // the other square root is an equally valid witness and must not
    // steer the output
    let flipped = bls381_hash::witness::g1_svdw::flip_first_sqrt(&witnesses);
    let mut alt = flipped;
    alt.extend_from_slice(MESSAGE);
    let same = run(&mollusk, 42, &alt);
    assert!(!same.program_result.is_err(), "flipped root rejected");
    assert_eq!(same.return_data, result.return_data, "flipped root changed the point");

    // corrupted witness must abort, not produce a different point
    let mut bad = payload.clone();
    let dx_last = witnesses.len() - 1;
    bad[dx_last] ^= 1;
    let rejected = run(&mollusk, 42, &bad);
    assert!(rejected.program_result.is_err(), "corrupt witness was accepted");
}

#[test]
fn bench_witness_hash_to_g2() {
    let mollusk = mollusk();

    let witnesses = bls381_hash::witness::g2::generate(MESSAGE);
    let mut payload = witnesses.clone();
    payload.extend_from_slice(MESSAGE);

    let result = run(&mollusk, 41, &payload);
    assert!(
        !result.program_result.is_err(),
        "witnessed hash_to_g2 failed: {:?}",
        result.program_result
    );

    let mut point = blst::blst_p2::default();
    let mut serialized = [0u8; 192];
    unsafe {
        blst::blst_hash_to_g2(
            &mut point,
            MESSAGE.as_ptr(),
            MESSAGE.len(),
            DST_G2.as_ptr(),
            DST_G2.len(),
            std::ptr::null(),
            0,
        );
        blst::blst_p2_serialize(serialized.as_mut_ptr(), &point);
    }
    assert_eq!(result.return_data, serialized.to_vec(), "differs from blst");

    println!(
        "witness-assisted hash_to_G2: {} CU ({} witness bytes)",
        result.compute_units_consumed,
        witnesses.len()
    );

    let mut bad = payload.clone();
    bad[120] ^= 1;
    let rejected = run(&mollusk, 41, &bad);
    assert!(rejected.program_result.is_err(), "corrupt witness was accepted");
}

// No witness byte can steer the output: every single-bit corruption of the G2
// witness must abort or reproduce the exact same point. Also covers the branch
// flag range, canonical-form parsing, and message binding.
#[test]
fn witness_g2_soundness() {
    let mollusk = mollusk();

    let witness = bls381_hash::witness::g2::generate(MESSAGE);
    let mut payload = witness.clone();
    payload.extend_from_slice(MESSAGE);

    let good = run(&mollusk, 41, &payload);
    assert!(!good.program_result.is_err(), "honest witness rejected");
    let truth = good.return_data.clone();

    for i in 0..witness.len() {
        let mut bad = payload.clone();
        bad[i] ^= 1;
        let r = run(&mollusk, 41, &bad);
        if !r.program_result.is_err() {
            assert_eq!(r.return_data, truth, "witness byte {i} steered the output");
        }
    }

    // a branch flag above 1 is non-canonical
    for flag in [0usize, 97] {
        let mut bad = payload.clone();
        bad[flag] = 2;
        assert!(run(&mollusk, 41, &bad).program_result.is_err(), "flag=2 at {flag} accepted");
    }

    // a witness limb at or above the modulus must be rejected by the parser
    let mut oob = payload.clone();
    for byte in oob[1..49].iter_mut() {
        *byte = 0xff;
    }
    assert!(run(&mollusk, 41, &oob).program_result.is_err(), "out-of-range witness accepted");

    // the witness is bound to the message it was generated for
    let mut replay = witness.clone();
    replay.extend_from_slice(b"a different snapshot vote payload");
    assert!(run(&mollusk, 41, &replay).program_result.is_err(), "cross-message replay accepted");

    // the other square root is an equally valid witness and must not steer
    let mut alt = bls381_hash::witness::g2::flip_first_sqrt(&witness);
    alt.extend_from_slice(MESSAGE);
    let same = run(&mollusk, 41, &alt);
    assert!(!same.program_result.is_err(), "flipped root rejected");
    assert_eq!(same.return_data, truth, "flipped root changed the point");
}

// The G1 (min-sig) counterpart of the G2 soundness sweep.
#[test]
fn witness_g1_soundness() {
    let mollusk = mollusk();

    let witness = bls381_hash::witness::g1::generate(MESSAGE);
    let mut payload = witness.clone();
    payload.extend_from_slice(MESSAGE);

    let good = run(&mollusk, 40, &payload);
    assert!(!good.program_result.is_err(), "honest witness rejected");
    let truth = good.return_data.clone();

    for i in 0..witness.len() {
        let mut bad = payload.clone();
        bad[i] ^= 1;
        let r = run(&mollusk, 40, &bad);
        if !r.program_result.is_err() {
            assert_eq!(r.return_data, truth, "witness byte {i} steered the output");
        }
    }

    for flag in [0usize, 97] {
        let mut bad = payload.clone();
        bad[flag] = 2;
        assert!(run(&mollusk, 40, &bad).program_result.is_err(), "flag=2 at {flag} accepted");
    }

    let mut oob = payload.clone();
    for byte in oob[1..49].iter_mut() {
        *byte = 0xff;
    }
    assert!(run(&mollusk, 40, &oob).program_result.is_err(), "out-of-range witness accepted");

    let mut replay = witness.clone();
    replay.extend_from_slice(b"a different snapshot vote payload");
    assert!(run(&mollusk, 40, &replay).program_result.is_err(), "cross-message replay accepted");

    // the other square root is an equally valid witness and must not steer
    let mut alt = bls381_hash::witness::g1::flip_first_sqrt(&witness);
    alt.extend_from_slice(MESSAGE);
    let same = run(&mollusk, 40, &alt);
    assert!(!same.program_result.is_err(), "flipped root rejected");
    assert_eq!(same.return_data, truth, "flipped root changed the point");
}

#[test]
fn bench_witness_svdw_hash_to_g2() {
    use bls12_381::hash_to_curve::MapToCurve;

    let mollusk = mollusk();

    let witnesses = bls381_hash::witness::g2_svdw::generate(MESSAGE);
    let mut payload = witnesses.clone();
    payload.extend_from_slice(MESSAGE);

    let result = run(&mollusk, 43, &payload);
    assert!(
        !result.program_result.is_err(),
        "witnessed svdw hash_to_g2 failed: {:?}",
        result.program_result
    );

    // Reference: host-side pre-clearing sum, cofactor cleared through
    // zkcrypto's clear_h (same Budroni-Pintore construction).
    let pre = bls381_hash::witness::g2_svdw::reference_preclear(MESSAGE);
    let aff = Option::<G2Affine>::from(G2Affine::from_uncompressed_unchecked(&pre))
        .expect("reference point parses");
    let expected = G2Affine::from(G2Projective::from(aff).clear_h()).to_uncompressed();
    assert_eq!(result.return_data, expected.to_vec(), "differs from host reference");

    // Different map, different suite: must NOT match the SSWU hash.
    let mut point = blst::blst_p2::default();
    let mut sswu = [0u8; 192];
    unsafe {
        blst::blst_hash_to_g2(
            &mut point,
            MESSAGE.as_ptr(),
            MESSAGE.len(),
            DST_G2.as_ptr(),
            DST_G2.len(),
            std::ptr::null(),
            0,
        );
        blst::blst_p2_serialize(sswu.as_mut_ptr(), &point);
    }
    assert_ne!(result.return_data, sswu.to_vec(), "svdw output cannot equal the sswu suite");

    println!(
        "witness-assisted SvdW hash_to_G2: {} CU ({} witness bytes)",
        result.compute_units_consumed,
        witnesses.len()
    );

    // the other square root is an equally valid witness and must not
    // steer the output
    let flipped = bls381_hash::witness::g2_svdw::flip_first_sqrt(&witnesses);
    let mut alt = flipped;
    alt.extend_from_slice(MESSAGE);
    let same = run(&mollusk, 43, &alt);
    assert!(!same.program_result.is_err(), "flipped root rejected");
    assert_eq!(same.return_data, result.return_data, "flipped root changed the point");

    // corrupted witness must abort, not produce a different point
    let mut bad = payload.clone();
    let dx_last = witnesses.len() - 1;
    bad[dx_last] ^= 1;
    let rejected = run(&mollusk, 43, &bad);
    assert!(rejected.program_result.is_err(), "corrupt witness was accepted");
}

#[test]
fn bench_witness_nu_encode() {
    let mollusk = mollusk();
    const DST_G1_NU: &[u8] = b"BLS_SIG_BLS12381G1_XMD:SHA-256_SSWU_NU_POP_";
    const DST_G2_NU: &[u8] = b"BLS_SIG_BLS12381G2_XMD:SHA-256_SSWU_NU_POP_";

    let payload = bls381_hash::witness::g1::generate_nu(MESSAGE);
    let r1 = run(&mollusk, 44, &payload);
    assert!(!r1.program_result.is_err(), "nu g1: {:?}", r1.program_result);
    let mut pt = blst::blst_p1::default();
    let mut ser = [0u8; 96];
    unsafe {
        blst::blst_encode_to_g1(&mut pt, MESSAGE.as_ptr(), MESSAGE.len(), DST_G1_NU.as_ptr(), DST_G1_NU.len(), std::ptr::null(), 0);
        blst::blst_p1_serialize(ser.as_mut_ptr(), &pt);
    }
    assert_eq!(r1.return_data, ser.to_vec(), "nu g1 differs from blst encode");
    println!("witness-assisted NU encode_to_G1: {} CU ({} witness bytes)", r1.compute_units_consumed, payload.len() - MESSAGE.len());
    let mut bad = payload.clone();
    bad[60] ^= 1;
    assert!(run(&mollusk, 44, &bad).program_result.is_err(), "corrupt nu g1 accepted");

    let payload = bls381_hash::witness::g2::generate_nu(MESSAGE);
    let r2 = run(&mollusk, 45, &payload);
    assert!(!r2.program_result.is_err(), "nu g2: {:?}", r2.program_result);
    let mut pt = blst::blst_p2::default();
    let mut ser = [0u8; 192];
    unsafe {
        blst::blst_encode_to_g2(&mut pt, MESSAGE.as_ptr(), MESSAGE.len(), DST_G2_NU.as_ptr(), DST_G2_NU.len(), std::ptr::null(), 0);
        blst::blst_p2_serialize(ser.as_mut_ptr(), &pt);
    }
    assert_eq!(r2.return_data, ser.to_vec(), "nu g2 differs from blst encode");
    println!("witness-assisted NU encode_to_G2: {} CU ({} witness bytes)", r2.compute_units_consumed, payload.len() - MESSAGE.len());
    let mut bad = payload.clone();
    bad[120] ^= 1;
    assert!(run(&mollusk, 45, &bad).program_result.is_err(), "corrupt nu g2 accepted");
}

// Correctness guard: the bare Montgomery reduction (from_mont) against the
// general multiply, plus the adapted iso-11 chain against Horner. Host-side.
#[test]
fn field_arithmetic_selftest() {
    bls381_hash::witness::g1::iso11_adapted_selftest();
    bls381_hash::witness::g1::redc_selftest();
}

#[test]
fn bench_min_pk_verify_end_to_end() {
    use blst::min_pk::{AggregatePublicKey, AggregateSignature, PublicKey, SecretKey, Signature};

    let mollusk = mollusk();

    let keys: Vec<SecretKey> = (0..20u8)
        .map(|i| {
            let ikm = [i + 1; 32];
            SecretKey::key_gen(&ikm, &[]).unwrap()
        })
        .collect();
    let pks: Vec<PublicKey> = keys.iter().map(|s| s.sk_to_pk()).collect();
    let all_refs: Vec<&PublicKey> = pks.iter().collect();
    let agg_all = AggregatePublicKey::aggregate(&all_refs, false)
        .unwrap()
        .to_public_key();

    let witness = bls381_hash::witness::g2::generate(MESSAGE);

    for k in [14usize, 20] {
        let sigs: Vec<Signature> = keys[..k]
            .iter()
            .map(|s| s.sign(MESSAGE, DST_G2, &[]))
            .collect();
        let sig_refs: Vec<&Signature> = sigs.iter().collect();
        let agg_sig = AggregateSignature::aggregate(&sig_refs, false)
            .unwrap()
            .to_signature();

        let mut payload = vec![(20 - k) as u8];
        payload.extend_from_slice(&agg_all.serialize());
        payload.extend_from_slice(&G1Affine::generator().to_uncompressed());
        payload.extend_from_slice(&agg_sig.compress());
        for pk in &pks[k..] {
            payload.extend_from_slice(&pk.compress());
        }
        payload.extend_from_slice(&witness);
        payload.extend_from_slice(MESSAGE);

        let result = run(&mollusk, 51, &payload);
        assert!(
            !result.program_result.is_err(),
            "min-pk verify failed at k={k}: {:?}",
            result.program_result
        );
        println!(
            "min-pk end-to-end verify k={k}: {} CU",
            result.compute_units_consumed
        );

        // tampered signature must fail
        let mut bad = payload.clone();
        bad[1 + 96 + 96 + 10] ^= 1;
        let rejected = run(&mollusk, 51, &bad);
        assert!(rejected.program_result.is_err(), "forged min-pk verify accepted at k={k}");
    }
}
