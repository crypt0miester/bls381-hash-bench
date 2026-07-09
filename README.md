# bls381-hash

Witness-assisted RFC 9380 hash-to-curve for BLS12-381, for Solana SBF programs.
`no_std`, allocation-free on-chain, using only syscalls active on mainnet today,
no `big_mod_exp` (SIMD-0529) and no map-to-curve syscall. Everything expensive
(inverses, Legendre symbols, square roots) rides in as witness data and is
verified with a multiplication or two. Host-side witness generation ships in the
same crate.

```rust
use bls381_hash::{dst, hash_to_g1};

// on-chain: DST is a runtime parameter, the payload is the witness bytes
// followed by the message
let point = hash_to_g1(dst::G1_RO, payload)?; // Vec<u8>, the 96-byte G1 point
```

```rust
// off-chain (host): build the witness for a message
let witness = bls381_hash::witness::g1::generate(message);
```

## Features

| feature | pulls in |
|---|---|
| `ro` (default) | `g1-ro` + `g2-ro`, the blst-compatible pair |
| `g1-ro`, `g2-ro` | standard `_SSWU_RO_POP_` pipelines |
| `g1-nu`, `g2-nu` | RFC 9380 encode_to_curve variants |
| `g1-svdw`, `g2-svdw` | custom-suite SvdW variants (no isogeny) |
| `modexp` | big_mod_exp-assisted G1 path, needs SIMD-0529 |
| `full` | everything above |

The `lib/` crate is the product; `program/` is an SBF tag-dispatch fixture and
`bench/` the mollusk benchmark.

## Benchmark

Measured against the fixture with blst cross-checks.

Measured with mollusk 0.13.4 on the agave 4.0 stack, SBF v3.

| pipeline | CU | witness | compatibility |
|---|---|---|---|
| hash_to_G1 (RO, min-sig) | ~230k | 338 B | `_SSWU_RO_POP_`, byte-equal to blst |
| hash_to_G2 (RO, min-pk) | ~453k | 482 B | `_SSWU_RO_POP_`, byte-equal to blst |
| encode_to_G1 (NU) | ~172k | 193 B | `_SSWU_NU_POP_`, byte-equal to blst encode |
| encode_to_G2 (NU) | ~306k | 289 B | `_SSWU_NU_POP_`, byte-equal to blst encode |
| hash_to_G1 (SvdW) | ~156k | 242 to 434 B | custom suite |
| hash_to_G2 (SvdW) | ~411k | 482 to 866 B | custom suite |

For scale, a naive port of zkcrypto `bls12_381` costs 11.3M CU for G1 and
46.5M CU for G2, and a single 381-bit field multiplication bottoms out around
3.3k CU on sbpf.

The NU suites hash with a single map (RFC 9380 encode_to_curve). Note that
the CFRG BLS signature draft registers only hash_to_curve (RO) ciphersuites,
and RFC 9380 limits encode_to_curve to applications whose security analysis
does not rely on a random oracle. Using NU for BLS rests on the argument in
section 5 of Wahby-Boneh ([eprint 2019/403](https://eprint.iacr.org/2019/403))
and the BCI+10 reference there: hashing onto a constant fraction of the group
suffices for unforgeability. That makes it a deliberate protocol choice, and
the hash must never be reused for anything that actually needs a random
oracle.

## Approach

Field inverses, Legendre symbols, and square roots are expensive to compute
but cheap to verify, so the caller supplies each as instruction data and the
program checks it:

- sqrt: witness `y`, check `y^2 == gx`
- inverse: witness `t` in Montgomery form, check `v*t == 1`
- non-square (SvdW branches): witness `s` with `s^2 = xi * f(x)` for a fixed
  non-residue `xi`

A wrong witness aborts the instruction, and no witness can steer the output,
so the hash stays a pure function of the message.

## Optimizations that landed

- One sqrt witness per SSWU branch. `g(x2) = (Z u^2)^3 g(x1)` with `Z` a
  non-residue (section 4.1 of the paper), so a single root proves its own
  branch.
- Montgomery pair inversion. One witness `w = (a*b)^-1` pins both inverses of
  a same-stage pair.
- Knuth-adapted polynomial evaluation (TAOCP 4.6.4 preprocessing). The iso-11
  runs in 27 multiplications instead of 51, the iso-3 in 5 plus a square
  instead of 11. Constants were derived and expansion-checked offline.
- Small constants (`Z = 11`, `Z2 = -(2+i)`, `A' = 240i`, the `2^256` fold)
  are multiplied with addition chains instead of field muls.
- Bare Montgomery reduction out of Montgomery form. `from_mont` skips the
  product loop of a multiply by one, about half the multiply-accumulates.

## SvdW variant (custom suite)

A direct Shallue-van de Woestijne map onto the curve (section 3 of the paper,
`u0 = -3` on E1 and `u0 = -1` on E2), which skips the isogeny entirely.
Branch selection is proven rather than trusted: SvdW takes the smallest j
with `f(x_j)` square, so claiming branch j requires a non-squareness proof
for each earlier branch plus the sqrt witness for its own. The output matches
no standard suite, and the NU suites beat it on G2 anyway, so this one stays
a finding.

## Further optimizations

Open knobs, in rough order of interest:

- The modexp path (tags 30 to 33) runs hash_to_G1 in ~273k CU with zero
  witness bytes, against ~230k plus 338 B for the witnessed path. It needs
  `big_mod_exp` (SIMD-0529, merged but not active). Once 0529 activates, a
  transaction that is byte-bound rather than CU-bound should prefer it.
- The min-pk verify transaction is closer to byte-bound than CU-bound (513k
  of the 1.4M CU ceiling, but witness plus keys eat real transaction space).
  Batching all five G2 inverses behind one witness would save 192 B for
  roughly 17k CU. Not implemented; the right trade depends on the consumer.
- G2 cofactor clearing costs ~65k CU: roughly ~45k across ~140 g2 add syscalls
  and ~20k of psi/psi2 Fp2 multiplication. The Budroni-Pintore chain is the best
  known construction, so the syscall share is pricing rather than structure, but
  the endomorphism field work is real. The verify path feeds the hash into the
  pairing uncompressed, so no decompression cost hides there.

## Layout

- `lib/src/fp.rs`, `fp2.rs`: Fp and Fp2 arithmetic (32-bit CIOS Montgomery)
- `lib/src/g1.rs`, `g2.rs`: RO and NU pipelines, cofactor clearing, host witness generation
- `lib/src/{g1,g2}_svdw.rs`: SvdW variants
- `lib/src/consts_g1.rs`, `consts_g2.rs`: SSWU, isogeny, psi, and adapted constants
- `lib/src/lib.rs`: public API, feature gates, `dst` module
- `program/src/lib.rs`: SBF tag-dispatch fixture
- `bench/tests/`: mollusk benchmarks, blst cross-checks, soundness tests

## Build and run

```
cd program && cargo build-sbf --arch v3
cd ../bench && cargo test -- --nocapture
```

Requires the Solana platform tools. The standard suites assert byte-equality
with blst at every stage. SvdW checks against a host-side reference. A
corrupted witness must abort, and supplying the other square root must not
change the output point.

## Status

Experimental. The witnessed hash is novel enough to warrant a hostile review
before it is used anywhere consensus depends on it.

## License

MIT. The SSWU, isogeny, and psi constants in `lib/src/consts_g{1,2}.rs` were
extracted from zkcrypto [`bls12_381`](https://github.com/zkcrypto/bls12_381)
(MIT/Apache-2.0); the map constructions follow Wahby-Boneh, eprint 2019/403.
