# Reviewed rank-three coefficient cutoff

Incorporated 2026-09-27 from [Grok's report](../reviews/rank-three-cutoff-review.md).
The submitted report is preserved verbatim. This note states the accepted
theorem, corrections and implementation. Historical packs and manuscript
inventory remain unchanged; their null cutoff fields describe the saved run.

## Precise theorem and proof

Let G be a positive definite integral 3x3 Gram matrix, a,b binary primal
and dual columns, b nonzero and a^t b even. Let K be any integral column
basis for L0={x:x^t b even}, G0=K^t G K, d0=det G0, and N0 the least
positive integer clearing the denominators of G0^-1. Then

    M=16N0,
    mu=[SL2(Z):Gamma0(M)]=product_(ell^e || M) ell^(e-1)(ell+1),
    B=mu/8.

If every signed norm4 coefficient c_N is exactly zero for 0<=N<=B,
then the characteristic theta function is identically zero. The constant
term and the upper endpoint are included. Missing occupied shells are
invalid evidence; unoccupied indices have coefficient zero once the full
sphere through B has been enumerated. No cuspidality is assumed.

The reviewed [coset modularity](coset-modularity.md) puts
F(z)=sum c_N exp(2 pi i N z) in weight 3/2 on Gamma0(M), character
chi_(4d0), and Theta(tau)=i^(a^t b) F(tau/8). This is the exact signed
doubled-coset difference, with a trivial coset based on 2L0. The Fourier
index in z is precisely the engine's norm4, not norm4/8 or norm4/4.

Under gamma=[[r,s],[u,v]], its transformation factor is
chi_(4d0)(v) (u/v) epsilon_v^(-3) (uz+v)^(3/2).
The character and symbol have values +/-1 on the group, epsilon_v^4=1,
and the fourth power of the analytic factor is (uz+v)^6. Thus F^4 has
integral weight six and trivial character on Gamma0(M), including -I.
Its holomorphy in the upper half-plane is immediate.

For completeness at cusps, the weight-six slash by a cusp frame is
periodic with the cusp width and has polynomial growth as Im z grows,
uniformly in a period strip. For a negative Fourier index -m, its
coefficient is bounded by a polynomial times exp(-2 pi m Im z/w),
using the Fourier coefficient integral at height Im z. Letting this
height grow forces the coefficient to be zero. Hence every negative
index vanishes and F^4 is holomorphic at all cusps. Periodicity alone
would not exclude an essential singularity.

The complex identity theorem in [Brunault, Theorem 1](https://perso.ens-lyon.fr/francois.brunault/recherche/Sturm-bound-general.pdf)
applies to holomorphic integral-weight forms. On Gamma0(M), the infinity
width is one and -I belongs to the group, so a nonzero weight-six form
has order at infinity at most mu/2. This is the only sourced inequality
needed; all power and rounding steps here are deductions.

Since M is divisible by 16, its index factor from the prime two is
divisible by eight; therefore B is an integer. If F were nonzero with
all c_N=0 for 0<=N<=B, its first nonzero index m would be at least B+1.
The first coefficient of F^4 is c_m^4 at index 4m, which is nonzero
over C. Its order would satisfy 4m>=4(B+1)=mu/2+4, contradicting the
inequality. Thus F=0. The result is an exact identity test, not a
modular congruence test or a proof of minimal level.

The source modular space and slash convention were reopened in
[Kane--Kim, sections 2.1--2.2 and Proposition 2.3](https://arxiv.org/html/2211.03987v2).
No spinor-genus theorem or unreviewed reduction in level is used.

## Corrections to the submitted explanation

1. The theorem must state a, b nonzero and even a^t b explicitly, as
   above. Grok's opening theorem relied on the brief for a and parity.
2. F^2 is a legitimate weight-three form with character chi_(-4):
   epsilon_v^(-6)=chi_(-4)(v). That odd character has value -1 at -1.
   The obstruction applies only to treating F^2 as having trivial
   character. F^4 is the correct trivial-character choice used here.
3. A nonzero element of a prime field has nonzero fourth power. The
   report's contrary sentence is false. No congruence theorem is
   derived here; it would require its own applicable hypotheses.
4. A smaller congruence subgroup increases the index. A proved larger
   modularity group, such as Gamma0 of a divisor of M, might reduce
   the cutoff. There is no such improvement asserted here.
5. At general width w, Brunault's order is the index n in q^(n/w),
   rather than the numerical q exponent n/w. Width one here removes
   that distinction. Periodicity alone also does not imply finitely
   many negative powers; the coefficient-integral argument above
   provides the needed cusp conclusion directly from growth.
6. Grok's sample-shell table heading says an earlier cancelled shell,
   but several listed cancellations occur after the first nonzero
   coefficient. The numerical coefficients themselves were checked.

These corrections preserve the bound and its fourth-power proof.

## Exact controls and certificate use

For b=(1,0,0), K=diag(2,1,1):

| G | N0 | M | mu | Inclusive norm4 cutoff B |
|---|---:|---:|---:|---:|
| I3 | 4 | 64 | 96 | 12 |
| A3 simple-root Gram | 16 | 256 | 384 | 48 |
| 2I3 | 8 | 128 | 192 | 24 |
| [[2,-1,-1],[-1,4,-1],[-1,-1,4]] | 80 | 1280 | 2304 | 288 |

The last example has coefficient zero at norm4=6 and coefficient two
at norm4=10. It cannot be certified zero from its first shell. Across
the completed census N0 ranges from 2 to 96; B(96)=384. The old census
already resolved every pair without this cutoff and its verdicts stand.

[code/rank3_cutoff.py](../code/rank3_cutoff.py) builds and independently replays a
separate work7-rank3-cutoff-v1 certificate. The implementation deliberately
admits only rank three, determinant<=24 and binary even pairs with b!=0.
It computes the exact kernel basis, inverse denominator lcm, integral
N0 G0^-1, factored index and inclusive cutoff. Construction uses LDL
spheres; default replay uses inverse-Gram boxes and a separate index
formula. A smaller complete zero search remains unresolved. Any resource
interruption discards all partial coefficient data and remains unresolved.
A complete nonzero shell proves nonvanishing even below the cutoff.

This certificate proves a theta identity only. It makes no full-group,
trivial-epsilon or symmetry-converse counterexample claim. Those require
the separate complete isometry/stabilizer evidence.

Example commands, choosing fresh output paths:

```powershell
python code/rank3_cutoff.py build data/a3.json --output reports/rebuilt-rank3/a3-cutoff.json
python code/rank3_cutoff.py verify reports/rebuilt-rank3/a3-cutoff.json
python scripts/run_tests.py --all
```

The initial independent control pack and reports are stored separately
under data/census-rank3/cutoff-controls.json and
reports/census-rank3/cutoff-controls-{box,ldl}.json. They exercise a
modular-zero proof, nonzero samples, insufficient bounds and interrupted
enumeration. Dedicated tests additionally check basis changes and rejected
metadata/coefficient tampering.

Validation: all five dedicated tests pass, and both replay backends check
the nine saved controls. The control pack intentionally contains three
unresolved runs (two insufficient bounds and one interrupted enumeration).
Those verdicts remain unresolved despite successful replay of their claims.
The accompanying [index audit](../reports/census-rank3/cutoff-index-audit.json)
checks all 840 stored kernel/level records and computes cutoffs ranging from
6 to 384; it does not rerun their coefficient enumeration. Basis changes
and tampered evidence are tested by test_rank3_cutoff.py.

The controls are rebuilt by [code/build_cutoff_controls.py](../code/build_cutoff_controls.py),
which replaces only its named new-stage control pack and reports. Its
pack SHA-256 is
a9dba5f73ae88cc7b17894686eea6a3d5b9a5ffad9808a54b31dc993f075487f.
The curated release has its own byte inventory in `checksums.sha256`.
