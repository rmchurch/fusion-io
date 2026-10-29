# Stellarator mapped-Hessian correction

`m3dc1_stell_field::eval` in `m3dc1_lib/m3dc1_stell_field.cpp` incorrectly
transforms scalar second derivatives from logical to physical coordinates.
This is present in upstream `nferraro/fusion-io` at
`ed100545760e8401e778c79776ba6a3e238fb60b` (the reviewed base).

Physical coordinate functions provide a simple reproducer: evaluated as
scalars, `R` and `Z` must have identity gradients and zero physical Hessians.
The original reader violates that identity on curved maps. Its mixed
`R,phi` and `Z,phi` expressions also reuse the same numerator incorrectly.
Selected second-derivative requests need logical first derivatives that
were not always requested.

For `X(q) = (R(q), Z(q), phi)`, `J = dX/dq`, and `K = inverse(J)`, the chain rule is

```text
Hq(v) = J^T HX(v) J + v_R Hq(R) + v_Z Hq(Z)
HX(v) = K^T [Hq(v) - v_R Hq(R) - v_Z Hq(Z)] K
```

The correction evaluates this identity for all six second derivatives and
requests the necessary logical derivatives. It retains the existing first
derivative arithmetic, inverse search, tolerances, interfaces and normalization.
It changes the FusionIO reader, not the M3D-C1 evolution solver.

Consumers of these Hessians include `m3dc1_magnetic_field::eval_deriv` in
`fusion_io/m3dc1_field.cpp`. Magnetic gradients can therefore change; they
must not be described as unchanged by this fix.

## Reproduction and evidence

The standalone test needs CMake, a C++17 compiler and existing HDF5 C libraries;
it requires no plasma data, MPI, tracing or installation:

```sh
cmake -S tests/mapped_hessian -B /tmp/fusionio-hessian -DCMAKE_BUILD_TYPE=Release
cmake --build /tmp/fusionio-hessian --parallel 1
ctest --test-dir /tmp/fusionio-hessian --output-on-failure
```

Use a new build directory and `-DFUSIONIO_SOURCE_DIR=/path/to/pristine-base`
to run the same fixtures against the base revision. The executable accepts
an optional CSV output path for individual observations.

- Affine, curved and toroidally sheared analytic fixtures: **336 checks;
  84 failures before, zero after**, including selective derivative requests.
- Paired analytic outputs: **12 values and 90 first derivatives bitwise unchanged**.
- A separate retained NCSX time-zero check at **512 cell interiors**:
  all **1,536 magnetic-vector components bitwise unchanged** for its modern
  `fp` reader path. Gradients change and agree with the respective original
  and corrected formulas within **3.24e-14 T/m** maximum Frobenius residual.
  This native-case evidence is separate; no raw case data are bundled here.

The analytic oracle uses `4096 * epsilon * max(1, abs(expected))` on its
well-conditioned fixtures. This is a software bound, not a physics threshold.

## Limits

Coverage is mapped scalar derivative order through two and the tested modern
`fp` case. Third/fourth mixed derivatives and legacy `f` magnetic-gradient
paths remain unqualified; legacy magnetic values may themselves depend on
second derivatives, so vector-value preservation is not claimed universally.
The 512-point check does not certify the full volume, axis, seams or exterior.
No installed library, scientific closure or acceptance threshold is changed.
The separate surface reconstruction still fails its original 1-mm
correspondence gate; its complete 0.1% flux criterion remains unmeasured.
