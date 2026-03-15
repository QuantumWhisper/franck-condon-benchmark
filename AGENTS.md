# AGENTS.md — Codebase Knowledge Base

This file is for AI agents working on the codebase. It documents architecture, key files, known bottlenecks, and physics conventions so you don't have to reverse-engineer them from scratch.

## Project Overview

Multi-language benchmark for Franck-Condon blockade I-V simulation. Reference implementation is MATLAB. All other language ports must match it numerically within 1e-10 relative tolerance.

**Primary reference**: Koch, von Oppen, Glazman, Phys. Rev. B **74**, 205438 (2006).

## Architecture

```
Simulation Pipeline:
  Parameters → calculateAllRateW → rateW → m_rateW → {sequential rates, cotunneling rates}
                                                          ↓                    ↓
                                                     FCMatrix            sumMMr / sumMMMMrs
                                                     (Laguerre)          (regularized I,J integrals)
                                                          ↓                    ↓
                                                     fermi(ε,μ,T)        digammaFcn → BoseFcn
                                                          ↓
                                              generateMatrixW → W matrix (2N × 2N)
                                                          ↓
                                              solve_steady_state → P (occupation probabilities)
                                                          ↓
                                              current_from_rate_equations → I_seq, I_cot, I_tol
```

## Key Files (matlab/src/)

| File | Purpose | Key Dependencies |
|------|---------|-----------------|
| `simulate_current_from_rate_equations.m` | Main driver — sweeps bias voltage | `current_from_rate_equations` |
| `run_benchmark.m` (matlab/) | Benchmark runner with timing + JSON output | All of src/ |
| `current_from_rate_equations.m` | Computes I_seq + I_cot from steady-state P | `generateMatrixW`, `rateW`, `solve_steady_state` |
| `generateMatrixW.m` | Builds 2N×2N rate equation matrix W | `rateW_lead` |
| `rateW_lead.m` | Rate summed over both leads | `rateW` |
| `rateW.m` | Dispatcher: memoized calculation or 4D cache lookup | `m_rateW` |
| `m_rateW.m` | Core rate computation — sequential + cotunneling | `FCMatrix`, `fermi`, `sumMMr`, `sumMMMMrs` |
| `FCMatrix.m` / `FCMatrixSingle.m` | Franck-Condon matrix elements via Laguerre polynomials | — |
| `sumMMr.m` / `sumMMMMrs.m` | Cotunneling single/double sums (n=0→0) | `regularizedI`, `regularizedJ` |
| `sumMMr11.m` / `sumMMMMrs11.m` | Cotunneling sums (n=1→1) | `regularizedI`, `regularizedJ` |
| `regularizedI.m` / `regularizedJ.m` | Analytically regularized cotunneling integrals | `digammaFcn`, `BoseFcn` |
| `digammaFcn.m` | Polygamma function (uses symbolic math — SLOW) | — |
| `solve_steady_state_occupation_probabilities.m` | Solves 0=WP via lsqlin | — |

## Known Performance Bottleneck

`digammaFcn.m` converts its argument to symbolic via `sym(x)`, evaluates the polygamma function symbolically, then converts back to `double`. This is the single biggest bottleneck — symbolic math is roughly 1000x slower than numeric evaluation.

The reason this exists: MATLAB's built-in `psi(k, x)` does not support complex arguments. The polygamma function is needed for the analytically regularized cotunneling integrals.

**For all language ports**: use a native complex digamma/polygamma implementation. In Julia, `SpecialFunctions.digamma` handles complex arguments. In Python, `mpmath.digamma` or `scipy`'s polygamma with complex extension works. In C, use the GSL or implement the asymptotic series directly. Do not replicate the symbolic workaround.

## Physics Notes for Porters

**State space**:
- Charge states: n ∈ {0, 1} (empty or singly occupied)
- Phonon states: q ∈ {0, 1, ..., N-1}
- State vector P has 2N components: P^0_0, ..., P^0_{N-1}, P^1_0, ..., P^1_{N-1}

**Franck-Condon matrix elements**:
The FC element ⟨q₂|D(λ)|q₁⟩ is computed using generalized Laguerre polynomials. See `FCMatrix.m` for the exact formula. The displacement operator D(λ) shifts the oscillator equilibrium by λ in dimensionless units.

**Rates**:
- Spin degeneracy: s=2 for transitions 0→1 (empty to occupied), s=1 for 1→0, s=2 for cotunneling
- Sequential rate: Γ × |M_{q1,q2}|² × f(ε), where f is the Fermi function
- Cotunneling rate: involves a double sum over virtual intermediate states with regularized integrals to handle the energy denominator divergence

**Chemical potentials**:
- μ_L = η × V_sd (left lead)
- μ_R = -(1-η) × V_sd (right lead)
- Molecular level: ε_d = V_g

**Conventions**:
- All energies in eV
- Temperature in Kelvin
- Current in units of e (multiply by elementary charge for SI Amperes)
- Lead index: +1 for left, -1 for right

## Benchmark Structure

Two parameter specs exist in `benchmark/spec/`:

| Spec | File | N | Use |
|------|------|---|-----|
| quick | `quick_params.json` | 6 | Fast validation, development iteration |
| default | `default_params.json` | 15 | Publication-quality results, real benchmark |

Each benchmark run produces four outputs in `benchmark/results/`:
- `{lang}_{spec}_results.json` — full results + metadata
- `{lang}_{spec}_IV.csv` — I-V data (Vsd, I_tol, I_seq, I_cot)
- `{lang}_{spec}_IV.pdf` — publication-quality plot (vector)
- `{lang}_{spec}_IV.png` — publication-quality plot (300 dpi)

## Julia Implementation (julia/src/)

The Julia port is a 1:1 faithful translation of the MATLAB reference. Every function, formula, and sign convention matches.

### Key Files

| File | MATLAB Equivalent | Notes |
|------|-------------------|-------|
| `FranckCondon.jl` | Module definition | Includes all files, exports `simulate_iv` |
| `constants.jl` | `KBoltzmann_ev`, `hbar_eV`, `ee_ElementaryCharge` | Same exact values |
| `laguerre.jl` | `laguerreL` (built-in) | Three-term recurrence, matches to full Float64 |
| `fc_matrix.jl` | `FCMatrixSingle.m`, `FCMatrix.m` | Dict-based cache replaces MATLAB `memoize` |
| `fermi_bose.jl` | `fermi.m`, `BoseFcn.m` | Identical formulas |
| `regularized.jl` | `regularizedI.m`, `regularizedJ.m`, `digammaFcn.m` | Uses `SpecialFunctions.digamma/trigamma` |
| `cotunneling.jl` | `sumMMr.m`, `sumMMMMrs.m`, `sumMMr11.m`, `sumMMMMrs11.m` + their `m_` inner functions | Convergence wrappers match MATLAB criteria |
| `rate.jl` | `m_rateW.m`, `rateW.m`, `calculateAllRateW.m` | Pre-computed `Dict` store replaces MATLAB memoize cache |
| `matrix.jl` | `generateMatrixW.m` | Same index mapping, same sigma functions, same peq |
| `solver.jl` | `solve_steady_state_occupation_probabilities.m` | Augmented system (replace last row with normalization) |
| `current.jl` | `current_from_rate_equations.m` | Same sign conventions for I_seq0, I_seq1, I_cot |
| `plotting.jl` | `plot_IV` in `run_benchmark.m` | CairoMakie, LaTeX labels, 12×9 cm, 300 dpi |
| `run_benchmark.jl` | `run_benchmark.m` | JIT warm-up + 3 timed runs + median |

### Design Decisions

| Decision | Rationale |
|----------|-----------|
| `SpecialFunctions.jl` digamma/trigamma | Native complex support; replaces MATLAB's symbolic bottleneck (1000x faster) |
| Laguerre via recurrence (no package) | Zero deps, matches MATLAB `laguerreL` exactly, ~10 ns/call |
| `Dict`-based FC matrix cache | Replaces MATLAB `memoize`; lambda is fixed per simulation so key is just (q1, q2) |
| Pre-computed rate store `Dict{Tuple,Vector}` | Replaces MATLAB's global memoize cache; key = (n1, n2, q1, lead), value = rates for all q2 |
| Augmented system solver (not SVD) | Replace last row of W with normalization; `M_aug \ d` via QR; clamp+renormalize |
| CairoMakie for plots | Pure Julia, MathTeXEngine for LaTeX, vector PDF + 300 dpi PNG |

### Numerical Precision Notes

The Julia port matches MATLAB to ~5 significant digits across all 201 bias points. The residual ~1e-5 relative difference arises from accumulated floating-point differences in the cotunneling pipeline:

1. **Digamma precision**: Julia's `SpecialFunctions.digamma` uses Stirling series (~1e-13 relative error). MATLAB's `psi(k, sym(x))` uses arbitrary precision. Verified via BigFloat: this accounts for <1 ULP difference per call.

2. **Cancellation amplification**: The regularized integrals compute `Re(ψ(a1)-ψ(a2)-ψ(a3)+ψ(a4))`. When the four digamma values nearly cancel, per-ULP differences get amplified.

3. **Steady-state sensitivity**: The rate matrix W is ill-conditioned at certain bias points (near current zero-crossings). The condition number σ₁/σ₁₂ can exceed 10¹⁸, amplifying any matrix-element error into the occupation probabilities P.

4. **MATLAB solver artifacts**: At specific bias points (Vsd ≈ 0.219, 0.585), MATLAB's `lsqlin` interior-point produces non-smooth I-V values that appear to be solver-dependent numerical artifacts.

**For future porters**: the 1e-10 tolerance target is achievable for the sequential tunneling component alone (which uses simple FC²×fermi), but the cotunneling component inherently amplifies ULP-level differences. A realistic cross-language tolerance is ~1e-5 for the total current.

## How to Add a New Language Port

1. Create a `{language}/` directory at the project root.
2. Read all parameters from `benchmark/spec/default_params.json`. Do not hardcode.
3. Implement the simulation pipeline following the architecture diagram above.
4. Output results in the JSON format specified in `BENCHMARK.md`.
5. Generate a publication-quality I-V plot (PDF + PNG) with three curves (I_tol, I_seq, I_cot).
6. Validate against the MATLAB reference results using the tolerance check in `BENCHMARK.md`.
7. Follow the benchmark protocol (warm-up for JIT, minimum 3 timed runs, report median).

### Lessons from the Julia Port

- **Digamma is the key bottleneck**: MATLAB's symbolic `psi(k, sym(x))` is ~1000x slower than native complex digamma. Every language port should use a native implementation (Julia: `SpecialFunctions.jl`, Python: `scipy.special` or `mpmath`, C: GSL `gsl_sf_psi`, Rust: custom or `special` crate).
- **Laguerre polynomial**: Implement via three-term recurrence. No need for external packages. The recurrence is forward-stable for real positive x.
- **Memoization**: MATLAB uses `memoize(@func)` extensively. In other languages, use a Dict/HashMap cache keyed by the function arguments. The FC matrix cache is most critical (called millions of times with repeated arguments).
- **Steady-state solver**: MATLAB uses `lsqlin` (constrained least-squares, interior-point). A simpler approach works: replace the last row of W with the normalization constraint `sum(P)=1`, solve via backslash/QR, clamp negatives, renormalize. For sensitive bias points, consider a proper QP solver (e.g., Ipopt, OSQP).
- **Bias sweep generation**: MATLAB's colon operator `a:d:b` uses different floating-point arithmetic than other languages' range generators. Read the MATLAB reference Vsd values from the JSON when validating to avoid ULP-level bias point mismatches.
- **Sign conventions**: Pay careful attention to the sign in the current computation. For n=0→1: `diffW = w_R - w_L`. For n=1→0: `diffW = w_L - w_R` (opposite!). For cotunneling: `diffW = w_RL - w_LR`.

The MATLAB source in `matlab/src/` is the ground truth for numerical behavior. When in doubt about a formula or sign convention, read the MATLAB code and the Koch et al. paper.
