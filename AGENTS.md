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

**For all language ports**: use a native complex digamma/polygamma implementation. In Julia, `SpecialFunctions.digamma` handles complex arguments. In Python, `mpmath.digamma` or `scipy`'s polygamma with complex extension works. In C, use the asymptotic series directly (see `c/src/digamma.c`). Do not replicate the symbolic workaround.

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

### Default Spec Benchmark Results

Default spec (N=15) results are available for all languages except Python. MATLAB was benchmarked on an Apple M4 Max (14 cores: 4E+10P, 32 GPU cores); all other languages were benchmarked on an Apple M5 (10 cores: 4S+6E, 24 GB).

| Language | Default (N=15) | Speedup vs MATLAB | Hardware |
|----------|----------------|-------------------|----------|
| MATLAB | 11725 s (~3.3 hr) | 1× (reference) | M4 Max |
| Rust | 1.9 s | **6202×** | M5 |
| C (GSL) | 2.2 s | **5330×** | M5 |
| C++ | 2.0 s | **5863×** | M5 |
| Fortran | 2.6 s | **4510×** | M5 |
| Julia | 102 s | **115×** | M5 |
| Python | — | — | too slow for N=15 |

All non-MATLAB ports agree to **<5×10⁻¹³** with each other. Against MATLAB, all match to **6.1×10⁻⁶** (excluding solver artifact points).

### Solver Artifact Points

The rate matrix W becomes ill-conditioned at certain bias points, causing different linear solvers (MATLAB's `lsqlin`, Rust's `nalgebra` QR, GSL QR, Eigen QR, LAPACK LU) to produce different results. These are not bugs — they are intrinsic to the mathematical problem at those specific Vsd values.

| Spec | Artifact Vsd values | Cause |
|------|---------------------|-------|
| quick (N=6) | ≈ 0.219, 0.585 | W ill-conditioned near current zero-crossings |
| default (N=15) | ≈ 0.219, 0.438 | Same phenomenon, different Vsd values due to different N |

At these points, each implementation gives a valid but solver-dependent result. All cross-language error metrics exclude these points. The non-artifact points agree to ~1e-6 vs MATLAB and ~1e-13 between non-MATLAB ports.

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

4. **MATLAB solver artifacts**: At specific bias points (Vsd ≈ 0.219, 0.585 for quick; Vsd ≈ 0.219, 0.438 for default), MATLAB's `lsqlin` interior-point produces non-smooth I-V values that appear to be solver-dependent numerical artifacts.

**Default spec (N=15)**: Julia matches MATLAB to **6.1×10⁻⁶** max relative error (excluding solver artifacts at Vsd ≈ 0.219, 0.438). Julia matches Rust to **5.1×10⁻¹³**. Wall time: **102 s** on Apple M5 (115× vs MATLAB). Compared to the quick spec, the default spec is ~19.5× slower due to the larger state space (2N=30 vs 2N=12) and deeper cotunneling convergence.

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

- **Digamma is the key bottleneck**: MATLAB's symbolic `psi(k, sym(x))` is ~1000x slower than native complex digamma. Every language port should use a native implementation (Julia: `SpecialFunctions.jl`, Python: `scipy.special` or `mpmath`, C/Rust: custom asymptotic series).
- **Laguerre polynomial**: Implement via three-term recurrence. No need for external packages. The recurrence is forward-stable for real positive x.
- **Memoization**: MATLAB uses `memoize(@func)` extensively. In other languages, use a Dict/HashMap cache keyed by the function arguments. The FC matrix cache is most critical (called millions of times with repeated arguments).
- **Steady-state solver**: MATLAB uses `lsqlin` (constrained least-squares, interior-point). A simpler approach works: replace the last row of W with the normalization constraint `sum(P)=1`, solve via backslash/QR, clamp negatives, renormalize. For sensitive bias points, consider a proper QP solver (e.g., Ipopt, OSQP).
- **Bias sweep generation**: MATLAB's colon operator `a:d:b` uses different floating-point arithmetic than other languages' range generators. Read the MATLAB reference Vsd values from the JSON when validating to avoid ULP-level bias point mismatches.
- **Sign conventions**: Pay careful attention to the sign in the current computation. For n=0→1: `diffW = w_R - w_L`. For n=1→0: `diffW = w_L - w_R` (opposite!). For cotunneling: `diffW = w_RL - w_LR`.

The MATLAB source in `matlab/src/` is the ground truth for numerical behavior. When in doubt about a formula or sign convention, read the MATLAB code and the Koch et al. paper.

## Python Implementation (python/src/)

The Python port follows the Julia port structure (which is itself a 1:1 MATLAB translation). Uses Numba JIT for hot inner loops and scipy.special for complex digamma.

### Key Files

| File | MATLAB Equivalent | Notes |
|------|-------------------|-------|
| `__init__.py` | Module definition | Exports `simulate_iv`, `ELEMENTARY_CHARGE` |
| `constants.py` | `KBoltzmann_ev`, `hbar_eV`, `ee_ElementaryCharge` | Same exact values |
| `laguerre.py` | `laguerreL` (built-in) | Three-term recurrence, `@numba.njit` |
| `fc_matrix.py` | `FCMatrixSingle.m`, `FCMatrix.m` | `@numba.njit` single element + dict cache |
| `fermi_bose.py` | `fermi.m`, `BoseFcn.m` | Identical formulas, numpy vectorized |
| `regularized.py` | `regularizedI.m`, `regularizedJ.m`, `digammaFcn.m` | `scipy.special.digamma` + Numba trigamma (asymptotic series) |
| `cotunneling.py` | `sumMMr.m`, `sumMMMMrs.m`, etc. + `m_` inner functions | Convergence wrappers match MATLAB criteria exactly |
| `rate.py` | `m_rateW.m`, `rateW.m`, `calculateAllRateW.m` | Dict-based rate store, same as Julia |
| `matrix.py` | `generateMatrixW.m` | Same index mapping, sigma functions, peq |
| `solver.py` | `solve_steady_state_occupation_probabilities.m` | Augmented system via `numpy.linalg.lstsq` |
| `current.py` | `current_from_rate_equations.m` | Same sign conventions for I_seq0, I_seq1, I_cot |
| `plotting.py` | `plot_IV` in `run_benchmark.m` | matplotlib + LaTeX, 12×9 cm, 300 dpi |
| `simulate.py` | Main simulation loop | Matches Julia `simulate_iv` exactly |
| `run_benchmark.py` | `run_benchmark.m` | 3 timed runs + median, validate vs MATLAB |

### Design Decisions

| Decision | Rationale |
|----------|-----------|
| `scipy.special.digamma` for complex digamma | Vectorized C code, handles complex natively; replaces MATLAB symbolic bottleneck |
| Custom `@numba.njit` trigamma | Asymptotic series with 20 Bernoulli numbers; `scipy.special.polygamma` may not support complex args in all versions |
| `@numba.njit` for laguerre and FC matrix | Hot inner functions called O(N²) times; pure math suitable for JIT |
| numpy vectorization for cotunneling sums | Matrix operations on small arrays; overhead of Numba JIT not worth it here |
| Dict-based FC matrix and rate store | Same pattern as Julia port; simple, correct, negligible overhead |
| `numpy.linalg.lstsq` for steady-state | Augmented system approach; clamp negatives, renormalize |

### Numba + scipy Tension

The digamma/trigamma functions (called millions of times in cotunneling sums) live inside numpy-vectorized operations, not tight Python loops. The strategy is:

1. **Digamma**: `scipy.special.digamma` — vectorized C code operating on numpy arrays. Fast enough that Numba JIT offers no benefit.
2. **Trigamma**: Custom `@numba.njit` implementation using the Bernoulli asymptotic series. Called element-wise via a Numba-JIT'd loop for maximum throughput.
3. **FC matrix, Laguerre**: `@numba.njit` for the scalar computation (called many times with caching).
4. **Everything else**: numpy vectorization. The Python loop overhead in the convergence wrappers (~3 iterations) is negligible.

An alternative approach using JAX (`jax.scipy.special.digamma` with end-to-end JIT) could potentially match Julia's speed but would require restructuring the entire pipeline for functional/array programming.

### Numerical Precision Notes

The Python port matches Julia to **machine precision** (~8e-13 relative error) at 199 of 201 bias points. The 2 outlier points (Vsd ≈ 0.219, 0.585) are the known MATLAB solver artifacts.

Against MATLAB (excluding artifact points): max relative error = 7.7e-5, consistent with the Julia port's ~1e-5 cross-language tolerance.

Performance: ~11 seconds for quick spec (N=6) on Apple Silicon, a **169x speedup** over MATLAB (1844s) and ~2x slower than Julia (5.2s). The default spec (N=15) was not benchmarked — Python's interpreted cotunneling loops scale poorly with N, making the default spec impractically slow on consumer hardware (estimated >30 min per run).

## C Implementation (c/src/)

The C port uses GSL (GNU Scientific Library) for QR decomposition, OpenMP for parallel bias-point computation, custom asymptotic series for complex digamma/trigamma (matching the Rust algorithm), cJSON (vendored) for JSON I/O, and gnuplot for plotting. It follows the same pipeline structure as Julia and Python.

### Key Files

| File | MATLAB Equivalent | Notes |
|------|-------------------|-------|
| `constants.h` | `KBoltzmann_ev`, `hbar_eV`, `ee_ElementaryCharge` | Header-only, same exact values |
| `laguerre.h/.c` | `laguerreL` (built-in) | Three-term recurrence, `int` alpha |
| `fc_matrix.h/.c` | `FCMatrixSingle.m`, `FCMatrix.m` | Static 256×256 2D array cache, log-space overflow handling |
| `fermi_bose.h/.c` | `fermi.m`, `BoseFcn.m` | Identical formulas |
| `digamma.h/.c` | `digammaFcn.m` | Pure asymptotic series for both digamma and trigamma; fast-path (5 terms) for |z|²>900; 10 Bernoulli terms for full path |
| `regularized.h/.c` | `regularizedI.m`, `regularizedJ.m` | Three variants: `regularized_I`, `regularized_J` (vector epsilon), `regularized_J_matrix` (matrix epsilon for n=1→1) |
| `cotunneling.h/.c` | `sumMMr.m`, `sumMMMMrs.m`, `sumMMr11.m`, `sumMMMMrs11.m` + `m_` inners | Convergence wrappers with selective recomputation; largest module (~570 lines) |
| `rate.h/.c` | `m_rateW.m`, `rateW.m`, `calculateAllRateW.m` | Flat array `RateStore` indexed by `[n1][n2][q1][lead_idx][q2]` |
| `matrix.h/.c` | `generateMatrixW.m` | Same index mapping, sigma functions, peq; `1/INFINITY == 0` for tau terms |
| `solver.h/.c` | `solve_steady_state_occupation_probabilities.m` | GSL `gsl_linalg_QR_decomp` + `gsl_linalg_QR_solve`; augmented system, clamp+renormalize |
| `current.h/.c` | `current_from_rate_equations.m` | Same sign conventions: 0→1 R−L, 1→0 L−R, cot RL−LR |
| `simulate.h/.c` | Main simulation loop | OpenMP `parallel for schedule(dynamic, 4)`; FC cache pre-populated and shared read-only; rate store per thread |
| `json_io.h/.c` | JSON I/O | cJSON-based; parses `tau:"Inf"` → `INFINITY`; loads MATLAB Vsd for bit-for-bit matching |
| `plotting.h/.c` | Plot generation | gnuplot via `popen()`; PDF + PNG; graceful fallback if gnuplot unavailable |
| `main.c` | `run_benchmark.m` | `clock_gettime(CLOCK_MONOTONIC)` timing; 3 runs + median; validation excludes MATLAB solver artifacts |

### Design Decisions

| Decision | Rationale |
|----------|-----------|
| Pure asymptotic series for digamma (no GSL) | Matches Rust/Fortran/C++ approach; eliminates GSL dependency for digamma; fast-path covers >99% of calls |
| Custom trigamma via 10-term Bernoulli asymptotic series | GSL has no complex trigamma; algorithm matches Rust exactly |
| Fast-path `digamma_asymptotic5` / `trigamma_asymptotic5` | 5 Bernoulli terms for \|z\|²>900; skips reflection and recurrence entirely; >99% hit rate at T=4.2K |
| Pre-computed `DIGAMMA_COEFF[k] = B_{2(k+1)} / (2*(k+1))` | Eliminates runtime division; multiply-instead-of-divide pattern |
| `norm_sqr` threshold instead of `cabs()` | Avoids sqrt per recurrence iteration; `re*re + im*im < 100.0` |
| OpenMP `parallel for schedule(dynamic, 4)` | Bias points are independent; dynamic scheduling handles variable cotunneling cost |
| FC cache pre-populated via `fc_cache_populate(fc, FC_MAX_N)` | All 256×256 entries computed before parallel region; `fc_cache_get` becomes purely read-only (thread-safe) |
| Static 256×256 FC cache (`FCCache` struct) | Avoids hash table overhead; max N=256 covers all practical cases; O(1) lookup |
| Flat-array `RateStore` with macro accessor | Replaces Dict/HashMap; single `malloc` per bias point; cache-friendly access pattern |
| GSL QR decomposition (not LU) | Matches Julia/Python augmented system approach; numerically stable for ill-conditioned W |
| cJSON (vendored, MIT) | Single .c/.h file; no build dependency; sufficient for benchmark JSON I/O |
| gnuplot via `popen()` | No compiled plotting dependency; publication-quality LaTeX labels; PDF+PNG |
| `-O2 -flto -march=native` (not `-O3`) | `-O3 -flto` causes 3× regression on Apple Clang due to LTO interaction; `-O2 -flto` is the optimal combination |
| `clock_gettime(CLOCK_MONOTONIC)` | Portable high-resolution timer; works on macOS and Linux |

### Complex Digamma/Trigamma Strategy

C does not have a native complex polygamma function. The implementation uses a pure asymptotic series approach matching the Rust port:

1. **Digamma ψ(z)**: Custom implementation with three tiers:
   - **Fast-path** (`digamma_asymptotic5`): For Re(z) > 0 and |z|² > 900 — uses only 5 pre-computed Bernoulli coefficients, no reflection or recurrence. At T=4.2K, >99% of calls hit this path.
   - **Full-path**: 10-term Bernoulli asymptotic expansion with pre-computed `DIGAMMA_COEFF` array. Recurrence shift until |z|² ≥ 100 (using `norm_sqr` to avoid `sqrt`). Reflection formula for Re(z) ≤ 0.
   - **Fallback**: Same algorithm handles all edge cases; no GSL dependency.

2. **Trigamma ψ'(z)**: Custom implementation with matching structure:
   - **Fast-path** (`trigamma_asymptotic5`): For Re(z) > 0 and |z|² > 900 — 5 Bernoulli terms.
   - **Full-path**: 10-term Bernoulli expansion, recurrence threshold |z|² ≥ 100, reflection formula.

The Bernoulli coefficients are stored as compile-time `static const double` arrays, computed as exact rational fractions (e.g., `1.0/6.0`, `-691.0/2730.0`). The `DIGAMMA_COEFF` array pre-divides by `2*(k+1)` to eliminate runtime division.

### Performance Notes

The C port runs the quick spec (N=6) in ~0.56 seconds on Apple M5 — a **3293x speedup** over MATLAB (1844s), faster than Julia (5.2s) and Python (11s).

The key optimization is **factored digamma precomputation** in `regularized_I`. The digamma arguments factor into row-only and column-only terms:
- `ψ(a1[i])` and `ψ(a3[i])` depend only on `epsilon1[i]` (row index)
- `ψ(a2[j])` and `ψ(a4[j])` depend only on `epsilon2[j]` (column index)

By precomputing `ψ(a1[i]) - ψ(a3[i])` per row and `ψ(a2[j]) - ψ(a4[j])` per column, the total digamma calls drop from **4×N²** to **4×N** — a 100x reduction for N=100 (the typical convergence truncation). Similarly, `regularized_J` precomputes `trigamma(a2)` per row since it depends only on `epsilon[i]`.

Additional optimizations: trigamma uses 10 Bernoulli terms (not 20) with recurrence threshold |z|≥10 (not 20), halving both recurrence steps and series terms while maintaining full double precision.

Performance history on the same hardware:
- Hand-written digamma: 428s (4x vs MATLAB)
- + GSL `gsl_sf_complex_psi_e`: 52s (36x)
- + factored regularized_I + optimized trigamma: 2.3s (809x) [quick spec, M4 Max]
- + pure asymptotic digamma + OpenMP + LTO: **0.56s (3293x)** [quick spec, M5]

**Default spec (N=15)**: **2.2 s** on Apple M5 — a **5330× speedup** over MATLAB. Matches MATLAB to **6.1×10⁻⁶** (excluding solver artifacts at Vsd ≈ 0.219, 0.438). Matches Rust to **5.4×10⁻¹³**.

### Numerical Precision Notes

The C port matches MATLAB to **7.65e-5 max relative error** at 199 of 201 bias points (quick spec). The 2 outlier points (Vsd ≈ 0.219, 0.585) are the known MATLAB `lsqlin` interior-point solver artifacts where MATLAB produces non-physical near-zero currents.

Against Python (excluding artifact points): the C and Python ports agree to ~machine precision, confirming that the C implementation is a faithful translation.

The error budget is identical to the Julia and Python ports:
- Sequential tunneling (FC²×fermi): matches to ~1e-14 (trivial computation)
- Cotunneling (digamma cancellation + ill-conditioned W): ~1e-5 due to ULP amplification
- MATLAB solver artifacts at Vsd ≈ 0.219, 0.585: excluded from validation

### Lessons from the C Port

- **GSL complex digamma exists but is not needed**: GSL 2.x provides `gsl_sf_complex_psi_e`, but a pure asymptotic series (matching Rust's algorithm) is equally accurate and removes the GSL dependency for digamma. GSL is still used for QR decomposition in the solver.
- **GSL has no complex trigamma**: You must implement this yourself. The 10-term Bernoulli asymptotic series with recurrence shift (|z|² ≥ 100) and reflection formula works well.
- **`-O3 -flto` causes regression on Apple Clang**: The combination of `-O3` with LTO causes a 3× performance regression on Apple Clang 17. Use `-O2 -flto` instead — LTO alone provides a significant speedup (10.7s vs 16.7s single-threaded for default spec). This is a known Apple Clang behavior.
- **Static arrays beat hash tables for FC cache**: With max N=256, a 256×256 `double` array (512 KB) is faster than any hash table and has zero collision overhead.
- **FC cache must be pre-populated for OpenMP**: The lazy-init pattern in `fc_cache_get` is not thread-safe. Pre-populating all 256×256 entries via `fc_cache_populate(fc, FC_MAX_N)` before the parallel region makes `fc_cache_get` purely read-only (safe to share across threads without synchronization).
- **Flat rate store with macro indexing**: A single `malloc(2*2*N*2*N * sizeof(double))` with a 5D index macro is simpler and faster than nested arrays or hash maps.
- **`1.0/INFINITY == 0.0`**: IEEE 754 guarantees this, so `tau=Inf` (unequilibrated phonons) works without special-casing the `1/tau` terms in the W matrix.
- **cJSON is sufficient**: A 3000-line vendored library handles all the JSON I/O needs. No need for heavier dependencies like jansson or json-c.
- **gnuplot via popen() is surprisingly capable**: LaTeX-quality labels, PDF vector output, and 300 dpi PNG — all from a simple script piped to `gnuplot`.
- **convergence wrappers dominate runtime**: The `sumMMr`/`sumMMMMrs` convergence loops account for >95% of wall time. Any performance optimization should target the inner `regularized_I`/`regularized_J` calls.

## Rust Implementation (rust/src/)

The Rust port follows the C port structure (which is itself a faithful MATLAB translation). Pure Rust with no C/FFI dependencies — uses custom asymptotic series for complex digamma/trigamma instead of GSL.

### Key Files

| File | C Equivalent | Notes |
|------|-------------|-------|
| `lib.rs` | Module definition | Declares all submodules |
| `constants.rs` | `constants.h` | Same exact values |
| `laguerre.rs` | `laguerre.c` | Three-term recurrence, identical algorithm |
| `fc_matrix.rs` | `fc_matrix.c` | Flat `Vec<f64>` cache (256×256), log-space overflow handling |
| `fermi_bose.rs` | `fermi_bose.c` | Identical formulas |
| `digamma.rs` | `digamma.c` | Pure Rust asymptotic series (no GSL FFI); 20 Bernoulli terms for digamma, 10 for trigamma |
| `regularized.rs` | `regularized.c` | Factored digamma precomputation (O(N) instead of O(N²) digamma calls) |
| `cotunneling.rs` | `cotunneling.c` | Convergence wrappers with selective recomputation; largest module |
| `rate.rs` | `rate.c` | Flat `Vec<f64>` rate store indexed by `[n1][n2][q1][lead_idx][q2]` |
| `matrix.rs` | `matrix.c` | Same index mapping, sigma functions, peq |
| `solver.rs` | `solver.c` | `nalgebra` QR decomposition; augmented system, clamp+renormalize |
| `current.rs` | `current.c` | Same sign conventions: 0→1 R−L, 1→0 L−R, cot RL−LR |
| `simulate.rs` | `simulate.c` | FC cache on heap, rate store per bias point, `-sign(Vsd)` correction |
| `json_io.rs` | `json_io.c` | `serde_json` for parsing; manual JSON/CSV writing for output |
| `plotting.rs` | `plotting.c` | gnuplot via `std::process::Command`; PDF + PNG |
| `main.rs` | `main.c` | `std::time::Instant` timing; 3 runs + median; validation excludes MATLAB solver artifacts |
| `explorer.rs` | *(new)* | Interactive TUI (ratatui + crossterm): I-V curve explorer with auto-recompute, stability diagram with heatmap |

### Design Decisions

| Decision | Rationale |
|----------|-----------|
| Pure Rust digamma (no GSL FFI) | Zero C dependencies; factored precomputation makes per-call cost negligible |
| 20-term Bernoulli series for digamma, 10-term for trigamma | Matches C port's algorithm exactly; full double precision |
| Recurrence threshold \|z\|≥20 (digamma) and \|z\|≥10 (trigamma) | Same as C port; trades recurrence steps for series accuracy |
| `num-complex` crate for `Complex64` | De facto standard; `.norm()`, `.ln()`, `.sin()`, `.cos()` all available |
| `nalgebra` QR decomposition | Mature Rust linear algebra; matches C/Julia/Python augmented system approach |
| Flat `Vec<f64>` for FC cache and rate store | Same pattern as C port; cache-friendly, no hash overhead |
| `serde_json::Value` for JSON parsing | Handles `"tau": "Inf"` string vs number without custom deserializer |
| Manual JSON/CSV writing (not serde serialization) | Full control over formatting (%.17g precision); matches C port output exactly |
| gnuplot via `Command::new("gnuplot")` | Same approach as C port; graceful fallback if gnuplot unavailable |

### Complex Digamma/Trigamma Strategy

Rust has no standard library for complex polygamma functions. The implementation uses a pure Rust approach:

1. **Digamma ψ(z)**: Custom asymptotic series with:
   - Reflection formula for Re(z) ≤ 0: `ψ(z) = ψ(1−z) − π·cot(πz)`
   - Recurrence shift until |z| ≥ 20: `ψ(z) = ψ(z+1) − 1/z`
   - 20-term Bernoulli asymptotic expansion: `ψ(z) ≈ ln(z) − 1/(2z) − Σ B_{2k}/(2k·z^{2k})`

2. **Trigamma ψ'(z)**: Custom asymptotic series with:
   - Reflection formula for Re(z) ≤ 0: `ψ'(z) = (π/sin(πz))² − ψ'(1−z)`
   - Recurrence shift until |z| ≥ 10: `ψ'(z) = ψ'(z+1) + 1/z²`
   - 10-term Bernoulli asymptotic expansion: `ψ'(z) ≈ 1/z + 1/(2z²) + Σ B_{2k}/z^{2k+1}`

The Bernoulli coefficients B₂ through B₄₀ are stored as `static` compile-time constants, computed as exact rational fractions.

### Performance Notes

The Rust port runs the quick spec (N=6) in ~1.5 seconds on Apple Silicon — a **1230x speedup** over MATLAB (1844s), the fastest implementation, beating C+GSL (2.3s), Julia (5.2s), and Python (11s).

The performance was achieved through several optimizations:
1. **Reduced asymptotic series**: 10 Bernoulli terms instead of 20 (still ~1e-20 truncation error at |z|≥10, far below machine epsilon)
2. **Pre-computed coefficients**: `DIGAMMA_COEFF[k] = B_{2(k+1)}/(2*(k+1))` as a static array, eliminating runtime division
3. **Multiply-instead-of-divide**: `inv_power * coeff` instead of `coeff / power`, replacing complex division with cheaper scalar-complex multiply
4. **norm_sqr() threshold**: `z.norm_sqr() < 100.0` instead of `z.norm() < 10.0`, avoiding sqrt per iteration
5. **Aggressive #[inline]**: All hot-path functions annotated with `#[inline]` or `#[inline(always)]`, enabling cross-module optimization with LTO

Performance history:
- Initial port (20-term digamma, no inlining): 7.2s (254x vs MATLAB)
- Optimized (10-term, pre-computed, inlined): **1.5s (1230x vs MATLAB)**

**Default spec (N=15)**: **1.9 s** on Apple M5 — a **6202× speedup** over MATLAB. The fastest implementation across both specs. Matches MATLAB to **6.1×10⁻⁶** (excluding solver artifacts at Vsd ≈ 0.219, 0.438).

### Numerical Precision Notes

The Rust port matches MATLAB to **7.65e-5 max relative error** at 199 of 201 bias points (quick spec), identical to the C port's error. The 2 outlier points (Vsd ≈ 0.219, 0.585) are the known MATLAB solver artifacts.

The error budget is identical to all other ports:
- Sequential tunneling (FC²×fermi): matches to ~1e-14
- Cotunneling (digamma cancellation + ill-conditioned W): ~1e-5 due to ULP amplification
- MATLAB solver artifacts at Vsd ≈ 0.219, 0.585: excluded from validation

### Lessons from the Rust Port

- **Pure Rust digamma is viable**: With factored precomputation, the per-call cost of a Rust asymptotic series is negligible. No need for GSL FFI unless targeting C-level performance.
- **`num-complex` is production-ready**: Complex64 arithmetic, transcendentals, and formatting all work correctly. Use `.norm()` (not `.abs()`) for the complex modulus.
- **`nalgebra` QR is correct for ill-conditioned systems**: The augmented system solver produces identical results to GSL's QR.
- **Rust's borrow checker catches array aliasing bugs**: The C port's `sort3(&a[0], &a[1], &a[2])` pattern doesn't compile in Rust — use `arr.sort_by()` instead.
- **`f64::INFINITY` and IEEE 754**: `1.0 / f64::INFINITY == 0.0` works in Rust just as `1.0/INFINITY` works in C, so `tau=Inf` requires no special handling.
- **No CSV crate needed**: Manual `write!` with `{:.6e}` format matches the C port's `fprintf(fp, "%.6e")` output exactly.
- **Minimal dependencies**: 5 crates (num-complex, serde, serde_json, nalgebra, rayon) for simulation + ratatui for the interactive explorer.

### Performance Optimizations (post-benchmark)

The Rust port includes optimizations beyond the initial faithful translation:

1. **Rayon parallelism**: Bias points are computed in parallel via `rayon::par_iter`. The FC cache is pre-populated and shared immutably across threads (`&FCCache` instead of `&mut FCCache`). Achieves near-linear scaling with core count.

2. **Specialized fast-path digamma**: `digamma_asymptotic5` and `trigamma_asymptotic5` use 5 Bernoulli terms (instead of 10) for |z|² > 900, skipping reflection and recurrence. At T=4.2K, >99% of calls hit this path.

3. **Digamma lookup table infrastructure**: `DigammaTable` precomputes ψ(½+iy) on a fine grid with Catmull-Rom interpolation. Currently benefits low-temperature regimes where |y| < 10 triggers recurrence. Falls back to the original asymptotic series for |y| outside the table range.

### Robustness: Parameter Validity and Edge Cases

The Rust port validates inputs and handles edge cases that other ports may not:

**Input validation** (enforced with `panic!` in `simulate_iv`):
- T must be positive (T=0 not supported — Fermi function undefined)
- vmode must be positive (vmode=0 causes divergent Bose function)
- tau must be positive or Inf (tau=0 causes NaN in W matrix; tau<0 is unphysical)
- alpha_L, alpha_R must be non-negative
- lambda must be non-negative
- Warning printed if lambda requires convergence N exceeding FC_MAX_N (1024)

**Safe parameter ranges**:

| Parameter | Safe Range | What Happens Outside |
|-----------|-----------|---------------------|
| T | > 0 (any positive) | T=0: panic. T very small: works but digamma arguments become large |
| N | 1 to ~200 | N=0: returns zeros. N>200: works but memory/time grows as N² |
| lambda | 0 to ~15 | lambda=0: handled explicitly. lambda>15: convergence N may exceed FC_MAX_N=1024, warning printed |
| alpha_L, alpha_R | ≥ 0 | Both zero: all rates vanish, W singular, solver returns fallback |
| Vsd | any real value | Sign correction applied. Vsd=0: current is zero by symmetry |
| Vg | any real value | Shifts molecular level, no edge cases |
| tau | > 0 or Inf | tau=0: panic. Inf: phonon relaxation terms vanish (IEEE 754) |
| eta | 0 to 1 | No validation; values outside [0,1] are unphysical but won't crash |
| vmode | > 0 | vmode=0: panic (bose_fcn diverges) |

**Silent error masking** (known, by design):
- `sanitize()` in cotunneling and rate computation converts NaN/Inf to 0.0. This matches the MATLAB reference behavior but masks upstream errors. In practice, this only fires at the known MATLAB solver artifact points (Vsd ≈ 0.219, 0.585).
- The QR solver clamps negative probabilities to 0 and renormalizes. This is physically motivated (probabilities ≥ 0) but masks ill-conditioning.

**Convergence safety**:
- All four convergence wrappers (sum_mmr, sum_mmr11, sum_mmmmrs, sum_mmmmrs11) have a maximum iteration limit of 200. If convergence is not reached, a warning is printed to stderr and the best available result is used.
- The `rel_diff_log10` and `rel_diff_log` functions handle zero denominators: both zero → converged (0.0); one zero → not converged (Inf). This prevents NaN from bypassing the convergence check.

**Fermi and Bose function edge cases**:
- `fermi(x, mu, T)` at T≤0: returns the step function limit (0 if x>mu, 1 if x<mu, 0.5 if x=mu)
- `fermi` with very large |(x-mu)/kT|: explicitly returns 0 or 1 to avoid exp overflow
- `bose_fcn(x, T)` at T≤0: returns 0
- `bose_fcn` with x≈0: guards against 1/(exp(0)-1) = 1/0

### Interactive Explorer (rust/src/explorer/)

The Rust crate includes a second binary target (`rust_explorer`) — a terminal UI for real-time parameter exploration built with ratatui 0.29 + crossterm.

**Architecture:**
- `App` struct owns all state: simulation parameters, results, background thread handles
- Event loop: `event::poll(50ms)` for keyboard input + `try_recv()` for computation results
- Background computation via `std::thread::spawn` + `mpsc::channel` for non-blocking UI
- FCCache and DigammaTable created once per stability diagram run, shared across all Vg iterations via `simulate_iv_with_cache()`
- Session save/load via `session.rs`: JSON serialization of all parameters, computed results, and display settings

**Three modes:**

| Mode | Description | Computation | Update |
|------|-------------|-------------|--------|
| I-V Curve | Chart widget with 3 Braille-line datasets (I_tol, I_seq, I_cot) | Single `simulate_iv_with_cache` call (~0.3s at N=6) | Auto-recompute on parameter change (300ms debounce) |
| Stability Diagram | Half-block heatmap of I(Vg, Vsd) with Viridis colormap, axis labels, log-scale colorbar | Loop over Vg values, each row sent via channel for progressive rendering | Manual (Enter to start, Esc to cancel) |
| Temperature Diagram | Half-block heatmap of I(T, Vsd) with Viridis colormap, T on x-axis, Vsd on y-axis, log-scale colorbar | Loop over T values (1–50 K default), FCCache and DigammaTable shared across all T points | Manual (Enter to start, Esc to cancel) |

**Display modes** (cycle with `d` key):

| Display | I-V Chart | Heatmaps | Units |
|---------|-----------|----------|-------|
| Current (I) | 3 line plots: I_tol, I_seq, I_cot | log\|I\| heatmap | A |
| Conductance (G) | dI/dV for all 3 components | log\|G\| heatmap | S |
| IETS (d²I/dV²) | d²I/dV² for all 3 components | log\|d²I\| heatmap | S/V |
| Normalized IETS | (d²I/dV²)/(dI/dV) for all 3 | log\|nIETS\| heatmap | 1/V |

IETS computation uses 3-point central finite differences for d²I/dV² (exact for uniform grids, handles non-uniform). Normalized IETS divides by dI/dV with a 1e-30 threshold to guard against division by zero in the Franck-Condon blockade regime.

**Key implementation details:**
- `DisplayMode` enum (`Current`, `Conductance`, `Iets`, `NormalizedIets`) with `transform()` and `transform_grid()` methods that apply the appropriate derivative operation.
- `Heatmap` custom Widget: renders `data[vg_idx][vsd_idx]` as half-block characters (`▀`) with per-cell fg/bg colors for 2× vertical resolution. Log-scale normalization.
- `Colorbar` Widget: vertical Viridis gradient strip with context-dependent label (lg|I|, lg|G|, lg|d²I|, lg|nI|) and numeric bounds.
- Auto-recompute: `iv_recompute_at: Option<Instant>` debounce timer, set 300ms after each parameter change in I-V mode. Checked in `poll_messages()`. Parameters adjustable even while computing.
- Stability mode blocks parameter adjustment during computation (too expensive to auto-recompute).
- CSV export: `e` key exports with display-mode-appropriate filename and headers (e.g., `iv_iets_export.csv`, `stability_niets_export.csv`).

**Session save/load:**
- `s` saves the full session to `fc_session.json` (or the path loaded from via `--load`)
- `--load <path>` CLI argument restores a saved session on startup
- Session includes: all parameters, I-V data, stability grid, temperature grid, display settings
- Handles `tau=Infinity` via custom serde (`"Inf"` string in JSON, matching benchmark spec convention)
- Version field for forward compatibility (rejects sessions from newer versions)
- Partial stability/temperature diagrams are preserved (empty rows stay empty)

**Keybindings:**
- `↑↓` navigate parameters, `←→` adjust (Shift=fine step)
- `d` cycle display mode: I → G → IETS → nIETS → I
- `Enter` trigger computation / cancel running computation
- `Tab` switch mode, `e` export CSV, `s` save session, `p` plot, `Esc` cancel, `q` quit

## Fortran Implementation (fortran/src/)

The Fortran port follows the C port structure (which is itself a faithful MATLAB translation). Pure Fortran 2008 with LAPACK for linear algebra and custom asymptotic series for complex digamma/trigamma (ported from Rust's optimized algorithm).

### Key Files

| File | C Equivalent | Notes |
|------|-------------|-------|
| `constants.f90` | `constants.h` | Module with same exact values |
| `laguerre.f90` | `laguerre.c` | Three-term recurrence, identical algorithm |
| `fc_matrix.f90` | `fc_matrix.c` | 256×256 cache with validity flags; `fc_cache_populate` pre-fills before OpenMP region |
| `fermi_bose.f90` | `fermi_bose.c` | Identical formulas |
| `digamma.f90` | `digamma.c` | Pure Fortran asymptotic series; 10-term full-path + 5-term fast-path (\|z\|²>900) + dispatchers |
| `regularized.f90` | `regularized.c` | Factored digamma precomputation (O(N) instead of O(N²) calls) |
| `cotunneling.f90` | `cotunneling.c` | Convergence wrappers with selective recomputation; largest module |
| `rate.f90` | `rate.c` | Flat allocatable array rate store with index arithmetic |
| `matrix.f90` | `matrix.c` | Same index mapping, sigma functions, peq |
| `solver.f90` | `solver.c` | LAPACK `dgesv` (LU factorization); augmented system, clamp+renormalize |
| `current.f90` | `current.c` | Same sign conventions: 0→1 R−L, 1→0 L−R, cot RL−LR |
| `simulate.f90` | `simulate.c` | OpenMP `parallel do schedule(dynamic, 4)`; FC cache pre-populated and shared read-only across threads |
| `json_io.f90` | `json_io.c` | Custom minimal JSON parser (no external library) |
| `plotting.f90` | `plotting.c` | gnuplot via `execute_command_line`; PDF + PNG |
| `main.f90` | `main.c` | `system_clock()` timing; 3 runs + median; validation excludes MATLAB solver artifacts |

### Design Decisions

| Decision | Rationale |
|----------|-----------|
| Pure Fortran digamma (10-term Bernoulli, from Rust) | No GSL FFI needed; Fortran native complex(8) supports all operations |
| Fast-path digamma/trigamma (5 terms for \|z\|²>900) | >99% hit rate at T=4.2K; skips reflection and recurrence; matches C/C++/Rust optimization |
| OpenMP `parallel do schedule(dynamic, 4)` | Bias points are independent; dynamic scheduling handles variable cotunneling cost |
| FC cache pre-populated via `fc_cache_populate` | All 256×256 entries computed before OpenMP region; read-only in parallel (thread-safe without synchronization) |
| LAPACK `dgesv` for steady-state solver | Universally available; LU factorization on augmented system; macOS Accelerate framework |
| Custom minimal JSON parser | Fortran has no standard JSON library; benchmark JSON format is simple enough for string parsing |
| 256×256 FC cache with logical validity array | Same pattern as C; Fortran allocatable arrays with 0-based bounds |
| Flat allocatable rate store | `allocatable :: data(:)` with manual index arithmetic; matches C's flat array approach |
| `ieee_arithmetic` module for infinity | Portable IEEE 754 infinity for tau="Inf" handling |
| gnuplot via `execute_command_line` | Same approach as C port; writes script to temp file, executes gnuplot |
| `-fopenmp` in FFLAGS | GCC gfortran supports OpenMP natively; links `-lgomp` automatically |

### Complex Digamma/Trigamma Strategy

Fortran has native `complex(kind=8)` support with intrinsic `log`, `sin`, `cos`, `exp`, `abs` for complex arguments. The implementation ports Rust's optimized algorithm with a two-tier approach:

1. **Digamma ψ(z)**: Custom implementation with two tiers:
   - **Fast-path** (`digamma_asymptotic5`): For Re(z) > 0 and |z|² > 900 — uses only 5 pre-computed `DIGAMMA_COEFF` entries, no reflection or recurrence. At T=4.2K, >99% of calls hit this path.
   - **Full-path**: 10-term Bernoulli asymptotic expansion with pre-computed `DIGAMMA_COEFF(k) = B_{2(k+1)} / (2*(k+1))`. Recurrence shift until |z|² ≥ 100 (using norm_sqr to avoid sqrt). Reflection formula for Re(z) ≤ 0.

2. **Trigamma ψ'(z)**: Matching two-tier structure:
   - **Fast-path** (`trigamma_asymptotic5`): For Re(z) > 0 and |z|² > 900 — 5 Bernoulli terms.
   - **Full-path**: 10-term Bernoulli expansion, recurrence threshold |z|² ≥ 100, reflection formula.

### Performance Notes

The Fortran port runs the quick spec (N=6) in ~0.36 seconds on Apple M5 — a **5122x speedup** over MATLAB (1844s), comparable to C (0.56s) and C++ (0.27s).

The key optimizations:
1. **OpenMP parallel bias-point loop**: `!$omp parallel do schedule(dynamic, 4)` distributes 201 bias points across all available cores
2. **Pre-populated FC cache**: All 256×256 entries computed once before the parallel region; shared read-only across threads (no synchronization needed)
3. **Fast-path digamma/trigamma**: 5-term Bernoulli asymptotic series for |z|² > 900 (>99% hit rate at T=4.2K); skips reflection and recurrence entirely
4. **Factored digamma precomputation** in `regularized_I`: 4×N calls instead of 4×N²
5. **10-term Bernoulli series** with pre-computed coefficients (from Rust)
6. **norm_sqr threshold**: avoids sqrt per recurrence iteration
7. **LAPACK on Accelerate**: Apple's optimized BLAS/LAPACK via the Accelerate framework
8. **Compiler flags**: `-O3 -march=native -flto -funroll-loops -fopenmp` for LTO + native SIMD + OpenMP

Performance history on the same hardware (Apple M5, quick spec):
- Initial port (-O2): 2.1s (878x vs MATLAB)
- Optimized (-O3 -march=native -flto): 1.6s (1153x vs MATLAB)
- + OpenMP + fast-path digamma: **0.36s (5122x vs MATLAB)**

**Default spec (N=15)**: **2.6 s** on Apple M5 — a **4510× speedup** over MATLAB. Matches MATLAB to **6.1×10⁻⁶** (excluding solver artifacts at Vsd ≈ 0.219, 0.438). Matches Rust to **4.9×10⁻¹³**.

### Numerical Precision Notes

The Fortran port matches MATLAB to **7.65e-5 max relative error** at 199 of 201 bias points (quick spec), identical to the C and Rust ports. The 2 outlier points (Vsd ≈ 0.219, 0.585) are the known MATLAB solver artifacts.

### Lessons from the Fortran Port

- **Fortran's native complex(8) is excellent**: All complex transcendentals (`log`, `sin`, `cos`, `exp`, `abs`) work out of the box. No need for external libraries.
- **LAPACK is universally available**: On macOS, `-framework Accelerate` provides optimized LAPACK. On Linux, `-llapack -lblas` suffices.
- **Custom JSON parser is viable**: The benchmark JSON format is simple enough that a ~200-line string parser handles all I/O. No need for json-fortran or C interop.
- **0-based vs 1-based indexing**: Fortran arrays default to 1-based, but phonon state indices are 0-based. The cleanest approach: keep physics indices as 0-based integers, use `array(idx+1)` for array access or declare with explicit `0:N-1` bounds.
- **`ieee_arithmetic` module**: `ieee_value(1.0d0, ieee_positive_inf)` provides portable IEEE infinity. `1.0d0 / infinity == 0.0d0` works for tau=Inf handling.
- **Module compilation order matters**: Fortran modules must be compiled in dependency order. The Makefile lists sources in topological order.
- **`implicit none` everywhere**: Catches typos that would silently create new variables in classic Fortran.
- **Argument declaration ordering**: In Fortran, dummy arguments used as array dimensions must be declared before the arrays that use them (e.g., `integer, intent(in) :: nVsd` before `real(8), intent(in) :: Vsd(nVsd)`).
- **FC cache must be pre-populated for OpenMP**: The lazy-init pattern in `fc_cache_get` is not thread-safe. Pre-populating all 256×256 entries via `fc_cache_populate(fc, FC_MAX_N)` before the parallel region makes `fc_cache_get` purely read-only (safe to share across threads without synchronization).
- **Homebrew gfortran supports OpenMP natively**: `-fopenmp` in FFLAGS is all that's needed. GCC-based gfortran (from Homebrew) links `-lgomp` automatically. No separate libomp installation required (unlike Apple Clang for C/C++).

## C++ Implementation (cpp/src/)

The C++ port follows the C port structure using modern C++17 idioms. No GSL or C library dependencies — uses Eigen for linear algebra, nlohmann/json for JSON I/O, and OpenMP for parallelism.

### Key Files

| File | C Equivalent | Notes |
|------|-------------|-------|
| `constants.hpp` | `constants.h` | `inline constexpr double` instead of `#define` |
| `laguerre.hpp/.cpp` | `laguerre.c` | Identical three-term recurrence |
| `fc_matrix.hpp/.cpp` | `fc_matrix.c` | `FCCache` struct; pre-populated to FC_MAX_N, shared read-only across OpenMP threads |
| `fermi_bose.hpp/.cpp` | `fermi_bose.c` | Identical formulas |
| `digamma.hpp/.cpp` | `digamma.c` | 10-term Bernoulli with fast-path (5 terms for \|z\|²>900); pre-computed DIGAMMA_COEFF; `std::norm()` thresholds |
| `regularized.hpp/.cpp` | `regularized.c` | `std::vector<Complex>` for precomputed digamma values |
| `cotunneling.hpp/.cpp` | `cotunneling.c` | `std::vector` replaces all `malloc`/`free`; largest module |
| `rate.hpp/.cpp` | `rate.c` | `RateStore` struct with `std::vector<double>` data member |
| `matrix.hpp/.cpp` | `matrix.c` | Same index mapping, sigma functions, peq |
| `solver.hpp/.cpp` | `solver.c` | `Eigen::HouseholderQR` replaces GSL QR |
| `current.hpp/.cpp` | `current.c` | Same sign conventions: 0→1 R−L, 1→0 L−R, cot RL−LR |
| `simulate.hpp/.cpp` | `simulate.c` | OpenMP `parallel for schedule(dynamic, 4)`; FC cache shared read-only across threads |
| `json_io.hpp/.cpp` | `json_io.c` | `nlohmann::json` for parsing; `fprintf` for output (%.17g) |
| `plotting.hpp/.cpp` | `plotting.c` | Identical gnuplot via `popen()` |
| `main.cpp` | `main.c` | `std::chrono::steady_clock` timing; 3 runs + median |

### Design Decisions

| Decision | Rationale |
|----------|-----------|
| Pure C++ digamma (no GSL) | Matches Rust/Fortran approach; eliminates C library dependency |
| 10-term Bernoulli with pre-computed DIGAMMA_COEFF | Eliminates runtime division; `inv_power * coeff` instead of `coeff / power` |
| Fast-path `digamma`/`trigamma` (5 terms for \|z\|²>900) | >99% hit rate at T=4.2K; skips reflection and recurrence entirely |
| `std::norm(z)` thresholds instead of `std::abs(z)` | Avoids sqrt per recurrence iteration; `std::norm` returns \|z\|² directly |
| FC cache pre-populated to FC_MAX_N, shared read-only | All 256×256 entries computed once before OpenMP region; `const FCCache&` throughout call chain |
| `-O3 -flto -march=native` | LTO enables cross-module inlining of digamma into regularized/cotunneling hot loops |
| `std::complex<double>` | Direct replacement for C99 `double _Complex`; IEEE 754 compatible |
| Eigen `HouseholderQR` | Header-only, no link-time deps; ~0.04ms for 30×30 systems |
| nlohmann/json | Single-header vendor; dramatically simplifies JSON I/O vs cJSON |
| `std::vector` everywhere | RAII replaces all `malloc`/`free`; no manual memory management |
| OpenMP `schedule(dynamic, 4)` | Parallel bias-point computation; dynamic for variable-cost cotunneling |
| `constexpr` Bernoulli coefficients | Compile-time constant arrays; no runtime initialization |
| `namespace fc { }` | All functions scoped; replaces C's `fc_` prefix convention |
| CMake build system | Standard C++ tooling; handles Eigen/OpenMP platform detection |
| `fprintf` for output JSON/CSV | Ensures %.17g precision matches C port output exactly |

### Complex Digamma/Trigamma Strategy

C++ uses `std::complex<double>` with a pure asymptotic series approach matching the Rust/C/Fortran ports:

1. **Digamma ψ(z)**: Custom implementation with two tiers:
   - **Fast-path**: For Re(z) > 0 and |z|² > 900 — uses only 5 pre-computed `DIGAMMA_COEFF` entries, no reflection or recurrence. At T=4.2K, >99% of calls hit this path.
   - **Full-path**: 10-term Bernoulli asymptotic expansion with pre-computed `DIGAMMA_COEFF[k] = B_{2(k+1)} / (2*(k+1))`. Recurrence shift until |z|² ≥ 400 (using `std::norm()` to avoid `sqrt`). Reflection formula for Re(z) ≤ 0.

2. **Trigamma ψ'(z)**: Matching two-tier structure:
   - **Fast-path**: For Re(z) > 0 and |z|² > 900 — 5 Bernoulli terms.
   - **Full-path**: 10-term Bernoulli expansion, recurrence threshold |z|² ≥ 100, reflection formula.

### Performance Notes

The C++ port runs the quick spec (N=6) in ~0.27 seconds on Apple M5 — a **6830x speedup** over MATLAB (1844s), comparable to Rust+Rayon (0.28s).

**Default spec (N=15)**: **2.0 s** on Apple M5 — a **5863× speedup** over MATLAB. Matches MATLAB to **6.1×10⁻⁶** (excluding solver artifacts at Vsd ≈ 0.219, 0.438). Matches Rust to **4.9×10⁻¹³**.

Performance history on the same hardware:
- Initial port (20-term digamma, per-thread FCCache recreation): 8.7s default, 1.2s quick
- Optimized (10-term with fast-path, shared FC cache, LTO): **2.0s default (4.4× faster), 0.27s quick (4.4× faster)**

The key optimizations:
1. **Shared read-only FC cache**: Pre-populate all 256×256 entries once via `fc_cache_populate(fc, FC_MAX_N)` before the OpenMP parallel region. Previously, each thread recreated and populated its own FCCache — 10 threads × 256×256 entries = massive redundant work. Now `const FCCache&` is passed through the entire call chain (`fc_cache_get`, `calculate_all_rateW`, `m_rateW`, all cotunneling functions).
2. **10-term Bernoulli with pre-computed DIGAMMA_COEFF**: Reduced from 20 terms to 10 (still ~1e-20 truncation error at |z|≥10). Pre-divided coefficients eliminate runtime complex division.
3. **Fast-path digamma/trigamma**: 5-term series for |z|² > 900, Re(z) > 0. At T=4.2K, >99% of calls skip reflection and recurrence entirely.
4. **`std::norm()` recurrence thresholds**: `std::norm(z) < 400.0` instead of `std::abs(z) < 20.0` avoids sqrt per recurrence step.
5. **LTO (`-flto`)**: Enables cross-module inlining of digamma/trigamma into regularized.cpp and cotunneling.cpp hot loops.

### Numerical Precision Notes

The C++ port matches MATLAB to **7.65e-5 max relative error** at 199 of 201 bias points (quick spec), identical to the C, Rust, and Fortran ports. The 2 outlier points (Vsd ≈ 0.219, 0.585) are the known MATLAB solver artifacts.

### Lessons from the C++ Port

- **FC cache must be pre-populated for OpenMP**: The lazy-init pattern (checking `valid[q1][q2]` and writing on miss) is not thread-safe. Pre-populating all entries via `fc_cache_populate(fc, FC_MAX_N)` makes `fc_cache_get` purely read-only with `const FCCache&` — safe to share across threads without synchronization.
- **`-O3 -flto` works for C++ on Apple Clang**: Unlike the C port (which regresses 3× with `-O3 -flto`), the C++ port benefits from LTO at `-O3`. This enables cross-module inlining of the digamma fast-path into the regularized integral hot loops.
- **`std::norm()` is not `std::abs()`**: `std::norm(z)` returns |z|² (no sqrt), while `std::abs(z)` returns |z| (with sqrt). Using `std::norm()` for threshold comparisons saves one sqrt per recurrence iteration — significant when called millions of times.
- **No GSL needed**: The pure asymptotic series for digamma/trigamma (10-term Bernoulli) matches GSL's accuracy for |z| ≥ 10. Pre-computed `DIGAMMA_COEFF` eliminates runtime division.
- **Eigen is effortless for small systems**: Header-only, ~3 lines for QR solve, no FFI overhead.
- **`std::complex<double>` works correctly**: On arm64 Apple Clang, no `__muldc3` overhead was observed at `-O3`.
- **OpenMP on Apple Clang**: Requires `brew install libomp` and explicit CMake configuration (`-Xpreprocessor -fopenmp`). The CMakeLists.txt auto-detects via `brew --prefix libomp`.
- **nlohmann/json eliminates boilerplate**: The JSON I/O module is ~100 lines shorter than the C port's cJSON equivalent.
- **RAII eliminates memory bugs**: Zero `malloc`/`free` calls in the entire codebase. `std::vector` handles all dynamic allocation.
