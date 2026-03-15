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

## How to Add a New Language Port

1. Create a `{language}/` directory at the project root.
2. Read all parameters from `benchmark/spec/default_params.json`. Do not hardcode.
3. Implement the simulation pipeline following the architecture diagram above.
4. Output results in the JSON format specified in `BENCHMARK.md`.
5. Validate against the MATLAB reference results using the tolerance check in `BENCHMARK.md`.
6. Follow the benchmark protocol (warm-up for JIT, minimum 3 timed runs, report median).

The MATLAB source in `matlab/src/` is the ground truth for numerical behavior. When in doubt about a formula or sign convention, read the MATLAB code and the Koch et al. paper.
