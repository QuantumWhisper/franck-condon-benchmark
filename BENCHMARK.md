# Benchmark Specification

This document is the contract all language implementations must follow. Deviating from it makes results incomparable.

## Objective

Compare wall-clock time and numerical accuracy of the Franck-Condon I-V simulation across languages and compilers. The goal is an apples-to-apples comparison of the core numerical work, not I/O or startup overhead.

## Input Parameters

All implementations read from `benchmark/spec/default_params.json`. Do not hardcode parameters.

## What Is Timed

The entire I-V sweep loop, from the first bias point to the last. This includes:

- Rate matrix construction (FC matrix elements, sequential + cotunneling rates)
- Steady-state probability solve (linear system 0 = WP)
- Current calculation (I_seq + I_cot)

Start the timer immediately before the first bias point. Stop it immediately after the last current value is computed.

## What Is NOT Timed

- File I/O, JSON parsing, result saving
- Path setup, imports, module loading
- Compilation (for JIT languages, a warm-up run must be completed before timing starts)
- Plotting or visualization

## Output Format

Save results as JSON with the following fields:

```json
{
  "language": "string",
  "version": "string",
  "wall_time_seconds": 0.0,
  "parameters": { ... },
  "Vsd": [...],
  "I_tol": [...],
  "I_seq": [...],
  "I_cot": [...]
}
```

Save to `benchmark/results/{language}_results.json`.

## Correctness Validation

All implementations must produce I-V curves matching the MATLAB reference within a relative tolerance of 1e-10. The comparison formula is:

```
max(abs(I_impl - I_ref) ./ max(abs(I_ref), 1e-30)) < 1e-10
```

Run this check before reporting timing results. A fast but wrong implementation is not a valid entry.

## Benchmark Protocol

1. **Cold start**: Fresh process or session, no cached results from a previous run.
2. **Warm-up** (JIT languages only): Complete one full run and discard the result before starting the timer. This accounts for JIT compilation overhead.
3. **Timed runs**: Minimum 3 runs. Report the median wall time.
4. **Hardware**: Report CPU model, RAM, and OS version alongside results.

## Fair Comparison Notes

- **Single-threaded** comparison is the primary metric. This is the apples-to-apples baseline.
- Parallel versions may be reported separately, with thread count noted explicitly.
- No GPU acceleration. The algorithm has sequential dependencies across the bias sweep that make GPU offloading impractical for this problem size.
- Use the same `default_params.json` for all runs. Do not tune N or other parameters to make your language look faster.
