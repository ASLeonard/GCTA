# Plan: Hutch++ Trace Estimator for calcu_tr_PA

## Context
- `calcu_tr_PA` computes tr(P·Aᵢ) for each variance component i using exact Frobenius inner product: O(n²) per component.
- P is materialised as full n×n symmetric matrix `_P` (8n² bytes). At n=20k → 3.2 GB; n=50k → 20 GB.
- `_Vi` (V⁻¹) is already a class member. `Vi_X` and `Xt_Vi_X_i` are currently local vars in `reml_iteration`.
- Goal: implement Hutch++ using implicit P matvecs (no need to store _P), opt-in via flag.
- Primary benefit: memory (avoid storing n×n P); throughput is traded for memory efficiency.

## Steps

### Phase 1 — Promote P-matvec state
1. Add `eigenMatrix _Vi_X` and `eigenMatrix _Xt_Vi_X_i` as class members in `main/gcta.h` near `_Vi` (line 598).
2. In `calcu_P` (est_hsq.cpp ~L1545): after computing `Vi_X` and `Xt_Vi_X_i`, cache: `_Vi_X = Vi_X; _Xt_Vi_X_i = Xt_Vi_X_i;`

### Phase 2 — Implicit P matvec helpers
3. Declare in gcta.h (private section near calcu_tr_PA at line 340):
   - `eigenVector applyP_vec(const eigenVector& v) const;`
   - `eigenVector applyPA_vec(const eigenVector& v, int comp_idx) const;`
4. Implement `applyP_vec` in est_hsq.cpp:
   - `w = _Vi.selfadjointView<Lower>() * v`  → BLAS dsymv
   - `a = _Vi_X.transpose() * v`             → BLAS dgemv
   - `b = _Xt_Vi_X_i.selfadjointView<Lower>() * a` → tiny p×p dsymv
   - `return w - _Vi_X * b`                  → BLAS dgemv
5. Implement `applyPA_vec(v, comp_idx)`:
   - `w = _A[_r_indx[comp_idx]].selfadjointView<Lower>() * v` → dsymv
   - `return applyP_vec(w)`

### Phase 3 — Hutch++ implementation
6. Declare `void calcu_tr_PA_hutchpp(eigenVector& tr_PA, int m_probes);` in gcta.h.
7. Implement `calcu_tr_PA_hutchpp` in est_hsq.cpp:
   For each variance component i:
   - k = m_probes / 3  (equal-third split)
   - Seed `std::mt19937_64 rng(std::random_device{}())`, `std::uniform_int_distribution<int> coin{0,1}` → map 0→-1, 1→+1
   - **Phase A (range sketch):** Draw S ∈ Rⁿˣᵏ Rademacher; K.col(j) = applyPA_vec(S.col(j), i)
   - **QR:** `HouseholderQR<eigenMatrix>(K)` → extract thin Q (n×k orthonormal)
   - **Phase B (low-rank trace):** MQ.col(j) = applyPA_vec(Q.col(j), i); t_lr = Q.cwiseProduct(MQ).sum()
   - **Phase C (Hutchinson residual):** Draw G ∈ Rⁿˣᵏ Rademacher; MG.col(j) = applyPA_vec(G.col(j), i)
   - Gamma = Q.transpose() * G  (k×k)
   - R = G - Q * Gamma; MR = MG - MQ * Gamma
   - t_hutch = (1.0/k) * R.cwiseProduct(MR).sum()
   - tr_PA(i) = t_lr + t_hutch
   - Total: 3k matvecs per component; each matvec = 2 dsymv(n×n) + O(np) corrections

### Phase 4 — Integration and options
8. Add flags to gcta.h private section:
   - `bool _reml_trace_approx = false;`
   - `int  _reml_trace_approx_nprobes = 60;`  (default k=20, ~5% error)
9. Declare and implement `set_reml_trace_approx(bool, int)` setter in gcta.h/est_hsq.cpp.
10. In `calcu_tr_PA` (est_hsq.cpp ~L1718): in the `else` (dense) branch, add:
    `if (_reml_trace_approx) { calcu_tr_PA_hutchpp(tr_PA, _reml_trace_approx_nprobes); return; }`
11. In `main/option.cpp` (~L725 near `--reml-alg`): parse `--reml-trace-approx` (bool) and `--reml-trace-nprobes <N>` (int with bounds check N >= 9).
12. In `main/option.h`: add corresponding local variable declarations.
13. In dispatch section of option.cpp (~L1378): call `pter_gcta->set_reml_trace_approx(...)`.

## Relevant Files
- `main/gcta.h` — add `_Vi_X`, `_Xt_Vi_X_i`, `_reml_trace_approx`, `_reml_trace_approx_nprobes`; declare 4 new methods near line 338/598
- `main/est_hsq.cpp` — modify `calcu_P` to cache state (L1545); implement `applyP_vec`, `applyPA_vec`, `calcu_tr_PA_hutchpp`; modify `calcu_tr_PA` (L1750, dense branch)
- `main/option.cpp` — add `--reml-trace-approx`, `--reml-trace-nprobes` near L725 and dispatch at L1378
- `main/option.h` — add option variable declarations

## Verification
1. For small n (n=500), run REML with and without `--reml-trace-approx`; verify variance component estimates agree within 1% across multiple seeds
2. Check convergence iteration counts are similar (trace noise can add ~1-3 extra AI-REML iterations)
3. For EM-REML (`--reml-alg 2`): note this is more sensitive to trace error; warn in LOGGER when used together
4. Memory profiling: at n=20k, confirm ~3 GB saving by skipping `_P` materialisation (Phase 4 extension: skip `_P` assembly when flag is set)
5. Test existing unit tests pass (tests/run_tests.py)

## Decisions
- Implicit P via `_Vi` + `_Vi_X` + `_Xt_Vi_X_i` (all already computed in reml_iteration cycle)
- M = P·Aᵢ applied directly (non-symmetric); generalized Hutch++ theorem applies
- Bivariate/within-family (`_bivar_reml`, `_within_family`) always uses existing exact sparse path — no change
- Default m_probes=60 (k=20) balances accuracy (~5% relative error) vs 20× extra bandwidth vs exact
- Log a warning when `--reml-trace-approx` is active so users know results are approximate

## Further Considerations
1. **P materialisation in ai_reml/em_reml**: The calls `P.selfadjointView<Lower>() * _y` and `P.selfadjointView<Lower>() * APy.col(i)` in ai_reml remain unchanged for now. Full memory savings require replacing these P matvecs with `applyP_vec` too, plus avoiding `calcu_P` writing to `_P`. This is a natural Phase 5.
2. **Parallelism**: outer `for (int i = 0; ...) _r_indx.size()` over components can be `#pragma omp parallel for` since components are independent; inner k-loop over probe vectors can also be parallelised.
3. **Convergence noise**: Hutch++ introduces ~5% trace noise with k=20. AI-REML's Fisher information step absorbs this; EM-REML is sensitive. Consider emitting a recommendation: prefer `--reml-alg 0` (AI) when using `--reml-trace-approx`.

## Reference
Meyer, Musco, Musco, Woodruff (2021). "Hutch++: Optimal Stochastic Trace Estimation."
SIAM SOSA 2021. arXiv:2010.09649.
