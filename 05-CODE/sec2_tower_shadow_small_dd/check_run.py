# -*- coding: utf-8 -*-
"""Full independent sweep -- see check_shadow_d6.py for the method. Writes
check_results.json + prints a progress log. Traces branch B+ by CONTINUITY
(value-cosine to previous step), not by "nearest |lambda| in a window"."""
import json
import time
import sys

from mpmath import mp, mpf

mp.dps = 40
sys.path.insert(0, r".")
from check_shadow_d6 import (
    load_seed, solve_fp, newton_continue, full_eig_window,
    basisA_vals, basisB_vals, values_monomial, values_chebB,
    cos_raw, cos_values, PHYS_GRID, log,
)

T0 = time.time()


def elapsed():
    return time.time() - T0


def get_fp(d_mpf, N, seed_free, seed_alpha):
    a_free, alpha, it, r = solve_fp(d_mpf, N, seed_free, seed_alpha, ladder_from=min(32, N))
    return a_free, alpha, r


def find_by_lambda(modes, target_lam, tol_rel=0.05):
    best = None
    bestd = None
    for m in modes:
        dd = abs(m["lam"] - target_lam)
        if bestd is None or dd < bestd:
            bestd = dd
            best = m
    return best, bestd


def trace_step(modes_new, prev_vec_values_fn, prev_values, grid):
    """pick the mode in modes_new whose function VALUES on the common grid
    have highest cos with prev_values (continuity tracing)."""
    best = None
    bestc = mpf(-1)
    for m in modes_new:
        if m["complex"]:
            continue
        vals = prev_vec_values_fn(m["vec"], grid)
        c = cos_values(prev_values, vals)
        if c > bestc:
            bestc = c
            best = m
    return best, bestc


def run_basis(basis_name, basis_vals_fn, values_fn, N):
    log(f"=== basis {basis_name}  N={N} ===")
    a48_4, alpha48_4 = load_seed("4.0")
    a48_6, alpha48_6 = load_seed("6.0")

    # --- base points d=4.0, d=6.0 ---
    seed4 = a48_4[1:N + 1] + [mpf(0)] * max(0, N - (len(a48_4) - 1))
    seed6 = a48_6[1:N + 1] + [mpf(0)] * max(0, N - (len(a48_6) - 1))

    af4, al4, r4 = get_fp(mpf("4.0"), N, seed4, alpha48_4)
    a4 = [mpf(1)] + af4
    log(f"fp d=4.0 N={N} resid={r4} alpha={al4}")
    modes4 = full_eig_window(a4, al4, mpf("4.0"), N, basis_vals_fn, basis_name)
    tower1_4 = al4 ** mpf(-4)
    tower2_4 = al4 ** mpf(-8)
    phi1_4, d1 = find_by_lambda(modes4, tower1_4)
    phi2_4, d2 = find_by_lambda(modes4, tower2_4)
    log(f"  Phi1(d=4): lam={phi1_4['lam']} |d-tower|={d1}")
    log(f"  Phi2(d=4): lam={phi2_4['lam']} |d-tower|={d2}")
    phi1_4_vals = values_fn(phi1_4["vec"], PHYS_GRID)
    phi2_4_vals = values_fn(phi2_4["vec"], PHYS_GRID)

    af6, al6, r6 = get_fp(mpf("6.0"), N, seed6, alpha48_6)
    a6 = [mpf(1)] + af6
    log(f"fp d=6.0 N={N} resid={r6} alpha={al6}")
    modes6 = full_eig_window(a6, al6, mpf("6.0"), N, basis_vals_fn, basis_name)
    tower1_6 = al6 ** mpf(-6)
    tower2_6 = al6 ** mpf(-12)
    phi1_6, e1 = find_by_lambda(modes6, tower1_6)
    phi2_6, e2 = find_by_lambda(modes6, tower2_6)
    log(f"  Phi1(d=6): lam={phi1_6['lam']} |d-tower|={e1}")
    log(f"  Phi2(d=6): lam={phi2_6['lam']} |d-tower|={e2}")
    phi1_6_vals = values_fn(phi1_6["vec"], PHYS_GRID)
    phi2_6_vals = values_fn(phi2_6["vec"], PHYS_GRID)

    results = {"N": N, "basis": basis_name, "d4": {}, "d6": {}}

    # --- ladder d=4: 4.0 -> 4.001 -> 4.01, continuity tracing ---
    cur_free, cur_alpha = af4, al4
    prev_values = phi1_4_vals
    for dstr in ["4.001", "4.01"]:
        dv = mpf(dstr)
        cur_free, cur_alpha, it, r = newton_continue(dv, N, cur_free, cur_alpha)
        a_cur = [mpf(1)] + cur_free
        modes = full_eig_window(a_cur, cur_alpha, dv, N, basis_vals_fn, basis_name)
        best, c_prev = trace_step(modes, values_fn, prev_values, PHYS_GRID)
        vals_best = values_fn(best["vec"], PHYS_GRID)
        cA_phi1 = cos_raw(best["vec"], phi1_4["vec"])
        cA_phi2 = cos_raw(best["vec"], phi2_4["vec"])
        cL_phi1 = cos_values(vals_best, phi1_4_vals)
        cL_phi2 = cos_values(vals_best, phi2_4_vals)
        log(f"  d={dstr} N={N} [{basis_name}] lam={best['lam']} cos_prev={c_prev} "
            f"cos_raw_Phi1={cA_phi1} cos_L2_Phi1={cL_phi1} cos_L2_Phi2={cL_phi2} resid_fp={r}")
        results["d4"][dstr] = {
            "lam": mp.nstr(best["lam"], 15),
            "cos_continuity_prev": mp.nstr(c_prev, 8),
            "cos_raw_Phi1": mp.nstr(cA_phi1, 8),
            "cos_raw_Phi2": mp.nstr(cA_phi2, 8),
            "cos_L2_Phi1": mp.nstr(cL_phi1, 8),
            "cos_L2_Phi2": mp.nstr(cL_phi2, 8),
        }
        prev_values = vals_best

    # --- ladder d=6: 6.0 -> 6.0005 -> 6.001 -> 6.002 -> 6.005 -> 6.01, continuity tracing ---
    cur_free, cur_alpha = af6, al6
    prev_values = phi1_6_vals
    for dstr in ["6.0005", "6.001", "6.002", "6.005", "6.01"]:
        dv = mpf(dstr)
        cur_free, cur_alpha, it, r = newton_continue(dv, N, cur_free, cur_alpha)
        a_cur = [mpf(1)] + cur_free
        modes = full_eig_window(a_cur, cur_alpha, dv, N, basis_vals_fn, basis_name)
        best, c_prev = trace_step(modes, values_fn, prev_values, PHYS_GRID)
        vals_best = values_fn(best["vec"], PHYS_GRID)
        cA_phi1 = cos_raw(best["vec"], phi1_6["vec"])
        cA_phi2 = cos_raw(best["vec"], phi2_6["vec"])
        cL_phi1 = cos_values(vals_best, phi1_6_vals)
        cL_phi2 = cos_values(vals_best, phi2_6_vals)
        log(f"  d={dstr} N={N} [{basis_name}] lam={best['lam']} cos_prev={c_prev} "
            f"cos_raw_Phi1={cA_phi1} cos_L2_Phi1={cL_phi1} cos_L2_Phi2={cL_phi2} resid_fp={r}")
        results["d6"][dstr] = {
            "lam": mp.nstr(best["lam"], 15),
            "cos_continuity_prev": mp.nstr(c_prev, 8),
            "cos_raw_Phi1": mp.nstr(cA_phi1, 8),
            "cos_raw_Phi2": mp.nstr(cA_phi2, 8),
            "cos_L2_Phi1": mp.nstr(cL_phi1, 8),
            "cos_L2_Phi2": mp.nstr(cL_phi2, 8),
        }
        prev_values = vals_best

    log(f"basis {basis_name} N={N} done, elapsed={elapsed():.1f}s")
    return results


if __name__ == "__main__":
    all_results = {}
    for N in [48, 64]:
        for basis_name, bvals, vvals in [("A_monomial", basisA_vals, values_monomial),
                                          ("B_chebshift", basisB_vals, values_chebB)]:
            key = f"{basis_name}_N{N}"
            all_results[key] = run_basis(basis_name, bvals, vvals, N)
            with open("check_results.json", "w", encoding="utf-8") as f:
                json.dump(all_results, f, indent=2, ensure_ascii=False)
            log(f"saved partial results after {key}, total elapsed={elapsed():.1f}s")
    log(f"ALL DONE, total elapsed={elapsed():.1f}s")
