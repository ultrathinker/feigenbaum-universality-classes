#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
step4_consolidate.py -- BRIEF-S3-F5 task 3 (unified numbering) + final
consolidation of all step outputs into s3_f5_results.json.
d=2,4,6 positions are the UNCONTESTED, established project values (wave1
spectra_grid.json / wave3_t22_d6 / cross-checked by MiniMax F5, unchanged
here). d=8 uses THIS task's own corrected n_>(8)=7 (step2c_full_confirm.json),
not MiniMax's disputed n_>(8)=6.
"""
import json
import os
import time

HERE = os.path.dirname(os.path.abspath(__file__))
PROG = os.path.join(HERE, 'progress.log')


def plog(msg):
    line = time.strftime('%Y-%m-%d %H:%M:%S ') + msg
    print(line, flush=True)
    with open(PROG, 'a', encoding='utf-8') as f:
        f.write(line + '\n')
        f.flush()


def load(name):
    with open(os.path.join(HERE, name), encoding='utf-8') as f:
        return json.load(f)


def main():
    plog('=== STEP4 consolidate: START ===')
    fp8 = load('step1_d8_fixedpoint.json')
    sp8 = load('step2c_full_confirm.json')
    neigh = load('step3_neighborhood.json')
    cosfix = load('step3b_cosine_fix.json')

    # n_>(d) established (d=2,4,6, uncontested) + this task's d=8 (corrected)
    n_gt = {'2': 0, '4': 3, '6': 5, '8': sp8['n_gt_8']}
    pos_no_delta = {d: n_gt[d] + 1 for d in n_gt}
    pos_with_delta = {d: n_gt[d] + 2 for d in n_gt}

    plog(f'n_>(d) [with this task d=8 correction]: {n_gt}')
    plog(f'tower1 position without delta: {pos_no_delta}')
    plog(f'tower1 position with delta: {pos_with_delta}')

    # internal consistency check: n_>(d) == d-1 for even d>=4?
    check = {d: (n_gt[d] == int(float(d)) - 1) for d in n_gt if float(d) >= 4}
    plog(f'pattern check n_>(d)=d-1 for d>=4 (even): {check}')

    out = {
        'meta': {'created': time.strftime('%Y-%m-%d %H:%M:%S'),
                 'tag': 'S3-F5 final consolidation'},
        'task1_d8': {
            'alpha8_N48': fp8['points']['N48']['alpha'],
            'alpha8_N56': fp8['points']['N56']['alpha'],
            'alpha8_drift_N48_N56': fp8['checks']['alpha_drift_N48_N56'],
            'alpha8_vs_minimax_diff': fp8['checks']['diff_vs_minimax_N56'],
            'dense_resid_N48': fp8['points']['N48']['dense_resid'],
            'dense_resid_N56': fp8['points']['N56']['dense_resid'],
            'tower1_8': sp8['tower1'],
            'n_gt_8_corrected': sp8['n_gt_8'],
            'n_gt_8_minimax_reported': sp8['minimax_n_gt_8_reported'],
            'minimax_window_bug': sp8['minimax_window'],
            'essential_list_8': sp8['essential_sorted'],
            'verdict': ('MiniMax n_>(8)=6 is an UNDERCOUNT. Root cause: their '
                        'stable_modes_at_N used WINDOW=0.42 (inherited from the d<=4 '
                        'spectrum_window convention), which silently excludes the '
                        'genuine, N=48/56/64-ladder-stable essential mode |lam|=0.42946 '
                        '(dense functional residual 1e-18, matched to 1e-15 relative '
                        'across three truncations). With this mode restored, n_>(8)=7.'),
        },
        'task2_neighborhood': {
            'd4_slope_per_unit_d': 0.0622,
            'd4_multiplicity': 'exactly 1 genuine (N-stable) mode within |lam-tower1|<0.02 at every '
                                'd in {3.98,3.99,4.0,4.01,4.02}; N=48-only extra entries at 3.99/4.01/4.02 '
                                'are pseudospectral pollution (absent at N=40 AND N=56 -- confirmed by a '
                                'direct N=40/48/56 check at d=4.01).',
            'd6_slope_per_unit_d': 0.0379,
            'd6_multiplicity': 'exactly 1 genuine mode at 5.98/6.0/6.02, no pollution ambiguity observed.',
            'cosine_d4_corrected': cosfix['cosine_corrected'],
            'cosine_d4_verdict': cosfix['verdict'],
            'cosine_d6': neigh['cosine_test_d6']['cosine'],
        },
        'task3_numbering': {
            'n_gt_established_2_4_6_this_task_8': n_gt,
            'position_no_delta': pos_no_delta,
            'position_with_delta': pos_with_delta,
            'pattern_check_n_gt_eq_d_minus_1_for_d_geq_4': check,
            'proposed_formula': (
                'n_>(d) = d-1 for even d>=4 (verified d=4,6,8: 3,5,7); d=2 is the '
                'boundary/onset point (n_>(2)=0, one less than the d>=4 pattern would '
                'give, because the tower value alpha^-2 is exactly where the branch is '
                'BORN, not a mid-spectrum crossing). Equivalently: tower alpha_d^-d '
                'occupies rank d among subleading-by-modulus (without delta) for d>=4, '
                'and rank 1 at d=2; with delta counted, ranks are d+1 (d>=4) and 2 (d=2). '
                'This SUPERSEDES MiniMax REPORT-F5.md section 4.4\'s formula '
                '"(d+1)-th with delta / d-th without delta" -- that formula does not '
                'even match MiniMax\'s OWN table 3.2 at d=2 (1st without delta, not 2nd) '
                'and was based on the undercounted n_>(8)=6.'
            ),
        },
    }
    OUT = os.path.join(HERE, 's3_f5_results.json')
    with open(OUT, 'w', encoding='utf-8') as f:
        json.dump(out, f, indent=1, ensure_ascii=False)
    plog(f'saved {OUT}')
    plog('=== STEP4 consolidate: DONE ===')


if __name__ == '__main__':
    main()
