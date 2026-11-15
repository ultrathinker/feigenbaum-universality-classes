#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
p2_build_zeros.py -- build the degree-2 zero-of-H catalog for agent A's
lambda2 retro-check.  [oid:lam2A]  NOT RIGOROUS.

*** REVISION 2 (post-discovery, see progress.log): a first run of Phase 2 with
the raw results.json "u_cat" (== zero_map.json "representative_u") showed
d_u(m) INCREASING and SATURATING to a nonzero constant (~2-4e-5) over
m=8..11 for every complex zero, instead of decaying.  Root cause found:
results.json's u_cat field is numerically IDENTICAL to zero_map.json's
representative_u, which is itself just the level_trajectories[m=7] value
(the initial float scan hit), NOT a converged/extrapolated H-zero.  The
archive's own level_trajectories (m=7..10) show u(m) still drifting by
~1e-5..1e-6 per level at m=10, geometrically at rate ~1/delta -- i.e. u_cat
carries a leftover truncation bias of the SAME order as the d_u(m) signal
this whole study is trying to measure.  Fix: Richardson-extrapolate EVERY
zero (real AND complex) from zero_map.json's m=8,9,10 level trajectories to
m=infinity using the precise delta, exactly as already planned for the real
zeros -- this is now done uniformly for the whole catalog. ***

READS (read-only, per BRIEF rules):
  ../../06-EXTERNAL-INPUTS/study1_d2_frozen_experiment/results.json
      -- catalog membership (which 22+1 disks are "accepted"/certified) and
         Z72's kilo_xcheck-corrected true position (stage3b_z72_corrected).
         The raw u_cat field itself is NOT used as ground truth anymore (see
         above) -- kept only for logging/comparison (bias_vs_raw_u_cat).
  D:\\Исследования\\И-1 — Нули функции H\\Feigenbaum-archive\\ideation_sweep1\\codex\\zero_map.json
      -- PRIMARY source for actual zero positions: level_trajectories at
         m=7,8,9,10 for every root (real AND complex), Richardson-extrapolated
         (complex Aitken-style, using the precise delta as the known
         geometric ratio) to project u(m) to m=infinity.

WRITES: zeros_deg2.json (into agentA_lambda2\\, the only writable directory).
"""
import json
import time

from mpmath import mp, mpf, mpc, nstr, fabs

mp.dps = 60

TAG = "[oid:lam2A] NOT RIGOROUS float pilot (agent A), phase 2 zero catalog, REVISION 2 (Richardson-corrected)"
T0 = time.time()
PROGRESS_LOG = 'progress.log'


def log(msg):
    line = "[%7.1fs] %s" % (time.time() - T0, msg)
    print(line, flush=True)
    with open(PROGRESS_LOG, 'a', encoding='utf-8') as f:
        f.write(time.strftime('[%Y-%m-%d %H:%M:%S] ') + '[phase2-build-v2] ' + msg + '\n')


RESULTS_JSON = r'../../06-EXTERNAL-INPUTS/study1_d2_frozen_experiment/results.json'
ZEROMAP_JSON = r'D:\Исследования\И-1 — Нули функции H\Feigenbaum-archive\ideation_sweep1\codex\zero_map.json'

# precise delta (Investigation-1 meta, same constant as calibration -- NOT recomputed)
DELTA = mpf('4.669201609102990671853203820466201617258185577475768632745651343004133555868548426390712381266351242')


def richardson_extrap_complex(root):
    """Complex Richardson (geometric-tail) extrapolation of level_trajectories[m]['u']
    to m=infinity, using the KNOWN precise delta as the geometric ratio (1/delta) --
    NOT the noisy empirical shift ratio, which is only float64-accurate."""
    lt = {e['m']: mpc(mpf(e['u'][0]), mpf(e['u'][1])) for e in root['level_trajectories']}
    ms = sorted(lt)
    if len(ms) < 3:
        return None
    m_hi = ms[-1]
    u_hi, u_mid, u_lo = lt[ms[-1]], lt[ms[-2]], lt[ms[-3]]
    shift_hi = u_hi - u_mid
    shift_mid = u_mid - u_lo
    ratio_obs = shift_hi / shift_mid if shift_mid != 0 else None
    rL = 1 / DELTA
    u_inf = u_hi + shift_hi * rL / (1 - rL)
    return {'u_inf': u_inf, 'u_levels': lt, 'shift_hi': shift_hi, 'shift_mid': shift_mid,
            'ratio_obs': ratio_obs, 'm_hi': m_hi,
            'residual_estimate': abs(shift_hi) * rL ** 2 / (1 - rL)}


def main():
    log(TAG)
    d = json.load(open(RESULTS_JSON, encoding='utf-8'))
    sc = d['stage3_catalog']
    raw_u_cat = {}
    for z in sc:
        if z['accepted'] or z['id'] == 'Z72':
            raw_u_cat[z['id']] = mpc(mpf(z['u_cat']['re_str']), mpf(z['u_cat']['im_str']))
    log('loaded %d raw u_cat reference values from results.json (for bias logging only)' % len(raw_u_cat))

    zm = json.load(open(ZEROMAP_JSON, encoding='utf-8'))
    roots = {r['id']: r for r in zm['roots']}
    log('loaded %d roots from zero_map.json (primary position source)' % len(roots))

    complex_chosen_ids = ['Z01', 'Z02', 'Z05', 'Z17', 'Z19', 'Z39', 'Z40', 'Z57', 'Z63', 'Z70', 'Z74']
    real_chosen_ids = ['Z06', 'Z44', 'Z50']

    def extrapolate_and_log(zid):
        r = roots[zid]
        ex = richardson_extrap_complex(r)
        bias = None
        if zid in raw_u_cat:
            bias = abs(ex['u_inf'] - raw_u_cat[zid])
        log('  %-6s u(m=%d)=%s+%sj  ratio_obs=%s  Richardson u_inf=%s+%sj  '
            'residual_est~%.2e  |u_inf - raw_u_cat|=%s'
            % (zid, ex['m_hi'], nstr(ex['u_levels'][ex['m_hi']].real, 12),
               nstr(ex['u_levels'][ex['m_hi']].imag, 12),
               nstr(ex['ratio_obs'], 6) if ex['ratio_obs'] else '-',
               nstr(ex['u_inf'].real, 14), nstr(ex['u_inf'].imag, 14),
               float(ex['residual_estimate']),
               ('%.4e' % float(bias)) if bias is not None else 'n/a'))
        return ex

    log('=== complex zeros: Richardson-extrapolating m->infinity (fixes the m=7-seed bias in raw u_cat) ===')
    complex_zeros = []
    for zid in complex_chosen_ids:
        ex = extrapolate_and_log(zid)
        complex_zeros.append({'id': zid, 're': ex['u_inf'].real, 'im': ex['u_inf'].imag,
                               'abs_u': float(fabs(ex['u_inf'])),
                               'arg_deg': float(mp.arg(ex['u_inf']) * 180 / mp.pi),
                               'richardson_residual_estimate': float(ex['residual_estimate']),
                               'bias_vs_raw_u_cat': float(abs(ex['u_inf'] - raw_u_cat[zid]))})

    log('=== real zeros: Richardson-extrapolating m->infinity ===')
    real_zeros = []
    for zid in real_chosen_ids:
        ex = extrapolate_and_log(zid)
        real_zeros.append({'id': zid, 're': ex['u_inf'].real, 'im': mpf(0),
                            'abs_u': float(ex['u_inf'].real),
                            'richardson_residual_estimate': float(ex['residual_estimate']),
                            'bias_vs_raw_u_cat': float(abs(ex['u_inf'] - raw_u_cat.get(zid, ex['u_inf'])))})

    # Z72: BRIEF-mandated true position (kilo_xcheck-corrected, NOT the certified disk
    # centre, NOT re-extrapolated -- its ground truth is defined differently, see note).
    z72c = d['stage3b_z72_corrected']
    z72_true = mpc(mpf(z72c['u_cat']['re_str']), mpf(z72c['u_cat']['im_str']))
    log('Z72 TRUE position (kilo_xcheck-corrected, per article Sec 7.6, BRIEF-mandated -- '
        'NOT the certified disk centre, NOT Richardson-extrapolated): u = %s + %sj'
        % (nstr(z72_true.real, 20), nstr(z72_true.imag, 20)))
    log('  NOTE: this position is itself the m=10 superstable-cascade match point (cross-validated '
        'by an independent implementation) -- used as ground truth per BRIEF, but its own m=10 '
        'residual is a construction artifact (~0); Z72 EXCLUDED from ratio statistics, kept only as '
        'a qualitative diagnostic at m=8,9,11.')

    # ---- assemble output catalog -------------------------------------------------------
    catalog = []
    for c in complex_zeros:
        catalog.append({'id': c['id'], 'kind': 'complex_richardson_extrapolated',
                         're': nstr(c['re'], 25), 'im': nstr(c['im'], 25),
                         'abs_u': c['abs_u'], 'arg_deg': c['arg_deg'],
                         'richardson_residual_estimate': c['richardson_residual_estimate'],
                         'bias_vs_raw_u_cat': c['bias_vs_raw_u_cat'],
                         'source': 'Feigenbaum-archive/ideation_sweep1/codex/zero_map.json (m=8,9,10 level trajectories, Richardson-extrapolated to m=inf using precise delta); raw results.json u_cat (=m=7 seed) kept only for bias comparison'})
    for r in real_zeros:
        catalog.append({'id': r['id'], 'kind': 'real_richardson_extrapolated',
                         're': nstr(r['re'], 25), 'im': '0',
                         'abs_u': r['abs_u'], 'arg_deg': 0.0,
                         'richardson_residual_estimate': r['richardson_residual_estimate'],
                         'bias_vs_raw_u_cat': r['bias_vs_raw_u_cat'],
                         'source': 'Feigenbaum-archive/ideation_sweep1/codex/zero_map.json (m=8,9,10 level trajectories, Richardson-extrapolated to m=inf using precise delta)'})
    catalog.append({'id': 'Z72_true', 'kind': 'z72_true_diagnostic_only',
                     're': nstr(z72_true.real, 25), 'im': nstr(z72_true.imag, 25),
                     'abs_u': float(fabs(z72_true)), 'arg_deg': None,
                     'source': 'Investigation-1/.../results.json stage3b_z72_corrected (kilo_xcheck-verified true position, article Sec 7.6)',
                     'note': 'EXCLUDED from ratio statistics (construction is m=10-circular); qualitative check only'})

    out = {
        'meta': {'tag': TAG, 'date': time.strftime('%Y-%m-%d %H:%M:%S'), 'not_rigorous': True,
                 'n_total': len(catalog), 'n_complex': len(complex_zeros), 'n_real': len(real_zeros),
                 'delta_used_for_richardson': nstr(DELTA, 50),
                 'revision_note': 'REVISION 2: all zeros (real+complex) now Richardson-extrapolated from zero_map.json level trajectories to m=infinity; raw results.json u_cat was found to equal the m=7 seed value, biased by ~1e-5..1e-4 relative to the converged H-zero -- see bias_vs_raw_u_cat per zero and progress.log for the discovery trail.',
                 'elapsed_sec': time.time() - T0},
        'zeros': catalog,
    }
    with open('zeros_deg2.json', 'w', encoding='utf-8') as f:
        json.dump(out, f, indent=2, ensure_ascii=False)
    log('wrote zeros_deg2.json (%d zeros: %d complex, %d real, +Z72 diagnostic) [REVISION 2]'
        % (len(catalog), len(complex_zeros), len(real_zeros)))


if __name__ == '__main__':
    main()
