#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""F8: Consolidate all results into f8_results.json."""
import os
import json
import time

ROOT_IN = r'D:/Исследования/И-3 — Классы универсальности/Исследование-3-Принцип'
OUT_DIR = r'.'
T0 = time.time()
LOG_PATH = os.path.join(OUT_DIR, 'progress.log')


def log(msg):
    line = '[%7.1fs] %s' % (time.time() - T0, msg)
    print(line, flush=True)
    with open(LOG_PATH, 'a', encoding='utf-8') as f:
        f.write(line + '\n')


def main():
    log('--- consolidate ---')
    R = json.load(open(os.path.join(OUT_DIR, 'rayleigh_results.json'), encoding='utf-8'))
    Q = json.load(open(os.path.join(OUT_DIR, 'quasiperiodicity_results.json'), encoding='utf-8'))
    Z = json.load(open(os.path.join(OUT_DIR, 'zeros_vs_m.json'), encoding='utf-8'))
    L = json.load(open(os.path.join(OUT_DIR, 'local_nodes.json'), encoding='utf-8'))

    out = dict(meta=R['meta'])
    out['rayleigh'] = dict(
        peaks=R['peaks'],
        coarse_peaks=dict(
            d6_re_q=float(R['coarse']['q'][max(0, max(range(len(R['coarse']['q'])),
                                                     key=lambda i: R['coarse']['R6_re'][i]))]),
            d6_re_R=float(max(R['coarse']['R6_re'])),
        ),
    )
    out['quasiperiodicity_summary'] = {}
    for k, v in Q['lines'].items():
        ap = v.get('autocorr_peaks', {})
        out['quasiperiodicity_summary'][k] = dict(
            n_valid=v.get('n_valid'),
            n_finite=v.get('n_finite'),
            logH_median=v.get('logH_median'),
            ac_logH_at_1_16=ap.get('lag_0.0625', {}).get('ac_logH'),
            ac_logH_at_1_24=ap.get('lag_0.0417', {}).get('ac_logH'),
            ac_logH_at_1_8=ap.get('lag_0.1250', {}).get('ac_logH'),
            ac_logH_at_1_4=ap.get('lag_0.2500', {}).get('ac_logH'),
            ratio_n=v.get('ratio_n'),
            ratio_abs_median=v.get('ratio_abs_median'),
            ratio_abs_std=v.get('ratio_abs_std'),
            ratio_arg_R=v.get('ratio_arg_R'),
            ratio_arg_mean=v.get('ratio_arg_mean'),
            fourier_power_q16=v.get('fourier_power_q16'),
            fourier_power_q8=v.get('fourier_power_q8'),
            fourier_power_q32=v.get('fourier_power_q32'),
            fourier_total=v.get('fourier_total'),
            top_freq=v.get('fourier_top_freqs', [{}])[0].get('freq'),
        )
    out['zeros_vs_m'] = Z
    out['local_nodes'] = L
    with open(os.path.join(OUT_DIR, 'f8_results.json'), 'w', encoding='utf-8') as f:
        json.dump(out, f, indent=2, ensure_ascii=False)
    log('wrote f8_results.json  elapsed=%.1fs' % (time.time() - T0))


if __name__ == '__main__':
    main()