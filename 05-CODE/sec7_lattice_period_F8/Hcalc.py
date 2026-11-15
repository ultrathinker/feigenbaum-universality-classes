#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""F8 — vectorized float64 numpy evaluator for H_d(u) at level m=8 (and any m).
   Uses Horner-style evaluation for the polynomial coeffs in t = z^d.
   Reuses the rule from agentB_deg6/sextic_fixedpoint.json (max|t|<1.15, |x_{2^m}|>1e-12)
   and the analogous rule from quartic_pilot.
   Validates by reproducing 5-10 catalog zeros to |H|<1e-8 and H(-1) ~ 0.
"""
import os
import sys
import json
import time
import math
import numpy as np

ROOT_IN = r'D:/Исследования/И-3 — Классы универсальности/Исследование-3-Принцип'
OUT_DIR = r'.'

T0 = time.time()
LOG_PATH = os.path.join(OUT_DIR, 'progress.log')


def log(msg):
    line = '[%7.1fs] %s' % (time.time() - T0, msg)
    print(line, flush=True)
    with open(LOG_PATH, 'a', encoding='utf-8') as f:
        f.write(line + '\n')


def horner_d(z, coeffs, d):
    """Evaluate sum_{k} coeffs[k] * z^(d*k) via Horner in t = z^d.
       coeffs: list/array length N (highest power last).
       z: complex ndarray (or scalar).
       Returns: same shape as z, dtype complex128.
    """
    z = np.asarray(z, dtype=np.complex128)
    t = z ** d
    r = np.full(z.shape, float(coeffs[-1]), dtype=np.float64).astype(np.complex128)
    for c in reversed(coeffs[:-1]):
        r = r * t + complex(c)
    return r


def H_eval(u_arr, g_coeffs, h_coeffs, AU0, delta, m=8, TMAX=1.15, DEN_MIN=1e-12):
    """Vectorized H_d at level m.  u_arr: complex ndarray (or broadcastable to complex).
       Returns: (H_abs, valid_mask, maxt, errtrunc) — same shape as u_arr.
    """
    u_arr = np.asarray(u_arr, dtype=np.complex128)
    s_coef = AU0 / (delta ** m)
    s = s_coef * u_arr  # broadcasts
    n1 = 2 ** m
    n2 = 2 ** (m + 1)

    x = np.zeros(u_arr.shape, dtype=np.complex128)
    maxt = np.zeros(u_arr.shape, dtype=np.float64)
    x2m = np.zeros(u_arr.shape, dtype=np.complex128)
    for k in range(n2):
        t_d = x ** 6 if False else (x ** 4 if False else x ** 4)  # placeholder
        # We need d from closure; we can't pass it cleanly through generic fn.
        raise NotImplementedError('Use H_eval_d(..., d=...)')


def make_H_eval(d, g_coeffs, h_coeffs, AU0, delta, TMAX=1.15):
    """Factory producing H_eval_d(u_arr, m=8) bound to degree d, coeffs etc.
       Returns a function (u_arr, m=8) -> (H, valid, maxt, x2m).
    """
    g_np = np.array(g_coeffs, dtype=np.float64)
    h_np = np.array(h_coeffs, dtype=np.float64)
    AU0_f = float(AU0)
    delta_f = float(delta)

    def H_eval_d(u_arr, m=8, DEN_MIN=1e-12):
        u_arr = np.asarray(u_arr, dtype=np.complex128)
        s_coef = AU0_f / (delta_f ** m)
        s = s_coef * u_arr
        n1 = 2 ** m
        n2 = 2 ** (m + 1)
        x = np.zeros(u_arr.shape, dtype=np.complex128)
        maxt = np.zeros(u_arr.shape, dtype=np.float64)
        x2m = np.zeros(u_arr.shape, dtype=np.complex128)
        for k in range(n2):
            t_d = x ** d
            gv = horner_d(x, g_np, d)
            hv = horner_d(x, h_np, d)
            x = gv + s * hv
            np.maximum(maxt, np.abs(t_d), out=maxt)
            if k + 1 == n1:
                x2m = x.copy()
        Habs = np.abs(x) / np.abs(x2m)
        valid = (np.isfinite(Habs) & (np.abs(x2m) > DEN_MIN)
                 & (maxt < TMAX) & (np.abs(x) < 1e100))
        Habs = np.where(valid, Habs, np.nan)
        return Habs, valid, maxt, x2m
    return H_eval_d


def load_inputs():
    fp6 = json.load(open(os.path.join(ROOT_IN, 'agentB_deg6', 'sextic_fixedpoint.json'),
                         encoding='utf-8'))
    hz6 = json.load(open(os.path.join(ROOT_IN, 'agentB_deg6', 'H6_zeros.json'),
                         encoding='utf-8'))
    fp4 = json.load(open(os.path.join(ROOT_IN, 'quartic_pilot', 'quartic_fixedpoint.json'),
                         encoding='utf-8'))
    hz4 = json.load(open(os.path.join(ROOT_IN, 'quartic_pilot', 'H4_zeros.json'),
                         encoding='utf-8'))
    return fp6, hz6, fp4, hz4


def main():
    log('--- Hcalc: load inputs ---')
    fp6, hz6, fp4, hz4 = load_inputs()
    log('delta6=%.20f  alpha6=%.20f' % (float(fp6['delta6']), float(fp6['alpha6'])))
    log('AU0_6=%.20f  ncat6=%d' % (float(hz6['AU0_6']), len(hz6['zeros_level8'])))
    log('delta4=%.20f  alpha4=%.20f' % (float(fp4['delta4']), float(fp4['alpha4'])))
    log('AU0_4=%.20f  ncat4=%d' % (float(hz4['AU0_4']), len(hz4['zeros_level8'])))

    # Build evaluators
    H6 = make_H_eval(6, fp6['g_coeffs'], fp6['h_coeffs'], hz6['AU0_6'], fp6['delta6'])
    H4 = make_H_eval(4, fp4['g_coeffs'], fp4['h_coeffs'], hz4['AU0_4'], fp4['delta4'])

    # Validate: H(-1) ~ 0 (because x_{2^{m+1}} = 0 for u = -1)
    u_test = np.array([-1.0 + 0j], dtype=np.complex128)
    h6, _, _, x2m = H6(u_test)
    h4, _, _, x2m4 = H4(u_test)
    log('check H6(-1): |H|=%.3e  |x_2m|=|x_{512}|=%.3e (should be 0)'
        % (h6[0], np.abs(x2m[0])))
    log('check H4(-1): |H|=%.3e  |x_2m|=|x_{512}|=%.3e (should be 0)'
        % (h4[0], np.abs(x2m4[0])))

    # Validate on 5-10 catalog zeros
    log('--- validate on catalog zeros (d=6) ---')
    z6 = hz6['zeros_level8']
    picks6 = [z for z in z6 if abs(z['im']) < 5][:8]
    u = np.array([complex(z['re'], z['im']) for z in picks6], dtype=np.complex128)
    H, valid, maxt, x2m = H6(u)
    for i, z in enumerate(picks6):
        log('  cat d=6 re=%.6f im=%.6f  absH8(f64)=%.3e  valid=%s  maxt=%.4f'
            % (z['re'], z['im'], H[i], bool(valid[i]), maxt[i]))

    log('--- validate on catalog zeros (d=4) ---')
    z4 = hz4['zeros_level8']
    picks4 = z4[:8]
    u = np.array([complex(z['re'], z['im']) for z in picks4], dtype=np.complex128)
    H, valid, maxt, x2m = H4(u)
    for i, z in enumerate(picks4):
        log('  cat d=4 re=%.6f im=%.6f  absH8(f64)=%.3e  valid=%s  maxt=%.4f'
            % (z['re'], z['im'], H[i], bool(valid[i]), maxt[i]))

    # Save meta for downstream
    meta = dict(
        delta6=float(fp6['delta6']), alpha6=float(fp6['alpha6']), AU0_6=float(hz6['AU0_6']),
        delta4=float(fp4['delta4']), alpha4=float(fp4['alpha4']), AU0_4=float(hz4['AU0_4']),
        n_coeffs_g6=len(fp6['g_coeffs']), n_coeffs_h6=len(fp6['h_coeffs']),
        n_coeffs_g4=len(fp4['g_coeffs']), n_coeffs_h4=len(fp4['h_coeffs']),
    )
    with open(os.path.join(OUT_DIR, 'hcalc_meta.json'), 'w', encoding='utf-8') as f:
        json.dump(meta, f, indent=2, ensure_ascii=False)
    log('wrote hcalc_meta.json  elapsed=%.1fs' % (time.time() - T0))


if __name__ == '__main__':
    os.environ.setdefault('OMP_NUM_THREADS', '2')
    os.environ.setdefault('MKL_NUM_THREADS', '2')
    os.environ.setdefault('OPENBLAS_NUM_THREADS', '2')
    main()