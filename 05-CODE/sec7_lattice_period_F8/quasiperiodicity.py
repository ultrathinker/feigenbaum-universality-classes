#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""F8 Task 3: Quasi-periodicity of H(u) along lines.

For each line (horizontal Im u = c, or vertical Re u = c), compute:
  - log|H(u)| on grid step 1/256, length >=20 (=> 5120+ points)
  - arg H(u) (unwrapped via numpy)
  - autocorrelations of log|H| and d(arg)/du vs shift Delta in (0, 1] step 1/256
  - ratio H(u+1/16)/H(u): distribution of |ratio| and arg(ratio) (Robust: drop near-zeros)
  - Fourier spectrum of log|H|: peak at freq 16?

Output: quasiperiodicity_results.json + figure.
"""
import os
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
    z = np.asarray(z, dtype=np.complex128)
    t = z ** d
    r = np.full(z.shape, float(coeffs[-1]), dtype=np.float64).astype(np.complex128)
    for c in reversed(coeffs[:-1]):
        r = r * t + complex(c)
    return r


def H_line_factory(d, g_coeffs, h_coeffs, AU0, delta, m=8, TMAX=1.15, DEN_MIN=1e-12):
    g_np = np.array(g_coeffs, dtype=np.float64)
    h_np = np.array(h_coeffs, dtype=np.float64)
    AU0_f = float(AU0); delta_f = float(delta)
    def H_at(u_arr):
        u_arr = np.asarray(u_arr, dtype=np.complex128)
        s_coef = AU0_f / (delta_f ** m)
        s = s_coef * u_arr
        n1 = 2 ** m; n2 = 2 ** (m + 1)
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
        return Habs, valid, maxt, x2m, x
    return H_at


def autocorr(x, maxlag):
    """Normalized autocorrelation xcorr(x, x) for lags 0..maxlag-1 (1-based for shift).
       Returns array length maxlag with xcorr[0] = 1 (lag 0)."""
    x = np.asarray(x, dtype=np.float64)
    x = x - x.mean()
    n = len(x)
    if n < 2 * maxlag:
        # use np.correlate
        f = np.fft.fft(x, n=2 ** int(math.ceil(math.log2(2 * n))))
        p = np.abs(f) ** 2
        ac = np.fft.ifft(p).real[:n]
        ac /= ac[0]
    else:
        ac = np.array([np.dot(x[:n - l], x[l:]) / (n - l)
                        / (np.dot(x[:n - l], x[:n - l]) / (n - l))
                        for l in range(maxlag)])
    return ac


def autocorr_fft(x, maxlag):
    x = np.asarray(x, dtype=np.float64)
    x = x - x.mean()
    n = len(x)
    N = 2 ** int(math.ceil(math.log2(2 * n)))
    f = np.fft.fft(x, n=N)
    p = np.abs(f) ** 2
    ac = np.fft.ifft(p).real[:maxlag]
    return ac / ac[0]


def analyse_line(H_abs, valid, maxt, x_full, sample_param, line_label, q_targets=(8, 16, 32),
                 shifts=(1/24, 1/16, 1/10, 1/12, 1/8, 1/4)):
    """H_abs: |H(u)| ndarray; valid: bool; x_full: complex H (for arg); sample_param: real grid.
       Returns dict of stats.
    """
    res = dict(line=line_label, n_pts=int(H_abs.size),
               n_valid=int(valid.sum()),
               min_maxt=float(maxt[valid].min()) if valid.any() else None,
               max_maxt=float(maxt[valid].max()) if valid.any() else None)
    if valid.sum() < 200:
        return res

    # Mask & filter outliers (near zeros) — use robust cutoff
    logH = np.full(H_abs.shape, np.nan, dtype=np.float64)
    logH[valid] = np.log10(H_abs[valid])
    # Robust: clip extremes (zeros cause logH -> -inf; we want finite)
    logH_clean = logH.copy()
    logH_clean[~valid] = np.nan
    finite_mask = np.isfinite(logH_clean)
    if finite_mask.sum() < 200:
        return res
    res['n_finite'] = int(finite_mask.sum())

    # Interpolate over nans (where zeros caused failures): use scipy-ish linear
    # For simplicity, fill nan with nearest valid neighbor via pandas? Not allowed.
    # Use numpy: forward fill then backward fill
    arr = logH_clean.copy()
    idx = np.arange(len(arr))
    valid_idx = idx[finite_mask]
    if valid_idx.size < 2:
        return res
    arr_filled = np.interp(idx, valid_idx, arr[valid_idx])
    res['logH_median'] = float(np.median(arr_filled))
    res['logH_p10'] = float(np.percentile(arr_filled, 10))
    res['logH_p90'] = float(np.percentile(arr_filled, 90))

    # arg H (unwrapped) — only where valid
    argH = np.full(H_abs.shape, np.nan, dtype=np.float64)
    argH[valid] = np.unwrap(np.angle(x_full[valid]))
    arg_filled = np.interp(idx, valid_idx, argH[valid_idx])
    darg = np.gradient(arg_filled, sample_param)
    res['arg_range'] = float(arg_filled.max() - arg_filled.min())

    # Autocorrelation of logH
    logH_centered = arr_filled - arr_filled.mean()
    ac_logH = autocorr_fft(logH_centered, maxlag=257)
    # Autocorrelation of darg
    darg_centered = darg - darg.mean()
    ac_darg = autocorr_fft(darg_centered, maxlag=257)
    # Find peaks at target shifts: shift in (0,1] with step 1/256 -> lag = shift*256
    lags = np.arange(257) / 256.0
    peaks = {}
    for s in shifts:
        lag = int(round(s * 256))
        peaks['lag_%.4f' % s] = dict(lag=lag, shift=s,
                                     ac_logH=float(ac_logH[lag]),
                                     ac_darg=float(ac_darg[lag]))
    # Control lags: peaks at *no*-shift (lag=1)
    peaks['lag_1_256=%.4f' % (1/256)] = dict(lag=1, shift=1/256,
                                              ac_logH=float(ac_logH[1]),
                                              ac_darg=float(ac_darg[1]))
    # Get the maximum near lag in [3, 256] excluding near-zero
    peaks['lag_3_256=%.4f' % (3/256)] = dict(lag=3, shift=3/256,
                                              ac_logH=float(ac_logH[3]),
                                              ac_darg=float(ac_darg[3]))
    peaks['lag_64_256=%.4f' % (64/256)] = dict(lag=64, shift=64/256,
                                                ac_logH=float(ac_logH[64]),
                                                ac_darg=float(ac_darg[64]))
    res['autocorr_peaks'] = peaks
    # max autocorr outside lag=0
    ac_max_idx = int(np.argmax(ac_logH[1:])) + 1
    res['ac_logH_max_idx'] = int(ac_max_idx)
    res['ac_logH_max_shift'] = float(ac_max_idx / 256)
    res['ac_logH_max_val'] = float(ac_logH[ac_max_idx])

    # Ratio H(u+1/16) / H(u)
    shift16 = 16  # points at step 1/256 -> 16/256 = 1/16
    H1 = H_abs[: -shift16]
    H2 = H_abs[shift16:]
    valid_r = valid[: -shift16] & valid[shift16:]
    # Also avoid zeros in ratio: require both > 1e-4
    eps_abs = 1e-3
    nz_mask = valid_r & (H1 > eps_abs) & (H2 > eps_abs)
    if nz_mask.sum() > 50:
        ratio = H2[nz_mask] / H1[nz_mask]
        x1 = x_full[: -shift16][nz_mask]
        x2 = x_full[shift16:][nz_mask]
        ratio_c = x2 / x1
        res['ratio_n'] = int(nz_mask.sum())
        res['ratio_abs_median'] = float(np.median(ratio))
        res['ratio_abs_p25'] = float(np.percentile(ratio, 25))
        res['ratio_abs_p75'] = float(np.percentile(ratio, 75))
        res['ratio_abs_std'] = float(np.std(ratio))
        # circular stats for arg
        arg_ratio = np.angle(ratio_c)
        res['ratio_arg_mean'] = float(np.angle(np.mean(np.exp(1j * arg_ratio))))
        res['ratio_arg_R'] = float(np.abs(np.mean(np.exp(1j * arg_ratio))))
        res['ratio_arg_std'] = float(np.sqrt(-2 * np.log(max(res['ratio_arg_R'], 1e-12))))
    # Control: ratio at shift 1/10 (25.6 points — use nearest int)
    for ctrl_shift, ctrl_label in [(10, '1/10'), (12, '1/12'), (8, '1/8')]:
        pts = int(round(256 / ctrl_shift))
        if pts < 1 or pts * 2 >= len(H_abs):
            continue
        H1c = H_abs[: -pts]; H2c = H_abs[pts:]
        valid_c = valid[: -pts] & valid[pts:]
        nz_c = valid_c & (H1c > eps_abs) & (H2c > eps_abs)
        if nz_c.sum() > 50:
            ratio_c = H2c[nz_c] / H1c[nz_c]
            x1c = x_full[: -pts][nz_c]; x2c = x_full[pts:][nz_c]
            arg_c = np.angle(x2c / x1c)
            key = 'ctrl_ratio_shift_%s' % ctrl_label
            res[key] = dict(shift=ctrl_label, n=int(nz_c.sum()),
                            abs_median=float(np.median(ratio_c)),
                            abs_std=float(np.std(ratio_c)),
                            arg_R=float(np.abs(np.mean(np.exp(1j * arg_c)))),
                            arg_mean=float(np.angle(np.mean(np.exp(1j * arg_c)))),
                            arg_std=float(np.sqrt(-2 * np.log(max(np.abs(np.mean(np.exp(1j * arg_c))), 1e-12)))))

    # Fourier spectrum of log|H|
    # Use arr_filled (length ~5120).  Detrend: subtract mean.
    F = np.fft.fft(logH_centered)
    P = np.abs(F) ** 2
    freqs = np.fft.fftfreq(len(logH_centered), d=1/256)
    pos = freqs > 0
    pos_freqs = freqs[pos]
    pos_power = P[pos]
    res['fourier_top_freqs'] = []
    top_idx = np.argsort(pos_power)[-10:][::-1]
    for ti in top_idx:
        res['fourier_top_freqs'].append(dict(freq=float(pos_freqs[ti]),
                                              power=float(pos_power[ti])))
    # Power at q=16, 8, 32 (exact)
    for q in q_targets:
        idx_q = int(np.argmin(np.abs(pos_freqs - q)))
        res['fourier_power_q%d' % q] = float(pos_power[idx_q])
    # Total power
    res['fourier_total'] = float(pos_power.sum())
    return res


def main():
    log('--- Quasi-periodicity: load ---')
    fp6 = json.load(open(os.path.join(ROOT_IN, 'agentB_deg6', 'sextic_fixedpoint.json'), encoding='utf-8'))
    hz6 = json.load(open(os.path.join(ROOT_IN, 'agentB_deg6', 'H6_zeros.json'), encoding='utf-8'))
    fp4 = json.load(open(os.path.join(ROOT_IN, 'quartic_pilot', 'quartic_fixedpoint.json'), encoding='utf-8'))
    hz4 = json.load(open(os.path.join(ROOT_IN, 'quartic_pilot', 'H4_zeros.json'), encoding='utf-8'))
    H6 = H_line_factory(6, fp6['g_coeffs'], fp6['h_coeffs'], hz6['AU0_6'], fp6['delta6'])
    H4 = H_line_factory(4, fp4['g_coeffs'], fp4['h_coeffs'], hz4['AU0_4'], fp4['delta4'])

    step = 1/256  # 0.00390625

    # d=6 horizontal lines: Im u = c, Re u in [-10, 10]
    horiz_c6 = [2.0, 6.0, 10.0, 15.0, 20.0, 25.0, 30.0, 40.0]
    vert_c6 = [2.0, 10.0, 20.0, 30.0]

    # d=4 lines: |u|<=20
    horiz_c4 = [1.0, 3.0, 6.0, 10.0, 14.0, 18.0]
    vert_c4 = [2.0, 6.0, 14.0]

    out = dict(meta=dict(delta6=fp6['delta6'], AU0_6=hz6['AU0_6'],
                         delta4=fp4['delta4'], AU0_4=hz4['AU0_4'],
                         step=step, n_pts_horiz=int(round(20 / step))))
    results = {}

    # ---- d=6 horizontal ----
    for c in horiz_c6:
        re_grid = np.arange(-10.0, 10.0 + step / 2, step)
        u = re_grid + 1j * c
        log('[d=6 horiz c=%.1f] n_pts=%d' % (c, u.size))
        t0 = time.time()
        Habs, valid, maxt, x2m, x = H6(u)
        log('  computed %.1fs valid=%d' % (time.time() - t0, valid.sum()))
        key = 'd6_horiz_%.1f' % c
        results[key] = analyse_line(Habs, valid, maxt, x, re_grid, key)
        log('  n_finite=%d  ac_logH_max at lag=%d (shift=%.4f) val=%.3f'
            % (results[key].get('n_finite', 0),
               results[key].get('ac_logH_max_idx', 0),
               results[key].get('ac_logH_max_shift', 0),
               results[key].get('ac_logH_max_val', 0)))
    # ---- d=6 vertical ----
    for c in vert_c6:
        im_grid = np.arange(0.5, 20.5 + step / 2, step)
        u = c + 1j * im_grid
        log('[d=6 vert c=%.1f] n_pts=%d' % (c, u.size))
        t0 = time.time()
        Habs, valid, maxt, x2m, x = H6(u)
        log('  computed %.1fs valid=%d' % (time.time() - t0, valid.sum()))
        key = 'd6_vert_%.1f' % c
        results[key] = analyse_line(Habs, valid, maxt, x, im_grid, key)
        log('  n_finite=%d  ac_logH_max at lag=%d (shift=%.4f) val=%.3f'
            % (results[key].get('n_finite', 0),
               results[key].get('ac_logH_max_idx', 0),
               results[key].get('ac_logH_max_shift', 0),
               results[key].get('ac_logH_max_val', 0)))
    # ---- d=4 horizontal ----
    for c in horiz_c4:
        re_grid = np.arange(-10.0, 10.0 + step / 2, step)
        u = re_grid + 1j * c
        log('[d=4 horiz c=%.1f] n_pts=%d' % (c, u.size))
        t0 = time.time()
        Habs, valid, maxt, x2m, x = H4(u)
        log('  computed %.1fs valid=%d' % (time.time() - t0, valid.sum()))
        key = 'd4_horiz_%.1f' % c
        results[key] = analyse_line(Habs, valid, maxt, x, re_grid, key)
        log('  n_finite=%d  ac_logH_max at lag=%d (shift=%.4f) val=%.3f'
            % (results[key].get('n_finite', 0),
               results[key].get('ac_logH_max_idx', 0),
               results[key].get('ac_logH_max_shift', 0),
               results[key].get('ac_logH_max_val', 0)))
    # ---- d=4 vertical ----
    for c in vert_c4:
        im_grid = np.arange(0.5, 20.5 + step / 2, step)
        u = c + 1j * im_grid
        log('[d=4 vert c=%.1f] n_pts=%d' % (c, u.size))
        t0 = time.time()
        Habs, valid, maxt, x2m, x = H4(u)
        log('  computed %.1fs valid=%d' % (time.time() - t0, valid.sum()))
        key = 'd4_vert_%.1f' % c
        results[key] = analyse_line(Habs, valid, maxt, x, im_grid, key)
        log('  n_finite=%d  ac_logH_max at lag=%d (shift=%.4f) val=%.3f'
            % (results[key].get('n_finite', 0),
               results[key].get('ac_logH_max_idx', 0),
               results[key].get('ac_logH_max_shift', 0),
               results[key].get('ac_logH_max_val', 0)))

    out['lines'] = results
    with open(os.path.join(OUT_DIR, 'quasiperiodicity_results.json'), 'w', encoding='utf-8') as f:
        json.dump(out, f, indent=2, ensure_ascii=False)
    log('wrote quasiperiodicity_results.json  elapsed=%.1fs' % (time.time() - T0))


if __name__ == '__main__':
    os.environ.setdefault('OMP_NUM_THREADS', '2')
    os.environ.setdefault('MKL_NUM_THREADS', '2')
    os.environ.setdefault('OPENBLAS_NUM_THREADS', '2')
    main()