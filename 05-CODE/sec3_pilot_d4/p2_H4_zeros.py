#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
p2_H4_zeros.py -- Phase 2 of the quartic float pilot: build H4, scan its zeros.
[oid:quartK5]  NOT RIGOROUS -- float exploration (numpy complex128 + mpmath polish),
no interval arithmetic.

H4 by the SAME recipe as the quadratic-article H (article 2.2), all parameters
replaced by their quartic analogues:

    universal perturbed family:  f_s(x) = g4(x) + s * h4(x),
        g4(x) = sum_j a_j x^{4j}  (a_0=1),  h4(x) = sum_j c_j x^{4j}  (h4(1)=1)
    orbit:    x_0 = 0,  x_{n+1} = g4(x_n) + s*h4(x_n)
    coupling: s_m(u) = u * AU0_4 / delta4^m
    H4_m(u)  = x_{2^{m+1}}(u) / x_{2^m}(u),   H4 = lim H4_m

AU0_4 calibration (analogue of the article's H(-1)=0 gauge): superstable events
of the universal family, s*_m (real principal cascade, x_{2^{m+1}}=0, exact
period 2^{m+1}), obey s*_m -> Ls * delta4^{-m}; AU0_4 := -Ls puts the limit
cascade at u = -1.

Zero scan: numpy complex128 grid over the upper half-plane + real axis; local
minima of |H4_m| -> Newton refinement (holomorphic derivative recursion) ->
mpmath polish at level m=8.  Domain validity (analogue of the article's |w|<1
domain condition): per-orbit (a) max |x_n|^4 < 2.2 (below the conservative
series-radius estimate) and (b) truncation error N=48 vs N=40 below 1e-11.

Output: H4_zeros.json
"""
import json
import time

import numpy as np
from mpmath import mp, mpf, mpc, fabs, nstr, re as mpre, im as mpim

mp.dps = 50

TAG = "[oid:quartK5] NOT RIGOROUS float pilot, phase 2"
T0 = time.time()


def log(msg):
    print("[%7.1fs] %s" % (time.time() - T0, msg), flush=True)


FP = json.load(open('quartic_fixedpoint.json', encoding='utf-8'))
A_MP = [mpf(s) for s in FP['g_coeffs']]           # a_0..a_48, a_0 = 1
C_MP = [mpf(s) for s in FP['h_coeffs']]           # c_0..c_48, c_0 = 0, sum c_j = 1
ALPHA4 = mpf(FP['alpha4'])
DELTA4 = mpf(FP['delta4'])
N_MOD = len(A_MP) - 1

A_PY = [float(x) for x in A_MP]
C_PY = [float(x) for x in C_MP]
A_np = np.array(A_PY)
C_np = np.array(C_PY)
A_np40 = A_np.copy()
A_np40[41:] = 0.0
C_np40 = C_np.copy()
C_np40[41:] = 0.0

# --------------------------------------------------------------------------- #
#  mpmath evaluators
# --------------------------------------------------------------------------- #
def g4(z):
    t = z ** 4
    r = A_MP[-1]
    for c in reversed(A_MP[:-1]):
        r = r * t + c
    return r


def h4(z):
    t = z ** 4
    r = C_MP[-1]
    for c in reversed(C_MP[:-1]):
        r = r * t + c
    return r


def g4p(z):
    t = z ** 4
    s = mpf(0)
    for j in range(len(A_MP) - 1, 0, -1):
        s = s * t + j * A_MP[j]
    return 4 * z ** 3 * s


def h4p(z):
    t = z ** 4
    s = mpf(0)
    for j in range(len(C_MP) - 1, 0, -1):
        s = s * t + j * C_MP[j]
    return 4 * z ** 3 * s


def orbit_H_mpmath(u, m, AU0, want_grad=False):
    """x_0=0, x_{n+1} = g4(x_n) + s*h4(x_n), s = u*AU0/delta4^m.
    Returns (x_{2^m}, x_{2^{m+1}}) (+ d/du of both if want_grad)."""
    s = u * AU0 / DELTA4 ** m
    n1, n2 = 2 ** m, 2 ** (m + 1)
    x = mpc(0)
    dx = mpc(0)
    dsdu = AU0 / DELTA4 ** m
    x2m = dx2m = None
    for k in range(n2):
        if want_grad:
            dx = (g4p(x) + s * h4p(x)) * dx + h4(x) * dsdu
        x = g4(x) + s * h4(x)
        if k + 1 == n1:
            x2m, dx2m = x, dx
    if want_grad:
        return x2m, x, dx2m, dx
    return x2m, x


# --------------------------------------------------------------------------- #
#  numpy scan
# --------------------------------------------------------------------------- #
def poly4(z, co):
    t = z * z
    t = t * t
    r = np.full_like(z, co[-1])
    for c in reversed(co[:-1]):
        r = r * t + c
    return r


def scan_level(m, AU0, grid):
    """Returns Habs, valid mask."""
    s = complex(float(AU0) / float(DELTA4) ** m) * grid
    n1, n2 = 2 ** m, 2 ** (m + 1)
    x = np.zeros(grid.shape, dtype=np.complex128)
    maxt = np.zeros(grid.shape, dtype=np.float64)
    errmax = np.zeros(grid.shape, dtype=np.float64)
    x2m = None
    for k in range(n2):
        gfull = poly4(x, A_np)
        gtrunc = poly4(x, A_np40)
        np.maximum(errmax, np.abs(gfull - gtrunc), out=errmax)
        x = gfull + s * poly4(x, C_np)
        t4 = x * x
        np.maximum(maxt, np.abs(t4 * t4), out=maxt)
        if k + 1 == n1:
            x2m = x.copy()
    num = np.abs(x)
    den = np.abs(x2m)
    with np.errstate(divide='ignore', invalid='ignore'):
        Habs = num / den
    valid = np.isfinite(Habs) & (den > 1e-12) & (maxt < 2.2) & (errmax < 1e-11)
    return Habs, valid, maxt, errmax


def eval_gh_py(x):
    """(g, h, g', h') at complex x, plain Python."""
    t = x ** 4
    g = A_PY[-1]
    for c in reversed(A_PY[:-1]):
        g = g * t + c
    h = C_PY[-1]
    for c in reversed(C_PY[:-1]):
        h = h * t + c
    dg = 0j
    for j in range(len(A_PY) - 1, 0, -1):
        dg = dg * t + j * A_PY[j]
    dh = 0j
    for j in range(len(C_PY) - 1, 0, -1):
        dh = dh * t + j * C_PY[j]
    x3 = x ** 3
    return g, h, 4.0 * x3 * dg, 4.0 * x3 * dh


def newton_refine(u0, m, AU0, iters=60, rmax=22.0, tol=1e-12):
    """Newton on H_m(u)=0, plain Python complex with overflow-guarded damping.
    Returns (u, |H|) or (None, inf)."""
    s_coef = float(AU0) / float(DELTA4) ** m
    n1, n2 = 2 ** m, 2 ** (m + 1)

    def H_at(u):
        s = u * s_coef
        x = 0j
        x2m = None
        for k in range(n2):
            g, h, _gp, _hp = eval_gh_py(x)
            x = g + s * h
            if k + 1 == n1:
                x2m = x
        if x2m == 0 or not (abs(x2m) > 1e-290):
            return None
        return x / x2m

    u = u0
    for _ in range(iters):
        s = u * s_coef
        x = 0j
        dx = 0j
        x2m = None
        dx2m = None
        for k in range(n2):
            g, h, gp, hp = eval_gh_py(x)
            dx = (gp + s * hp) * dx + h * s_coef
            x = g + s * h
            if k + 1 == n1:
                x2m, dx2m = x, dx
        if not (abs(x2m) > 1e-290) or not np.isfinite(abs(x) + abs(dx)):
            return None, float('inf')
        H = x / x2m
        dH = (dx * x2m - x * dx2m) / (x2m * x2m)
        if abs(H) < tol:
            break
        if dH == 0 or not np.isfinite(H.real + H.imag + dH.real + dH.imag):
            return None, float('inf')
        step = H / dH
        # damped, overflow-guarded step
        moved = False
        for _h in range(20):
            u_try = u - step
            if abs(u_try) > rmax:
                step = step / 2
                continue
            Ht = H_at(u_try)
            if Ht is not None and abs(Ht) < abs(H):
                u = u_try
                moved = True
                break
            step = step / 2
        if not moved:
            break
        if abs(step) < 1e-15:
            break
    Hf = H_at(u)
    if Hf is None:
        return None, float('inf')
    return u, abs(Hf)


def domain_check(u, m, AU0, tmax=2.2, errtol=1e-11):
    """Validity of the refined zero: orbit stays inside |x|^4 < tmax and the
    series truncation error (N=48 vs N=40) stays below errtol."""
    s = u * float(AU0) / float(DELTA4) ** m
    n1, n2 = 2 ** m, 2 ** (m + 1)
    x = 0j
    maxt4 = 0.0
    maxerr = 0.0
    den = None
    for k in range(n2):
        t = x ** 4
        g48 = A_PY[-1]
        for c in reversed(A_PY[:-1]):
            g48 = g48 * t + c
        g40 = A_PY[40]
        for c in reversed(A_PY[:40]):
            g40 = g40 * t + c
        err = abs(g48 - g40)
        if err > maxerr:
            maxerr = err
        if abs(t) > maxt4:
            maxt4 = abs(t)
        h = C_PY[-1]
        for c in reversed(C_PY[:-1]):
            h = h * t + c
        x = g48 + s * h
        if k + 1 == n1:
            den = x
    return maxt4, maxerr, abs(den)


def main():
    log(TAG)

    # ---- AU0_4 calibration -------------------------------------------------------
    log('--- AU0_4 calibration: superstable cascade of the universal family ---')
    smax = 12
    s_stars = {}
    s_prev = None
    for m in range(2, smax + 1):
        n = 2 ** (m + 1)

        def fres(s):
            x = mpf(0)
            ds = mpf(0)
            for _ in range(n):
                ds = (g4p(x) + s * h4p(x)) * ds + h4(x)
                x = g4(x) + s * h4(x)
            return x, ds

        s0 = s_prev / DELTA4 if s_prev is not None else mpf('-1e-3')
        s = s0
        ok_conv = False
        for _it in range(80):
            x, ds = fres(s)
            if fabs(x) < mpf('1e-40'):
                ok_conv = True
                break
            if ds == 0:
                break
            ss = x / ds
            improved = False
            for _h in range(30):
                x2, _ = fres(s - ss)
                if fabs(x2) < fabs(x):
                    s = s - ss
                    improved = True
                    break
                ss = ss / 2
            if not improved:
                break
        if not ok_conv:
            raise RuntimeError('AU0 calibration Newton failed at m=%d' % m)
        for j in range(1, m + 1):
            x = mpf(0)
            for _ in range(2 ** j):
                x = g4(x) + s * h4(x)
            if fabs(x) < mpf('1e-30'):
                raise RuntimeError('inherited root at m=%d j=%d' % (m, j))
        s_stars[m] = s
        s_prev = s
        log('m=%2d  s*_m = %s   Ls^(m) = %s'
            % (m, nstr(s, 20), nstr(s * DELTA4 ** m, 20)))
    rL = 1 / DELTA4
    Ls_seq = {m: s_stars[m] * DELTA4 ** m for m in s_stars}
    Ls = (Ls_seq[smax] - rL * Ls_seq[smax - 1]) / (1 - rL)
    AU0 = -Ls
    log('Ls = %s   AU0_4 = -Ls = %s' % (nstr(Ls, 30), nstr(AU0, 30)))
    for m in (6, 8, 10):
        x2m, x2m1 = orbit_H_mpmath(mpf(-1), m, AU0)
        log('gauge check H_%d(-1) = %s' % (m, nstr(x2m1 / x2m, 12)))

    # ---- float zero scan ----------------------------------------------------------
    R_MAX = 20.0
    STEP = 0.0625
    xs = np.arange(-R_MAX, R_MAX + STEP / 2, STEP)
    ys = np.arange(0.0, R_MAX + STEP / 2, STEP)
    UU, VV = np.meshgrid(xs, ys, indexing='ij')
    grid = UU + 1j * VV

    results = {}
    for m in (7, 8):
        Habs, valid, maxt, errmax = scan_level(m, AU0, grid)
        n_valid = int(valid.sum())
        log('scan m=%d: grid %dx%d, valid %d (%.1f%%), min|H| on valid = %.3e'
            % (m, grid.shape[0], grid.shape[1], n_valid,
               100.0 * n_valid / grid.size,
               float(Habs[valid].min()) if n_valid else float('nan')))
        Hm = np.where(valid, Habs, np.inf)
        # pad with inf so the real-axis row and edges can host minima
        P = np.full((Hm.shape[0] + 2, Hm.shape[1] + 2), np.inf)
        P[1:-1, 1:-1] = Hm
        core = P[1:-1, 1:-1]
        neighbours = np.stack([
            P[:-2, :-2], P[:-2, 1:-1], P[:-2, 2:],
            P[1:-1, :-2], P[1:-1, 2:],
            P[2:, :-2], P[2:, 1:-1], P[2:, 2:]])
        ismin = (core <= neighbours).all(axis=0) & (core < 0.2)
        ii, jj = np.nonzero(ismin)
        cand = [complex(grid[i, j]) for i, j in zip(ii, jj)]
        # also seed from the previous level's zeros (same u, tiny drift)
        if m == 8:
            cand = cand + [z['u'] for z in results[7]]
        log('  candidates: %d' % len(cand))
        zeros = []
        for u0 in cand:
            u, res = newton_refine(u0, m, AU0, tol=1e-9)
            if u is None or res > 1e-5:
                continue
            mt4, merr, den = domain_check(u, m, AU0)
            zeros.append({'u': u, 'absH': res, 'maxt4': mt4, 'maxerr': merr,
                          'den': den,
                          'domain_ok': bool(mt4 < 2.2 and merr < 1e-11 and den > 1e-12)})
        dedup = []
        for z in zeros:
            if all(abs(z['u'] - v['u']) > 1e-4 for v in dedup):
                dedup.append(z)
        log('  refined zeros (|H|<1e-5): %d (domain-ok: %d)'
            % (len(dedup), sum(1 for z in dedup if z['domain_ok'])))
        results[m] = dedup

    # ---- mpmath polish at m=8 of all domain-ok candidates --------------------------
    cat7 = [z for z in results[7] if z['domain_ok']]
    cat8 = [z for z in results[8] if z['domain_ok']]

    def mp_polish(u_start, m, iters=14, tol_step=mpf('1e-24')):
        u_mp = mpc(mpf(u_start.real), mpf(u_start.imag))
        for _ in range(iters):
            x2m, x2m1, dx2m, dx2m1 = orbit_H_mpmath(u_mp, m, AU0, want_grad=True)
            H = x2m1 / x2m
            dH = (dx2m1 * x2m - x2m1 * dx2m) / x2m ** 2
            step = H / dH
            u_mp = u_mp - step
            if abs(step) < tol_step:
                break
        x2m, x2m1 = orbit_H_mpmath(u_mp, m, AU0)
        return u_mp, abs(x2m1 / x2m)

    cat8_mp = []
    for z in cat8:
        try:
            u_mp, r = mp_polish(z['u'], 8)
        except Exception as exc:  # noqa: BLE001
            log('  polish failed near %.4f%+.4fi: %s' % (z['u'].real, z['u'].imag, exc))
            continue
        if float(r) < 1e-20:
            cat8_mp.append({'u8_mp': u_mp, 'absH8': float(r),
                            'u8_seed': z['u'], 'maxt4': z['maxt4'], 'maxerr': z['maxerr']})
    log('polished level-8 zeros (|H8|<1e-20): %d' % len(cat8_mp))

    # ---- pair with level-7 zeros (survival filter) ----------------------------------
    paired = []
    for z8 in cat8_mp:
        u8c = complex(mpre(z8['u8_mp']), mpim(z8['u8_mp']))
        best = min(cat7, key=lambda z: abs(z['u'] - u8c))
        if abs(best['u'] - u8c) < 5e-3:
            paired.append({**z8, 'u7': best['u'], 'u8': u8c})
        else:
            log('  level-8 zero at %+.6f%+.6fi has NO level-7 partner (nearest %.4f)'
                % (u8c.real, u8c.imag, abs(best['u'] - u8c)))
    log('paired zeros (survive m=7 and m=8): %d' % len(paired))

    # ---- level-9 refinement, drift-ratio check, limit extrapolation ------------------
    for z in paired:
        try:
            u9_mp, r9 = mp_polish(z['u8_mp'], 9)
        except Exception:  # noqa: BLE001
            u9_mp, r9 = None, mpf('inf')
        z['u9_mp'] = u9_mp
        z['absH9'] = float(r9)
        if float(r9) < 1e-20:
            u9c = complex(mpre(u9_mp), mpim(u9_mp))
            z['u9'] = u9c
            z['drift78'] = abs(z['u8'] - z['u7'])
            z['drift89'] = abs(u9c - z['u8'])
            z['drift_ratio'] = z['drift78'] / z['drift89'] if z['drift89'] > 0 else None
            z['uinf'] = u9c + (u9c - z['u8']) / float(DELTA4 - 1)
        else:
            z['drift78'] = abs(z['u8'] - z['u7'])
            z['drift89'] = None
            z['drift_ratio'] = None
            z['uinf'] = z['u8'] + (z['u8'] - z['u7']) / float(DELTA4 - 1)

    paired.sort(key=lambda z: (abs(z['uinf']), z['uinf'].real, z['uinf'].imag))
    log('=== CATALOG (limit-extrapolated zeros of H4) ===')
    for z in paired:
        dr = ('%.3f' % z['drift_ratio']) if z['drift_ratio'] else 'n/a'
        m9 = ('%.1e' % z['absH9']) if z['absH9'] is not None else 'fail'
        log('   uinf = %+.8f %+.8f i  |u|=%6.4f  |H8|=%.1e |H9|=%s  d78=%.2e ratio=%s'
            % (z['uinf'].real, z['uinf'].imag, abs(z['uinf']), z['absH8'], m9,
               z['drift78'], dr))

    out = {
        'meta': {'tag': TAG, 'date': time.strftime('%Y-%m-%d %H:%M:%S'),
                 'recipe': 'article 2.2 analogue: x_{n+1}=g4(x_n)+s*h4(x_n), s=u*AU0_4/delta4^m, H4_m=x_{2^{m+1}}/x_{2^m}',
                 'scan': {'R_max': R_MAX, 'step': STEP, 'levels': [7, 8],
                          'domain_rule': 'max|x_n|^4 < 2.2 and |g_N48-g_N40| < 1e-11 along orbit, |x_{2^m}|>1e-12'},
                 'elapsed_sec': time.time() - T0},
        'AU0_4': nstr(AU0, 40),
        'Ls': nstr(Ls, 40),
        's_star': {str(m): nstr(s_stars[m], 30) for m in s_stars},
        'zeros_level8': [{'re': c['u8'].real, 'im': c['u8'].imag,
                          'absH': c['absH8'], 'str': nstr(c['u8_mp'], 25),
                          'maxt4': c['maxt4'], 'maxerr': c['maxerr']} for c in paired],
        'zeros_level7': [{'re': c['u7'].real, 'im': c['u7'].imag} for c in paired],
        'zeros_limit_extrapolated': [
            {'re': c['uinf'].real, 'im': c['uinf'].imag,
             'u7_re': c['u7'].real, 'u7_im': c['u7'].imag,
             'u8_re': c['u8'].real, 'u8_im': c['u8'].imag,
             'u9_re': (c['u9'].real if c.get('u9') else None),
             'u9_im': (c['u9'].imag if c.get('u9') else None),
             'absH9': c.get('absH9'),
             'drift78': c['drift78'], 'drift89': c.get('drift89'),
             'drift_ratio': c.get('drift_ratio')}
            for c in paired],
        'counts': {'level8_polished': len(cat8_mp), 'paired': len(paired),
                   'domain_ok_7': len(cat7), 'domain_ok_8': len(cat8),
                   'all_refined_7': len(results[7]), 'all_refined_8': len(results[8])},
    }
    with open('H4_zeros.json', 'w', encoding='utf-8') as f:
        json.dump(out, f, indent=2, ensure_ascii=False)
    log('wrote H4_zeros.json')


if __name__ == '__main__':
    main()
