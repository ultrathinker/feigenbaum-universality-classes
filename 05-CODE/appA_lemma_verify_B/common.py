# -*- coding: utf-8 -*-
"""Shared loaders/evaluators for the S1-VERIFY task. Own code (not orchestrator's)."""
import os
import time
os.environ.setdefault('OMP_NUM_THREADS', '2')
os.environ.setdefault('MKL_NUM_THREADS', '2')
os.environ.setdefault('OPENBLAS_NUM_THREADS', '2')

import json
from mpmath import mp, mpf, mpc, nint

mp.dps = 50

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
PROG = os.path.join(HERE, 'progress.log')


def plog(msg):
    line = time.strftime('%Y-%m-%d %H:%M:%S ') + msg
    print(line, flush=True)
    with open(PROG, 'a', encoding='utf-8') as f:
        f.write(line + '\n')
        f.flush()


def G_and_Gp(a, t):
    """Simultaneous Horner evaluation of G(t)=sum a_j t^j and G'(t)."""
    n = len(a) - 1
    val = a[n]
    der = mpf(0)
    for j in range(n - 1, -1, -1):
        der = der * t + val
        val = val * t + a[j]
    return val, der


def load_spectra_grid():
    """d(string) -> dict(a=[mpf,...], alpha=mpf) using a48/alpha48."""
    path = os.path.join(ROOT, 'wave1_spectra_d', 'spectra_grid.json')
    with open(path, encoding='utf-8') as f:
        data = json.load(f)
    out = {}
    for p in data['points']:
        d = p['d']
        a = [mpf(x) for x in p['_a48']]
        alpha = mpf(p['alpha48'])
        out[d] = {'a': a, 'alpha': alpha, 'raw': p}
    return out


def load_sextic_n72():
    path = os.path.join(ROOT, 'agentB_deg6', 'sextic_fixedpoint.json')
    with open(path, encoding='utf-8') as f:
        data = json.load(f)
    a = [mpf(x) for x in data['g_coeffs']]
    alpha = mpf(data['alpha6'])
    return {'a': a, 'alpha': alpha, 'raw': data}


def load_quartic_n48():
    path = os.path.join(ROOT, 'quartic_pilot', 'quartic_fixedpoint.json')
    with open(path, encoding='utf-8') as f:
        data = json.load(f)
    a = [mpf(x) for x in data['g_coeffs']]
    alpha = mpf(data['alpha4'])
    return {'a': a, 'alpha': alpha, 'raw': data}


def load_spectrum_window():
    path = os.path.join(ROOT, 'wave2_vector_homotopy', 'spectrum_window.json')
    with open(path, encoding='utf-8') as f:
        data = json.load(f)
    return data


def load_spectrum_d6():
    path = os.path.join(ROOT, 'wave3_t22_d6', 'spectrum_d6.json')
    with open(path, encoding='utf-8') as f:
        data = json.load(f)
    return data


def Phi_k_real(a, d, k, t):
    """Phi_k(t) = G(t)^(dk+1) - d*t^(k+1)*G'(t) - K(t), real path.
    Uses an exact python-int exponent when dk is (numerically) an integer,
    else falls back to a complex principal-value power (returns mpc)."""
    Gt, Gpt = G_and_Gp(a, t)
    dk = d * k
    n_round = nint(dk)
    if abs(dk - n_round) < mpf('1e-20'):
        n_int = int(n_round)
        Gn = Gt ** (n_int + 1)
    else:
        Gn = mpc(Gt) ** (dk + 1)
    term2 = d * (t ** (k + 1)) * Gpt
    Kt = Gt - d * t * Gpt
    return Gn - term2 - Kt


def K_of(a, d, t):
    Gt, Gpt = G_and_Gp(a, t)
    return Gt - d * t * Gpt


def DT_of_Phi(a, d, alpha, Phi_func, t):
    """(DT Phi)(t) using the brief's closed formula, Phi_func: t -> value (mpf or mpc)."""
    ad = alpha ** d
    ta = t / ad
    Gt1, _ = G_and_Gp(a, mpf(1))
    Kt = K_of(a, d, t)
    U, _ = G_and_Gp(a, ta)
    aU = abs(U)
    sU = 1 if U >= 0 else -1
    W = aU ** d
    _, Gpw = G_and_Gp(a, W)
    Phi1 = Phi_func(mpf(1))
    PhiW = Phi_func(W)
    Phita = Phi_func(ta)
    return -(Phi1 / Gt1) * Kt - alpha * PhiW - d * alpha * sU * (aU ** (d - 1)) * Gpw * Phita


def sample_points(m=7):
    """m Chebyshev-like interior points of (0,1), avoiding the endpoints/degenerate t=0."""
    from mpmath import cos, pi
    return [(1 - cos(pi * (i + mpf('0.5')) / m)) / 2 for i in range(m)]
