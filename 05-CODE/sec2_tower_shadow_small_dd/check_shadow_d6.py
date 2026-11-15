# -*- coding: utf-8 -*-
"""
check_shadow_d6.py -- independent re-check of the r9_main.py shadow test (own discretisation, two bases, basis-independent cosine).

Does NOT import/copy engine.py. Own implementation:
  - own Newton continuation for the fixed point G_d(t), t=|x|^d, using
    Chebyshev-Gauss-LOBATTO-type interior nodes t_i = (1-cos(i*pi/(N+1)))/2,
    i=1..N+1  (DIFFERENT formula from engine.py's Gauss-type half-integer
    nodes (1-cos(pi*(i+0.5)/npts))/2 used for colloc_points there).
  - TWO representations of the linearised doubling-operator eigenproblem on
    the same N-dim polynomial space span{t,...,t^N}:
      (A) monomial basis {t^j}, j=1..N            (same family R9 used,
          reproduced independently to sanity-check the qualitative trend)
      (B) shifted-Chebyshev basis e_j(t) = T_j(2t-1) - (-1)^j, j=1..N
          (genuinely different, well-conditioned basis for the SAME
          function space -- e_j(0)=0 automatically)
  - cosine computed THREE ways for every (d, N) eigenvector found:
      cos_raw_A  -- raw coefficient dot product in basis (A) [R9's metric]
      cos_raw_B  -- raw coefficient dot product in basis (B)
      cos_L2     -- values of Phi at a FIXED common physical grid of M
                    points in (0,1], dot product in that value space
                    (basis-independent functional/L2-style metric)
  - mode identified by CONTINUITY tracing across a short Delta-d ladder
    (nearest by |Delta lambda| AND high value-cosine to the previous step),
    not just "nearest |lambda| in a fixed window at the target d" (the trap
    named in the audit brief: window-nearest can jump onto a different
    branch at an avoided crossing).

Seeds: monomial coefficients a_j (N=48 truncation) and alpha at d=4.0, d=6.0
taken from ../sec2_spectra_grid_d/spectra_grid.json
(fields '_a48', 'alpha48') -- read-only input data, not engine.py code.
"""
import json
import os
import sys
import time

from mpmath import mp, mpf, mpc, pi, cos as mcos, matrix, lu_solve
from mpmath import eig as mp_eig

mp.dps = 40

HERE = os.path.dirname(os.path.abspath(__file__))
LOG_PATH = os.path.join(HERE, "check_progress.log")


def log(msg):
    line = f"{time.strftime('%Y-%m-%dT%H:%M:%S')}  check_shadow_d6  {msg}"
    with open(LOG_PATH, "a", encoding="utf-8") as f:
        f.write(line + "\n")
    print(line, flush=True)


WAVE1_JSON = r"../sec2_spectra_grid_d/spectra_grid.json"


def load_seed(d_str):
    with open(WAVE1_JSON, encoding="utf-8") as f:
        data = json.load(f)
    for p in data["points"]:
        if p["d"] == d_str:
            a48 = [mpf(x) for x in p["_a48"]]
            alpha48 = mpf(p["alpha48"])
            return a48, alpha48
    raise KeyError(d_str)


# ---------------------------------------------------------------------
# G(t), G'(t) evaluation (monomial power series a_0=1, a_1..a_N)
# ---------------------------------------------------------------------
def G_eval(a, t):
    r = a[-1]
    for c in reversed(a[:-1]):
        r = r * t + c
    return r


def Gp_eval(a, t):
    s = mpf(0)
    for j in range(1, len(a)):
        s += j * a[j] * t ** (j - 1)
    return s


# ---------------------------------------------------------------------
# OWN collocation nodes -- Chebyshev-Gauss-LOBATTO interior points of an
# (N+2)-point Lobatto grid on [0,1], i.e. DIFFERENT node family from
# engine.py's colloc_points((1-cos(pi*(i+0.5)/npts))/2.
# ---------------------------------------------------------------------
def lobatto_interior(npts):
    return [(1 - mcos(pi * (i + 1) / (npts + 1))) / 2 for i in range(npts)]


# ---------------------------------------------------------------------
# OWN Newton continuation for the fixed point (power-series / collocation
# match of G(t) + alpha*G(W)=0, W=|G(t/alpha^d)|^d, plus G(1)=-1/alpha from
# t=0 limit). Independent code, different node set from engine.py.
# ---------------------------------------------------------------------
def fp_residual_jac(a_free, alpha, d, N, pts):
    """a_free = [a_1,...,a_N] (a_0=1 fixed). Returns F (N+1,), J (N+1 x N+1)
    where last column/unknown is alpha."""
    a = [mpf(1)] + list(a_free)
    ad = alpha ** d
    F = []
    J = []
    for t in pts:
        ta = t / ad
        U = G_eval(a, ta)
        aU = abs(U)
        sU = 1 if U >= 0 else -1
        dWdU = d * sU * aU ** (d - 1)
        W = aU ** d
        Gt = G_eval(a, t)
        GW = G_eval(a, W)
        GpW = Gp_eval(a, W)
        Up = Gp_eval(a, ta)
        F.append(Gt + alpha * GW)
        row = [mpf(0)] * (N + 1)
        pow_t = mpf(1)
        pow_W = mpf(1)
        pow_ta = mpf(1)
        base = alpha * GpW * dWdU
        for j in range(1, N + 1):
            pow_t *= t
            pow_W *= W
            pow_ta *= ta
            row[j - 1] = pow_t + alpha * pow_W + base * pow_ta
        row[N] = GW - d * t * dWdU * GpW * Up / ad
        J.append(row)
    # extra equation from t->0 limit: G(1) = -1/alpha  <=>  1 + alpha*G(1) = 0
    G1 = G_eval(a, mpf(1))
    F.append(mpf(1) + alpha * G1)
    row = [mpf(0)] * (N + 1)
    pow_1 = mpf(1)
    for j in range(1, N + 1):
        row[j - 1] = alpha  # d/da_j (1 + alpha*G(1)) = alpha * 1^j = alpha
    row[N] = G1
    J.append(row)
    return F, J


def newton_continue(d, N, a_free0, alpha0, tol=mpf("1e-30"), max_it=50):
    pts = lobatto_interior(N)  # N interior points; +1 extra eq from t->0 => N+1 eqs, N+1 unknowns
    a_free = list(a_free0)
    alpha = mpf(alpha0)
    for it in range(max_it):
        F, J = fp_residual_jac(a_free, alpha, d, N, pts)
        nrm = max(abs(x) for x in F)
        if nrm < tol:
            return a_free, alpha, it, nrm
        Jm = matrix(J)
        rhs = matrix([-x for x in F])
        dx = lu_solve(Jm, rhs)
        lam = mpf(1)
        for _h in range(40):
            a_new = [a_free[j] + lam * dx[j] for j in range(N)]
            al_new = alpha + lam * dx[N]
            try:
                F2, _ = fp_residual_jac(a_new, al_new, d, N, pts)
                nrm2 = max(abs(x) for x in F2)
            except Exception:
                nrm2 = mpf("inf")
            if nrm2 < nrm:
                a_free, alpha, nrm = a_new, al_new, nrm2
                break
            lam /= 2
        else:
            if nrm < mpf("1e-20"):
                return a_free, alpha, it, nrm
            raise RuntimeError(f"newton stalled d={d} res={nrm}")
    if nrm < mpf("1e-20"):
        return a_free, alpha, max_it, nrm
    raise RuntimeError(f"newton no convergence d={d} res={nrm}")


def solve_fp(d, N, seed_a_free, seed_alpha, ladder_from=20):
    """Continuation in N first (cheap, from a coarse truncation of the
    seed) then in d via a short step if |seed_d - d| is not tiny -- here we
    always start from an on-target-d seed (Delta-d already applied by
    caller via a short internal Newton chain), so this mainly upsamples N."""
    Ns = sorted(set([n for n in (16, 24, ladder_from, N) if n <= N]))
    a_free = seed_a_free[: max(Ns[0], 1)]
    if len(a_free) < Ns[0]:
        a_free = a_free + [mpf(0)] * (Ns[0] - len(a_free))
    alpha = seed_alpha
    it_tot = 0
    for Nl in Ns:
        cur = (a_free[:Nl] + [mpf(0)] * max(0, Nl - len(a_free)))
        a_free, alpha, it, r = newton_continue(d, Nl, cur, alpha)
        it_tot += it
    return a_free, alpha, it_tot, r


# ---------------------------------------------------------------------
# Basis (A): monomial t^j.  Basis (B): shifted Chebyshev e_j(t) = T_j(2t-1)
# - (-1)^j, j=1..N (own recurrence evaluator, valid for any real argument).
# ---------------------------------------------------------------------
def cheb_T_all(s, N):
    """Return [T_0(s),...,T_N(s)] via the standard 3-term recurrence
    (valid for any real s, not just [-1,1])."""
    Ts = [mpf(1), s]
    for k in range(2, N + 1):
        Ts.append(2 * s * Ts[-1] - Ts[-2])
    return Ts


def basisA_vals(t, N):
    """monomial: returns [t^1,...,t^N]"""
    out = []
    p = mpf(1)
    for j in range(1, N + 1):
        p *= t
        out.append(p)
    return out


def basisB_vals(t, N):
    """shifted Chebyshev e_j(t) = T_j(2t-1) - (-1)^j, j=1..N"""
    s = 2 * t - 1
    Ts = cheb_T_all(s, N)
    out = []
    sign = -1  # (-1)^1
    for j in range(1, N + 1):
        out.append(Ts[j] - (mpf(-1) ** j))
    return out


def build_DT_matrix(a, alpha, d, N, basis_vals, basis_name):
    """Generic builder: A = V^{-1} F on span{basis_1,...,basis_N}."""
    pts = lobatto_interior(N)
    ad = alpha ** d
    G1 = G_eval(a, mpf(1))
    # basis values at t=1, needed for the gauge term's j-dependence
    b1 = basis_vals(mpf(1), N)
    V = matrix(N, N)
    F = matrix(N, N)
    for i, t in enumerate(pts):
        ta = t / ad
        U = G_eval(a, ta)
        aU = abs(U)
        sU = 1 if U >= 0 else -1
        W = aU ** d
        Kt = G_eval(a, t) - d * t * Gp_eval(a, t)
        c3 = -alpha * d * sU * aU ** (d - 1) * Gp_eval(a, W)
        bt = basis_vals(t, N)
        bW = basis_vals(W, N)
        bta = basis_vals(ta, N)
        for j in range(N):
            V[i, j] = bt[j]
            F[i, j] = -(b1[j] / G1) * Kt - alpha * bW[j] + c3 * bta[j]
    LU, p = mp.LU_decomp(V)
    A = matrix(N, N)
    for jcol in range(N):
        bcol = matrix(N, 1)
        for i in range(N):
            bcol[i] = F[i, jcol]
        y = mp.L_solve(LU, bcol, p)
        x = mp.U_solve(LU, y)
        for i in range(N):
            A[i, jcol] = x[i]
    return A


def full_eig_window(a, alpha, d, N, basis_vals, basis_name, window=mpf("0.5")):
    A = build_DT_matrix(a, alpha, d, N, basis_vals, basis_name)
    E, ER = mp_eig(A, left=False, right=True)
    out = []
    for k in range(N):
        lam = E[k]
        is_cplx = abs(lam.imag) > abs(lam) * mpf("1e-9") + mpf("1e-28")
        lamval = lam if is_cplx else lam.real
        if is_cplx:
            if abs(lam) > window:
                continue
        else:
            if abs(lamval) > window:
                continue
        rvec = [ER[i, k] for i in range(N)]
        out.append({"lam": lamval, "vec": rvec, "complex": is_cplx})
    out.sort(key=lambda r: -abs(r["lam"]))
    return out


# ---------------------------------------------------------------------
# common physical grid for the L2/functional cosine (independent of basis)
# ---------------------------------------------------------------------
M_GRID = 120
PHYS_GRID = [mpf(i + 1) / (M_GRID + 1) for i in range(M_GRID)]  # (0,1) interior


def values_monomial(vec, grid):
    N = len(vec)
    out = []
    for t in grid:
        p = mpf(1)
        s = mpf(0)
        for j in range(N):
            p *= t
            s += vec[j] * p
        out.append(s)
    return out


def values_chebB(vec, grid):
    N = len(vec)
    out = []
    for t in grid:
        s2t1 = 2 * t - 1
        Ts = cheb_T_all(s2t1, N)
        s = mpf(0)
        for j in range(1, N + 1):
            ej = Ts[j] - (mpf(-1) ** j)
            s += vec[j - 1] * ej
        out.append(s)
    return out


def cos_raw(v1, v2):
    n = max(len(v1), len(v2))
    x1 = list(v1) + [mpf(0)] * (n - len(v1))
    x2 = list(v2) + [mpf(0)] * (n - len(v2))
    n1 = sum((abs(x) ** 2 for x in x1)) ** mpf("0.5")
    n2 = sum((abs(x) ** 2 for x in x2)) ** mpf("0.5")
    if n1 == 0 or n2 == 0:
        return mpf(0)
    dot = sum(x1[i].conjugate() * x2[i] if hasattr(x1[i], "conjugate") else x1[i] * x2[i] for i in range(n))
    return abs(dot) / (n1 * n2)


def cos_values(vals1, vals2):
    n1 = sum(abs(x) ** 2 for x in vals1) ** mpf("0.5")
    n2 = sum(abs(x) ** 2 for x in vals2) ** mpf("0.5")
    if n1 == 0 or n2 == 0:
        return mpf(0)
    dot = sum((vals1[i].conjugate() if hasattr(vals1[i], "conjugate") else vals1[i]) * vals2[i] for i in range(len(vals1)))
    return abs(dot) / (n1 * n2)


if __name__ == "__main__":
    log("start")
    print("smoke test placeholder -- see check_shadow_d6_run.py for the full sweep")
