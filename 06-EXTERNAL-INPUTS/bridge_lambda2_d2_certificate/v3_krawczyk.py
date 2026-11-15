#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
v3_krawczyk.py -- BONUS interval stage (python-flint 0.9, arb balls), REWRITE 3.

Interval Chebyshev layer in TRIGONOMETRIC form: T_k(cos(theta)) = cos(k*theta),
theta = arccos(2t-1) -- ball radius growth LINEAR in k (Clenshaw/monomial forms
amplify radii by (2|s|)^k resp. 2^{2k} -- both found unusable on 1e-30 boxes).

  K1: Krawczyk for the CF fixed point (b_1..b_48, alpha) -> alpha ~1e-30 cert.
  K2: Krawczyk for lambda2 eigenpair (v, lambda2), (b,alpha) from K1 box (chained)
      -> lambda2 ~1e-24 cert. K2b: same for delta.

Rigorous: K(X) subset int(X) componentwise => unique zero in X (Krawczyk/Moore).
"""
import json
import time

from flint import arb, arb_mat, ctx
from mpmath import mp, mpf, mpmathify, matrix as mpmatrix, lu_solve, cos, pi as mppi, nstr

mp.dps = 140
ctx.dps = 130

TAG = '[kid:verify] interval Krawczyk stage, rewrite 3 (trig-form Chebyshev balls)'
T0 = time.time()
LOG_PATH = 'v3_run.log'


def log(msg):
    line = '[%7.1fs] %s' % (time.time() - T0, msg)
    print(line, flush=True)
    with open(LOG_PATH, 'a', encoding='utf-8') as f:
        f.write(line + '\n')


log(TAG)
_v1 = json.load(open('v1_cheb.json', encoding='utf-8'))
N = 48
ALPHA_STR = _v1['alpha']
X48_STRS = _v1['x48_cheb_b1_to_b48']
LAM2_STR = _v1['lambda2']
DELTA_STR = _v1['delta']
V2_STRS = _v1['v_lambda2_B']
VD_STRS = _v1['v_delta_B']

sgn_int = [(-1) ** k for k in range(N + 1)]


def ball(rstr):
    return arb('-' + rstr).union(arb(rstr))


def arb_eye(n):
    M = arb_mat(n, n)
    for i in range(n):
        M[i, i] = arb(1)
    return M


def gauss_nodes_arb(npts):
    ts = []
    for i in range(npts):
        t = (cos(mppi * (mpf(i) + mpf('0.5')) / npts) + 1) / 2
        ts.append(arb(nstr(t, 120)))
    return ts


class TrigNode:
    """Chebyshev data at ball argument t in (0,1): T_k(2t-1)=cos(k*th), etc."""

    def __init__(self, t, b):
        s = 2 * t - 1
        self.t = t
        th = s.acos()
        self.th = th
        self.T = [arb(1)]
        for k in range(1, N + 1):
            self.T.append((arb(k) * th).cos())
        st = th.sin()
        self.sinth = st
        # G(t) = sum b_k T_k ; G'(t) = (2/sin th) * sum k b_k sin(k th)
        self.G = b[0]
        dsum = arb(0)
        for k in range(1, N + 1):
            self.G = self.G + b[k] * self.T[k]
            dsum = dsum + arb(k) * b[k] * (arb(k) * th).sin()
        self.Gp = 2 * dsum / st


def b0_of(bs):
    s = arb(1)
    for k in range(1, len(bs) + 1):
        s = s - arb(sgn_int[k]) * bs[k - 1]
    return s


# =============================== K1 =============================================
log('=== K1: fixed point Krawczyk (box radius 1e-44) ===')
bhat = [arb(s) for s in X48_STRS]
ahat = arb(ALPHA_STR)
R = '1e-44'
Xb = [b + ball(R) for b in bhat]
Xa = ahat + ball(R)
X_vec = Xb + [Xa]
xhat_vec = bhat + [ahat]

bX = [b0_of(Xb)] + Xb          # ball Chebyshev coeffs on the box
bH = [b0_of(bhat)] + bhat      # at midpoint

pts49 = gauss_nodes_arb(N + 1)


def K1_F_and_J(bs, a, pts, want_J):
    a2 = a * a
    nodes = [TrigNode(t, bs) for t in pts]
    rows = []
    for idx, t in enumerate(pts):
        ta = t / a2
        nta = TrigNode(ta, bs)
        U = nta.G
        W = U * U
        nW = TrigNode(W, bs)
        F_i = nodes[idx].G + a * nW.G
        if want_J:
            row = []
            GpW = nW.Gp
            Gp_ta = nta.Gp
            for j in range(1, N + 1):
                dG_t = nodes[idx].T[j] - arb(sgn_int[j])
                dG_W = nW.T[j] - arb(sgn_int[j])
                dG_ta = nta.T[j] - arb(sgn_int[j])
                row.append(dG_t + a * (dG_W + 2 * U * GpW * dG_ta))
            row.append(nW.G - 4 * t * U * GpW * Gp_ta / a2)
            rows.append((F_i, row))
        else:
            rows.append((F_i, None))
    return rows


log('F(x_hat) (trig balls)...')
resH = K1_F_and_J(bH, ahat, pts49, False)
log('  F(x_hat) max radius %.3e ; max |mid| %.3e'
    % (float(max(f.rad() for f, _ in resH)), float(max(abs(f.mid()) for f, _ in resH))))

log('F(X), J(X) on the box...')
resX = K1_F_and_J(bX, Xa, pts49, True)
FX = [f for f, _ in resX]
log('  F(X) max radius %.3e ; max |mid| %.3e'
    % (float(max(f.rad() for f in FX)), float(max(abs(f.mid()) for f in FX))))
JX = arb_mat(N + 1, N + 1)
for i, (_f, row) in enumerate(resX):
    for j in range(N + 1):
        JX[i, j] = row[j]
log('  J(X) max entry radius %.3e'
    % float(max(JX[i, j].rad() for i in range(N + 1) for j in range(N + 1))))

# mid-Jacobian inverse in mpmath
Jm = mpmatrix(N + 1, N + 1)
for i in range(N + 1):
    for j in range(N + 1):
        Jm[i, j] = mpf(JX[i, j].mid().str(140, radius=False))


def mp_inverse(M, n):
    cols = []
    for jcol in range(n):
        e = mpmatrix(n, 1)
        e[jcol, 0] = mpf(1)
        x = lu_solve(M, e)
        cols.append([x[i, 0] for i in range(n)])
    Y = mpmatrix(n, n)
    for jcol in range(n):
        for i in range(n):
            Y[i, jcol] = cols[jcol][i]
    return Y


log('inverting mid Jacobian (dps 140)...')
Ym = mp_inverse(Jm, N + 1)
Y = arb_mat([[arb(nstr(Ym[i, j], 120)) for j in range(N + 1)] for i in range(N + 1)])
log('  max |Y| = %.3e' % float(max(abs(Ym[i, j]) for i in range(N + 1) for j in range(N + 1))))

Fh = arb_mat(N + 1, 1)
for i in range(N + 1):
    Fh[i, 0] = resH[i][0]
Xm = arb_mat(N + 1, 1)
for i in range(N + 1):
    Xm[i, 0] = X_vec[i]
xhm = arb_mat(N + 1, 1)
for i in range(N + 1):
    xhm[i, 0] = xhat_vec[i]

log('Krawczyk step K1...')
K1m = xhm - Y * Fh + (arb_eye(N + 1) - Y * JX) * (Xm - xhm)
n_contained = sum(1 for i in range(N + 1) if X_vec[i].contains_interior(K1m[i, 0]))
log('K1: contained %d/%d' % (n_contained, N + 1))
k1_ok = (n_contained == N + 1)
if k1_ok:
    ka = K1m[N, 0]
    half = float((ka.rad() + abs(ka.mid() - ahat.mid())).mid())
    log('K1 CERTIFICATE OK: alpha certified to +-%.3e' % half)
else:
    log('K1 FAILED')

# =============================== K2 =============================================
log('=== K2: lambda2 eigenpair Krawczyk (box r_v=1e-18, r_lam=1e-21) ===')
vh = [mpmathify(s).real for s in V2_STRS]  # imag parts are eig-noise ~1e-60
s_odd = sum(vh[j - 1] for j in range(1, N + 1, 2))
vh = [v / (2 * s_odd) for v in vh]
lam2_hat = mpf(LAM2_STR)
pts48 = gauss_nodes_arb(N)

# ball trig nodes with the K1-box coefficients (bX, Xa)
a2X = Xa * Xa
node_t = [TrigNode(t, bX) for t in pts48]
node_ta = []
node_W = []
for t in pts48:
    nta = TrigNode(t / a2X, bX)
    node_ta.append(nta)
    node_W.append(TrigNode(nta.G * nta.G, bX))

DTcol_cache = {}


def DTcol_ball(j, ti):
    key = (j, ti)
    if key in DTcol_cache:
        return DTcol_cache[key]
    nt, nta, nW = node_t[ti], node_ta[ti], node_W[ti]
    t = pts48[ti]
    U = nta.G
    GpW = nW.Gp
    Gp_ta = nta.Gp
    phi1 = arb(1 - sgn_int[j])  # B_j(1) = 1 - (-1)^j
    val = (Xa * phi1 * nt.G
           + 4 * t * phi1 * U * GpW * Gp_ta
           - Xa * (nW.T[j] - arb(sgn_int[j]))
           - 2 * Xa * U * GpW * (nta.T[j] - arb(sgn_int[j])))
    DTcol_cache[key] = val
    return val


def eigen_F_ball(vb, lam):
    out = []
    for ti in range(N):
        t = pts48[ti]
        nt = node_t[ti]
        dt = arb(0)
        ph = arb(0)
        for j in range(1, N + 1):
            c = vb[j - 1]
            dt = dt + c * DTcol_ball(j, ti)
            ph = ph + c * (nt.T[j] - arb(sgn_int[j]))
        out.append(dt - lam * ph)
    norm = arb(-1)
    for j in range(1, N + 1, 2):
        norm = norm + 2 * vb[j - 1]
    out.append(norm)
    return out


# mpmath mirrors for mid-matrices
b_mp = [mpf(1) - sum(((-1) ** k) * mpf(X48_STRS[k - 1]) for k in range(1, N + 1))] + [mpf(s) for s in X48_STRS]
a_mp = mpf(ALPHA_STR)
a2_mp = a_mp ** 2
t48_mp = [(cos(mppi * (mpf(i) + mpf('0.5')) / N) + 1) / 2 for i in range(N)]


class TrigNodeMp:
    def __init__(self, t, b):
        from mpmath import acos, cos as mcos, sin as msin
        s = 2 * t - 1
        th = acos(s)
        self.T = [mpf(1)]
        for k in range(1, N + 1):
            self.T.append(mcos(k * th))
        st = msin(th)
        self.G = b[0]
        dsum = mpf(0)
        for k in range(1, N + 1):
            self.G = self.G + b[k] * self.T[k]
            dsum = dsum + k * b[k] * msin(k * th)
        self.Gp = 2 * dsum / st


node_mp = []
for t in t48_mp:
    nt = TrigNodeMp(t, b_mp)
    nta = TrigNodeMp(t / a2_mp, b_mp)
    W = nta.G * nta.G
    nW = TrigNodeMp(W, b_mp)
    node_mp.append((nt, nta, nW))


def DTcol_mp(j, t):
    nt, nta, nW = node_mp[t48_mp.index(t)]
    U = nta.G
    phi1 = mpf(1 - sgn_int[j])
    return (a_mp * phi1 * nt.G
            + 4 * t * phi1 * U * nW.Gp * nta.Gp
            - a_mp * (nW.T[j] - sgn_int[j])
            - 2 * a_mp * U * nW.Gp * (nta.T[j] - sgn_int[j]))


def eigen_F_mp(v, lam):
    out = []
    for t in t48_mp:
        nt = node_mp[t48_mp.index(t)][0]
        dt = sum(v[j - 1] * DTcol_mp(j, t) for j in range(1, N + 1))
        ph = sum(v[j - 1] * (nt.T[j] - sgn_int[j]) for j in range(1, N + 1))
        out.append(dt - lam * ph)
    out.append(mpf(-1) + sum(2 * v[j - 1] for j in range(1, N + 1, 2)))
    return out


log('mid eigen residual (mpmath dps140): %.3e' % float(max(abs(x) for x in eigen_F_mp(vh, lam2_hat))))


def eigen_J_mp(v, lam):
    M = mpmatrix(N + 1, N + 1)
    for i, t in enumerate(t48_mp):
        nt = node_mp[i][0]
        for j in range(1, N + 1):
            M[i, j - 1] = DTcol_mp(j, t) - lam * (nt.T[j] - sgn_int[j])
        ph = sum(v[j - 1] * (nt.T[j] - sgn_int[j]) for j in range(1, N + 1))
        M[i, N] = -ph
    for j in range(1, N + 1):
        M[N, j - 1] = mpf(2) if j % 2 == 1 else mpf(0)
    M[N, N] = mpf(0)
    return M


J2 = eigen_J_mp(vh, lam2_hat)
Y2m = mp_inverse(J2, N + 1)
Y2 = arb_mat([[arb(nstr(Y2m[i, j], 120)) for j in range(N + 1)] for i in range(N + 1)])
log('Y2 built; max |Y2| = %.3e' % float(max(abs(Y2m[i, j]) for i in range(N + 1) for j in range(N + 1))))

RV = '1e-18'
RL = '1e-21'
vhat_arb = [arb(nstr(v, 105)) for v in vh]   # spectrum real; imag noise dropped
lamhat_arb = arb(nstr(lam2_hat, 110))
Xv = [v + ball(RV) for v in vhat_arb]
Xlam = lamhat_arb + ball(RL)
X2_vec = Xv + [Xlam]
x2hat_vec = vhat_arb + [lamhat_arb]

Fh2 = eigen_F_ball(vhat_arb, lamhat_arb)
log('F(x2_hat) max radius %.3e ; max |mid| %.3e'
    % (float(max(f.rad() for f in Fh2)), float(max(abs(f.mid()) for f in Fh2))))
FX2 = eigen_F_ball(Xv, Xlam)
log('F(X2) max radius %.3e ; max |mid| %.3e'
    % (float(max(f.rad() for f in FX2)), float(max(abs(f.mid()) for f in FX2))))

J2X = arb_mat(N + 1, N + 1)
for i in range(N):
    nt = node_t[i]
    for j in range(1, N + 1):
        J2X[i, j - 1] = DTcol_ball(j, i) - Xlam * (nt.T[j] - arb(sgn_int[j]))
    ph = arb(0)
    for j in range(1, N + 1):
        ph = ph + Xv[j - 1] * (nt.T[j] - arb(sgn_int[j]))
    J2X[i, N] = -ph
for j in range(1, N + 1):
    J2X[N, j - 1] = arb(2) if j % 2 == 1 else arb(0)
J2X[N, N] = arb(0)
log('J2(X2) max entry radius %.3e' % float(max(J2X[i, j].rad() for i in range(N + 1) for j in range(N + 1))))

Fh2m = arb_mat(N + 1, 1)
for i in range(N + 1):
    Fh2m[i, 0] = Fh2[i]
X2m = arb_mat(N + 1, 1)
for i in range(N + 1):
    X2m[i, 0] = X2_vec[i]
x2hm = arb_mat(N + 1, 1)
for i in range(N + 1):
    x2hm[i, 0] = x2hat_vec[i]

log('Krawczyk step K2...')
K2m = x2hm - Y2 * Fh2m + (arb_eye(N + 1) - Y2 * J2X) * (X2m - x2hm)
n2 = 0
fail2 = []
for i in range(N + 1):
    if X2_vec[i].contains_interior(K2m[i, 0]):
        n2 += 1
    else:
        fail2.append((i, float((K2m[i, 0].rad() + abs(K2m[i, 0].mid() - x2hat_vec[i].mid())).mid())))
log('K2: contained %d/%d ; failing (idx, K-radius): %s' % (n2, N + 1, fail2[:4]))
k2_ok = (n2 == N + 1)
lam2_half = None
if k2_ok:
    kl = K2m[N, 0]
    lam2_half = float((kl.rad() + abs(kl.mid() - lamhat_arb.mid())).mid())
    log('K2 CERTIFICATE OK: lambda2 certified to +-%.3e' % lam2_half)
else:
    log('K2 FAILED')

# =============================== K2b: delta ======================================
log('=== K2b: delta eigenpair ===')
vd = [mpmathify(s).real for s in VD_STRS]
sd = sum(vd[j - 1] for j in range(1, N + 1, 2))
vd = [v / (2 * sd) for v in vd]
dhat = mpf(DELTA_STR)
J3 = eigen_J_mp(vd, dhat)
Y3m = mp_inverse(J3, N + 1)
Y3 = arb_mat([[arb(nstr(Y3m[i, j], 120)) for j in range(N + 1)] for i in range(N + 1)])
vd_arb = [arb(nstr(v, 105)) for v in vd]
dhat_arb = arb(nstr(dhat, 110))
Xv3 = [v + ball(RV) for v in vd_arb]
Xd3 = dhat_arb + ball(RL)
X3_vec = Xv3 + [Xd3]
x3hat = vd_arb + [dhat_arb]
Fh3 = eigen_F_ball(vd_arb, dhat_arb)
J3X = arb_mat(N + 1, N + 1)
for i in range(N):
    nt = node_t[i]
    for j in range(1, N + 1):
        J3X[i, j - 1] = DTcol_ball(j, i) - Xd3 * (nt.T[j] - arb(sgn_int[j]))
    ph = arb(0)
    for j in range(1, N + 1):
        ph = ph + Xv3[j - 1] * (nt.T[j] - arb(sgn_int[j]))
    J3X[i, N] = -ph
for j in range(1, N + 1):
    J3X[N, j - 1] = arb(2) if j % 2 == 1 else arb(0)
J3X[N, N] = arb(0)
Fh3m = arb_mat(N + 1, 1)
for i in range(N + 1):
    Fh3m[i, 0] = Fh3[i]
X3m = arb_mat(N + 1, 1)
for i in range(N + 1):
    X3m[i, 0] = X3_vec[i]
x3hm = arb_mat(N + 1, 1)
for i in range(N + 1):
    x3hm[i, 0] = x3hat[i]
K3m = x3hm - Y3 * Fh3m + (arb_eye(N + 1) - Y3 * J3X) * (X3m - x3hm)
n3 = sum(1 for i in range(N + 1) if X3_vec[i].contains_interior(K3m[i, 0]))
log('K2b: contained %d/%d' % (n3, N + 1))
k3_ok = (n3 == N + 1)
delta_half = None
if k3_ok:
    kd = K3m[N, 0]
    delta_half = float((kd.rad() + abs(kd.mid() - dhat_arb.mid())).mid())
    log('K2b CERTIFICATE OK: delta certified to +-%.3e' % delta_half)

out = {
    'meta': {'tag': TAG, 'date': time.strftime('%Y-%m-%d %H:%M:%S'),
             'arb_dps': ctx.dps, 'mpmath_dps': mp.dps,
             'k1_ok': bool(k1_ok), 'k2_lambda2_ok': bool(k2_ok), 'k2b_delta_ok': bool(k3_ok),
             'note': 'phase-1 strings have b_41..b_48 = 0.0 (below Newton tol); K1 box 1e-30 covers true values',
             'elapsed_sec': time.time() - T0},
    'K1_fixed_point': {'box_radius': R, 'contained': '%d/%d' % (n_contained, N + 1),
                       'alpha_mid': ALPHA_STR},
    'K2_lambda2': {'box_r_v': RV, 'box_r_lambda': RL, 'contained': '%d/%d' % (n2, N + 1),
                   'lambda2_mid': LAM2_STR},
    'K2b_delta': {'delta_mid': DELTA_STR},
}
with open('v3_krawczyk.json', 'w', encoding='utf-8') as f:
    json.dump(out, f, indent=2, ensure_ascii=False)
log('wrote v3_krawczyk.json')
