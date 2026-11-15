# [oid:quartK5] NOT RIGOROUS — аудит оркестратора, третий кодовый путь.
# Никаких заимствований кода kilo: свой цикл H4, своя орбита M-стороны,
# своя сопровождающая производная. Данные — только числа из JSON.
import json
from mpmath import mp, mpf, mpc, fabs

mp.dps = 50

fx = json.load(open('quartic_fixedpoint.json', encoding='utf-8'))
mr = json.load(open('match_results.json', encoding='utf-8'))
hz = json.load(open('H4_zeros.json', encoding='utf-8'))

G = [mpf(str(c)) for c in fx['g_coeffs']]
H = [mpf(str(c)) for c in fx['h_coeffs']]
delta4 = mpf(fx['delta4'])
AU04 = mpf(hz['AU0_4'])
a_F = mpf(mr['a_F'])
A4 = mpf(mr['A4'])
M = mr['meta']['frozen']['M_MATCH']          # 9
PER = mr['meta']['frozen']['PERIOD']         # 1024

def series(coeffs, x):
    # sum coeffs[k] * x^(4k) — Horner по t = x^4
    t = x**4
    acc = coeffs[-1]
    for c in reversed(coeffs[:-1]):
        acc = acc * t + c
    return acc

def H4(u, m):
    s = u * AU04 / delta4**m
    x = mpc(0)
    n2m = 2**m
    xm = None
    for n in range(1, 2**(m+1) + 1):
        x = series(G, x) + s * series(H, x)
        if n == n2m:
            xm = x
    return x / xm

def orbit_mside(a, nsteps):
    # f_a(x) = 1 - a*x^4, x0 = 0 (критическая точка), сопровождающая производная по a
    x = mpc(0); da = mpc(0)
    snaps = {}
    for n in range(1, nsteps + 1):
        da = -4 * a * x**3 * da - x**4
        x = 1 - a * x**4
        if n & (n - 1) == 0:
            snaps[n] = x
    return x, da, snaps

def parse_c(s):
    return mpc(complex(s.replace(' ', '').replace('(', '').replace(')', '')))

by_id = {r['id']: r for r in mr['catalog_results']}
PICK = ['Q01', 'Q05', 'Q10', 'Q28']

print("=== АУДИТ 1: |H4_m(u_cat)| моим циклом (нуль ли это) ===")
for qid in PICK + ['Q57']:
    r = by_id[qid]
    u = mpc(r['u_cat_re'], r['u_cat_im'])
    h8 = H4(u, 8)
    h9 = H4(u, 9)
    print(f"{qid}: |H4_8|={float(fabs(h8)):.3e}  |H4_9|={float(fabs(h9)):.3e}  u=({r['u_cat_re']:.4f},{r['u_cat_im']:.4f})")

print("\n=== АУДИТ 2: M-сторона — a_ss суперстабилен с ТОЧНЫМ периодом? (моя орбита) ===")
for qid in PICK:
    r = by_id[qid]
    a = parse_c(r['a_ss'])
    x_end, da, snaps = orbit_mside(a, PER)
    min_early = min(float(fabs(snaps[2**j])) for j in range(0, M + 1) if 2**j < PER)
    # два моих шага Ньютона — корень настоящий или рядом?
    a1 = a - x_end / da
    x2, da2, _ = orbit_mside(a1, PER)
    a2 = a1 - x2 / da2
    print(f"{qid}: |x_1024|={float(fabs(x_end)):.3e}  min|x_2^j, j<10|={min_early:.3e}  "
          f"|Ньютон-сдвиг1|={float(fabs(a1 - a)):.3e}  |Ньютон-сдвиг2|={float(fabs(a2 - a1)):.3e}")

print("\n=== АУДИТ 3: калибровка u = A4*(a_ss - a_F)*delta4^m против u_cat ===")
for qid in PICK:
    r = by_id[qid]
    a = parse_c(r['a_ss'])
    u_sc = A4 * (a - a_F) * delta4**M
    u_cat = mpc(r['u_cat_re'], r['u_cat_im'])
    d = float(fabs(u_sc - u_cat))
    print(f"{qid}: мой d_u={d:.6e}  kilo d_u={r['d_u']:.6e}  совпадение={'ДА' if abs(d - r['d_u'])/max(r['d_u'],1e-300) < 1e-3 else 'НЕТ'}")

print("\n=== АУДИТ 4: закон дрейфа — 1/lambda2 ===")
lam2 = mpf(str(fx['checks']['eigen2_deflated']))
print(f"lambda2={float(lam2):.6f}  1/lambda2={float(1/lam2):.4f}  Q01 r89={mr['rate_study'][0]['ratio_89']:.4f}  r910={mr['rate_study'][0]['ratio_910']:.4f}")

print("\n=== АУДИТ 5: порча — сдвиг a_ss на 1e-9 должен убить суперстабильность ===")
r = by_id['Q05']
a_bad = parse_c(r['a_ss']) + mpf('1e-9')
x_end, _, _ = orbit_mside(a_bad, PER)
print(f"Q05+1e-9: |x_1024|={float(fabs(x_end)):.3e} (должно быть НЕ мало)")
