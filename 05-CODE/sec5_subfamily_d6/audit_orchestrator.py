# -*- coding: utf-8 -*-
# АУДИТ ОРКЕСТРАТОРА (ворота 2, задача G): третий кодовый путь.
# Независим от common.py: свой орбитный код (x^6 = (x^2)^3 иной ассоциацией),
# недемпфированный Ньютон, свои формулы u<->a. Константы — из тех же первичных
# JSON пилота B (a_F, A6, delta6), что и у агента G.
import json, sys
from mpmath import mp, mpf, mpc, fabs

mp.dps = 100

FP = json.load(open('../sec3_pilot_d6/sextic_fixedpoint.json', encoding='utf-8'))
MS = json.load(open('../sec3_pilot_d6/mside_cascade.json', encoding='utf-8'))
TG = json.load(open('targets.json', encoding='utf-8'))
RF = json.load(open('rate_full.json', encoding='utf-8'))
D6 = mpf(FP['delta6']); AF = mpf(MS['a_F']); A6 = mpf(MS['A6'])
print(f"[конст] delta6={mp.nstr(D6,17)} a_F={mp.nstr(AF,17)} A6={mp.nstr(A6,17)}", flush=True)

def get_ucat(sid):
    it = TG['targets'][sid]
    return mpc(mpf(repr(it['u_re'])), mpf(repr(it['u_im'])))

def orbit(a, n):
    x = mpc(0); dx = mpc(0)
    for _ in range(n):
        x2 = x*x; x3 = x2*x; x5 = x2*x3; x6 = x3*x3
        dx = -x6 - 6*a*x5*dx
        x = 1 - a*x6
        if fabs(x) > 100: return x, dx, True
    return x, dx, False

def check(sid, m, du_lo, du_hi):
    u = get_ucat(sid)
    P = 2**(m+1)
    a = AF + u/(A6 * D6**m)          # прямой сид (tier-1), без поиска
    ok_conv = False
    for it in range(60):
        F, dF, esc = orbit(a, P)
        if esc: break
        if fabs(F) < mpf('1e-50'): ok_conv = True; break
        a = a - F/dF
    F, dF, _ = orbit(a, P)
    # точный период: min по j=1..m |x_{2^j}| отделён от нуля
    minper = min(fabs(orbit(a, 2**j)[0]) for j in range(1, m+1))
    du = fabs(A6*(a-AF)*D6**m - u)
    ok = ok_conv and minper > mpf('0.005') and du_lo <= du <= du_hi
    print(f"[{sid} m={m} P={P}] conv={ok_conv} |F|={mp.nstr(fabs(F),3)} "
          f"min|x_2^j|={mp.nstr(minper,4)} d_u={mp.nstr(du,4)} "
          f"(заявлено {mp.nstr(mpf(str(RF[sid]['d_u'][str(m)])),4)}, tag={RF[sid]['tag'][str(m)]}) -> {'OK' if ok else 'FAIL'}", flush=True)
    return ok

ok1 = check('S036', 12, mpf('1e-5'), mpf('1e-3'))     # заявлено 2.244e-4, direct
ok2 = check('S070', 10, mpf('1e-4'), mpf('3e-3'))     # заявлен первый чистый хит на m=10
ok3 = check('S036', 9, mpf('0.02'), mpf('0.06'))      # контроль: на m=9 d_u крупный (~0.037)
print('\nВЕРДИКТ АУДИТА G:', 'ЗЕЛЁНЫЙ' if (ok1 and ok2 and ok3) else 'ЕСТЬ ПРОБЛЕМЫ')
sys.exit(0 if (ok1 and ok2 and ok3) else 1)
