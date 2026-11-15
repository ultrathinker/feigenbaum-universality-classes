# -*- coding: utf-8 -*-
# АУДИТ ОРКЕСТРАТОРА (задача L): третий кодовый путь.
# (1) подлинность каскада: свой орбитный код, остаток+точный период у выборки уровней;
# (2) свой анализ слоя поправок: r(m)=R(m)-δ_ext, знаки и скорость БЕЗ их скриптов.
import json, sys
from mpmath import mp, mpf, fabs, log

mp.dps = 150

def orbit(a, n, d):
    x = mpf(0)
    for _ in range(n):
        x = 1 - a*(abs(x)**d)
    return x

def audit(fn, dstr, delta_ext, muA, muB):
    D = json.load(open(fn, encoding='utf-8'))
    lv = {int(k): mpf(v['a']) for k, v in D['levels'].items()}
    d = mpf(dstr)
    ks = sorted(lv)
    print(f"\n=== d={dstr}: уровни {ks[0]}..{ks[-1]} ===", flush=True)
    # (1) подлинность: уровни 10 и 13 — остаток и точный период своим кодом
    for k in (10, 13):
        a = lv[k]; P = 2**k
        res = fabs(orbit(a, P, d))
        minper = min(fabs(orbit(a, 2**j, d)) for j in range(1, k))
        print(f"  [подлинность] m={k}: |x_{{2^{k}}}|={mp.nstr(res,3)}  min_j|x_2^j|={mp.nstr(minper,3)}", flush=True)
        assert res < mpf('1e-40') and minper > mpf('1e-6')  # min|x_2^j| ~ alpha^-(m-1), сжимается с уровнем
    # (2) свой анализ: R(m), r(m)=R(m)-δ_ext, знаки, скорость
    R = {m: (lv[m]-lv[m-1])/(lv[m+1]-lv[m]) for m in ks[1:-1]}
    r = {m: R[m] - delta_ext for m in R}
    ms = sorted(r)
    signs = ['+' if r[m] > 0 else '-' for m in ms]
    print("  m:", ms)
    print("  sign(r):", ''.join(signs))
    # офсет-иммунный анализ: s(m)=r(m+1)-r(m) гасит постоянное смещение δ_ext;
    # s(m) = c·μ^m·(μ-1) → s(m+1)/s(m) = μ, знак s не чередуется ⇔ μ>0
    s = {m: r[m+1]-r[m] for m in ms[:-1]}
    sm = sorted(s)
    ssigns = ['+' if s[m] > 0 else '-' for m in sm]
    print("  sign(s):", ''.join(ssigns))
    print("  s(m+1)/s(m) хвост:")
    for m in sm[-6:-1]:
        q = s[m+1]/s[m]
        print(f"    m={m}->{m+1}: μ={mp.nstr(q,6)}  (1/μ={mp.nstr(1/q,5)})", flush=True)
    tail = [s[m+1]/s[m] for m in sm[-5:-1]]
    med = sorted(tail)[len(tail)//2]
    print(f"  медиана хвоста μ={mp.nstr(med,5)} → скорость 1/μ={mp.nstr(1/med,5)}; "
          f"кандидаты: A μ={muA} (1/μ={mp.nstr(1/mpf(muA),4)}), B μ={muB} (алътерн.)")
    flips = sum(1 for i in range(len(ssigns)-1) if ssigns[i] != ssigns[i+1])
    print(f"  смен знака s(m): {flips} из {len(ssigns)-1}")
    return med, flips

mA, fA = audit('cascade_d300.json', '3.0', mpf('6.08468907426212189604612106626'), '0.23714611165', '-0.19849370741')
mB, fB = audit('cascade_d275.json', '2.75', mpf('5.75581828327934'), '0.22065', '-0.18190')
okA = fA <= 1 and mpf('0.18') < mA < mpf('0.30')
okB = fB <= 1 and mpf('0.17') < mB < mpf('0.30')
print("\nВЕРДИКТ АУДИТА L:", "ЗЕЛЁНЫЙ — монотонность и скорость слоя подтверждены моим анализом" if (okA and okB) else "ПРОВЕРИТЬ ВРУЧНУЮ")
