# -*- coding: utf-8 -*-
# АУДИТ ОРКЕСТРАТОРА (ворота 2, третий кодовый путь): mpmath.mpc, dps 70.
# Независим от генератора (flint acb) и верификатора (ручные кортежи re/im/r).
import json, sys, re
from mpmath import mp, mpf, mpc, fabs

mp.dps = 70

certs = json.load(open('krawczyk4_certificates.json', encoding='utf-8'))
match = json.load(open('../sec3_pilot_d4/match_results.json', encoding='utf-8'))
roots = certs['roots']
print(f"[1] инвентарь: {len(roots)} корней; метки:", flush=True)
labels = [r['label'] for r in roots]
assert len(roots) == 61, "не 61"
assert 'Q57' not in labels, "Q57 должен отсутствовать"
assert len(set(labels)) == 61, "дубли меток"
print("    OK: 61 уникальная метка, Q57 отвергнут", flush=True)

# [2] центры = сиды пилота (строково). Ищем сиды в match_results.
def find_seeds(obj, out):
    if isinstance(obj, dict):
        if 'a_ss' in obj:
            out.append(obj)
        for v in obj.values(): find_seeds(v, out)
    elif isinstance(obj, list):
        for v in obj: find_seeds(v, out)
seeds = []
find_seeds(match, seeds)
print(f"[2] найдено {len(seeds)} записей с a_ss в match_results.json", flush=True)
def keyf(s):
    for k in ('label','id','name'):
        if k in s: return s[k]
    return None
seed_by_label = {}
for s in seeds:
    k = keyf(s)
    if k: seed_by_label[str(k)] = s
mismatch = 0
checked = 0
for r in roots:
    lab = r['label']
    s = seed_by_label.get(lab)
    if s is None:
        continue
    a = s['a_ss']
    if isinstance(a, str):
        # сид может быть комплексной строкой или парой
        cre, cim = r['box_center_re'], r['box_center_im']
        if cre.rstrip('0').rstrip('.') not in a and cre[:25] not in a:
            mismatch += 1
            print(f"    MISMATCH {lab}: center_re={cre[:30]} vs a_ss={str(a)[:60]}")
        checked += 1
    elif isinstance(a, dict):
        cre = r['box_center_re']
        are = str(a.get('re', a.get('real','')))
        if are[:25] != cre[:25]:
            mismatch += 1
            print(f"    MISMATCH {lab}: {cre[:28]} vs {are[:28]}")
        checked += 1
print(f"    сверено {checked}, несовпадений {mismatch}", flush=True)

# [3] попарная изоляция независимо (float достаточно: разделение ~1e-9 >> r=1e-26)
import math
cs = [(complex(float(mpf(r['box_center_re'])), float(mpf(r['box_center_im']))), r['label']) for r in roots]
worst = (1e99, '', '')
for i in range(61):
    for j in range(i+1, 61):
        d = abs(cs[i][0]-cs[j][0])
        if d < worst[0]: worst = (d, cs[i][1], cs[j][1])
print(f"[3] изоляция: худшая пара {worst[1]}/{worst[2]}, |Δa| = {worst[0]:.3e} (r=1e-26 → запас {worst[0]/2e-26:.1e}x)", flush=True)
assert worst[0] > 1e-10

# [4] выборочная независимая проверка корней: Ньютон в mpc
def orbit_F_dF(a, n=1024, want_mins=False):
    x = mpc(0); dx = mpc(0)
    mins = {}
    nxt = 2
    for k in range(1, n+1):
        # x_{k+1} = 1 - a x_k^4 ; dx_{k+1} = -x_k^4 - 4 a x_k^3 dx
        x2 = x*x; x3 = x2*x; x4 = x2*x2
        dx = -x4 - 4*a*x3*dx
        x = 1 - a*x4
        if want_mins and k == nxt and k < n:
            mins[k] = fabs(x)
            nxt *= 2
    return x, dx, mins

SAMPLE = ['Q01','Q10','Q28','Q39','Q51','Q52','Q61']
rby = {r['label']: r for r in roots}
allok = True
for lab in SAMPLE:
    r = rby[lab]
    c = mpc(mpf(r['box_center_re']), mpf(r['box_center_im']))
    a = c
    for it in range(12):
        F, dF, _ = orbit_F_dF(a)
        step = F/dF
        a = a - step
        if fabs(step) < mpf('1e-55'): break
    F, dF, mins = orbit_F_dF(a, want_mins=True)
    resid = fabs(F); dist = fabs(a-c); minper = min(mins.values())
    ok = (resid < mpf('1e-45')) and (dist < mpf('1e-26')) and (minper > mpf('0.005'))
    allok &= ok
    print(f"[4] {lab}: |F(a*)|={mp.nstr(resid,3)} |a*-c|={mp.nstr(dist,3)} "
          f"min|x_2^j|={mp.nstr(minper,4)} |F'|={mp.nstr(fabs(dF),4)} -> {'OK' if ok else 'FAIL'}", flush=True)

print(f"\nВЕРДИКТ АУДИТА (третий путь): {'ВСЕ ПРОВЕРКИ ЗЕЛЁНЫЕ' if (allok and mismatch==0) else 'ЕСТЬ ПРОБЛЕМЫ'}")
sys.exit(0 if (allok and mismatch == 0) else 1)
