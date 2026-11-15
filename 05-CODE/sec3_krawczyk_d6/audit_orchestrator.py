# -*- coding: utf-8 -*-
# АУДИТ ОРКЕСТРАТОРА (ворота 2/3, задача H «Кравчик-6»): третий кодовый путь,
# mpmath.mpc — независим от flint-генератора и кортежного верификатора.
import json, sys
from mpmath import mp, mpf, mpc, fabs

certs = json.load(open('krawczyk6_certificates.json', encoding='utf-8'))
match = json.load(open('../sec3_pilot_d6/match_results.json', encoding='utf-8'))
roots = certs['roots']
labels = [r['label'] for r in roots]
print(f"[1] инвентарь: {len(roots)} корней, слой2:", [r['label'] for r in roots if 'e-3' in str(r.get('box_radius',''))][:5], flush=True)
assert len(roots) == 287 and len(set(labels)) == 287

# сиды слепого прогона
def find_seeds(obj, out):
    if isinstance(obj, dict):
        if 'a_ss_re' in obj: out.append(obj)
        for v in obj.values(): find_seeds(v, out)
    elif isinstance(obj, list):
        for v in obj: find_seeds(v, out)
seeds = []
find_seeds(match, seeds)
sby = {}
for s in seeds:
    for k in ('label','id','name'):
        if k in s: sby[str(s[k])] = s; break
print(f"[2] сидов в match_results: {len(seeds)}, с метками: {len(sby)}", flush=True)

mism = 0; l2checked = 0
for r in roots:
    s = sby.get(r['label'])
    if s is None: continue
    layer2 = float(mpf(r['box_radius'].split('[')[-1].split(' ')[0]) if '[' in str(r['box_radius']) else mpf(str(r['box_radius']))) < 1e-28
    if not layer2:
        if str(s['a_ss_re'])[:28] != str(r['box_center_re'])[:28]:
            mism += 1; print(f"    MISMATCH {r['label']}")
    else:
        # слой 2: центр = уточнённый сид; пин тождественности |c_ref - сид| < 5e-24
        c = mpc(mpf(r['box_center_re']), mpf(r['box_center_im']))
        st = mpc(mpf(str(s['a_ss_re'])), mpf(repr(float(s['a_ss_im']))))
        dpin = fabs(c - st)
        print(f"    слой2 {r['label']}: |c_ref-сид|={mp.nstr(dpin,3)} (< 5e-24: {float(dpin)<5e-24})", flush=True)
        if float(dpin) >= 5e-24: mism += 1
        l2checked += 1
print(f"[2] строковых несовпадений: {mism}; слой2 проверен: {l2checked}", flush=True)

# [3] изоляция независимо
cs = [(complex(float(mpf(r['box_center_re'])), float(mpf(r['box_center_im']))), r['label']) for r in roots]
worst = (1e99, '', '')
for i in range(287):
    for j in range(i+1, 287):
        d = abs(cs[i][0]-cs[j][0])
        if d < worst[0]: worst = (d, cs[i][1], cs[j][1])
print(f"[3] изоляция: худшая пара {worst[1]}/{worst[2]}, |Δa| = {worst[0]:.3e}", flush=True)
assert worst[0] > 1e-11

# [4] независимый Ньютон (свой орбитный код: x6 = (x^3)^2 иной ассоциацией)
def orbit(a, n=1024, mins=False):
    x = mpc(0); dx = mpc(0); mm = {}; nxt = 2
    for k in range(1, n+1):
        x2 = x*x; x3 = x2*x; x6 = x3*x3
        dx = -x6 - 6*a*(x3*x2)*dx
        x = 1 - a*x6
        if mins and k == nxt and k < n: mm[k] = fabs(x); nxt *= 2
    return x, dx, mm

def probe(lab, dps, rmax):
    mp.dps = dps
    r = next(rr for rr in roots if rr['label'] == lab)
    c = mpc(mpf(r['box_center_re']), mpf(r['box_center_im']))
    a = c
    for it in range(15):
        F, dF, _ = orbit(a)
        if fabs(F) < mpf('1e-'+str(dps-15)): break
        a = a - F/dF
    F, dF, mm = orbit(a, mins=True)
    dist = fabs(a-c); minper = min(mm.values())
    ok = fabs(F) < mpf('1e-40') and dist < mpf(rmax) and minper > mpf('0.01')
    print(f"[4] {lab} (dps{dps}): |F|={mp.nstr(fabs(F),3)} |a*-c|={mp.nstr(dist,3)} (r_max {rmax}) "
          f"min|x_2^j|={mp.nstr(minper,3)} |F'|={mp.nstr(fabs(dF),3)} -> {'OK' if ok else 'FAIL'}", flush=True)
    return ok

ok = True
for lab in ('S196','S203','S076'):
    ok &= probe(lab, 70, '1e-22')
for lab in ('S078','S090'):
    ok &= probe(lab, 140, '1e-29')

print('\nВЕРДИКТ АУДИТА H:', 'ЗЕЛЁНЫЙ' if (ok and mism == 0) else 'ЕСТЬ ПРОБЛЕМЫ')
sys.exit(0 if (ok and mism == 0) else 1)
