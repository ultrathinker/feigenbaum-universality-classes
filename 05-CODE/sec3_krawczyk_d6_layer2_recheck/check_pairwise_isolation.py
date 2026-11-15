# -*- coding: utf-8 -*-
# Independent check (Sonnet reviewer): after shrinking S078/S090 to r=1e-34
# (same center), does pairwise isolation among all 287 certified boxes still
# hold, and does the margin for S078/S090 vs their nearest neighbors improve?
import json
import flint

flint.ctx.prec = 512
A = flint.arb

orig = json.load(open(r"../sec3_krawczyk_d6/krawczyk6_certificates.json", encoding="utf-8"))
recert = json.load(open("S090_S078_recert.json", encoding="utf-8"))
recert_by = {r["label"]: r for r in recert["roots"]}

records = []
for r in orig["roots"]:
    rec = dict(r)
    if r["label"] in ("S078", "S090"):
        rec["box_center_re"] = recert_by[r["label"]]["box_center_re"]
        rec["box_center_im"] = recert_by[r["label"]]["box_center_im"]
        rec["box_radius"] = recert_by[r["label"]]["box_radius_str"]  # "1e-34"
    records.append(rec)

def worst_and_neighbors(records, labels):
    worst = None
    worst_pair = None
    nn = {lab: (None, None) for lab in labels}
    all_disjoint = True
    for i in range(len(records)):
        for j in range(i + 1, len(records)):
            a, b = records[i], records[j]
            dre = A(a["box_center_re"]) - A(b["box_center_re"])
            dim = A(a["box_center_im"]) - A(b["box_center_im"])
            rsum = A(a["box_radius"]) + A(b["box_radius"])
            sep_re = abs(dre).lower()
            sep_im = abs(dim).lower()
            sep = sep_re if bool(sep_re > sep_im) else sep_im
            margin = sep / rsum.upper()
            disjoint = bool(sep > rsum.upper())
            all_disjoint = all_disjoint and disjoint
            if worst is None or bool(margin < worst):
                worst, worst_pair = margin, (a["label"], b["label"])
            for lab in labels:
                if a["label"] == lab or b["label"] == lab:
                    other = b["label"] if a["label"] == lab else a["label"]
                    cur = nn[lab]
                    if cur[0] is None or bool(margin < cur[0]):
                        nn[lab] = (margin, other)
    return all_disjoint, worst, worst_pair, nn

all_disjoint, worst, worst_pair, nn = worst_and_neighbors(records, ("S078", "S090"))
print("all_disjoint (287 boxes, S078/S090 shrunk to r=1e-34):", all_disjoint)
print("worst pair overall:", worst_pair, "margin_lb=", worst)
for lab in ("S078", "S090"):
    print(lab, "nearest neighbor margin:", nn[lab][1], nn[lab][0])

# Compare against ORIGINAL (unshrunk) margins for S078/S090 to show
# isolation only improves.
all_disjoint0, worst0, worst_pair0, nn0 = worst_and_neighbors(orig["roots"], ("S078", "S090"))
print("--- original (unshrunk) ---")
print("all_disjoint:", all_disjoint0, "worst pair:", worst_pair0, worst0)
for lab in ("S078", "S090"):
    print(lab, "nearest neighbor margin (orig):", nn0[lab][1], nn0[lab][0])
