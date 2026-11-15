import os
BASE = r"../sec7_principle_status_F1F2"
OUT = os.path.dirname(os.path.abspath(__file__))
text = open(os.path.join(BASE, "oc_run.log"), encoding="utf-8").read()
for needle in ["d6 m=8:", "d6 m=9:", "d4 m=8:", "d4 m=9:"]:
    idx = text.find(needle)
    print(needle, "idx=", idx)
    if idx != -1:
        chunk = text[idx:idx+400]
        with open(os.path.join(OUT, "f2c_hits.txt"), "a", encoding="utf-8") as f:
            f.write("=== %s ===\n%s\n\n" % (needle, chunk))
print("done")
