import os
BASE = r"../sec7_principle_status_F1F2"
OUT = os.path.dirname(os.path.abspath(__file__))
text = open(os.path.join(BASE, "oc_run.log"), encoding="utf-8").read()
idx = text.find("S036   d(m)=")
if idx == -1:
    idx = text.find("S036")
chunk = text[max(0,idx-200):idx+1600]
with open(os.path.join(OUT, "rate_full_dump.txt"), "w", encoding="utf-8") as f:
    f.write(chunk)
print("wrote", len(chunk), "chars")
