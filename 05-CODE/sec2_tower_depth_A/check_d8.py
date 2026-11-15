import json
with open('f5_depth_table.json') as f:
    data = json.load(f)
pts = data['points']
print('All d values:')
for p in pts:
    print(repr(p['d']), type(p['d']))
