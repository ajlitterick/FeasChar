"""Check that every parsed memoir row has the right module dimensions (sum of deg * mult).

Usage: python3 check_dims.py tables_gen.json
"""
import json, re, sys

def deg(lab):
    return int(re.match(r'\d+', lab).group(0))

def count(lab):
    return 1 + lab.count(',')

tables = json.load(open(sys.argv[1]))
bad = 0
for i, t in enumerate(tables, 1):
    dims = [int(re.sub(r'\D', '', m['name'])) for m in t['mods']]
    for r in t['rows']:
        for g, md in enumerate(t['mods']):
            s = sum(deg(l) * count(l) * m for l, m in zip(md['labels'], r['mult'][g]))
            if s != dims[g]:
                bad += 1
                print('#%d %s row %d, %s: sum %d != %d' % (i, t['caption'][:45], r['n'], md['name'], s, dims[g]))
print(bad, 'bad rows')
