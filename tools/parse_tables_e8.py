"""Parse the E8 feasible-character tables of AJL-thetables.tex into JSON.

Usage: python -I parse_tables.py <AJL-thetables.tex> <out.json>
"""
import sys, re, json, io

src = io.open(sys.argv[1], encoding='utf-8').read()
start = src.index(r'\section{$E_8$}')
body = src[start:]

def clean_label(s):
    s = s.strip()
    s = s.replace('$', '').replace('{', '').replace('}', '').replace(' ', '')
    s = s.replace('\\ast', '*').replace('^*', '*')
    s = s.replace('(', '').replace(')', '')
    return s

ENV = r'\\begin\{(table|longtable)\}(.*?)\\end\{\1\}'
tables = []
for m in re.finditer(ENV, body, re.S):
    kind, block = m.group(1), m.group(2)
    cap = re.search(r'\\caption\{(.*?)\}\s*(\\\\)?\s*\n', block)
    if not cap:
        continue
    caption = cap.group(1).strip()
    if kind == 'table':
        tab = re.search(r'\\begin\{tabular\}\{(.*?)\}\s*\n(.*?)\\end\{tabular\}', block, re.S)
        if not tab:
            continue
        spec, rows_tex = tab.group(1), tab.group(2)
    else:
        spec = re.match(r'\{(.*?)\}\s*\n', block).group(1)
        rows_tex = block[cap.end():]
    has_flag = spec.startswith('rc')
    header, rows = None, []
    for l in rows_tex.split('\n'):
        l2 = l.replace('\\hline', '').replace('\\endhead', '').replace('\\endfirsthead', '').strip()
        l2 = re.sub(r'\\\\\s*$', '', l2).strip()
        if not l2 or l2.startswith('\\label') or l2.startswith('%'):
            continue
        cells = [c.strip() for c in l2.split('&')]
        if re.match(r'^\d+\)$', cells[0]):
            i = 0
            while i < len(cells):
                rm = re.match(r'^(\d+)\)$', cells[i])
                if not rm:
                    raise ValueError('bad row in %s: %s' % (caption, l2))
                j = i + 1
                flags = ''
                # flag column present unless the next cell is already a multiplicity
                # (some side-by-side tables have a flag column only in the left half)
                if has_flag and j < len(cells) and not re.match(r'^\d+$', cells[j]):
                    flags = cells[j]; j += 1
                k = j
                while k < len(cells) and not re.match(r'^\d+\)$', cells[k]):
                    k += 1
                vals = [int(v) for v in cells[j:k] if v != '']
                rows.append({'n': int(rm.group(1)), 'P': '\\possprim' in flags,
                             'N': '\\nongcr' in flags, 'mult': vals})
                i = k
        elif header is None and cells[0] == '' and 'multicolumn' not in l2:
            labs = [clean_label(c) for c in cells if c != '']
            half = len(labs) // 2
            if '||' in spec or (len(labs) % 2 == 0 and labs[:half] == labs[half:]):
                labs = labs[:half]
            header = labs
    after = body[m.end():m.end() + 300]
    perm = re.search(r'Permutations?:(.*?)\n', block + '\n' + after.split('\\begin')[0])
    tables.append({'caption': caption, 'labels': header, 'rows': sorted(rows, key=lambda r: r['n']),
                   'permutations': perm.group(1).strip() if perm else None})

bad = [t['caption'] for t in tables
       if t['labels'] is None or any(len(r['mult']) != len(t['labels']) for r in t['rows'])]
io.open(sys.argv[2], 'w', encoding='utf-8').write(json.dumps(tables, indent=1))
print(len(tables), 'tables;', sum(len(t['rows']) for t in tables), 'rows;',
      sum(r['P'] for t in tables for r in t['rows']), 'P-rows')
print('malformed:', bad)
