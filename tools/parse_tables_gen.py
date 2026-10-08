"""Parse the F4, E6 and E7 feasible-character tables of AJL-thetables.tex (two column groups,
L(G) | V_min) into a GAP file.

Usage: python3 parse_tables_gen.py <AJL-thetables.tex> <out.g> [<out.json>]

Each table becomes
  rec(type, caption, group (GAP name), p, mods := [rec(name, labels), rec(name, labels)],
      rows := [rec(n, P, N, mult := [adjoint multiplicities, V_min multiplicities])],
      perms, notes)
Labels are as printed (e.g. "8_a", "10*").  When a caption says "the duals of X, Y also occur,
with the same multiplicity", a column "X*" is added for each, with X's multiplicities.
The one-line "Irreducible on V and V_min" entries are added by hand.
"""
import sys, re, io, json

src = io.open(sys.argv[1], encoding='utf-8').read()
sec = {}
for typ, nxt in [('F_{4}', 'E_6'), ('E_6', 'E_7'), ('E_7', 'E_8')]:
    a = src.index(r'\section{$%s$}' % typ)
    b = src.index(r'\section{$%s$}' % nxt)
    sec[typ.replace('_', '').replace('{', '').replace('}', '')] = src[a:b]

def clean_label(s):
    s = s.strip().replace('$', '').replace('{', '').replace('}', '').replace(' ', '')
    s = s.replace('\\ast', '*').replace('^*', '*').replace('(', '').replace(')', '')
    return s

def gap_name(cap):
    g = cap.split('<')[0].strip().strip('$').strip()
    g = g.replace('$', '').replace('\\cdot', '.').replace(' ', '').replace('{', '').replace('}', '')
    g = re.sub(r'Alt_(\d+)', r'A\1', g)
    m = re.match(r'SL_2\((\d+)\)', g)
    if m:
        return '2.L2(%s)' % m.group(1)
    g = g.replace('_', '')
    rep = {'^2B2(8)': 'Sz(8)', "^2F4(2)'": "2F4(2)'", '^3D4(2)': '3D4(2)', 'Sp6(2)': 'S6(2)',
           '\\Omega8^+(2)': 'O8+(2)', 'U4(2)': 'U4(2)'}
    return rep.get(g, g)

def char_p(cap):
    m = re.search(r'p = (\d+)\$', cap)
    if m and 'p = 0' not in cap:
        return int(m.group(1))
    return 0

def parse_block(block, has_flag, context=''):
    """One tabular body: returns (groups [(name, count)], labels, rows {n: rec})"""
    groups, labels, rows = None, None, {}
    for l in block.split('\n'):
        l2 = l.replace('\\hline', '').replace('\\endhead', '').replace('\\endfirsthead', '').strip()
        l2 = re.sub(r'\\\\\s*$', '', l2).strip()
        if not l2 or l2.startswith('\\label') or l2.startswith('%') or l2.startswith('\\caption'):
            continue
        if 'multicolumn' in l2 and groups is None:
            groups = [(re.sub(r'[\$\{\}_\\]', '', name).replace('V', 'V'), int(k))
                      for k, name in re.findall(r'\\multicolumn\{(\d+)\}\{[^}]*\}\{([^}]*\}?\$?)\}', l2)]
            groups = [('V' + re.sub(r'\D', '', g), k) for g, k in groups]
            continue
        cells = [re.sub(r'\\+$', '', c.strip()).strip() for c in l2.split('&')]   # stray trailing backslashes
        if re.match(r'^\d+\)$', cells[0]):
            if sum(k for _, k in groups) != len(labels):    # \multicolumn count typo: trust labels
                FIXES.append('%d labels but \\multicolumn counts %s' % (len(labels), groups))
                groups[-1] = (groups[-1][0], groups[-1][1] + len(labels) - sum(k for _, k in groups))
            ncols = len(labels)
            i = 0
            while i < len(cells):
                rm = re.match(r'^(\d+)\)$', cells[i])
                if not rm:
                    raise ValueError('bad row: %s' % l2)
                j = i + 1
                flags = ''
                if has_flag and j < len(cells) and (cells[j] == '' or '\\' in cells[j]):
                    flags = cells[j]; j += 1
                k = j
                while k < len(cells) and not re.match(r'^\d+\)$', cells[k]):
                    k += 1
                svals = [re.sub(r'[${} ]', '', v) for v in cells[j:k]]
                svals = [v for v in svals if v != '']
                if len(svals) != ncols:
                    raise ValueError('row length %d != %d: %s' % (len(svals), ncols, l2))
                n = int(rm.group(1))
                if all(re.match(r'^-?\d+$', v) for v in svals):
                    rows[n] = {'P': '\\possprim' in flags, 'N': '\\nongcr' in flags, 'vals': [int(v) for v in svals]}
                else:   # parametrised row, "where $r = a,\ldots,b$"
                    pm = re.search(r'where \$r = (\d+),\s*\\ldots,\s*(\d+)\$', context)
                    FIXES.append('row %d expanded for r = %s..%s' % (n, pm.group(1), pm.group(2)))
                    for r in range(int(pm.group(1)), int(pm.group(2)) + 1):
                        rows[n * 1000 + r] = {'P': '\\possprim' in flags, 'N': '\\nongcr' in flags,
                                              'vals': [int(eval(v, {'r': r})) for v in svals]}
                i = k
        elif labels is None and groups is not None and cells[0] == '':
            labels = [clean_label(c) for c in cells if c != '']
    return groups, labels, rows

ENV = r'\\begin\{(table|longtable)\}(.*?)\\end\{\1\}'
tables = []
for typ, body in sec.items():
    for m in re.finditer(ENV, body, re.S):
        kind, block = m.group(1), m.group(2)
        cap = re.search(r'\\caption\{(.*?)\}\s*(\\\\)?\s*\n', block)
        if not cap:
            continue
        caption = cap.group(1).strip()
        if kind == 'table':
            blocks = re.findall(r'\\begin\{tabular\}\{(.*?)\}\s*\n(.*?)\\end\{tabular\}', block, re.S)
        else:
            spec = re.match(r'\{(.*?)\}\s*\n', block).group(1)
            blocks = [(spec, block[cap.end():])]
        mods, rows = [], {}
        FIXES = []
        for spec, rows_tex in blocks:
            try:
                groups, labels, brows = parse_block(rows_tex, spec.startswith('rc'), block)
            except Exception as ex:
                raise ValueError('%s: %s' % (caption, ex))
            pos = 0
            gstart = len(mods)
            for name, k in groups:
                mods.append({'name': name, 'labels': labels[pos:pos + k]})
                pos += k
            if pos != len(labels):
                raise ValueError('labels/groups mismatch in %s' % caption)
            for n, r in brows.items():
                rr = rows.setdefault(n, {'n': n, 'P': False, 'N': False, 'mult': [None] * 0})
                rr['P'] = rr['P'] or r['P']; rr['N'] = rr['N'] or r['N']
                off = 0
                for gi, (name, k) in enumerate(groups):
                    rr['mult'].append(r['vals'][off:off + k]); off += k
        if len(mods) != 2:
            raise ValueError('expected 2 column groups in %s: %s' % (caption, mods))
        # implicit duals
        dm = re.search(r'duals of (.*?) also occur', block)
        dm2 = re.search(r'\n([^\n]*?) occur with the same multiplicities as their duals', block)
        notes = ['parser fix: ' + f for f in FIXES]
        if dm2 and not dm:   # "16_a^*, 16_b^* and 26_b^* occur ...": the starred ones are implicit
            dm = dm2
        if dm:
            labs = [clean_label(x).rstrip('*') for x in re.split(r',|\band\b', dm.group(1)) if x.strip()]
            notes.append('duals added: ' + ', '.join(labs))
            for gi, md in enumerate(mods):
                for lab in labs:
                    if lab in md['labels']:
                        idx = md['labels'].index(lab)
                        md['labels'].append(lab + '*')
                        for r in rows.values():
                            r['mult'][gi].append(r['mult'][gi][idx])
        after = body[m.end():m.end() + 300]
        perm = re.search(r'Permutations?:(.*?)\n', block + '\n' + after.split('\\begin')[0])
        tables.append({'type': typ, 'caption': caption, 'group': gap_name(caption), 'p': char_p(caption),
                       'mods': mods, 'rows': [rows[n] for n in sorted(rows)],
                       'perms': perm.group(1).strip() if perm else '', 'notes': notes})

# one-line tables: irreducible on both modules
for typ, grp, p, cap, a, v in [
        ('F4', '3D4(2)', 0, '$^{3}D_4(2) < F_4$, $p \\neq 2, 3$', 52, 26),
        ('E6', '3.Fi22', 2, '$3 \\cdot Fi_{22} < E_6$, $p = 2$', 78, 27),
        ('E6', '3.O7(3)', 2, '$3 \\cdot \\Omega_7(3) < E_6$, $p = 2$', 78, 27),
        ('E6', '3.G2(3)', 2, '$3 \\cdot G_2(3) < E_6$, $p = 2$', 78, 27),
        ('E6', "2F4(2)'", 0, "$^{2}F_4(2)' < E_6$, $p \\neq 2$, $3$", 78, 27)]:
    tables.append({'type': typ, 'caption': cap + ' (irreducible on both)', 'group': grp, 'p': p,
                   'mods': [{'name': 'V%d' % a, 'labels': [str(a)]}, {'name': 'V%d' % v, 'labels': [str(v)]}],
                   'rows': [{'n': 1, 'P': True, 'N': False, 'mult': [[1], [1]]}], 'perms': '', 'notes': ['hand-entered']})

def gstr(s):
    return '"%s"' % s.replace('\\', '\\\\').replace('"', '\\"')

out = ['TABLESG := [];']
for t in tables:
    mods = ', '.join('rec(name := %s, labels := [%s])' % (gstr(md['name']), ', '.join(gstr(l) for l in md['labels']))
                     for md in t['mods'])
    rows = ',\n    '.join('rec(n := %d, P := %s, N := %s, mult := %s)' % (
        r['n'], 'true' if r['P'] else 'false', 'true' if r['N'] else 'false', r['mult']) for r in t['rows'])
    out.append('Add(TABLESG, rec(type := %s, caption := %s, group := %s, p := %d,\n  mods := [%s],\n  perms := %s, notes := [%s],\n  rows := [\n    %s ]));'
               % (gstr(t['type']), gstr(t['caption']), gstr(t['group']), t['p'], mods, gstr(t['perms']),
                  ', '.join(gstr(x) for x in t['notes']), rows))
io.open(sys.argv[2], 'w', encoding='ascii').write('\n'.join(out) + '\n')
if len(sys.argv) > 3:
    io.open(sys.argv[3], 'w', encoding='utf-8').write(json.dumps(tables, indent=1))
from collections import Counter
print(len(tables), 'tables', dict(Counter(t['type'] for t in tables)), sum(len(t['rows']) for t in tables), 'rows')
print('groups:', sorted(set(t['group'] for t in tables)))
