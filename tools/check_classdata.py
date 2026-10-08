"""Validation (a): compare the GAP class data of feasgen.g (classdata.py, from fg_dump.g) with
FeasChar's precomputed element data n.M (n = 2..17), as sets of power tuples.

Usage: python3 check_classdata.py <dir with n.M files> <classdata.py>

Each n.M file defines <TYPE>_ELTS<n> as a Magma set of [V_min tuple, adjoint tuple] (adjoint
tuple only for E8), each tuple the traces of t^(n/j), j | n, j > 1.  Magma cyclotomics are
written as [ CyclotomicField(k) | [ RationalField() | a_0, ..., a_{phi(k)-1} ] ], i.e.
coefficients on the power basis of zeta_k.
"""
import sys, re, cmath, os

def cyc(k, coeffs):
    return sum(a * cmath.exp(2j * cmath.pi * i / k) for i, a in enumerate(coeffs))

def parse_seq(s, i):
    """Parse a Magma literal '[ Type | e1, e2, ... ]' or '{ Type | ... }' starting at s[i].
    Returns (value, next index).  Sequences over RationalField() are coefficient lists; a
    sequence over CyclotomicField(k) has coefficient-list elements, converted to complex."""
    assert s[i] in '[{'
    close = ']' if s[i] == '[' else '}'
    bar = s.index('|', i)
    typ = s[i + 1:bar].strip()
    i = bar + 1
    elems = []
    while True:
        while s[i] in ' \n\t\r,':
            i += 1
        if s[i] == close:
            i += 1
            break
        if s[i] in '[{':
            v, i = parse_seq(s, i)
        else:
            m = re.match(r'-?\d+(/\d+)?', s[i:])
            tok = m.group(0)
            v = (int(tok.split('/')[0]) / int(tok.split('/')[1])) if '/' in tok else int(tok)
            i += len(tok)
        elems.append(v)
    mk = re.match(r'CyclotomicField\((\d+)\)$', typ)
    if mk:
        k = int(mk.group(1))
        elems = [cyc(k, co) for co in elems]
    return elems, i

def parse_magma(txt):
    out = {}
    for m in re.finditer(r'(\w+)_ELTS(\d+)\s*:=\s*', txt):
        typ, n = m.group(1), int(m.group(2))
        if typ == 'E8_3875':
            continue
        val, _ = parse_seq(txt, m.end())
        out[(typ, n)] = val
    return out

def key(z):
    z = complex(z)
    return (round(z.real, 4) + 0.0, round(z.imag, 4) + 0.0)

def canon_magma(elt, typ):
    # elt = [vmin tuple, adjoint tuple], or the adjoint tuple alone (E8)
    if typ == 'E8':
        elt = [elt]
    return tuple(tuple(key(z) for z in tup) for tup in elt)

def canon_gap(elt):
    return tuple(tuple(key(complex(*z)) for z in tup) for tup in elt)

ns = {}
exec(open(sys.argv[2]).read().replace(chr(92) + chr(10), ''), ns)   # GAP's line continuations
gap = ns['data']
allok = True
for n in range(2, 18):
    path = os.path.join(sys.argv[1], '%d.M' % n)
    mag = parse_magma(open(path).read())
    for typ in ['G2', 'F4', 'E6', 'E7', 'E8']:
        if (typ, n) not in mag or (typ, n) not in gap:
            print('n = %2d %s: missing (magma %s, gap %s)' % (n, typ, (typ, n) in mag, (typ, n) in gap))
            allok = False
            continue
        M = set(canon_magma(e, typ) for e in mag[(typ, n)])
        G = set(canon_gap(e) for e in gap[(typ, n)])
        # values were rounded twice (GAP to 6 digits, here to 4): match the leftovers with a tolerance
        flat = lambda t: [x for tup in t for z in tup for x in z]
        close = lambda a, b: len(a) == len(b) and all(abs(x - y) < 1e-3 for x, y in zip(flat(a), flat(b)))
        nM, nG = len(M), len(G)
        onlyM, onlyG = list(M - G), list(G - M)
        for a in list(onlyM):
            b = next((b for b in onlyG if close(a, b)), None)
            if b is not None:
                onlyM.remove(a); onlyG.remove(b)
        M, G = set(onlyM), set(onlyG)
        ok = (M == G)
        allok &= ok
        print('n = %2d %s: Magma %4d distinct tuples, GAP %4d distinct (%4d elements): %s'
              % (n, typ, nM, nG, len(gap[(typ, n)]), 'AGREE' if ok else 'DIFFER'))
        if not ok:
            print('   only Magma:', sorted(M - G)[:3])
            print('   only GAP:  ', sorted(G - M)[:3])
print('ALL AGREE' if allok else 'MISMATCH')
