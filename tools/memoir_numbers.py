r"""Number the tables of Chapter 6 of the Memoir as LaTeX does (\numberwithin{table}{chapter}).

Each table or longtable environment with a \caption, and each \addtocounter{table}{1} (used for
the one-line "irreducible on ..." entries), advances the counter.

Usage: python3 memoir_numbers.py AJL-thetables.tex [chapter=6] > numbers.json
Output: a list of {"number": "6.N", "caption": ..., "type": "F4"|"E6"|"E7"|"E8", "line": ...}
"""
import sys, re, json

src = open(sys.argv[1], encoding='utf-8').read()
chap = sys.argv[2] if len(sys.argv) > 2 else '6'
tok = re.compile(r'\\section\{\$([^$]*)\$\}|\\begin\{(table|longtable)\}|\\end\{(table|longtable)\}'
                 r'|\\caption\{|\\addtocounter\{table\}\{1\}')
out, n, typ, inenv, captioned, envstart = [], 0, None, None, False, 0
for m in tok.finditer(src):
    line = src.count('\n', 0, m.start()) + 1
    if m.group(1):
        typ = m.group(1).replace('_', '').replace('{', '').replace('}', '')
    elif m.group(2):
        inenv, captioned, envstart = m.group(2), False, m.start()
    elif m.group(3):
        inenv = None
    elif m.group(0).startswith('\\caption') and inenv and not captioned:
        captioned = True
        n += 1
        cap = re.match(r'\\caption\{(.*?)\}\s*(\\\\)?\s*\n', src[m.start():]).group(1)
        out.append({'number': '%s.%d' % (chap, n), 'caption': cap.strip(), 'type': typ, 'line': line})
    elif m.group(0).startswith('\\addtocounter'):
        n += 1
        rest = src[m.end():m.end() + 400]
        cap = re.search(r'\\arabic\{table\}:\s*(.*?)\\end\{center\}|\\thetable:\s*(.*?)\\end\{center\}', rest, re.S)
        text = (cap.group(1) or cap.group(2)).strip() if cap else '(one-line entry)'
        out.append({'number': '%s.%d' % (chap, n), 'caption': text, 'type': typ, 'line': line})
json.dump(out, sys.stdout, indent=1)
sys.stderr.write('%d tables numbered; last %s\n' % (len(out), out[-1]['number']))
for k in out:
    if k['number'] == '%s.298' % chap:
        sys.stderr.write('6.298: %s (%s, line %d)\n' % (k['caption'], k['type'], k['line']))
