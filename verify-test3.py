#!/usr/bin/env python3
"""Checks the point lists of testbase3 (or of the file named on the command
line) independently of the sieve.

The file is what test3.sh prints: before each invocation a line
"# <arg> <arg> ...", then the output.  For every invocation that prints
points this script computes the set of points it should print:

 - up to the height bound 60, by brute force with Magma's Kummer package: for
   every coprime triple (a, b, c) up to the bound, Points(K, [a, b, c]) gives
   the points of the Kummer surface above it and Points(J, P) says whether
   they lift to the Jacobian.  This shares nothing with j-points (Magma's own
   search for points, Points(J : Bound), is a port of j-points and is not
   used).  With the bound on all four coordinates the points beyond it are
   dropped; with -a, which bounds only the triple, every point counts, printed
   with all four coordinates coprime, so its first three can exceed the bound;
 - above 60, where the brute force is too slow, by the unsieved run of
   j-points itself (-n 0 -N 0: every coprime triple checked exactly), which
   is independent of the sieve but not of the exact check -- that is what
   test4 compares with Magma.

A custom -f format is turned into a pattern to read the points back; with
-1 the printed point must be one of the set, and there must be one exactly
when the set is not empty.  Invocations that print no list (the errors, the
report of a run without -q) are not checked and are listed as such.

One Magma script is written (verify-test3.m) and run once; the unsieved runs
take about 25 minutes, most of it the invocation at height 4500.  Prints one
line per invocation and exits 1 if any check fails.
"""
import re, subprocess, sys

SMALL = 60
JP = './j-points'
ref = sys.argv[1] if len(sys.argv) > 1 else 'testbase3'
text = open(ref).read().split('\n')

blocks = []
for line in text:
    if line.startswith('#'):
        blocks.append({'header': line, 'args': re.findall(r'<([^>]*)>', line), 'out': []})
    elif blocks:
        blocks[-1]['out'].append(line)

magma_jobs = {}    # (coeffs, h) -> job number: brute force, all points above the triples
for b in blocks:
    a = b['args']
    b['check'] = False
    if len(a) < 2 or not re.fullmatch(r'-?\d+( -?\d+)*', a[0]) or not a[1].isdigit():
        b['why'] = 'not a point list (error or usage)'; continue
    coeffs = [int(x) for x in a[0].split()]
    while len(coeffs) > 1 and coeffs[-1] == 0: coeffs.pop()
    if len(coeffs) < 6:
        b['why'] = 'not a point list (degree below 5)'; continue
    opts = a[2:]
    if '-q' not in opts:
        b['why'] = 'no -q: a report or an error, not checked'; continue
    b['coeffs'] = coeffs; b['h'] = int(a[1])
    b['all'] = '-a' in opts; b['one'] = '-1' in opts
    b['w'] = int(opts[opts.index('-w') + 1]) if '-w' in opts else None
    fmt = opts[opts.index('-f') + 1] if '-f' in opts else None
    pat = re.compile(r'\((-?\d+), (-?\d+), (-?\d+), (-?\d+)\)' if fmt is None
                     else re.escape(fmt).replace(r'%ld', r'(-?\d+)'))
    found = pat.findall('\n'.join(b['out']))
    b['found'] = set('(%s, %s, %s, %s)' % tuple(p) for p in found)
    b['nfound'] = len(found)
    b['check'] = True
    if b['h'] <= SMALL:
        b['how'] = 'Magma'
        magma_jobs.setdefault((tuple(coeffs), b['h']), len(magma_jobs))
    else:
        b['how'] = 'unsieved'

# the brute force in Magma, once per (curve, bound)
results = {}
if magma_jobs:
    m = ['P<x> := PolynomialRing(Rationals());',
         'function normalise(s)',
         '  d := LCM([Denominator(c) : c in s]); s := [Integers() | c*d : c in s];',
         '  g := GCD(s); s := [c div g : c in s];',
         '  i := 1; while s[i] eq 0 do i +:= 1; end while;',
         '  if s[i] lt 0 then s := [-c : c in s]; end if; return s;',
         'end function;',
         'function allpoints(f, h)',
         '  J := Jacobian(HyperellipticCurve(f)); K := KummerSurface(J); S := {};',
         '  for a in [0..h] do for b in [-h..h] do for c in [-h..h] do',
         '    if a eq 0 and (b lt 0 or (b eq 0 and c le 0)) then continue; end if;',
         '    if GCD([a, b, c]) ne 1 then continue; end if;',
         '    for P in Points(K, [a, b, c]) do',
         '      if #Points(J, P) eq 0 then continue; end if;',
         '      s := normalise(Eltseq(P));',
         '      Include(~S, Sprintf("(%o, %o, %o, %o)", s[1], s[2], s[3], s[4]));',
         '    end for;',
         '  end for; end for; end for;',
         '  return S;',
         'end function;']
    for (coeffs, h), k in magma_jobs.items():
        m.append('print "#JOB %d";' % k)
        m.append('f := P!%s;' % list(coeffs))
        m.append('for s in allpoints(f, %d) do print s; end for;' % h)
    m.append('quit;')
    open('verify-test3.m', 'w').write('\n'.join(m) + '\n')
    out = subprocess.run(['magma', '-b', 'verify-test3.m'], capture_output=True, text=True)
    if out.returncode != 0 or 'error' in out.stdout.lower():
        print(out.stdout[-2000:]); sys.exit(2)
    cur = None
    for line in out.stdout.split('\n'):
        if line.startswith('#JOB'):
            cur = int(line.split()[1]); results[cur] = set()
        elif line.startswith('(') and cur is not None:
            results[cur].add(line.strip())

def coords(s):
    return [int(x) for x in s[1:-1].split(', ')]

unsieved = {}
def unsieved_set(coeffs, h, allp, w=None):
    key = (tuple(coeffs), h, allp, w)
    if key not in unsieved:
        cmd = ([JP, ' '.join(str(c) for c in coeffs), str(h), '-q', '-n', '0', '-N', '0']
               + (['-a'] if allp else []) + (['-w', str(w)] if w is not None else []))
        out = subprocess.run(cmd, capture_output=True, text=True).stdout
        unsieved[key] = set(l.strip() for l in out.split('\n') if l.startswith('('))
    return unsieved[key]

status = 0
for b in blocks:
    label = b['header'][:70]
    if not b['check']:
        print('  --   %-70s %s' % (label, b['why'])); continue
    if b['how'] == 'Magma':
        want = results[magma_jobs[(tuple(b['coeffs']), b['h'])]]
        # the brute force has every point above the triples up to h: the
        # plain run wants all four coordinates within h, -w the first three
        # within h and the fourth within its bound, -a everything (with -w,
        # the fourth within its bound)
        if not b['all']:
            want = set(s for s in want if max(abs(c) for c in coords(s)[:3]) <= b['h'])
        dmax = b['w'] if b['w'] is not None else (None if b['all'] else b['h'])
        if dmax is not None:
            want = set(s for s in want if abs(coords(s)[3]) <= dmax)
    else:
        want = unsieved_set(b['coeffs'], b['h'], b['all'], b['w'])
    if b['one']:
        ok = (b['nfound'] == min(1, len(want))) and b['found'] <= want
        what = '-1: %d printed, %d exist' % (b['nfound'], len(want))
    else:
        ok = (b['found'] == want) and (b['nfound'] == len(b['found']))
        what = '%d points, %s %d' % (b['nfound'], b['how'], len(want))
        if not ok:
            what += '; missing %s; extra %s' % (sorted(want - b['found']), sorted(b['found'] - want))
    print('  %s %-70s %s' % ('ok  ' if ok else 'FAIL', label, what))
    if not ok: status = 1
sys.exit(status)
