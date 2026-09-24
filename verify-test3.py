#!/usr/bin/env python3
"""Checks the point lists of testbase3 (or of the file named on the command
line) against Magma, independently of j-points.

The file is what test3.sh prints: before each invocation a line
"# <arg> <arg> ...", then the output.  For every invocation that prints
points this script computes the set of points it should print:

 - with the height bound on all four coordinates, Magma's own search,
   Points(J : Bound := h);
 - with -a, which bounds only the image in P^2 (the coprime triple of the
   first three coordinates), a search over every coprime triple (a, b, c) up
   to the bound, asking Magma for the points of the Kummer surface above it
   (Points(K, [a, b, c])) and for their lifts to the Jacobian (Points(J, P));
   the point is printed with all four coordinates coprime, so its first
   three can exceed the bound.

A custom -f format is turned into a pattern to read the points back; with
-1 the printed point must be one of the set, and there must be one exactly
when the set is not empty.  Invocations that print no list (the errors, the
report of a run without -q) are not checked and are listed as such.

One Magma script is written (verify-test3.m) and run once; Magma must be on
the path.  Prints one line per invocation and exits 1 if any check fails.
"""
import re, subprocess, sys

ref = sys.argv[1] if len(sys.argv) > 1 else 'testbase3'
text = open(ref).read().split('\n')

# split into blocks: header line and output lines
blocks = []
for line in text:
    if line.startswith('#'):
        args = re.findall(r'<([^>]*)>', line)
        blocks.append({'header': line, 'args': args, 'out': []})
    elif blocks:
        blocks[-1]['out'].append(line)

def normalise_points(strs):
    return set(strs)

# what each block needs from Magma
jobs = []
for n, b in enumerate(blocks):
    a = b['args']
    b['check'] = None
    if len(a) < 2 or not re.fullmatch(r'-?\d+( -?\d+)*', a[0]) or not a[1].isdigit():
        b['why'] = 'not a point list (error or usage)'; continue
    coeffs = [int(x) for x in a[0].split()]
    while len(coeffs) > 1 and coeffs[-1] == 0: coeffs.pop()
    if len(coeffs) < 6:
        b['why'] = 'not a point list (degree below 5)'; continue
    h = int(a[1])
    opts = a[2:]
    if '-q' not in opts:
        b['why'] = 'no -q: a report or an error, not checked'; continue
    fmt = None
    if '-f' in opts: fmt = opts[opts.index('-f') + 1]
    b['all'] = '-a' in opts
    b['one'] = '-1' in opts
    # the points printed, read back with the format
    if fmt is None:
        pat = re.compile(r'\((-?\d+), (-?\d+), (-?\d+), (-?\d+)\)')
    else:
        pat = re.compile(re.escape(fmt).replace(r'%ld', r'(-?\d+)'))
    found = pat.findall('\n'.join(b['out']))
    b['found'] = set('(%s, %s, %s, %s)' % tuple(p) for p in found)
    b['nfound'] = len(found)
    b['check'] = len(jobs)
    jobs.append((n, coeffs, h, b['all']))

# the Magma script
m = ['P<x> := PolynomialRing(Rationals());',
     'function normalise(s)',
     '  d := LCM([Denominator(c) : c in s]); s := [Integers() | c*d : c in s];',
     '  g := GCD(s); s := [c div g : c in s];',
     '  i := 1; while s[i] eq 0 do i +:= 1; end while;',
     '  if s[i] lt 0 then s := [-c : c in s]; end if; return s;',
     'end function;',
     'function str(s) return Sprintf("(%o, %o, %o, %o)", s[1], s[2], s[3], s[4]); end function;',
     'function bounded(f, h)',
     '  J := Jacobian(HyperellipticCurve(f)); K := KummerSurface(J); S := {};',
     '  for pt in Points(J : Bound := h) do s := normalise(Eltseq(K!pt));',
     '    if s[1] eq 0 and s[2] eq 0 and s[3] eq 0 then continue; end if;',
     '    Include(~S, str(s)); end for;',
     '  return S;',
     'end function;',
     'function allpoints(f, h)',
     '  J := Jacobian(HyperellipticCurve(f)); K := KummerSurface(J); S := {};',
     '  for a in [0..h] do for b in [-h..h] do for c in [-h..h] do',
     '    if a eq 0 and (b lt 0 or (b eq 0 and c le 0)) then continue; end if;',
     '    if GCD([a, b, c]) ne 1 then continue; end if;',
     '    for P in Points(K, [a, b, c]) do',
     '      s := normalise(Eltseq(P));',
     '      if #Points(J, P) eq 0 then continue; end if;',
     '      Include(~S, str(s));',
     '    end for;',
     '  end for; end for; end for;',
     '  return S;',
     'end function;']
for k, (n, coeffs, h, allp) in enumerate(jobs):
    m.append('print "#JOB %d";' % k)
    m.append('f := P!%s;' % coeffs)
    m.append('for s in %s(f, %d) do print s; end for;' % ('allpoints' if allp else 'bounded', h))
m.append('quit;')
open('verify-test3.m', 'w').write('\n'.join(m) + '\n')

out = subprocess.run(['magma', '-b', 'verify-test3.m'], capture_output=True, text=True)
if out.returncode != 0 or 'error' in out.stdout.lower():
    print(out.stdout[-2000:]); sys.exit(2)
results = {}
cur = None
for line in out.stdout.split('\n'):
    if line.startswith('#JOB'):
        cur = int(line.split()[1]); results[cur] = set()
    elif line.startswith('(') and cur is not None:
        results[cur].add(line.strip())

status = 0
for n, b in enumerate(blocks):
    label = b['header'][:70]
    if b['check'] is None:
        print('  --   %-70s %s' % (label, b['why'])); continue
    want = results[b['check']]
    if b['one']:
        ok = (b['nfound'] == min(1, len(want))) and b['found'] <= want
        what = '-1: %d printed, %d exist' % (b['nfound'], len(want))
    else:
        ok = (b['found'] == want) and (b['nfound'] == len(b['found']))
        what = '%d points, Magma %d' % (b['nfound'], len(want))
        if not ok:
            what += '; missing %s; extra %s' % (sorted(want - b['found']), sorted(b['found'] - want))
    print('  %s %-70s %s' % ('ok  ' if ok else 'FAIL', label, what))
    if not ok: status = 1
sys.exit(status)
