#!/bin/sh
# The options, the messages and the errors of j-points, at small height
# bounds: a list of invocations whose output "make test3" compares with
# testbase3.  Before each invocation the script prints its arguments, one
# per <...>, so that a failing test can be found from the diff of test3.out
# against testbase3, and so that verify-test3.py knows what it is checking.
#
# The reference was checked independently: verify-test3.py reads testbase3
# and, for every invocation whose output is a list of points, compares it
# with Magma -- with Magma's own search for the points when all four
# coordinates are bounded, and with a search over every coprime triple
# (a, b, c) up to the bound that asks Magma for the points of the Kummer
# surface above it and for their lifts to the Jacobian when -a bounds only
# the first three (see the script).
#
#  1-6   -a: no bound on the fourth coordinate, coordinates beyond a machine
#        word printed through gmp (6, with the height bound 60, has them)
#  7-11  the same curves with the bound on all four coordinates
# 12-13  -1: the first point found, or nothing when there is none
# 14-16  -f: the format, also with -a and coordinates beyond a machine word
# 17-24  -n, -N, -p: pinned numbers of primes, no sieving at all (-n 0 -N 0
#        and -N 0 alone, which lost the sieved points before 2.1), a
#        first-stage count above the total, more primes than -p allows, -p
#        after -n
# 25-26  -M: the number of primes of the first two stages pinned, with and
#        without a third stage
# 27     -s 1: a bit array of 1 KB, three chunks at this height bound
# 28-29  a zero leading coefficient is dropped (before 2.1 it printed the
#        origin); a zero polynomial is an error
# 30-38  the errors: no argument, too few or too many coefficients, a height
#        bound of 0 or above the maximum, an unknown option, an option
#        without its value, a number of primes and a size that are not
#        positive
# 39     the report of a run that is not -q, with the primes of the three
#        stages pinned and the line naming the version dropped
# 40     a polynomial that is not squarefree ((x^2 + 1)^2 (x^2 + 2)): an
#        error
# 41     -f with the escapes \t, \n and \\ interpreted
# 42-44  -w: the fourth coordinate bounded separately (its bound below and
#        above the height bound), and together with -a (the point sets are
#        those of the -a runs of the same curve at the same height bound,
#        cut to g m <= h and |d| <= bound4)
# 45     the report of a run with -w, the line naming the version dropped
# 46-48  a sextic with f0 = 0 and f1 = -1, searched as a monic quintic
#        (x -> 1/x and x -> -x) with the points transformed back: the
#        report of the run, the points with -a (the same set as before the
#        transformation, checked against the unsieved run and Magma), and
#        -1 (the first point of the transformed search)
# 49     -c with a name the model does not have: an error
# 50     -c: two constants of the cost model set (the cell of the tube's
#        analysis priced out, so the box is sieved), listed in the report
# The reports (39, 45, 46, 50) pin the four constants of the cost model
# whose defaults depend on the width of the sieve's passes (VECTOR) to the
# scalar build's values with -c, since they name the primes chosen: the
# reference is then the same for every build.
JP=${JP:-./j-points}
# no program to run: exit 2 (the comparison with the reference is make's,
# whose target fails with 1 when they differ)
[ -x "$JP" ] || { echo "$JP: not an executable" >&2; exit 2; }
# options added to every invocation of t and f, not to those of e
JPOPTS=${JPOPTS:-}
# run j-points, the arguments announced first, the points sorted
t() { printf '#'; for a; do printf ' <%s>' "$a"; done; echo; "$JP" "$@" $JPOPTS | LC_ALL=C sort; }
# the same, with the output run through a filter (the first argument)
f() { filter=$1; shift; printf '#'; for a; do printf ' <%s>' "$a"; done; printf ' | %s\n' "$filter"; "$JP" "$@" $JPOPTS | eval "$filter"; }
# an invocation that is expected to fail: no $JPOPTS, output unsorted
e() { printf '#'; for a; do printf ' <%s>' "$a"; done; echo; "$JP" "$@"; }
t '1 178 817 -274 16 1' 20 -q -a
t '21 116 171 128 55 12 1' 20 -q -a
t '3 -1 2 5 -4 1 4' 20 -q -a
t '0 1 2 -1 3 1 1' 20 -q -a
t '1 2 -1 3 1 2' 20 -q -a
t '1 178 817 -274 16 1' 60 -q -a
t '1 178 817 -274 16 1' 20 -q
t '21 116 171 128 55 12 1' 20 -q
t '3 -1 2 5 -4 1 4' 20 -q
t '0 1 2 -1 3 1 1' 20 -q
t '1 2 -1 3 1 2' 20 -q
t '1 178 817 -274 16 1' 200 -q -1
t '2 0 0 0 0 0 3' 300 -q -1
f 'cat; echo' '1 6 5 22 22 8 1' 200 -q -f '[%ld:%ld:%ld:%ld]'
f 'cat; echo' '1 178 817 -274 16 1' 60 -q -a -f '<%ld|%ld|%ld|%ld>'
f 'cat; echo' '1 178 817 -274 16 1' 100 -q -f '%ld,%ld,%ld,%ld;'
t '21 116 171 128 55 12 1' 100 -q -n 0 -N 0
t '21 116 171 128 55 12 1' 100 -q -N 0
t '21 116 171 128 55 12 1' 300 -q -n 0
t '21 116 171 128 55 12 1' 300 -q -n 2
t '21 116 171 128 55 12 1' 300 -q -n 8 -N 5
t '21 116 171 128 55 12 1' 300 -q -n 3 -N 30 -p 30
t '21 116 171 128 55 12 1' 100 -q -p 3
t '21 116 171 128 55 12 1' 300 -q -n 25 -p 10
t '21 116 171 128 55 12 1' 300 -q -M 6
t '21 116 171 128 55 12 1' 300 -q -n 2 -M 4 -N 4
t '1 178 817 -274 16 1' 4500 -q -s 1
t '1 6 5 22 22 8 0' 100 -q
e '0 0 0 0 0 0 0' 100
e
e '1 2 3 4 5' 100
e '1 2 3 4 5 6 7 8' 100
e '1 2 3 4 5 6 7' 0
e '1 2 3 4 5 6 7' 3037000500
e '1 2 3 4 5 6 7' 100 -x
e '1 2 3 4 5 6 7' 100 -n
e '1 2 3 4 5 6 7' 100 -p 0
e '1 2 3 4 5 6 7' 100 -s 0
f "grep -v '^This is j-points'" '21 116 171 128 55 12 1' 100 -n 3 -M 6 -N 9 -c AND=0.63,SIZE=756,TEST2=5.6,TTEST=2
e '2 0 5 0 4 0 1' 100
f 'cat; echo' '1 178 817 -274 16 1' 60 -q -f '<%ld:%ld:%ld:%ld>\t\\\n'
t '21 116 171 128 55 12 1' 60 -q -w 20
t '21 116 171 128 55 12 1' 20 -q -w 500
t '1 178 817 -274 16 1' 60 -q -a -w 100000
f "grep -v '^This is j-points'" '21 116 171 128 55 12 1' 100 -w 50 -n 3 -M 6 -N 9 -c AND=0.63,SIZE=756,TEST2=5.6,TTEST=2
f "grep -v '^This is j-points'" '0 -1 2 -1 3 1 1' 100 -n 3 -M 6 -N 9 -c AND=0.63,SIZE=756,TEST2=5.6,TTEST=2
t '0 -1 2 -1 3 1 1' 60 -q -a
t '0 -1 2 -1 3 1 1' 300 -q -1
e '21 116 171 128 55 12 1' 100 -c FOO=1
f "grep -v '^This is j-points'" '21 116 171 128 55 12 1' 100 -c CELL=1e6,AND=0.45,SIZE=756,TEST2=5.6,TTEST=2
