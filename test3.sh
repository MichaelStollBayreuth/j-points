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
# 51     -a on a curve with real roots of f near 10^10 and 2 10^10, between
#        which f is negative: before 3.0.1 the real region cast the far ends
#        of the excluded intervals of c to longs and dropped whole rows, six
#        of the eleven points
# 52     -a with the cost constants set so that the real region is sieved, on
#        lambda (x - 10^17)(x^4 - 1): before 3.0.1 the Cauchy bound, in
#        doubles, fell below the root 10^17, the signs of f were taken the
#        wrong way round, and three points were lost
# 53     -t 12, the points unsorted: the point over (0, 0, 1) and the double
#        points first, as with one thread (before 3.0.1 they could come
#        after the points of later units)
# 54     -t far above the number of units of work, which bounds the threads
# 55     -w with a bound that is not a long: an error (before 3.0.1 it was
#        taken as 2^63 - 1)
# 56     -s with a size whose conversion to words overflowed (a crash before
#        3.0.1): the size is cut to the bit array
# 57     -f with flags, a width and %% in the gmp path of -a (before 3.0.1
#        only a bare %ld was passed to gmp, and a width printed a pointer)
# 58     a height bound above 2^26: an error (before 3.0.2 heights up to
#        3037000499 were accepted, at which the real region could lose points)
# 59-60  -1 with the first point found before the units of work, in the
#        report of a run that is not -q: the point over (0, 0, 1), a double
#        point (before 3.0.2 the report printed it and then "Found 0
#        rational points")
# 61-62  -N 0, no sieving at all, in the report of a run that is not -q:
#        neither a mode of the sieve nor the condition at 2 is announced, as
#        the box is searched (before 3.0.2 the report announced the mode the
#        cost model chose, here the tube with the condition at 2 and the tube
#        cut to the real region); 61 stops at the point over (0, 0, 1) (-1)
# 63-64  -f with flags that gmp_printf reads otherwise than printf (+ with a
#        space, 0 with a precision, - with 0, a lone . as the precision) in
#        the gmp path of -a, with one thread (gmp_printf) and with two
#        (gmp_asprintf): every point in the style of printf (before 3.0.2 the
#        point beyond a machine word printed as [ 0|00000000|10000|...], its
#        third coordinate 1 padded with zeros on the right)
# 65     the report of a run on a curve where the modulus 16 admits every
#        class of the condition at 2 and the modulus 64 only 6/7 of them:
#        the condition at 2 is used (before 3.0.2 it was given up whenever
#        the modulus 16 admitted every class)
# 66-67  an option word that is not two characters long: an error (before
#        3.0.2 an option was recognised by its second character, so that
#        -quiet was -q and -n5 7 was -n 7)
# 68-70  -n, -M and -N with a negative number: an error (before 3.0.2 a
#        negative number of primes was taken as 0)
# 71-75  -c with a value that is infinite, overflows or is subnormal, with
#        the empty string and with a final comma: an error (before 3.0.2 each
#        was accepted, 1e-320 as 9.99989e-321)
# The reports (39, 45, 46, 50, 59-62, 65) pin the four constants of the cost
# model whose defaults depend on the width of the sieve's passes (VECTOR) to
# the scalar build's values with -c, since they name the primes chosen: the
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
# an invocation that is expected to fail: no $JPOPTS, output unsorted, its
# exit status not that of the script
e() { printf '#'; for a; do printf ' <%s>' "$a"; done; echo; "$JP" "$@"; return 0; }
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
t '0 400000000000000000000 -200000000060000000000 -399999999969999999998 200000000059999999999 -30000000002 1' 10 -q -a
t '-15999999989999968000000020000000000000000000 159999999899999680000000200 0 0 15999999989999968000000020000000000000000000 -159999999899999680000000200' 256 -q -a -c ROW=10000000,RROW=1,RANGE=1,PASS=1,REXCL=1
f cat '21 116 171 128 55 12 1' 300 -q -t 12
t '21 116 171 128 55 12 1' 20 -q -t 1000000000
e '21 116 171 128 55 12 1' 100 -w 1000000000000000000000000000000
t '21 116 171 128 55 12 1' 100 -q -s 72057594037927936
f 'cat; echo' '0 400000000000000000000 -200000000060000000000 -399999999969999999998 200000000059999999999 -30000000002 1' 10 -q -a -f '[%24ld|%-3ld|%+ld|%ld] %%;'
e '1 2 3 4 5 6 7' 67108865
f "grep -v '^This is j-points'" '1 0 0 0 0 0 1' 10 -1 -c AND=0.63,SIZE=756,TEST2=5.6,TTEST=2
f "grep -v '^This is j-points'" '1 0 0 0 0 0 2' 10 -1 -c AND=0.63,SIZE=756,TEST2=5.6,TTEST=2
f "grep -v '^This is j-points'" '1 0 0 0 0 0 1' 500 -N 0 -1 -c AND=0.63,SIZE=756,TEST2=5.6,TTEST=2
f "grep -v '^This is j-points'" '21 116 171 128 55 12 1' 60 -n 0 -N 0 -c AND=0.63,SIZE=756,TEST2=5.6,TTEST=2
f 'cat; echo' '1 1 1 1 1 30000000000 1' 20 -a -q -t 1 -f '[%+ ld|%08.3ld|%-05ld|%.ld]'
f 'cat; echo' '1 1 1 1 1 30000000000 1' 20 -a -q -t 2 -f '[%+ ld|%08.3ld|%-05ld|%.ld]'
f "grep -v '^This is j-points'" '1 6 5 22 22 8 1' 500 -c AND=0.63,SIZE=756,TEST2=5.6,TTEST=2
e '1 2 3 4 5 6 7' 100 -quiet
e '1 2 3 4 5 6 7' 100 -n5 7
e '1 2 3 4 5 6 7' 100 -n -1
e '1 2 3 4 5 6 7' 100 -M -2
e '1 2 3 4 5 6 7' 100 -N -1
e '21 116 171 128 55 12 1' 100 -c AND=inf
e '21 116 171 128 55 12 1' 100 -c AND=1e999
e '21 116 171 128 55 12 1' 100 -c AND=1e-320
e '21 116 171 128 55 12 1' 100 -c ''
e '21 116 171 128 55 12 1' 100 -c AND=0.63,
