# j-points

A program that finds rational points on the Jacobian of a curve of genus 2,
given as y^2 = f(x) with f of degree 5 or 6 with integral coefficients, up to
a bound on the naive height of their images on the Kummer surface. It sieves
the triples of the first three Kummer coordinates modulo small primes, in the
manner of [ratpoints](https://github.com/MichaelStollBayreuth/ratpoints), from
which it descends, and lifts the survivors to the Jacobian.

The program needs the GNU gmp library. Build and test it with

    make j-points
    make test

and install the executable in /usr/local/bin with `make install-bin`.
The default build is portable; `make VECTOR=avx2 j-points` (or `sse2`)
builds the sieve's passes over the bit array with vectors of four (two)
words, which can be faster on a machine with that instruction set (the
change log has the measurements) and does not run on one without.
`make test` runs seven suites and fails if any of them does: test1, one curve
with many points at height 2000; test2, 61 curves of every kind (testcurves2)
at heights up to 4500, whose reference is the program's output checked by an
unsieved run of every curve (verify-test2.sh, an hour of CPU time); test3,
the options, messages and errors, checked likewise (verify-test3.py); test4,
the same 61 curves at height 30, with and without `-a`, against a brute force
over every coprime triple that uses Magma's Kummer package for the points
above a triple and their lifts (mkref4.m); and testbrute, the sieve against
the unsieved search at height 50; and testrich, fourteen searches with `-a`
on point-rich curves (testcurves-rich: from the ratpoints suites, and the
curve with 642 known rational points of Müller and Stoll, ANT 10 (2016)),
the regime of the enumeration of points of bounded canonical height, whose
reference was checked by the unsieved run as well. Magma is only needed to
regenerate testbase4; its own search for points is a port of j-points and
is not used as a reference.

Please acknowledge use of the program in published work.

## Usage

    j-points 'a_0 a_1 ... a_d' h [-n num_primes1] [-M num_primes2] [-N num_primes3]
             [-p num_primes] [-s size] [-t threads] [-f format] [-w bound4]
             [-1] [-q] [-a]

where

+ f(x) = a_d x^d + ... + a_1 x + a_0 with d = 5 or 6 (zero leading
  coefficients are dropped); f must be squarefree, which the program checks;
+ h is the bound on the naive height: the program finds the points
  (a : b : c : d) on the Kummer surface, with coprime integral coordinates,
  that lift to the Jacobian and satisfy |a|, |b|, |c|, |d| <= h. On a machine
  with 64-bit longs, h can be at most 3037000499.

The optional arguments, which can be given in any order, are these.

+ `-n num_primes1`: the number of primes used for the first sieving stage.
  By default it is chosen automatically, like the primes themselves (see
  below). A plain run that sieves in the tube (see below) has no first
  stage, and `-n` does nothing then.
+ `-M num_primes2`: the number of primes used for the first two sieving
  stages together. By default it is chosen automatically.
+ `-N num_primes3`: the number of primes used for all three sieving stages
  together. By default it is chosen automatically. With `-N 0` there is no
  sieving at all: every coprime triple is checked exactly, which is slow but
  independent of the sieve.
+ `-p num_primes`: the number of primes, from the beginning of the table
  3, 5, 7, ..., 251, among which the sieving primes are chosen (default and
  maximum 53). Sieve tables, for the first two stages, are made only for the
  primes up to 127; the larger ones can serve the third stage.
+ `-s size`: the size in kilobytes of the bit array in which the sieving is
  done (default 10). This determines the amount of memory the program uses.
  A smaller value can give better performance if the array then fits into the
  cache in its entirety, so the effect of this parameter depends heavily on
  the size and speed of the cache.
+ `-t threads`: the number of threads (default 1). The values of the first
  coordinate are the units of work; the points are printed in the order of
  the search whatever the number of threads, so the output does not change.
  With `-1` the search runs in one thread.
+ `-f format`: a format string for printf to print the points with. It should
  take four long integers, the coordinates on the Kummer surface. The default
  is `"(%ld, %ld, %ld, %ld)\n"`. The escapes `\n`, `\t` and `\\` are
  interpreted (so `-f '%ld %ld %ld %ld\n'` prints one point per line); any
  other backslash is kept. With `-a`, coordinates that exceed the machine
  word are printed by gmp through the same format.
+ `-w bound4`: the bound on the fourth coordinate, in place of h: the program
  finds the points with |a|, |b|, |c| <= h and |d| <= bound4. This is the
  enumeration of the points of bounded modified naive height (Mueller-Stoll,
  Canonical heights on genus-2 Jacobians, Section 17), in which the fourth
  coordinate is divided by the maximum ||f|| of the |a_i|: the points with
  modified height at most log N are exactly those of `-w W` with h = N and
  W = ||f|| N. Together with `-a`, the bound on the first three coordinates
  is that of `-a` and the fourth is bounded by bound4.
+ `-1`: stop as soon as one point has been found.
+ `-q`: suppress all messages other than the points found.
+ `-a`: find all points whose image in P^2 under the first three coordinates
  (a triple of coprime integers) has naive height at most h, with no bound on
  the fourth coordinate. The points are printed with all four coordinates
  coprime, so their first three can exceed h; coordinates that exceed the
  machine word are printed through gmp.

## How the program works

Basically the program implements a loop over all triples (a, b, c) of first
Kummer coordinates with |a|, |b|, |c| <= h, where h is the height bound
given. For each of these triples all possible rational fourth coordinates d
are found such that

+ (a : b : c : d) is on the Kummer surface K, and
+ the point (a : b : c : d) in K(Q) lifts to J(Q).

The coordinates are scaled by the denominator of d, so that they are coprime
integers.

This test is preceded by a sieve. First it is checked whether (a, b, c)
mod p are the first three coordinates of a point in K(F_p) that lifts to
J(F_p). This is done for a number of primes p. Only the surviving triples are
then used to compute the possible d's and to check whether the points thus
obtained lift to J(Q).

There are a number of improvements to this basic scheme.

1. We use bits to represent the individual numerators. In this way we can
   sieve as many numerators as bits fit into a long word (usually 64) at the
   same time, using bit-wise "and" operations: for fixed a and b, a word of
   the bit array over c is ANDed with a word of a table that holds the
   admissible c modulo p (the first stage).
2. When this kind of sieving has reduced the candidates considerably, we
   continue the sieving with more primes for each surviving word separately
   (the second stage), and then with more primes for each surviving bit
   separately, by a table of the points of K(F_p) that lift to J(F_p) (the
   third stage, which can use primes too large for the tables of the first
   two stages).
3. We take more primes than necessary for the sieving procedure and determine
   in a first step which are the best ones. This is measured by the ratio of
   numbers surviving the corresponding step of the sieve (essentially the
   number of points mod p, divided by p^2), estimated from a sample of the
   residue classes for the larger primes. The primes of each stage and their
   numbers are then chosen by a cost model of the run: the passes of the
   first stage (a table of a large prime costs more, being large), the tests
   of the surviving words and bits, the exact checks, and the time to set up
   the tables, which matters at small height bounds. `-n`, `-M` and `-N`
   pin the numbers.
4. Rows (a, b) that the sieve would empty outright are skipped: when the
   leading coefficient of f is not a square modulo a sieving prime p, the
   curve has no point at infinity over F_p, and no point of J reduces to a
   Kummer point with a divisible by p unless b is too.
5. The bound on the fourth coordinate d cuts the box of the first three to
   a thin tube: d is a root of a quadratic whose coefficients are
   polynomials in (a, b, c), and a root in [-h, h] exists only close to
   the plane section d = 0 of the Kummer surface (or near the origin). A
   plain run (not `-a`, whose points have no such bound) therefore does
   not sieve whole rows: for each a, the plane of (b, c) is cut into cells
   of rows times words, a cell is dropped when rigorous bounds show that no
   root in [-h, h] exists in it, and only the words of the cells left are
   sieved, one by one. The gain depends on the curve and grows with the
   height bound (from 1.3x to 14x at 2000 on the test curves, 4.5x to 20x
   at 10000); a cost model decides between the tube and the sieve of whole
   rows, which stays better where the tube is fat.
6. The prime 2 contributes through a table modulo 64: a point's integer
   coordinates satisfy the Kummer equation, and the quantity A^2 of the
   lifting test is a square, so a class of (a, b, c) modulo 64 must admit
   some fourth coordinate (of odd or even denominator) with both. The
   admitted c of a class of (a, b) form the word that fills the bit array
   of the row, and rows whose class admits nothing are not sieved.

## Examples

User times on an Intel Core i7-6600U at 2.60 GHz under Linux; the last column
gives the numbers of primes used in the first stage and in both stages of the
sieve together.

| command                                   | time   | points | primes |
|-------------------------------------------|-------:|-------:|-------:|
| `j-points '1 6 5 22 22 8 1' 200`          | 0.08 s |     11 |   5+9  |
| `j-points '1 6 5 22 22 8 1' 500`          | 0.35 s |     12 |   5+9  |
| `j-points '1 6 5 22 22 8 1' 1000`         | 1.8 s  |     15 |   5+9  |
| `j-points '1 178 817 -274 16 1' 200`      | 0.05 s |     35 |   6+4  |
| `j-points '1 178 817 -274 16 1' 500`      | 0.06 s |     63 |   6+4  |
| `j-points '1 178 817 -274 16 1' 1000`     | 0.1 s  |    101 |   6+4  |
| `j-points '1 178 817 -274 16 1' 2000`     | 0.36 s |    124 |   6+4  |
| `j-points '21 116 171 128 55 12 1' 200`   | 0.07 s |     15 |   5+6  |
| `j-points '21 116 171 128 55 12 1' 500`   | 0.4 s  |     19 |   5+6  |
| `j-points '21 116 171 128 55 12 1' 1000`  | 2.2 s  |     23 |   5+6  |

## Change log

+ 25-Oct-1998: started writing version 0.1, taking ratpoints-1.2 as a
  starting point.
+ 12-Nov-1998:
  + Included some enhancements taken from ratpoints-1.4.
  + Avoid dealing with multiples of reduced triples.
  + This allows taking only squares for the first coordinate, when f is
    monic of degree 5.
+ 24-Nov-1998: eliminated a bug (the program didn't take into account the
  possibility that the Kummer equation gives a constant).
+ 23-Apr-2001, version 0.6: changed `!a&1` into `!(a&1)` (no real bug, but
  it prevented an optimisation from being used).
+ 02-May-2001, version 1.0:
  + Split the code into two files: j-points.c and j-sift.c.
  + Hand-optimised assembler code for j-sift for the i386 architecture.
  + Makefile.
  + A few minor bugfixes.
+ 08-Aug-2006, version 1.1: removed a call `mpz_mul_ui(&fff, &fff,
  (unsigned long)b)` in check_lifts that let the program miss points of the
  form (0 : b : c : d) when b was not a square (it was computing b^7 f(c/b)
  instead of b^6 f(c/b)).
+ 14-Apr-2016, version 1.2: added code to handle the case when `-a` is
  specified and the coordinates of the resulting point(s) exceed machine size.
+ 17-Jan-2022, version 2.0:
  + The table of primes extended to the 30 primes below 128 (it had the
    17 primes below 64), and the option `-p` to limit the number of primes
    considered.
  + The sieve tests coprimality by a binary gcd instead of the Euclidean one.
  + kummer_init rewritten.
  + A Makefile with the targets test, dist and install-bin, and the reference
    output testbase for the test.
+ 24-Sep-2026, version 2.1: an audit of the code and a test suite. The
  mathematics and the sieve were found sound; the points of 2.0 are
  unchanged wherever it printed the right ones. Fixed:
  + `-N 0` alone silently lost every point of the main search (the first
    stage was given one prime, the tables none); it is now an unsieved
    search.
  + A zero leading coefficient made the program print the origin as a point;
    zero leading coefficients are now dropped.
  + `-p` given after `-n` or `-N` left the first stage with more primes than
    the table had.
  + With `-a`, a `-f` format other than the default broke the printing of
    coordinates beyond the machine word; gmp prints them now.
  + Heights above 2^31 overflowed a difference in the lifting test.
  + The usage message names every option; the version banner said 1.2.
+ 25-Sep-2026, version 3.0: the sieve made faster, in steps; the points
  printed are the same.
  + The bookkeeping of a row without integer divisions: what the sieve
    needs to know about a chunk of the bit array (where each prime's table
    row starts in it, the residue of its first word) is computed once at
    set-up instead of per row, and the residues of the first two
    coordinates are kept and advanced from one row to the next. 1.2-1.4x.
  + The sieve tables allocated per prime in use, with rows of p+1 words,
    instead of a static array of 480 MB (of which a run touched about 100
    MB; now 15-20 MB), the row pointers advanced from row to row, and the
    fill of the bit array merged into the first prime's pass. Another
    1.1-1.15x.
  + The rows (a, b) that a sieving prime p empties outright are no longer
    sieved: when the leading coefficient is not a square modulo p, the
    curve has no point at infinity over F_p, so no point reduces to a
    Kummer point with first coordinate 0 and second coordinate nonzero
    modulo p; for an a divisible by such primes only the b divisible by
    them are searched. Nothing on curves with a square leading
    coefficient, 1.2x on the record curve (a quarter of its rows).
  + The table of primes extended to the 53 primes below 256 (it had the 30
    below 128), a third sieving stage per surviving bit for the primes
    without a sieve table, and the primes of each stage and their numbers
    chosen by a cost model of the run instead of the two ratios, whose
    options `-r` and `-R` are gone; `-M` pins the number of primes of the
    first two stages, `-N` now counts all three. The rates of the larger
    primes are estimated from a sample of the residue classes, and the
    tables are made for the primes chosen only, so that the set-up is
    cheaper than before at small height bounds. 1.6-1.8x on the record
    curve (which had been starved of primes), 1.0-1.03x on random curves.
  + The tube: a plain run analyses, for each a, the plane of (b, c) in
    cells and sieves only the words in which a fourth coordinate of height
    at most h can exist (see "How the program works"); rigorous bounds, so
    the points are the same. 1.3x to 14x at 2000 on the test curves, more
    at larger height bounds; nothing for `-a`.
  + The passes of the sieve over the bit array can be built with vectors of
    two or four words (`make VECTOR=sse2` or `avx2`), behind a compile-time
    switch: 1.2x on `-a` runs with AVX2 on a recent laptop; SSE2 gave 1.1x
    there and lost 0.1-0.2x on an older desktop (Ivy Bridge), so measure
    before choosing it.
  + A condition at 2: a point with coprime integer coordinates satisfies
    the Kummer equation and has a certain square among its coordinates
    (the quantity A^2 of the lifting test), so a class of (a, b, c) modulo
    64 can carry a point only if some fourth coordinate, of odd or of even
    denominator, satisfies both modulo 64. The admitted c of a class of
    (a, b) form one word, which fills the bit array of the row, and rows
    whose class admits nothing are skipped. Measured at height bounds of
    500 and more, and only when it admits less than everything modulo 16.
    1.5x on `-a` runs of small-coefficient curves, 1.1x on #28, nothing on
    the record curve (f is a square modulo 4 there).
  + A sextic with f0 = 0 and f1 = 1 or -1 is searched as the monic quintic
    it becomes under x -> 1/x (and x -> -x), a quintic with leading
    coefficient -1 under x -> -x: the first coordinate is a square then, so
    the search runs over sqrt(h) values of it instead of h. The Kummer
    coordinates go to (c : b : a : d) and (a : -b : c : d), which keeps the
    height, so the points of the curve given are exactly those found,
    transformed back. 17x on such a curve at 2000.
  + The polynomial is checked to be squarefree. The escapes `\n`, `\t` and
    `\\` in the format of `-f` are interpreted. `-w bound4` bounds the
    fourth coordinate separately from the first three, which gives the
    enumeration of the points of bounded modified naive height without a
    post-filter.
  + Threads (`-t`): the values of the first coordinate are handed out to the
    threads in order, each sieves its rows with its own state, and the
    points of each value are printed in order, so the output is that of a
    single thread. Four times faster with eight threads on the record curve
    on a laptop with two performance and eight efficiency cores, 2.6-2.8x
    with three threads on a four-core desktop, plain runs included.

Michael Stoll, October 1998 - September 2026.

Copyright (C) 1998-2026 Michael Stoll. Distributed under the GNU GPL, version 2
or (at your option) any later version; see gpl-2.0.txt.
