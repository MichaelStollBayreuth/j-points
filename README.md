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
words, which can be faster on a machine with that instruction set
(CHANGE_LOG.md has the measurements) and does not run on one without. The
choice sticks for the following `make`, `make test`, `make tune` and `make
install-bin` in that tree (vector.mk records it) until `make VECTOR=
j-points` returns to the scalar build.
`make tune` measures the constants of the cost model that chooses the
sieving primes, their stages and the mode of the sieve on the machine at
hand and writes them to tuning.mk, which the next `make` compiles in
(tune.sh; about two hours on an idle machine). The defaults were measured
this way on a laptop and a desktop of 2023, for each vector width, and a
tune there moves them by a few per cent at most; an older core may want
more.

`make test` runs seven suites and fails if any of them does: test1, one curve
with many points at height 2000; test2, 61 curves of every kind (testcurves2)
at heights up to 4500, whose reference is the program's output checked by an
unsieved run of every curve (verify-test2.sh, an hour of CPU time); test3,
the options, messages and errors, checked likewise (verify-test3.py); test4,
the same 61 curves at height 30, with and without `-a`, against a brute force
over every coprime triple that uses Magma's Kummer package for the points
above a triple and their lifts (mkref4.m); testbrute, the sieve against the
unsieved search at height 50; testrich, fourteen searches with `-a` on
point-rich curves (testcurves-rich: from the ratpoints suites, and the curve
with 642 known rational points of Müller and Stoll, ANT 10 (2016)), the
regime of the enumeration of points of bounded canonical height, whose
reference was checked by the unsieved run as well; and testthreads, test1,
test2 and testrich run with two and three threads against the same
references. Magma is only needed to regenerate testbase4; its own search for
points is a port of j-points and is not used as a reference.

Please acknowledge use of the program in published work.

## Usage

    j-points 'a_0 a_1 ... a_d' h [-n num_primes1] [-M num_primes2] [-N num_primes3]
             [-p num_primes] [-s size] [-t threads] [-f format] [-w bound4]
             [-c constants] [-1] [-q] [-a]

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
  enumeration of the points of bounded modified naive height (Müller-Stoll,
  Canonical heights on genus-2 Jacobians, Section 17), in which the fourth
  coordinate is divided by the maximum ||f|| of the |a_i|: the points with
  modified height at most log N are exactly those of `-w W` with h = N and
  W = ||f|| N. Together with `-a`, the bound on the first three coordinates
  is that of `-a` and the fourth is bounded by bound4.
+ `-c constants`: constants of the cost model that chooses the sieving
  primes, their stages and the mode of the sieve, as `NAME=value`, several
  separated by commas, the names those of the `COST_` constants of j-points.c
  without the prefix (`-c CELL=640,AND=0.45`); the run's report lists them.
  For experiments with the model and for `make tune`. `j-points -c list`
  alone prints the constants in force in the build.
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

1. We use bits to represent the individual coordinate values. In this way we can
   sieve as many values as bits fit into a long word (usually 64) at the
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
   Independently of any bound, a point needs f to be nonnegative at the
   Mumford roots x, u of a t^2 - b t + c when they are real (the two points
   of the divisor are real then; when they are complex there is nothing to
   ask, and a real fourth coordinate exists in any case). Along a row
   (a, b), x falls and u rises as c grows, so each interval where f is
   negative (between consecutive real roots of f, isolated once by Sturm's
   theorem in exact arithmetic) excludes one interval of c on each branch,
   in closed form. This real region is half to nine tenths of the box: an
   `-a` run sieves only the ranges of words of a row that lie in it, by the
   passes of the whole-row sieve, and a plain run cuts the tube's words to
   it, each where the cost model finds it worth the small cost per row.
6. The prime 2 contributes through a table modulo 64: a point's integer
   coordinates satisfy the Kummer equation, and the quantity A^2 of the
   lifting test is a square, so a class of (a, b, c) modulo 64 must admit
   some fourth coordinate (of odd or even denominator) with both. The
   admitted c of a class of (a, b) form the word that fills the bit array
   of the row, and rows whose class admits nothing are not sieved.

## Examples

Real times of single-threaded runs (`-t 1`, the default) of the scalar build
on a laptop of 2023 (Intel Core i7-1355U), the medians of three; a laptop's
clock drifts under load, so take them as indicative. The last column gives
the mode of the sieve -- whole rows, the tube of the bound on the fourth
coordinate, the tube cut to the real region, or the real region (see "How
the program works") -- and the numbers of primes of its three stages.

| command                                   | time    | points | sieve |
|-------------------------------------------|--------:|-------:|-------|
| `j-points '1 6 5 22 22 8 1' 200`          | 0.02 s  |     11 | tube, 0+10+13 |
| `j-points '1 6 5 22 22 8 1' 500`          | 0.05 s  |     12 | tube, 0+9+13 |
| `j-points '1 6 5 22 22 8 1' 1000`         | 0.22 s  |     15 | tube, 0+9+13 |
| `j-points '1 6 5 22 22 8 1' 2000`         | 0.92 s  |     16 | tube cut to the region, 0+9+13 |
| `j-points '1 6 5 22 22 8 1' 5000`         | 9.9 s   |     16 | tube cut to the region, 0+9+12 |
| `j-points '1 6 5 22 22 8 1' 10000`        | 70 s    |     17 | tube cut to the region, 0+8+12 |
| `j-points '1 178 817 -274 16 1' 200`      | 0.01 s  |     35 | tube, 0+10+14 |
| `j-points '1 178 817 -274 16 1' 500`      | 0.03 s  |     63 | tube, 0+9+13 |
| `j-points '1 178 817 -274 16 1' 1000`     | 0.04 s  |    101 | tube, 0+8+12 |
| `j-points '1 178 817 -274 16 1' 2000`     | 0.06 s  |    124 | tube cut to the region, 0+9+13 |
| `j-points '1 178 817 -274 16 1' 5000`     | 0.26 s  |    183 | tube, 0+8+12 |
| `j-points '1 178 817 -274 16 1' 10000`    | 0.90 s  |    250 | tube cut to the region, 0+8+12 |
| `j-points '21 116 171 128 55 12 1' 200`   | 0.02 s  |     15 | tube, 0+10+14 |
| `j-points '21 116 171 128 55 12 1' 500`   | 0.06 s  |     19 | tube, 0+9+13 |
| `j-points '21 116 171 128 55 12 1' 1000`  | 0.17 s  |     23 | tube, 0+8+11 |
| `j-points '21 116 171 128 55 12 1' 2000`  | 0.78 s  |     26 | tube, 0+8+11 |
| `j-points '21 116 171 128 55 12 1' 5000`  | 5.7 s   |     31 | tube, 0+8+11 |
| `j-points '21 116 171 128 55 12 1' 10000` | 28 s    |     35 | tube, 0+8+11 |
| `j-points '1 6 5 22 22 8 1' 1000 -a`      | 0.27 s  |     16 | rows, 9+17+19 |
| `j-points '1 6 5 22 22 8 1' 2000 -a`      | 1.3 s   |     17 | region, 8+14+18 |
| `j-points '1 6 5 22 22 8 1' 5000 -a`      | 15 s    |     18 | region, 8+15+19 |
| `j-points '1 6 5 22 22 8 1' 10000 -a`     | 93 s    |     18 | region, 8+20+24 |
| `j-points '21 116 171 128 55 12 1' 1000 -a` | 0.22 s |     39 | rows, 8+16+18 |
| `j-points '21 116 171 128 55 12 1' 2000 -a` | 1.3 s  |     43 | rows, 8+15+19 |
| `j-points '21 116 171 128 55 12 1' 5000 -a` | 13 s   |     47 | region, 8+16+20 |
| `j-points '21 116 171 128 55 12 1' 10000 -a` | 84 s  |     51 | region, 8+20+24 |
| `j-points "$C" 1000 -a`                   | 0.29 s  |   4923 | rows, 12+18+22 |
| `j-points "$C" 2000 -a`                   | 1.7 s   |   7650 | rows, 11+20+23 |
| `j-points "$C" 5000 -a`                   | 20 s    |  13395 | rows, 11+25+28 |
| `j-points "$C" 10000 -a`                  | 137 s   |  19516 | region, 11+25+29 |

The last four rows are for the
[curve with 642 known rational points](https://www.mathe2.uni-bayreuth.de/stoll/recordcurve.html),
whose Jacobian has many points of small height:
`C='247747600 -985905640 567207969 2396040466 52485681 -470135160 82342800'`.
With `-a` the bound is on the first three coordinates only, the
enumeration of the points of bounded canonical height (see the option); its
runs grow with the cube of the height bound, the plain runs in the tube
less. The threads (`-t`) divide these times by up to the number of cores, as
the change log records.

## Change log

The history of the versions, with their measurements, is in
[CHANGE_LOG.md](CHANGE_LOG.md).

Michael Stoll, October 1998 - September 2026.

Copyright (C) 1998-2026 Michael Stoll. Distributed under the GNU GPL, version 2
or (at your option) any later version; see gpl-2.0.txt.
