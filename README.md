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

Please acknowledge use of the program in published work.

## Usage

    j-points 'a_0 a_1 ... a_d' h [-n num_primes1] [-N num_primes2] [-p num_primes]
             [-r ratio1] [-R ratio2] [-s size] [-f format] [-1] [-q] [-a]

where

+ f(x) = a_d x^d + ... + a_1 x + a_0 with d = 5 or 6;
+ h is the bound on the naive height: the program finds the points
  (a : b : c : d) on the Kummer surface, with coprime integral coordinates,
  that lift to the Jacobian and satisfy |a|, |b|, |c|, |d| <= h. On a machine
  with 64-bit longs, h can be at most 3037000499.

The optional arguments, which can be given in any order, are these.

+ `-n num_primes1`: the number of primes used for the first sieving stage.
  By default it is chosen automatically, using `ratio1` below.
+ `-N num_primes2`: the number of primes used for the two sieving stages
  together. By default it is chosen automatically, using `ratio2` below.
+ `-p num_primes`: the number of primes, from the beginning of the table
  3, 5, 7, ..., 127, among which the sieving primes are chosen (default 20,
  at most 30).
+ `-r ratio1`: the ratio of the running time of the second versus the first
  stage of sieving (per bit). It is used to choose the number of sieving primes
  for the first stage automatically (default 5000) and is ignored when `-n` is
  given.
+ `-R ratio2`: the ratio of the running time needed for checking whether a
  surviving triple gives rise to points versus one step of the second sieving
  stage. It is used to choose the number of sieving primes for the second stage
  automatically (default 2.5) and is ignored when `-N` is given.
+ `-s size`: the size in kilobytes of the bit array in which the sieving is
  done (default 10). This determines the amount of memory the program uses.
  A smaller value can give better performance if the array then fits into the
  cache in its entirety, so the effect of this parameter depends heavily on
  the size and speed of the cache.
+ `-f format`: a format string for printf to print the points with. It should
  take four long integers, the coordinates on the Kummer surface. The default
  is `"(%ld, %ld, %ld, %ld)\n"`.
+ `-1`: stop as soon as one point has been found.
+ `-q`: suppress all messages other than the points found.
+ `-a`: find all points whose first three coordinates have naive height at
  most h, that is, without any bound on the fourth coordinate. Coordinates
  that exceed the machine word are handled through gmp.

## How the program works

Basically the program implements a loop over all triples (a, b, c) of first
Kummer coordinates with |a|, |b|, |c| <= h, where h is the height bound
given. For each of these triples all possible rational fourth coordinates d
are found such that

+ (a : b : c : d) is on the Kummer surface K, and
+ the point (a : b : c : d) in K(Q) lifts to J(Q).

The coordinates are scaled by the denominator of d, so that they are coprime
integers.

This test is split into two stages. First it is checked whether (a, b, c)
mod p are the first three coordinates of a point in K(F_p) that lifts to
J(F_p). This is done for a number of primes p. Only the surviving triples are
then used to compute the possible d's and to check whether the points thus
obtained lift to J(Q).

There are a number of improvements to this basic scheme.

1. We use bits to represent the individual numerators. In this way we can
   sieve as many numerators as bits fit into a long word (usually 64) at the
   same time, using bit-wise "and" operations.
2. When this kind of sieving has reduced the candidates considerably, we
   continue the sieving with more primes for each candidate separately.
3. We take more primes than necessary for the sieving procedure and determine
   in a first step which are the best ones. This is measured by the ratio of
   numbers surviving the corresponding step of the sieve (essentially the
   number of points mod p, divided by p^2).

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

Michael Stoll, October 1998 - January 2022.

Copyright (C) 1998-2022 Michael Stoll. Distributed under the GNU GPL, version 2
or (at your option) any later version; see gpl-2.0.txt.
