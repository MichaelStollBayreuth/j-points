# j-points

A program that finds rational points on the Jacobian of a curve of genus 2,
given as y^2 = f(x) with f of degree 5 or 6 with integral coefficients, up to
a bound on the naive height of their images on the Kummer surface. It sieves
the triples of the first three Kummer coordinates modulo small primes, in the
manner of [ratpoints](https://github.com/MichaelStollBayreuth/ratpoints), from
which it descends, and lifts the survivors to the Jacobian.

The full description, the options and the change log up to version 1.2 are in
[readme](readme). Version 2.0 (January 2022), which the readme does not describe
yet, extended the table of primes to those below 128 and added the option `-p`
to limit their number, replaced the Euclidean gcd of the sieve by a binary
coprimality test, reworked `kummer_init`, and gained a Makefile with `make test`.

The program needs the GNU gmp library. Build and test it with

    make j-points
    make test

Copyright (C) 1998-2022 Michael Stoll. Distributed under the GNU GPL, version 2
or (at your option) any later version; see gpl-2.0.txt.
