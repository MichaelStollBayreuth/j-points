#   j-points-2.1
#    - A program to find rational points on Jacobians of genus 2 curves
#   Copyright (C) 1998, 2006, 2016, 2022, 2026  Michael Stoll
#
#   This program is free software: you can redistribute it and/or
#   modify it under the terms of the GNU General Public License
#   as published by the Free Software Foundation, either version 2 of
#   the License, or (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of version 2 of the GNU General
#   Public License along with this program.
#   If not, see <http://www.gnu.org/licenses/>.
#
#
#   Makefile

CC = gcc
# the test target times its run with the shell's time, which /bin/sh lacks
SHELL = /bin/bash
RM = rm -f
INSTALL = cp

INSTALL_DIR = /usr/local

CCFLAGS0 = -Wall -O2 -fomit-frame-pointer
CCFLAGS =
LFLAGS = -lgmp -lgcc -lc -lm

VERSION = 2.1

# Files that make up the distribution
DISTFILES = Makefile j-points.h j-points.c j-sift.c README.md gpl-2.0.txt \
            testbase testcurves2 testbase2 testbase3 \
            test2.sh test3.sh testbrute.sh mkref2.m verify-test3.py

# Temporary files that are generated during build
# and can be removed afterwards
TEMPFILES = j-points.o j-sift.o j-sift.s test.out test2.out test3.out testbrute.out \
            testbrute-sieved.out testbrute-exact.out testbrute-failed.out verify-test3.m

# Executables produced when building
TARGETFILES = j-points

all: j-points

FAILED = "Test failed!"
# what a test does when its output differs from the reference: print the
# message and fail the recipe, so that the exit status of make reports it
FAIL = { echo ${FAILED}; false; }

# The suites "make test" runs, each a target below.  A test whose output
# differs from its reference prints "Test failed!" and fails its target;
# "make test" runs every suite whatever the earlier ones did and fails at
# the end if any of them failed.
TESTS = test1 test2 test3 testbrute

.PHONY: test
test:
	@status=0; for t in ${TESTS}; do \
	   ${MAKE} --no-print-directory $$t || status=1; done; \
	 exit $$status

# One curve with many points at the height bound 2000, timed: the run that
# measures the sieve.  The 26 points of testbase are those Magma finds.
test1: j-points testbase
	time ./j-points '21 116 171 128 55 12 1' 2000 -q > test.out
	cmp -s testbase test.out || ${FAIL}

# The 61 curves of testcurves2 (random ones at every height bound that
# matters for the bit arrays, square and negative leading coefficients,
# rational roots, bad reduction at the small primes, degree 5 monic and
# not, large coefficients) against testbase2, which mkref2.m made with
# Magma's own search for the points.
test2: j-points testbase2 testcurves2 test2.sh
	./test2.sh > test2.out 2>&1
	cmp -s testbase2 test2.out || ${FAIL}

# The options, the messages and the errors (see test3.sh), against
# testbase3.  The point lists of the reference were checked against Magma
# by verify-test3.py, which can be run again whenever testbase3 changes.
test3: j-points testbase3 test3.sh
	./test3.sh > test3.out 2>&1
	cmp -s testbase3 test3.out || ${FAIL}

# The sieve against no sieve, on the curves of testcurves2 at a small
# height bound (see testbrute.sh); the script itself reports a difference.
testbrute: j-points testcurves2 testbrute.sh
	./testbrute.sh > testbrute.out || ${FAIL}

install-bin: j-points
	${INSTALL} j-points ${INSTALL_DIR}/bin/
	chmod 755 ${INSTALL_DIR}/bin/j-points

dist: ${DISTFILES}
	mkdir -p j-points-${VERSION}
	cp ${DISTFILES} j-points-${VERSION}/
	tar --create --file=j-points-${VERSION}-`date --rfc-3339=date`.tar.gz --gzip --dereference j-points-${VERSION}
	rm -r j-points-${VERSION}

clean:
	${RM} ${TEMPFILES}

distclean: clean
	${RM} ${TARGETFILES}

j-points: j-points.o j-sift.o
	${CC} j-points.o j-sift.o -o j-points ${LFLAGS} ${CCFLAGS}

j-points.o: j-points.c j-points.h
	${CC} j-points.c -c -o j-points.o ${CCFLAGS0} ${CCFLAGS}

j-sift.o: j-sift.c j-points.h
	${CC} j-sift.c -c -o j-sift.o ${CCFLAGS0} -funroll-loops ${CCFLAGS}

j-sift.s: j-sift.c j-points.h
	${CC} j-sift.c -S -o j-sift.s ${CCFLAGS0} -funroll-loops ${CCFLAGS}


