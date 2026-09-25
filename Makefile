#   j-points-3.0
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

CCFLAGS0 = -Wall -O2 -fomit-frame-pointer -pthread
CCFLAGS =
# The width of the sieve's passes, in words (see VW in j-points.h): the
# default is the scalar walk; "make VECTOR=avx2" builds the passes with
# vectors of four words (AVX2), "make VECTOR=sse2" with two (SSE2, any
# x86-64).  The build is then tied to that instruction set.
VECTOR =
VFLAGS_avx2 = -DVW=4 -mavx2
VFLAGS_sse2 = -DVW=2 -msse2
VFLAGS = ${VFLAGS_${VECTOR}}
LFLAGS = -pthread -lgmp -lgcc -lc -lm

# Machine-dependent tuning of the cost model that chooses the sieving primes,
# their stages and the mode of the sieve (the COST_ constants of j-points.c).
# "make tune" measures them and writes tuning.mk; without it the values
# compiled into j-points.c are used.  This needs GNU make for the
# conditionals; if yours is not GNU make, delete the block and either leave
# tuning.mk out or trust it unconditionally.
-include tuning.mk

TUNE_CONFIG = VECTOR=${VECTOR}
ifdef TUNED_FOR
ifneq (${TUNED_FOR},${TUNE_CONFIG})
$(warning tuning.mk was measured for "${TUNED_FOR}", but this build is)
$(warning "${TUNE_CONFIG}" -- ignoring it; run "make tune" again)
TUNEFLAGS =
endif
endif

VERSION = 3.0

# Files that make up the distribution
DISTFILES = Makefile j-points.h j-points.c j-sift.c README.md gpl-2.0.txt \
            testbase testcurves2 testbase2 testbase3 testbase4 \
            testcurves-rich testbase-rich \
            test2.sh test3.sh test4.sh testbrute.sh tune.sh \
            verify-test2.sh verify-test3.py mkref4.m

# Temporary files that are generated during build
# and can be removed afterwards
TEMPFILES = j-points.o j-sift.o j-sift.s build.stamp build.stamp.tmp test.out test2.out test3.out test4.out testrich.out testthreads.out \
            testbrute.out testbrute-sieved.out testbrute-exact.out testbrute-failed.out \
            verify-test2-failed.out verify-test3.m

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
TESTS = test1 test2 test3 test4 testbrute testrich testthreads

.PHONY: test
test:
	@status=0; for t in ${TESTS}; do \
	   ${MAKE} --no-print-directory $$t || status=1; done; \
	 exit $$status

# One curve with many points at the height bound 2000, timed: the run that
# measures the sieve.
test1: j-points testbase
	time ./j-points '21 116 171 128 55 12 1' 2000 -q > test.out
	cmp -s testbase test.out || ${FAIL}

# The 61 curves of testcurves2 (random ones at every height bound that
# matters for the bit arrays, square and negative leading coefficients,
# rational roots, bad reduction at the small primes, degree 5 monic and
# not, large coefficients) against testbase2, the program's output, which
# verify-test2.sh checks against the unsieved run of every curve (an hour
# of CPU time; run it again whenever testbase2 changes).
test2: j-points testbase2 testcurves2 test2.sh
	./test2.sh > test2.out 2>&1
	cmp -s testbase2 test2.out || ${FAIL}

# The options, the messages and the errors (see test3.sh), against
# testbase3.  The point lists of the reference were checked by
# verify-test3.py -- by brute force with Magma's Kummer package up to height
# 60, by the unsieved run above -- which can be run again whenever testbase3
# changes.
test3: j-points testbase3 test3.sh
	./test3.sh > test3.out 2>&1
	cmp -s testbase3 test3.out || ${FAIL}

# The exact check against Magma: every curve of testcurves2 at height 30,
# with and without -a (see test4.sh), against testbase4, which mkref4.m made
# by brute force over every coprime triple with Magma's Kummer package.
test4: j-points testbase4 testcurves2 test4.sh
	./test4.sh > test4.out 2>&1
	cmp -s testbase4 test4.out || ${FAIL}

# The sieve against no sieve, on the curves of testcurves2 at a small
# height bound (see testbrute.sh); the script itself reports a difference.
testbrute: j-points testcurves2 testbrute.sh
	./testbrute.sh > testbrute.out || ${FAIL}

# Fourteen runs with -a on point-rich curves (testcurves-rich: from the
# ratpoints suites, and the record curve of Mueller-Stoll) against
# testbase-rich, the output of 2.1 checked by verify-test2.sh against the
# unsieved run (JPOPTS=-a CURVES=testcurves-rich REF=testbase-rich
# ./verify-test2.sh) -- the regime of the enumeration of points of bounded
# canonical height, and the timing suite for the work on the sieve (about
# 8 s).
testrich: j-points testcurves-rich testbase-rich test2.sh
	JPOPTS=-a CURVES=testcurves-rich ./test2.sh > testrich.out 2>&1
	cmp -s testbase-rich testrich.out || ${FAIL}

# The threads: test1 with two, test2 and testrich with three, against the
# same references (the points are printed in the order of the search
# whatever the number of threads).
testthreads: j-points testbase testbase2 testbase-rich test2.sh
	./j-points '21 116 171 128 55 12 1' 2000 -q -t 2 > testthreads.out
	cmp -s testbase testthreads.out || ${FAIL}
	JPOPTS=-t\ 3 ./test2.sh > testthreads.out 2>&1
	cmp -s testbase2 testthreads.out || ${FAIL}
	JPOPTS=-a\ -t\ 3 CURVES=testcurves-rich ./test2.sh > testthreads.out 2>&1
	cmp -s testbase-rich testthreads.out || ${FAIL}

# Measure the constants of the cost model on this machine and write them to
# tuning.mk (see tune.sh).  Deliberately not part of "make all": it takes
# about 40 minutes a pass, wants an otherwise idle machine, and must not run
# under "make -j".  Run "make" afterwards to rebuild with the result.
.PHONY: tune
tune: j-points
	@TUNE_CONFIG='${TUNE_CONFIG}' ./tune.sh

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
	${RM} ${TARGETFILES} tuning.mk

j-points: j-points.o j-sift.o
	${CC} j-points.o j-sift.o -o j-points ${LFLAGS} ${CCFLAGS}

# What is compiled depends on VECTOR and on tuning.mk, which make cannot see
# by itself: build.stamp records the flags in use, so that a change to
# either rebuilds the objects.
.PHONY: FORCE
build.stamp: FORCE
	@echo '${VFLAGS} ${TUNEFLAGS}' > $@.tmp; \
	 if cmp -s $@.tmp $@; then rm $@.tmp; else mv $@.tmp $@; fi

j-points.o: j-points.c j-points.h build.stamp
	${CC} j-points.c -c -o j-points.o ${CCFLAGS0} ${VFLAGS} ${TUNEFLAGS} ${CCFLAGS}

j-sift.o: j-sift.c j-points.h build.stamp
	${CC} j-sift.c -c -o j-sift.o ${CCFLAGS0} -funroll-loops ${VFLAGS} ${CCFLAGS}

j-sift.s: j-sift.c j-points.h build.stamp
	${CC} j-sift.c -S -o j-sift.s ${CCFLAGS0} -funroll-loops ${VFLAGS} ${CCFLAGS}


