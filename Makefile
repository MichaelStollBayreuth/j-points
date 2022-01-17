#   j-points-2.0
#    - A program to find rational points on Jacobians of genus 2 curves
#   Copyright (C) 1998, 2006, 2016, 2022  Michael Stoll
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
RM = rm -f
INSTALL = cp

INSTALL_DIR = /usr/local

CCFLAGS0 = -Wall -O2 -fomit-frame-pointer
CCFLAGS =
LFLAGS = -lgmp -lgcc -lc -lm

VERSION = 2.0

# Files that make up the distribution
DISTFILES = Makefile j-points.h j.points-${VERSION}.c j-sift-${VERSION}.c readme testbase

# Temporary files that are generated during build
# and can be removed afterwards
TEMPFILES = j-points-${VERSION}.o j-sift-${VERSION}.o j-sift-${VERSION}.s test.out

# Executables produced when building
TARGETFILES = j-points

all: j-points

test: j-points
	time ./j-points '21 116 171 128 55 12 1' 2000 -q > test.out
	cmp -s testbase test.out || echo "Test failed!"

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

j-points: j-points-${VERSION}.o j-sift-${VERSION}.o
	${CC} j-points-${VERSION}.o j-sift-${VERSION}.o -o j-points ${LFLAGS} ${CCFLAGS}

j-points-${VERSION}.o: j-points-${VERSION}.c j-points.h
	${CC} j-points-${VERSION}.c -c -o j-points-${VERSION}.o ${CCFLAGS0} ${CCFLAGS}

j-sift-${VERSION}.o: j-sift-${VERSION}.c j-points.h
	${CC} j-sift-${VERSION}.c -c -o j-sift-${VERSION}.o ${CCFLAGS0} -funroll-loops ${CCFLAGS}

j-sift-${VERSION}.s: j-sift-${VERSION}.c j-points.h
	${CC} j-sift-${VERSION}.c -S -o j-sift-${VERSION}.s ${CCFLAGS0} -funroll-loops ${CCFLAGS}


