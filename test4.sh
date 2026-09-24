#!/bin/sh
# The exact check against Magma's Kummer package: for every curve of
# testcurves2 the points at the height bound 30, with the bound on all four
# coordinates and with -a (the bound on the triple of the first three only),
# sorted.  "make test4" compares the output with testbase4, which mkref4.m
# made by brute force over every coprime triple (a, b, c) up to 30, asking
# Magma for the points of the Kummer surface above each triple and for their
# lifts to the Jacobian; that shares nothing with j-points, unlike Magma's
# own search for points, which is a port of it.  Height 30 is small, but it
# is the lifting test and the Kummer equation that are on trial here, on
# every kind of curve; the sieve is tested by test2 and testbrute.
# Before each invocation the script prints its arguments.
JP=${JP:-./j-points}
[ -x "$JP" ] || { echo "$JP: not an executable" >&2; exit 2; }
JPOPTS=${JPOPTS:-}
grep -v '^#' testcurves2 | grep -v '^$' | while IFS=';' read -r coeffs h comment; do
  echo "# j-points '$coeffs' 30"
  "$JP" "$coeffs" 30 -q $JPOPTS | LC_ALL=C sort
  echo "# j-points '$coeffs' 30 -a"
  "$JP" "$coeffs" 30 -q -a $JPOPTS | LC_ALL=C sort
done
