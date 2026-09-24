#!/bin/sh
# Runs j-points on every curve of testcurves2 at its height bound and prints
# the points sorted; "make test2" compares the output with testbase2, which
# mkref2.m made with Magma's own search for the points (an independent
# implementation of the Kummer surface and of the lifting to the Jacobian).
# Before each invocation the script prints its arguments, so that a failing
# curve can be found from the diff of test2.out against testbase2.
#
# The script takes the program from $JP, ./j-points by default, and adds
# $JPOPTS to every invocation, so that a build with other flags, or other
# options, can be run against the same reference.
JP=${JP:-./j-points}
# no program to run: exit 2 (the comparison with the reference is make's,
# whose target fails with 1 when they differ)
[ -x "$JP" ] || { echo "$JP: not an executable" >&2; exit 2; }
JPOPTS=${JPOPTS:-}
grep -v '^#' testcurves2 | grep -v '^$' | while IFS=';' read -r coeffs h comment; do
  echo "# j-points '$coeffs' $h"
  "$JP" "$coeffs" "$h" -q $JPOPTS | LC_ALL=C sort
done
