#!/bin/sh
# Runs j-points on every curve of testcurves2 at its height bound and prints
# the points sorted; "make test2" compares the output with testbase2, the
# program's own output, which verify-test2.sh checks against the unsieved run
# of every curve (-n 0 -N 0, every coprime triple tested exactly).  Before
# each invocation the script prints its arguments, so that a failing curve
# can be found from the diff of test2.out against testbase2.
#
# The script takes the program from $JP, ./j-points by default, and adds
# $JPOPTS to every invocation, so that a build with other flags, or other
# options, can be run against the same reference.
JP=${JP:-./j-points}
# no program to run: exit 2 (the comparison with the reference is make's,
# whose target fails with 1 when they differ)
[ -x "$JP" ] || { echo "$JP: not an executable" >&2; exit 2; }
JPOPTS=${JPOPTS:-}
# the list of curves, testcurves2 by default (testcurves-rich for "make testrich")
CURVES=${CURVES:-testcurves2}
grep -v '^#' "$CURVES" | grep -v '^$' | while IFS=';' read -r coeffs h comment; do
  echo "# j-points '$coeffs' $h"
  "$JP" "$coeffs" "$h" -q $JPOPTS | LC_ALL=C sort
done
