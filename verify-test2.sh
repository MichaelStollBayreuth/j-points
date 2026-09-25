#!/bin/sh
# Verifies testbase2 independently of the sieve: for every curve of
# testcurves2, the unsieved run at its height bound (-n 0 -N 0: every coprime
# triple checked exactly, no sieve table involved) must print the points that
# testbase2 lists for it.  This takes about an hour of CPU time (the second
# README curve at 4500 alone takes twenty minutes), so it is not part of
# "make test"; it runs $P curves at a time, 4 by default.  Prints one line per
# curve and exits 1 if any curve differs.  The exact check the unsieved run
# relies on is what test4 compares with Magma.  Another list and reference
# can be given in $CURVES and $REF (testcurves2 and testbase2 by default).
JP=${JP:-./j-points}
P=${P:-4}
CURVES=${CURVES:-testcurves2}
REF=${REF:-testbase2}
# options added to every run, as test2.sh adds them (-a for the rich list)
JPOPTS=${JPOPTS:-}
[ -x "$JP" ] || { echo "$JP: not an executable" >&2; exit 2; }
if [ "$1" = "--one" ]; then
  # one curve: the block of testbase2 against the unsieved run
  coeffs=$2; h=$3
  ref=$(awk -v hdr="# j-points '$coeffs' $h" '$0 == hdr {f = 1; next} /^#/ {f = 0} f' "$REF")
  got=$("$JP" "$coeffs" "$h" -q -n 0 -N 0 $JPOPTS | LC_ALL=C sort)
  if [ "$ref" = "$got" ]; then
    echo "# j-points '$coeffs' $h: same $(echo "$got" | grep -c '(') points"
  else
    echo "# j-points '$coeffs' $h: DIFFERENT"
    echo "$ref" > "verify-test2-ref-$$.out"; echo "$got" > "verify-test2-got-$$.out"
    diff "verify-test2-ref-$$.out" "verify-test2-got-$$.out"
    rm -f "verify-test2-ref-$$.out" "verify-test2-got-$$.out"
    exit 1
  fi
  exit 0
fi
rm -f verify-test2-failed.out
grep -v '^#' "$CURVES" | grep -v '^$' | cut -d';' -f1,2 | tr ';' '\n' \
  | xargs -d '\n' -n 2 -P "$P" sh -c 'CURVES="$CURVES" REF="$REF" JPOPTS="$JPOPTS" "$0" --one "$1" "$2" || echo failed >> verify-test2-failed.out' "$0"
[ -f verify-test2-failed.out ] && { rm -f verify-test2-failed.out; exit 1; }
exit 0
