#!/bin/sh
# The sieve against no sieve: for every curve of testcurves2, the points at a
# small height bound found with the sieving primes chosen as usual must be
# those found with no sieving at all (-n 0 -N 0 checks every coprime triple
# exactly), with and without -a.  This tests the sieve tables and the sieving
# loops against the exact checks of the same program; the exact checks
# themselves are tested against Magma by test2 and verify-test3.py.
# Prints one line per curve and exits 1 if any pair of runs differs.
JP=${JP:-./j-points}
[ -x "$JP" ] || { echo "$JP: not an executable" >&2; exit 2; }
JPOPTS=${JPOPTS:-}
H=${H:-50}
CURVES=${CURVES:-testcurves2}
status=0
grep -v '^#' "$CURVES" | grep -v '^$' | while IFS=';' read -r coeffs h comment; do
  for opt in "" "-a"; do
    sieved=$("$JP" "$coeffs" $H -q $opt $JPOPTS | LC_ALL=C sort)
    exact=$("$JP" "$coeffs" $H -q $opt -n 0 -N 0 $JPOPTS | LC_ALL=C sort)
    if [ "$sieved" = "$exact" ]; then
      echo "# j-points '$coeffs' $H $opt: same $(echo "$sieved" | grep -c '(') points"
    else
      echo "# j-points '$coeffs' $H $opt: DIFFERENT"
      echo "$sieved" > testbrute-sieved.out; echo "$exact" > testbrute-exact.out
      diff testbrute-sieved.out testbrute-exact.out
      echo DIFFERENT > testbrute-failed.out
    fi
  done
done
[ -f testbrute-failed.out ] && { rm -f testbrute-failed.out; exit 1; }
exit 0
