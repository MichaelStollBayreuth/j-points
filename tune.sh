#!/bin/bash
# Find good values for the constants of the cost model that chooses the
# sieving primes, their stages and the mode of the sieve -- the COST_
# constants of j-points.c -- on this machine, and write them to tuning.mk,
# which the Makefile includes.  Run it as "make tune"; it needs ./j-points
# and modifies no source file.
#
# The constants are the costs in cycles of the pieces of work of a run (a word
# of a first-stage pass, a word or a bit that survives a stage, an exact
# check, a row, a cell of the tube's analysis, ...), and only their ratios
# matter: they decide how many primes each stage gets and which ones, and
# whether a run sieves the box (whole rows), the tube of the bound on the
# fourth coordinate, the tube cut to the real region, or the ranges of the
# region.  Five runs at the height bound 2000 cover the modes: the curve of
# test1 (the tube), the record curve of Mueller-Stoll with -a (the box), c1
# of ratpoints with -a (the region), a monic quintic (the tube cut to the
# region) and #28 with -a (the box, on a curve whose rows the condition at 2
# thins).  A setting is judged on the sum of the five times, so that a value
# that suits one mode and ruins another loses.
#
# The sweep is a coordinate descent: each of the constants that decide the
# choices in turn, at half, 0.7, 1.4 and twice its value, the others at the
# best values found so far, the best of the stage winning however small its
# lead; then a second pass from the result, unless the first moved nothing.
# The cost surface is flat in most directions, so most of the winners are
# drift within the noise -- but a drifted constant can be what lets the next
# one gain (the cost of a surviving word and of its second-stage test only
# pay together), which a descent that moves nothing for less than the noise
# never finds.  So the descent runs free, and a pruning pass then puts every
# constant that moved back to its value in turn and leaves it back unless the
# set is then more than STEP worse than the full set of winners (against the
# full set, not against the last pruned set: the losses of the steps add up,
# and this bounds their sum): tuning.mk lists what matters.  The constants
# are passed to every run with -c, so no rebuild is needed; the winners go to
# tuning.mk as -D flags.
#
# Measuring is the delicate part, and naive timing does not work.  The cost
# surface is flat -- a factor of two in most constants costs under 3% -- while
# a laptop under sustained load drifts by 25% as it heats up, so a straight
# sweep measures the clock, not the settings.  Therefore:
#
#   * a warm-up brings the machine to a steady state before anything counts;
#   * every candidate is timed back to back with the current settings and
#     only the ratio of the two is kept, so drift cancels;
#   * the median over rounds is used, not a single reading;
#   * the current settings are themselves among the candidates of every
#     stage, so their measured ratio to themselves says how much noise is
#     left; if that is more than NOISE in any stage, the run is declared
#     inconclusive;
#   * nothing is written unless the winner beats the current settings by
#     more than MARGIN, which keeps a noisy run from making the build slower.
#
# The runs are pinned to one core (CPU, default 0; CPU= for no pinning) and
# timed by the shell.  About 40 minutes a pass on an idle machine.  Do not
# run it under "make -j" or on a busy machine.
set -u
export LC_ALL=C

ROUNDS=${ROUNDS:-3}
PASSES=${PASSES:-2}
MARGIN=${MARGIN:-0.98}     # accept only if the best ratio is below this
STEP=${STEP:-0.01}         # the pruning leaves the winners' gain within this
NOISE=${NOISE:-0.02}       # ... and the baseline's self-ratio is within this
WARMUP=${WARMUP:-20}       # seconds of load before measuring
CPU=${CPU-0}               # the core the runs are pinned to; empty: none
JP=${JP:-./j-points}
# the constants swept, in this order, and the factors tried on each
CONSTS=${CONSTS:-"AND SIZE WORD TEST2 BIT TEST3 EXACT ROW CELL TWORD TTEST RROW RANGE"}
FACTORS=${FACTORS:-"0.5 0.7 1.4 2"}
# the runs timed: the arguments of j-points but -q and -c
RUNS=(
  "'21 116 171 128 55 12 1' 2000"
  "'247747600 -985905640 567207969 2396040466 52485681 -470135160 82342800' 2000 -a"
  "'1 6 5 22 22 8 1' 2000 -a"
  "'1 178 817 -274 16 1' 2000"
  "'1 1 -2 -15 13 2 1' 2000 -a"
)

[ -x "$JP" ] || { echo "tune.sh: $JP is missing; build it first" >&2; exit 1; }

TMP=$(mktemp -d) || exit 1
trap 'rm -rf "$TMP"' EXIT INT TERM

# The settings to measure against: the constants in force in the binary
# ("j-points -c list": the defaults of j-points.c for its width of the
# passes, or those of tuning.mk when it is in effect and the binary was
# built with them).  They are passed explicitly to every run, the baseline
# included, so that the comparison never depends on what happens to be
# compiled in.
declare -A cur win trial
$JP -c list > "$TMP/list" 2>/dev/null || { echo "tune.sh: $JP -c list failed" >&2; exit 1; }
for k in $CONSTS; do
  v=$(awk -v k="$k" '$1 == k { print $2 }' "$TMP/list")
  [ -n "$v" ] || { echo "tune.sh: $JP -c list does not name COST_$k" >&2; exit 1; }
  cur[$k]=$v
done
if [ -f tuning.mk ] && [ "$(sed -n 's/^TUNED_FOR *= *//p' tuning.mk)" = "${TUNE_CONFIG:-}" ]
then echo "tuning.mk is already in effect; measuring against its values"; fi
# the -c argument for an assignment of the constants
settings() { local -n arr=$1; local s="" k
             for k in $CONSTS; do s+="${s:+,}$k=${arr[$k]}"; done; echo "$s"; }
BASE=$(settings cur)
for k in $CONSTS; do win[$k]=${cur[$k]}; done
echo "measuring against $BASE"
echo "on $(uname -n), ${TUNE_CONFIG:-no configuration given}, $ROUNDS rounds, up to $PASSES passes"

if [ -n "$CPU" ] && command -v taskset >/dev/null 2>&1; then PIN="taskset -c $CPU"; else PIN=""; fi

# the reference output of every run, which every candidate must reproduce
for i in "${!RUNS[@]}"; do
  eval "$JP ${RUNS[$i]} -q -c '$BASE'" > "$TMP/ref$i" || exit 1
done
if [ -f testbase ]; then
  cmp -s "$TMP/ref0" testbase || { echo "tune.sh: $JP disagrees with testbase" >&2; exit 1; }
fi

# one timing = the sum over the runs of the elapsed time with the settings
# $1; empty output when a run fails or gives other points than the baseline
time_runs() {
  local tot=0 i t
  for i in "${!RUNS[@]}"; do
    t=$( { TIMEFORMAT=%R; time eval "$PIN $JP ${RUNS[$i]} -q -c '$1'" > "$TMP/out"; } 2>&1 ) || return 1
    cmp -s "$TMP/out" "$TMP/ref$i" \
      || { echo "tune.sh: ${RUNS[$i]} with -c $1 gave other points" >&2; return 1; }
    tot=$(awk -v a="$tot" -v b="$t" 'BEGIN{printf "%.3f", a+b}')
  done
  echo "$tot"
}

printf 'warming up (%ss) ' "$WARMUP"
end=$(( $(date +%s) + WARMUP ))
while [ "$(date +%s)" -lt $end ]; do time_runs "$BASE" > /dev/null; printf '.'; done
echo

# measure(): "label<TAB>settings" lines from $1 -> "label median_ratio" in $2
measure() {
  : > "$TMP/raw"
  local n r off label args tc tb
  n=$(wc -l < "$1")
  for ((r = 1; r <= ROUNDS; r++)); do
    printf '  round %d/%d:' $r $ROUNDS
    off=$(( (r - 1) % n ))
    { tail -n +$((off + 1)) "$1"; [ $off -gt 0 ] && head -n $off "$1"; true; } > "$TMP/order"
    while IFS='	' read -r label args; do
      tc=$(time_runs "$args"); [ -n "$tc" ] || exit 1
      tb=$(time_runs "$BASE"); [ -n "$tb" ] || exit 1   # the current settings, right next to it
      echo "$label $(awk -v c="$tc" -v b="$tb" 'BEGIN{printf "%.5f", c/b}') $tc $tb" >> "$TMP/raw"
      printf ' .'
    done < "$TMP/order"
    echo
  done
  sort "$TMP/raw" | awk '
    { v[$1] = v[$1] " " $2; k[$1] = 1 }
    END { for (a in k) { m = split(v[a], x, " ")
                         for (i = 1; i < m; i++) for (j = i+1; j <= m; j++)
                           if (x[j]+0 < x[i]+0) { t = x[i]; x[i] = x[j]; x[j] = t }
                         # the median; for an even number of rounds the mean
                         # of the two central values, not the lower one
                         med = (m % 2) ? x[int((m+1)/2)] \
                                       : (x[m/2] + x[m/2+1]) / 2
                         printf "%s %.5f\n", a, med } }' | sort > "$2"
}

report() { awk '{ printf "    %-14s %+7.1f%%\n", $1, 100*($2-1) }' "$1"; }
# the best candidate of a stage (not the current settings, not the value in
# place, not the winners), and the ratio a label measured
best_of() { awk '$1 != "current" && $1 != "unchanged" && $1 != "winners" { if (m == "" || $2 < m) { m = $2; k = $1 } }
                 END { print k }' "$1"; }
ratio_of() { awk -v k="$2" '$1 == k { print $2 }' "$1"; }

: > "$TMP/selfs"
LASTBEST=1
for ((pass = 1; pass <= PASSES; pass++)); do
  moved=0
  for k in $CONSTS; do
    # the candidates: the current settings (the noise check), the winners
    # so far with this constant unchanged (unless they are the current
    # settings), and the winners with it at each factor
    : > "$TMP/cand"
    printf 'current\t%s\n' "$BASE" >> "$TMP/cand"
    if [ "$(settings win)" != "$BASE" ]; then
      printf 'unchanged\t%s\n' "$(settings win)" >> "$TMP/cand"; skip=current
    else skip=; fi
    for f in $FACTORS; do
      v=$(awk -v x="${win[$k]}" -v f="$f" 'BEGIN{printf "%.4g", x*f}')
      for kk in $CONSTS; do trial[$kk]=${win[$kk]}; done
      trial[$k]=$v
      printf '%s=%s\t%s\n' "$k" "$v" "$(settings trial)" >> "$TMP/cand"
    done
    echo
    echo "pass $pass: COST_$k, now ${win[$k]}$( [ -n "$skip" ] && echo ", the other winners so far in place" )"
    measure "$TMP/cand" "$TMP/res"
    report "$TMP/res"
    ratio_of "$TMP/res" current >> "$TMP/selfs"
    # the value in place is the winners' (unchanged), or the current
    # settings' in a first stage; the best candidate replaces it when it
    # measured better at all
    ref=$( [ -n "$skip" ] && echo unchanged || echo current )
    refratio=$(ratio_of "$TMP/res" $ref)
    best=$(best_of "$TMP/res")
    bestratio=$(ratio_of "$TMP/res" "$best")
    if awk -v b="$bestratio" -v r="$refratio" 'BEGIN { exit !(b < r) }'; then
      win[$k]=${best#*=}; moved=1
      echo "    -> COST_$k = ${win[$k]}"
      LASTBEST=$bestratio
    else
      LASTBEST=$refratio
    fi
  done
  [ $moved = 0 ] && break
done

# the pruning: each constant that moved goes back to its current value in
# turn, the moves left so far in place, and stays back unless the set is
# then more than STEP worse than the full set of winners (measured again in
# every stage: "winners"; "unchanged" is the set before this step)
declare -A full
for k in $CONSTS; do full[$k]=${win[$k]}; done
for k in $CONSTS; do
  [ "${win[$k]}" = "${cur[$k]}" ] && continue
  : > "$TMP/cand"
  printf 'current\t%s\n' "$BASE" >> "$TMP/cand"
  printf 'winners\t%s\n' "$(settings full)" >> "$TMP/cand"
  [ "$(settings win)" != "$(settings full)" ] \
    && printf 'unchanged\t%s\n' "$(settings win)" >> "$TMP/cand"
  for kk in $CONSTS; do trial[$kk]=${win[$kk]}; done
  trial[$k]=${cur[$k]}
  printf '%s=%s\t%s\n' "$k" "${cur[$k]}" "$(settings trial)" >> "$TMP/cand"
  echo
  echo "pruning: COST_$k back from ${win[$k]} to ${cur[$k]}, the moves left so far in place"
  measure "$TMP/cand" "$TMP/res"
  report "$TMP/res"
  ratio_of "$TMP/res" current >> "$TMP/selfs"
  fullratio=$(ratio_of "$TMP/res" winners)
  backratio=$(ratio_of "$TMP/res" "$k=${cur[$k]}")
  if awk -v b="$backratio" -v f="$fullratio" -v s="$STEP" 'BEGIN { exit !(b < f + s) }'; then
    win[$k]=${cur[$k]}
    echo "    -> COST_$k back to ${cur[$k]}: within $(awk -v s=$STEP 'BEGIN{printf "%.0f", 100*s}')% of the winners without it"
    LASTBEST=$backratio
  else
    echo "    -> COST_$k stays at ${win[$k]}"
    LASTBEST=$( [ "$(settings win)" != "$(settings full)" ] && ratio_of "$TMP/res" unchanged || echo "$fullratio" )
  fi
done

# how far the current settings measured from themselves, in any stage: all
# are the same comparison of the binary with itself, so any one being far
# from 1 means the machine could not be measured on
SELF=$(awk '{ d = $1 - 1; if (d < 0) d = -d; if (d > w) w = d }
            END { printf "%.5f", 1 + w }' "$TMP/selfs")
WIN=$(settings win)
flags=""; moves=""
for k in $CONSTS; do
  flags+=" -DCOST_$k=${win[$k]}"
  [ "${win[$k]}" != "${cur[$k]}" ] && moves+="${moves:+, }COST_$k ${cur[$k]} -> ${win[$k]}"
done
echo
verdict=$(awk -v b="$LASTBEST" -v s="$SELF" -v m="$MARGIN" -v z="$NOISE" \
  'BEGIN { d = s - 1; if (d < 0) d = -d
           if (d > z) print "noisy"; else if (b < m) print "accept"; else print "keep" }')
[ "$verdict" = accept ] && [ "$WIN" = "$BASE" ] && verdict=keep

case $verdict in
  noisy)
    echo "The current settings measured $(awk -v s=$SELF 'BEGIN{printf "%.1f", 100*(s-1)}')% away from themselves, so this"
    echo "machine is too noisy right now to tell the settings apart.  Nothing"
    echo "written.  Try again when it is idle, or with 'ROUNDS=6 make tune'."
    ;;
  keep)
    if [ "$WIN" = "$BASE" ]; then
      echo "No constant moved: the current settings are kept and nothing is written."
    else
      echo "The settings found ($moves) measured $(awk -v b=$LASTBEST 'BEGIN{printf "%.1f", 100*(1-b)}')% better than the"
      echo "current ones, less than the required $(awk -v m=$MARGIN 'BEGIN{printf "%.0f", 100*(1-m)}')%: kept, nothing written."
    fi
    ;;
  accept)
    cat > tuning.mk <<END
# Machine-dependent tuning, written by "make tune" on $(date +%Y-%m-%d) on $(uname -n).
# Delete this file to go back to the values compiled into j-points.c.
# TUNED_FOR records the configuration it was measured for; the Makefile
# ignores this file if the configuration has changed since.
# Measured against $BASE
# and moved: $moves
# (the constants not listed keep their values in j-points.c).
TUNED_FOR = ${TUNE_CONFIG:-unknown}
TUNEFLAGS =$flags
END
    echo "Wrote tuning.mk with the moves $moves"
    echo "$(awk -v b=$LASTBEST 'BEGIN{printf "%.1f", 100*(1-b)}')% better than the current settings."
    echo "Run 'make' to rebuild with it; delete tuning.mk to discard it."
    ;;
esac
