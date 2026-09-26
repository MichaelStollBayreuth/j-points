# j-points: change log

The versions of the program, oldest first, with what changed and what it
measured; the README describes the current version.

+ 25-Oct-1998: started writing version 0.1, taking ratpoints-1.2 as a
  starting point.
+ 12-Nov-1998:
  + Included some enhancements taken from ratpoints-1.4.
  + Avoid dealing with multiples of reduced triples.
  + This allows taking only squares for the first coordinate, when f is
    monic of degree 5.
+ 24-Nov-1998: eliminated a bug (the program didn't take into account the
  possibility that the Kummer equation gives a constant).
+ 23-Apr-2001, version 0.6: changed `!a&1` into `!(a&1)` (no real bug, but
  it prevented an optimisation from being used).
+ 02-May-2001, version 1.0:
  + Split the code into two files: j-points.c and j-sift.c.
  + Hand-optimised assembler code for j-sift for the i386 architecture.
  + Makefile.
  + A few minor bugfixes.
+ 08-Aug-2006, version 1.1: removed a call `mpz_mul_ui(&fff, &fff,
  (unsigned long)b)` in check_lifts that let the program miss points of the
  form (0 : b : c : d) when b was not a square (it was computing b^7 f(c/b)
  instead of b^6 f(c/b)).
+ 14-Apr-2016, version 1.2: added code to handle the case when `-a` is
  specified and the coordinates of the resulting point(s) exceed machine size.
+ 17-Jan-2022, version 2.0:
  + The table of primes extended to the 30 primes below 128 (it had the
    17 primes below 64), and the option `-p` to limit the number of primes
    considered.
  + The sieve tests coprimality by a binary gcd instead of the Euclidean one.
  + kummer_init rewritten.
  + A Makefile with the targets test, dist and install-bin, and the reference
    output testbase for the test.
+ From September 2026, work on improving the program was done with the
  help of **Claude Code** (mostly Fable 5.1).
+ 24-Sep-2026, version 2.1: an audit of the code and a test suite. The
  mathematics and the sieve were found sound; the points of 2.0 are
  unchanged wherever it printed the right ones. Fixed:
  + `-N 0` alone silently lost every point of the main search (the first
    stage was given one prime, the tables none); it is now an unsieved
    search.
  + A zero leading coefficient made the program print the origin as a point;
    zero leading coefficients are now dropped.
  + `-p` given after `-n` or `-N` left the first stage with more primes than
    the table had.
  + With `-a`, a `-f` format other than the default broke the printing of
    coordinates beyond the machine word; gmp prints them now.
  + Heights above 2^31 overflowed a difference in the lifting test.
  + The usage message names every option; the version banner said 1.2.
+ 26-Sep-2026, version 3.0: the sieve made faster, in steps; the points
  printed are the same.
  + The bookkeeping of a row without integer divisions: what the sieve
    needs to know about a chunk of the bit array (where each prime's table
    row starts in it, the residue of its first word) is computed once at
    set-up instead of per row, and the residues of the first two
    coordinates are kept and advanced from one row to the next. 1.2-1.4x.
  + The sieve tables allocated per prime in use, with rows of p+1 words,
    instead of a static array of 480 MB (of which a run touched about 100
    MB; now 15-20 MB), the row pointers advanced from row to row, and the
    fill of the bit array merged into the first prime's pass. Another
    1.1-1.15x.
  + The rows (a, b) that a sieving prime p empties outright are no longer
    sieved: when the leading coefficient is not a square modulo p, the
    curve has no point at infinity over F_p, so no point reduces to a
    Kummer point with first coordinate 0 and second coordinate nonzero
    modulo p; for an a divisible by such primes only the b divisible by
    them are searched. Nothing on curves with a square leading
    coefficient, 1.2x on the record curve (a quarter of its rows).
  + The table of primes extended to the 53 primes below 256 (it had the 30
    below 128), a third sieving stage per surviving bit for the primes
    without a sieve table, and the primes of each stage and their numbers
    chosen by a cost model of the run instead of the two ratios, whose
    options `-r` and `-R` are gone; `-M` pins the number of primes of the
    first two stages, `-N` now counts all three. The rates of the larger
    primes are estimated from a sample of the residue classes, and the
    tables are made for the primes chosen only, so that the set-up is
    cheaper than before at small height bounds. 1.6-1.8x on the record
    curve (which had been starved of primes), 1.0-1.03x on random curves.
  + The tube: a plain run analyses, for each a, the plane of (b, c) in
    cells and sieves only the words in which a fourth coordinate of height
    at most h can exist (see "How the program works"); rigorous bounds, so
    the points are the same. 1.3x to 14x at 2000 on the test curves, more
    at larger height bounds; nothing for `-a`.
  + The real region: f must be nonnegative at the real Mumford roots of
    a t^2 - b t + c, which excludes intervals of c from a row in closed form
    (see "How the program works"). An `-a` run sieves only the words of a
    row inside the region, a plain run cuts the tube's words to it, where
    the cost model finds it worth it; the model also counts the rows of a
    monic quintic correctly (only the square first coordinates), which
    improves its choice of primes there. 1.1x to 1.4x on `-a` runs of
    curves on which f is negative on a good part of the line (c1 at 2000
    and 4000, c3 at 2000 and 4000), nothing on the record curve (f > 0 nearly
    everywhere); 1.1x to 1.25x on plain runs of such curves. The range
    limit of the tube's rows raised from 8 to 24, and cells in which the
    condition holds throughout taken whole: 1.1x on c4's plain run.
  + The passes of the sieve over the bit array can be built with vectors of
    two or four words (`make VECTOR=sse2` or `avx2`), behind a compile-time
    switch; the vector walk runs through each period of a table row without
    a test per vector. 1.25x on `-a` runs with AVX2 and 1.2x with SSE2 on a
    recent laptop, 1.3-1.5x with AVX2 and 1.2-1.3x with SSE2 on two desktops
    (an eight-core Xeon E-2288G, a 16-core i7-13700K), 1.3x with SSE2 on an
    old Ivy Bridge; nothing for a plain run in the tube, which has no
    passes.
  + A condition at 2: a point with coprime integer coordinates satisfies
    the Kummer equation and has a certain square among its coordinates
    (the quantity A^2 of the lifting test), so a class of (a, b, c) modulo
    64 can carry a point only if some fourth coordinate, of odd or of even
    denominator, satisfies both modulo 64. The admitted c of a class of
    (a, b) form one word, which fills the bit array of the row, and rows
    whose class admits nothing are skipped. Measured at height bounds of
    500 and more, and only when it admits less than everything modulo 16.
    1.5x on `-a` runs of small-coefficient curves, 1.1x on #28, nothing on
    the record curve (f is a square modulo 4 there).
  + A sextic with f0 = 0 and f1 = 1 or -1 is searched as the monic quintic
    it becomes under x -> 1/x (and x -> -x), a quintic with leading
    coefficient -1 under x -> -x: the first coordinate is a square then, so
    the search runs over sqrt(h) values of it instead of h. The Kummer
    coordinates go to (c : b : a : d) and (a : -b : c : d), which keeps the
    height, so the points of the curve given are exactly those found,
    transformed back. 13x to 18x on such curves at 2000, 26x at 4000.
  + The polynomial is checked to be squarefree. The escapes `\n`, `\t` and
    `\\` in the format of `-f` are interpreted. `-w bound4` bounds the
    fourth coordinate separately from the first three, which gives the
    enumeration of the points of bounded modified naive height without a
    post-filter.
  + The exact check of a surviving triple with gmp, which the three stages
    and the condition at 2 leave 8-19 thousand of at height 2000 where the
    first stage alone left millions, takes at most 1.1% of a run that lasts
    longer than a blink (0.04-0.7% on the `-a` runs), so the pre-filter the
    plan held in reserve for it was not built.
  + `make tune` measures the constants of the cost model on the machine at
    hand and writes them to tuning.mk (tune.sh, on the model of ratpoints':
    a coordinate descent over thirteen constants on five runs at 2000 that
    cover the box, the tube, the tube cut to the region and the region,
    every candidate timed back to back with the current settings, then a
    pruning pass that keeps the moves that matter); `-c` sets the constants
    for one run. On a laptop and a desktop of 2023 the tune moved four of
    them the same way -- the price of a first-stage pass down, and with it
    the size penalty of the larger primes and the tube's test, the
    second-stage test cheaper -- and by more the wider the vectors of the
    passes; those are the defaults now, per vector width. Against the old
    ones: 3-9% on most runs of the scalar build, 5-12% with SSE2 and AVX2,
    the region runs unmoved, and a tune from the new defaults moves nothing
    on the desktop; an Ivy Bridge of 2012 still gains up to 9% on most runs,
    though its own tune moves other constants. A Xeon of 2019 points the
    same way but was too noisy to measure.
  + Threads (`-t`): the values of the first coordinate are handed out to the
    threads in order, each sieves its rows with its own state, and the
    points of each value are printed in order, so the output is that of a
    single thread. Four times faster with eight threads on the record curve
    on a laptop with two performance and eight efficiency cores; 6.4x with
    eight threads on an eight-core Xeon (7.4x with its hyperthreads), 7.4x
    with eight and 11.5x with 24 threads on a 16-core i7 with eight
    performance and eight efficiency cores; plain runs in the tube scale
    less (5-8x there), the sieve of whole rows and of the real region alike.
