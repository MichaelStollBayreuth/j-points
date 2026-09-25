/***********************************************************************
 * j-points-3.0                                                        *
 *  - A program to find rational points on Jacobians of genus 2 curves *
 * Copyright (C) 1998, 2006, 2016, 2022, 2026  Michael Stoll           *
 *                                                                     *
 * This program is free software: you can redistribute it and/or       *
 * modify it under the terms of the GNU General Public License         *
 * as published by the Free Software Foundation, either version 2 of   *
 * the License, or (at your option) any later version.                 *
 *                                                                     *
 * This program is distributed in the hope that it will be useful,     *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of      *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the       *
 * GNU General Public License for more details.                        *
 *                                                                     *
 * You should have received a copy of version 2 of the GNU General     *
 * Public License along with this program.                             *
 * If not, see <http://www.gnu.org/licenses/>.                         *
 ***********************************************************************/

/***********************************************************************
 * j-sift.c                                                            *
 *  - the sieve: for fixed first two coordinates (a, b), the bit array *
 *    over the third coordinate c, sieved with the first-stage primes, *
 *    its surviving words tested against the second-stage primes, and  *
 *    the surviving bits tested against the third-stage primes and    *
 *    handed to the exact check; the bookkeeping of a row done       *
 *    without integer divisions                                        *
 ***********************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gmp.h>

#include "j-points.h"


/* The row of the third coordinate is sieved in chunks of at most
   array_size words, and the chunks are the same for every row (they do
   not depend on a and b).  So what the sieve needs to know about a chunk
   is computed once, by init_sift(), and not per row: for each first-stage
   prime, where in its table row (period p words) the chunk begins and how
   the walk over the chunk splits into a head up to the end of the table
   row, full periods of p words and a tail; for each second-stage prime,
   the residue of the chunk's first word. */
typedef struct
{ long first;  /* the table word the chunk begins with: in [1, p] for the
                  scalar walk, in [0, period) for the vector walk */
  long head;   /* the words up to the end of the table row (or the chunk) */
  long nper;   /* the full periods of p words after the head */
  long tail;   /* the words after the last full period */
} walk_spec;

typedef struct
{ long w_low, w_high;  /* the words of the chunk: w_low <= i < w_high */
  long nvec;           /* the vectors of VW words that cover the chunk */
  walk_spec *walk;     /* one per first-stage prime */
  long *res;           /* per sieving prime: w_low mod p */
} chunk_spec;

/* The walk of a prime p over the bit array reads its table row VW words
   at a time from any offset, so the row holds the pattern over a period
   that is a multiple of p and at least VW words, plus VW words of its
   continuation (with VW = 1: p + 1 words, word p repeating word 0). */
long sieve_period(long p)
{ long P = p;
  while(P < VW) { P += p; }
  return(P);
}
long sieve_rowlen(long p) { return(sieve_period(p) + VW); }
#if VW > 1
typedef bit_array vec __attribute__((vector_size(VW*sizeof(bit_array))));
typedef bit_array vecu __attribute__((vector_size(VW*sizeof(bit_array)),
                                      aligned(sizeof(bit_array))));
#endif

static chunk_spec *chunks;
static long num_chunks;

static long npr;    /* the number of sieving primes, all stages */
static unsigned short *res0;
/* res0[k*npr + n] = k mod pr[n] for 0 <= k <= array_size: with the residue
   of the chunk's first word, the residue of any word of the chunk costs an
   addition and a comparison instead of a division */

/* The third stage tests a bit (a, b, c) that survived the second stage
   against its primes one by one: (a : b : c) mod p is a point of the
   Kummer surface with a point of J above it iff, for a nonzero a, the
   entry (b/a, c/a) of is_point_on_j is set, for a = 0 and b nonzero the
   entry c/b of is_f_square (and the curve has a point at infinity mod p;
   without one such rows are never sieved, see row_step), and always when
   a = b = 0.  So a row (a, b) needs per prime a table row and a factor
   (1/a or 1/b mod p), computed with divisions when a bit of the row first
   gets here (most rows never do), by row3_setup. */
static unsigned char all_ones[MAX_PRIME_EVEN];  /* the row for a = b = 0 */
static _Thread_local unsigned char *tab3[NUM_PRIMES];  /* per third-stage prime: the row */
static _Thread_local long mult3[NUM_PRIMES];           /* ... and the factor */
static unsigned long recip3[NUM_PRIMES];
/* recip3[m] = ceil(2^32 / p) for the m-th third-stage prime p: for
   0 <= x < 2^16, (x * recip3[m]) >> 32 is exactly x / p, since x / p is
   either an integer or at least 1/p away from one, and the error of the
   product is below 2^-16.  So x mod p costs two multiplications. */
#define MOD3(x, m) ((x) - pr[sieve_primes2 + (m)] * (long)(((x) * recip3[m]) >> 32))
static _Thread_local long row3_a = -1, row3_b = 0;     /* the row the data is for */

/* the sieving primes, the lengths of their table rows and of the block
   of rows of one residue class of a, and the table rows of the current
   (a, b): rowptr[n] for the n-th prime */
static long pr[NUM_PRIMES], rowlen[NUM_PRIMES], blocklen[NUM_PRIMES];
static long period[NUM_PRIMES];      /* the period of the walk, see above */
static _Thread_local bit_array *rowptr[NUM_PRIMES];
static long w_low_all, w_high_all;   /* the words of a whole row */
static double th, th2;               /* the bound on the fourth coordinate and its square */
static double invp[NUM_PRIMES];     /* 1/p, for the parts of a walk over a range */

/* the tube's bands (see there): the words collected for the rows of the
   band a thread analyses, per row up to TUBE_RANGES ranges [lo, hi) */
#define TUBE_BAND 4096    /* rows analysed and then sieved together */
#define TUBE_RANGES 24    /* word ranges of the tube kept per row; more are merged */
static _Thread_local long band_b0;                /* the first row */
static _Thread_local long *band_nr = NULL;        /* ranges per row */
static _Thread_local long (*band_lo)[TUBE_RANGES], (*band_hi)[TUBE_RANGES];


/**************************************************************************
 * Set-up: the chunks and their tables                                    *
 **************************************************************************/

void init_sift(void)
{
  long w_low = (-height)>>LONG_SHIFT;    /* FLOOR(-height, LONG_LENGTH) */
  long w_high = (height>>LONG_SHIFT)+1;  /* CEIL(height+1, LONG_LENGTH) */
  long np1 = sieve_primes1, k, n;
  walk_spec *walk;
  long *res;

  npr = sieve_primes3;
  w_low_all = w_low; w_high_all = w_high;
  for(n = 0; n < sieve_primes3; n++)
  { pr[n] = prime[pnn[n]];
    period[n] = sieve_period(pr[n]);
    rowlen[n] = sieve_rowlen(pr[n]);
    blocklen[n] = pr[n]*rowlen[n];
  }
  for(n = 0; n < MAX_PRIME_EVEN; n++) { all_ones[n] = 1; }
  for(n = 0; n < sieve_primes3 - sieve_primes2; n++)
  { unsigned long p = pr[sieve_primes2 + n];
    recip3[n] = ((1UL << 32) + p - 1) / p;
  }
  num_chunks = CEIL(w_high - w_low, array_size);
  th = (double)dbound; th2 = th*th;
  init_thread_sieve();
  chunks = (chunk_spec *)malloc(num_chunks*sizeof(chunk_spec));
  walk = (walk_spec *)malloc((num_chunks*np1 + 1)*sizeof(walk_spec));
  res = (long *)malloc((num_chunks*npr + 1)*sizeof(long));
  res0 = (unsigned short *)malloc(((array_size + 1)*npr + 1)*sizeof(unsigned short));
  if(chunks == NULL || walk == NULL || res == NULL || res0 == NULL)
  { error(7); }
  for(k = 0; k < num_chunks; k++)
  { chunk_spec *ch = &chunks[k];
    long range;
    ch->w_low = w_low + k*array_size;
    ch->w_high = ch->w_low + array_size;
    if(ch->w_high > w_high) { ch->w_high = w_high; }
    range = ch->w_high - ch->w_low;
    ch->nvec = (range + VW - 1)/VW;
    ch->walk = &walk[k*np1];
    ch->res = &res[k*npr];
    for(n = 0; n < np1; n++)
    { long p = prime[pnn[n]], start;
      walk_spec *w = &ch->walk[n];
      if(VW == 1)
      { start = ch->w_low % p;
        if(start < 0) { start += p; }
        w->first = (start == 0) ? p : start;
        w->head = p - w->first;
        if(w->head > range) { w->head = range; }
        w->nper = (range - w->head)/p;
        w->tail = range - w->head - w->nper*p;
      }
      else
      { start = ch->w_low % period[n];
        if(start < 0) { start += period[n]; }
        w->first = start;
        w->head = w->nper = w->tail = 0;
      }
    }
    for(n = 0; n < npr; n++)
    { long p = pr[n], start = ch->w_low % p;
      if(start < 0) { start += p; }
      ch->res[n] = start;
    }
  }
  for(n = 0; n < npr; n++)
  { long p = pr[n];
    for(k = 0; k <= array_size; k++) { res0[k*npr + n] = k % p; }
    invp[n] = 1.0 / (double)p;
  }
  return;
}

/* the state of the calling thread: the bit array (aligned for the
   vector walk and padded by a vector) and the lists of the tube's bands */
void init_thread_sieve(void)
{
  survivors = (bit_array *)aligned_alloc(64, ((array_size + VW)*sizeof(bit_array) + 63)/64*64);
  band_nr = (long *)malloc(TUBE_BAND*sizeof(long));
  band_lo = malloc(TUBE_BAND*TUBE_RANGES*sizeof(long));
  band_hi = malloc(TUBE_BAND*TUBE_RANGES*sizeof(long));
  if(survivors == NULL || band_nr == NULL || band_lo == NULL || band_hi == NULL)
  { error(7); }
  return;
}

/* the data of the third stage for the row (a, b), see above; the part
   that depends on a alone (its inverse mod p, or that p | a) is kept
   from row to row while a stays */
static _Thread_local long inva3[NUM_PRIMES];  /* 1/a mod p, or 0 when p | a */
static void row3_setup(long a, long b)
{
  long m;
  if(a != row3_a)
  { for(m = 0; m < sieve_primes3 - sieve_primes2; m++)
    { long n = sieve_primes2 + m, ap = a % pr[n];
      inva3[m] = (ap == 0) ? 0 : inverses[pnn[n]][ap];
    }
  }
  for(m = 0; m < sieve_primes3 - sieve_primes2; m++)
  { long n = sieve_primes2 + m, pn = pnn[n], p = pr[n];
    long bp = b % p;
    if(bp < 0) { bp += p; }
    if(inva3[m] != 0)
    { mult3[m] = inva3[m];
      tab3[m] = &is_point_on_j[pn][MOD3(bp * mult3[m], m)][0];
    }
    else if(bp != 0)
    { mult3[m] = inverses[pn][bp];
      tab3[m] = &is_f_square[pn][0];
    }
    else
    { mult3[m] = 0;
      tab3[m] = &all_ones[0];
    }
  }
  row3_a = a; row3_b = b;
  return;
}

/* the third stage for the bit (a, b, c): c = LONG_LENGTH * i + j with
   the word index i, whose residues are r0 + res (r0 the row of res0 for
   the word's position in its chunk, res the chunk's); returns 1 if the
   bit survives every prime */
static inline int check_bit(long a, long b, long c, unsigned short *r0,
                            long *res)
{
  long m, j = c & LONG_MASK;
  if(a != row3_a || b != row3_b) { row3_setup(a, b); }
  for(m = 0; m < sieve_primes3 - sieve_primes2; m++)
  { long n = sieve_primes2 + m, p = pr[n], r = r0[n] + res[n], x;
    if(r >= p) { r -= p; }
    x = MOD3(r * LONG_LENGTH + j, m);    /* c mod p */
    x = MOD3(x * mult3[m], m);          /* c/a or c/b mod p */
    if(!tab3[m][x]) { return(0); }
  }
  return(1);
}

static inline int relprime3(long, long, long);

/* the bits of a word that survived the sieve's tables: the coprimality,
   the third stage, the exact check; returns 1 when one point is enough
   and one was found */
static inline int check_bits(bit_array nums, long a, long b, long i,
                             unsigned short *r0, long *res)
{
  long c;
  /* c will be the coordinate corresponding to the selected bit */
  for(c = i<<LONG_SHIFT; nums; nums >>= 1, c++)
  {/* test one bit */
    if(nums & 1)
    { num_surv2++;
      if(relprime3(a, b, c) && check_bit(a, b, c, r0, res))
      { num_surv3++;
        if(check_one_point(a, b, c) && one_point) { return(1); }
      }
    }
  }
  return(0);
}


/**************************************************************************
 * Some helper functions                                                  *
 **************************************************************************/

/* Determine if n and m are coprime.
 * n, m are >= 0, and n and m are not both even. */
static inline int relprime2(long n, long m)
{
  if(n == 0) { return(m == 1); }
  if(m == 0) { return(n == 1); }
  /* remove pwers of 2 */
  while(!(m & 1)) { m >>= 1; }
  while(!(n & 1)) { n >>= 1; }
  /* successively subtract the smaller from the larger
   * and replace the result by its odd part,
   * until both are equal (to their gcd) */
  while(n != m)
  { if(n > m)
    { n -= m; do { n >>= 1; } while(!(n & 1)); }
    else
    { m -= n; do { m >>= 1; } while(!(m & 1)); }
  }
  return(m == 1);
}

/* Determine if gcd(a,b,c) == 1. Here a is >= 0. */
static inline int relprime3(long a, long b, long c)
{ /* all even --> no */
  if(((a | b | c) & 1) == 0) { return(0); }
  /* make b and c nonnegative */
  if(b < 0) { b = -b; }
  if(c < 0) { c = -c; }
  if(a == 0) { return(relprime2(b, c)); }
  if(b == 0) { return(relprime2(a, c)); }
  if(c == 0) { return(relprime2(a, b)); }
  /* now we have a, b, c > 0 */
  /* remove powers of 2 */
  while(!(a & 1)) { a >>= 1; }
  while(!(b & 1)) { b >>= 1; }
  while(!(c & 1)) { c >>= 1; }
  /* successively subtract the smaller of a, b from the larger
   * and replace the result by its odd part,
   * until both are equal (to their gcd) */
  while(a != b)
  { if(a > b)
    { a -= b; do { a >>= 1; } while(!(a & 1)); }
    else
    { b -= a; do { b >>= 1; } while(!(b & 1)); }
  }
  /* shortcut return if gcd(a,b) == 1 */
  if(a == 1) { return(1); }
  /* repeat procedure with a and c */
  while(a != c)
  { if(a > c)
    { a -= c; do { a >>= 1; } while(!(a & 1)); }
    else
    { c -= a; do { c >>= 1; } while(!(c & 1)); }
  }
  return(a == 1);
}

/* Test a surviving word of the chunk against the second-stage primes and
   hand the surviving bits to the exact check.  i is the index of the word,
   k = i - w_low its position in the chunk, res the residues of w_low. */
static inline int check_point(bit_array nums, long a, long b, long i,
                              long k, long *res)
{
  unsigned short *r0 = &res0[k*npr];
  long n;
  num_surv1++;
  for(n = sieve_primes1; n < sieve_primes2 && nums; n++)
  { long r = r0[n] + res[n];
    if(r >= pr[n]) { r -= pr[n]; }
    nums &= rowptr[n][r];
  }
  if(nums) { return(check_bits(nums, a, b, i, r0, res)); }
  return(0);
}


/**************************************************************************
 * The analysis of the plane of (b, c): the tube and the real region      *
 **************************************************************************/

/* For a fixed first coordinate a, the third coordinates c of a row (a, b)
   that can carry a point are limited by two conditions.

   The real region: for a > 0, when the Mumford roots x, u of a t^2 - b t
   + c are real (k2 = b^2 - 4 a c > 0; k2 = 0 is the double points, found
   elsewhere), the two points of the divisor are real, so f(x) >= 0 and
   f(u) >= 0 (for a = 0 the one point is (c/b, y), so f(c/b) >= 0); when
   they are complex there is nothing to ask, and a real fourth coordinate
   exists in any case.  The real roots of f are isolated once (see
   j-points.c), and x = (b + sqrt(k2)) / (2 a) and u = (b - sqrt(k2)) /
   (2 a) must avoid the intervals where f is negative.  Along a row,
   x runs down from infinity to b / (2 a) and u up from minus infinity to
   it as c grows, so each such interval excludes at most one interval of
   c on each branch, c = b x - a x^2 at its ends: row_region computes the
   words of the row that remain, some tens of operations, no cells.
   This cuts the box of an -a run to half or nine tenths.

   The tube: when the fourth coordinate d is bounded by W (a plain run:
   W = h; -w: W as given; not with -a), a root of Q(d) = k2 d^2 + k1 d +
   k0 lies in [-W, W]: either Q(W) and Q(-W) differ in sign (or one is
   0), or both roots lie inside, which needs |k1| <= 2 W |k2| and Q(W),
   Q(-W) on the side of k2.  For W = h this holds only in a thin tube
   around the plane section d = 0 of the Kummer surface (or near the
   origin): at h = 2000 on a few of the 63 words of a row, in many rows
   on none.  So the plane of (b, c) is cut into cells of rows times
   words, and a cell is dropped when the bounds show that the condition
   fails everywhere in it, or taken whole when they show that it holds
   everywhere: Q(W), Q(-W), k1, k2 are polynomials in (b, c), and with
   the Taylor expansion sum p_jk s^j t^k at the centre (bm, cm), |s| <=
   rb, |t| <= rc, the value keeps the sign of p_00 throughout when
   |p_00| > sum' |p_jk| rb^j rc^k; the rounding errors of the
   coefficients, which may cancel, are covered by TUBE_EPS times the
   absolute polynomial at (|bm| + rb, |cm| + rc), so that no word that
   can hold a point is lost.  An undecided cell is halved, in rows or in
   words, until it is at most TUBE_ROWS rows by one word; that word is
   then collected for each of its rows, as the words of a cell taken
   whole are.

   The rows of a band are then sieved in order, in one of two ways (the
   model chooses, see choose_primes): the tube -- the words of the cells,
   cut to the real region, each tested against all the table primes one
   by one, and its surviving bits against the third stage and the exact
   check -- or the region -- the passes of the box sieve run over the
   ranges of words of the real region. */

#define TUBE_EPS 1.5e-5   /* each coefficient with a margin of 1e-6 */
#define TUBE_ROWS 64      /* rows of a leaf cell (16 to 128 measured much the same) */
#define L 5               /* the row length of a coefficient array (degree at most 4) */
#define NC (L*L)
#define P(j, k) [L*(j) + (k)]   /* the coefficient of b^j c^k (or s^j t^k) */

/* for the current a, as polynomials in (b, c): k0 (degree 4), k1 (3) and
   k2 (2), and the absolute values of the coefficients of |k0| + W |k1| +
   W^2 |k2|, the scale of the errors */
static _Thread_local double K0[NC], K1[NC], K2[NC], KA[NC];
static _Thread_local long tube_a;   /* the current a */
static long tube_ncells;

static void tube_a_init(long a)
{
  double A = (double)a, A2 = A*A, A3 = A2*A, A4 = A3*A;
  long i;
  for(i = 0; i < NC; i++) { K0[i] = 0.0; K1[i] = 0.0; K2[i] = 0.0; }
  K0 P(0,0) = kd400*A4;
  K0 P(1,0) = kd310*A3; K0 P(0,1) = kd301*A3;
  K0 P(2,0) = kd220*A2; K0 P(1,1) = kd211*A2; K0 P(0,2) = kd202*A2;
  K0 P(3,0) = kd130*A;  K0 P(2,1) = kd121*A;  K0 P(1,2) = kd112*A;
  K0 P(0,3) = kd103*A;
  K0 P(4,0) = kd040; K0 P(3,1) = kd031; K0 P(2,2) = kd022; K0 P(1,3) = kd013;
  K0 P(0,4) = kd004;
  K1 P(0,0) = -4.0*A3*fd[0]; K1 P(1,0) = -2.0*A2*fd[1]; K1 P(0,1) = -4.0*A2*fd[2];
  K1 P(1,1) = -2.0*A*fd[3];  K1 P(0,2) = -4.0*A*fd[4];  K1 P(1,2) = -2.0*fd[5];
  K1 P(0,3) = -4.0*fd[6];
  K2 P(2,0) = 1.0; K2 P(0,1) = -4.0*A;
  for(i = 0; i < NC; i++)
  { KA[i] = fabs(K0[i]) + th*fabs(K1[i]) + th2*fabs(K2[i]); }
  tube_a = a;
  return;
}

/* The real region of the row (a, b): the ranges of words [lo[i], hi[i])
   of c that it meets, in order, their number returned.  On the branch
   x >= b / (2 a), c = b x - a x^2 falls as x grows, so an interval (nl,
   nh) where f is negative excludes the c strictly between C(nh) and
   C(nl) (from C(max(nl, m)) with m = b / (2 a), and everything below
   when nh is infinite); on the branch u <= m, c rises, and the interval
   excludes the c between C(nl) and C(min(nh, m)).  For a = 0 the root
   is c / b (b > 0), so the interval excludes the c between b nl and b
   nh.  The intervals are the shrunk ones (f is negative on them for
   sure), each further shrunk by a margin for the rounding, so that a c
   is excluded only when it is inside for certain; a word is excluded
   when all its c are. */
#define ROW_RANGES 12   /* at most 2 * 4 exclusions leave 9 ranges */
static long row_region(long a, long b, long *lo, long *hi)
{
  double xlo[8], xhi[8];   /* the excluded intervals of c, open */
  long n = 0, i, j, nr = 0, w = w_low_all;
  if(num_neg == 0)
  { lo[0] = w_low_all; hi[0] = w_high_all; return(1); }
  if(a > 0)
  { double A = (double)a, B = (double)b, m = 0.5*B/A;
    for(i = 0; i < num_neg; i++)
    { double nl = neg_shrunk_lo[i], nh = neg_shrunk_hi[i], l, h;
      /* the branch x >= m */
      l = (nl > m) ? nl : m;
      if(nh > l)
      { xlo[n] = (nh == HUGE_VAL) ? -HUGE_VAL : B*nh - A*nh*nh;
        xhi[n] = B*l - A*l*l;
        n++;
      }
      /* the branch u <= m */
      h = (nh < m) ? nh : m;
      if(h > nl)
      { xlo[n] = (nl == -HUGE_VAL) ? -HUGE_VAL : B*nl - A*nl*nl;
        xhi[n] = B*h - A*h*h;
        n++;
      }
    }
  }
  else
  { double B = (double)b;
    for(i = 0; i < num_neg; i++)
    { double nl = neg_shrunk_lo[i], nh = neg_shrunk_hi[i];
      xlo[n] = (nl == -HUGE_VAL) ? -HUGE_VAL : B*nl;
      xhi[n] = (nh == HUGE_VAL) ? HUGE_VAL : B*nh;
      n++;
    }
  }
  /* the margins (an infinite end stays), then the intervals in order of
     their lower ends */
  for(i = 0; i < n; i++)
  { if(xlo[i] > -HUGE_VAL) { xlo[i] += 1.0e-9*(fabs(xlo[i]) + 1.0) + 1.0; }
    if(xhi[i] < HUGE_VAL) { xhi[i] -= 1.0e-9*(fabs(xhi[i]) + 1.0) + 1.0; }
  }
  for(i = 1; i < n; i++)
  { double l = xlo[i], h = xhi[i];
    for(j = i; j > 0 && xlo[j-1] > l; j--) { xlo[j] = xlo[j-1]; xhi[j] = xhi[j-1]; }
    xlo[j] = l; xhi[j] = h;
  }
  /* the words left: w is the first word not yet placed */
  for(i = 0; i < n; i++)
  { long wl, wh;   /* the words wholly inside (xlo, xhi): [wl, wh) */
    if(xhi[i] - xlo[i] < (double)LONG_LENGTH) { continue; }
    /* the first c > xlo is at most (long) xlo + 2, the last c < xhi at
       least (long) xhi - 2 (the cast truncates towards 0; a c more or
       less at the ends is kept, which is safe); the words wholly inside
       those c, by the shifts, which round down */
    { if(xlo[i] < -3.0e9) { wl = w_low_all; }
      else { long c = (long)xlo[i] + 2; wl = (c + LONG_MASK) >> LONG_SHIFT; }
      if(xhi[i] > 3.0e9) { wh = w_high_all; }
      else { long c = (long)xhi[i] - 2; wh = (c + 1) >> LONG_SHIFT; }
    }
    if(wl < w) { wl = w; }
    if(wh > w_high_all) { wh = w_high_all; }
    if(wh <= wl) { continue; }
    if(wl > w) { lo[nr] = w; hi[nr] = wl; nr++; }
    w = wh;
  }
  if(w < w_high_all) { lo[nr] = w; hi[nr] = w_high_all; nr++; }
  return(nr);
}

/* the Taylor shift of the polynomial p[0], p[s], ..., p[deg*s] (the
   coefficients deg apart by the stride s) to the point x, in place */
static inline void shift1(double *p, long deg, long s, double x)
{
  long i, j;
  for(i = 0; i < deg; i++)
  { for(j = deg - 1; j >= i; j--) { p[j*s] += x * p[(j+1)*s]; } }
  return;
}

/* the shift of P[j][k] (b^j c^k, j + k <= deg) to (bm, cm) */
static inline void shift2(double *Q, long deg, double bm, double cm)
{
  long j, k;
  for(k = 0; k <= deg; k++) { shift1(Q + k, deg - k, L, bm); }
  for(j = 0; j <= deg; j++) { shift1(Q + L*j, deg - j, 1, cm); }
  return;
}

/* the sum over (j, k) != (0, 0), j + k <= deg, of |Q[j][k]| rp[j][k], rp
   holding the products rb^j rc^k */
static inline double tube_sum(double *Q, long deg, double *rp)
{
  double sum = 0.0;
  long j, k;
  for(j = 0; j <= deg; j++)
    for(k = (j == 0) ? 1 : 0; k <= deg - j; k++)
    { sum += fabs(Q[L*j + k]) * rp[L*j + k]; }
  return(sum);
}

/* the value of the polynomial with the coefficients QA at (x, y) */
static inline double tube_eval(double *QA, long deg, double x, double y)
{
  double sum = 0.0, xj = 1.0;
  long j, k;
  for(j = 0; j <= deg; j++)
  { double yk = 1.0;
    for(k = 0; k <= deg - j; k++) { sum += QA[L*j + k] * xj * yk; yk *= y; }
    xj *= x;
  }
  return(sum);
}

/* the sign of a polynomial throughout the cell: 1 or -1 when the bounds
   show it, else 0; *hi bounds its absolute value, err is the margin for
   the rounding errors */
static inline int poly_sign(double *Q, long deg, double *rp, double err,
                            double *lo, double *hi)
{
  double s = tube_sum(Q, deg, rp) + err;
  *lo = fabs(Q[0]) - s; *hi = fabs(Q[0]) + s;
  if(*lo > 0.0) { return((Q[0] > 0.0) ? 1 : -1); }
  return(0);
}

/* the words [i_lo, i_hi) can hold a point of the row b: collected into
   the row's ranges (merged with the last one when they touch; beyond
   TUBE_RANGES ranges the row keeps their hull) */
static inline void tube_emit(long b, long i_lo, long i_hi)
{
  long row = b - band_b0, nr = band_nr[row];
  if(nr > 0 && band_hi[row][nr-1] >= i_lo)
  { if(band_hi[row][nr-1] < i_hi) { band_hi[row][nr-1] = i_hi; }
    return;
  }
  if(nr == TUBE_RANGES)
  { band_hi[row][0] = i_hi; band_nr[row] = 1; return; }
  band_lo[row][nr] = i_lo; band_hi[row][nr] = i_hi;
  band_nr[row] = nr + 1;
  return;
}

/* the cell of the rows b_lo..b_hi and the words [i_lo, i_hi) */
static void tube_cell(long b_lo, long b_hi, long i_lo, long i_hi)
{
  static _Thread_local double rp[NC], last_rb = -1.0, last_rc = -1.0;
  long rows = b_hi - b_lo + 1, words = i_hi - i_lo;
  double p0[NC], p1[NC], p2[NC], q[NC];
  double bm = 0.5*(double)(b_lo + b_hi), rb = 0.5*(double)(b_hi - b_lo);
  double cm = 32.0*(double)(i_lo + i_hi) - 0.5;
  double rc = 32.0*(double)(i_hi - i_lo) - 0.5;
  double fb = fabs(bm) + rb, fc = fabs(cm) + rc;   /* the far corner */
  double ma, lo2, hi2, lop, hip, lom, him, err, m[NC];
  int s2, sp, sm, inside;
  long j, k;
  tube_ncells++;
  /* the products rb^j rc^k, the same for cells of one size */
  if(rb != last_rb || rc != last_rc)
  { double rbj;
    for(j = 0, rbj = 1.0; j < L; j++, rbj *= rb)
    { double rck = rbj;
      for(k = 0; k < L; k++, rck *= rc) { rp P(j,k) = rck; }
    }
    last_rb = rb; last_rc = rc;
  }
  for(j = 0; j < NC; j++) { p0[j] = K0[j]; p1[j] = K1[j]; p2[j] = K2[j]; }
  shift2(p0, 4, bm, cm); shift2(p1, 3, bm, cm); shift2(p2, 2, bm, cm);
  /* the scale of the rounding errors: |k0| + W |k1| + W^2 |k2| over the
     cell is bounded by its absolute polynomial at the far corner, and
     that bounds each of them */
  ma = tube_eval(KA, 4, fb, fc);
  err = TUBE_EPS*ma;
  s2 = poly_sign(p2, 2, rp, err/th2, &lo2, &hi2);   /* the sign of k2 */
  for(j = 0; j < NC; j++)
  { q[j] = p0[j] + th*p1[j] + th2*p2[j];    /* Q(W) */
    m[j] = p0[j] - th*p1[j] + th2*p2[j];    /* Q(-W) */
  }
  sp = poly_sign(q, 4, rp, err, &lop, &hip);
  sm = poly_sign(m, 4, rp, err, &lom, &him);
  /* a sign change of Q between -W and W throughout: a root inside */
  inside = (sp != 0 && sm != 0 && sp != sm);
  if(sp != 0 && sm != 0 && sp == sm)
  { /* the same sign throughout: both roots inside, or none; out when k2
       has the other sign, or |k1| > 2 W |k2|, throughout */
    double lo1, hi1;
    if(s2 != 0 && s2 != sp) { return; }
    poly_sign(p1, 3, rp, err/th, &lo1, &hi1);
    if(lo1 > 2.0*th*hi2) { return; }
  }
  /* the condition holds throughout: the cell whole */
  if(inside)
  { long b;
    for(b = b_lo; b <= b_hi; b++) { tube_emit(b, i_lo, i_hi); }
    return;
  }
  /* halve the cell: in words while it has more than one and is not much
     longer in rows (TUBE_ROWS rows to a word), else in rows while it has
     more than TUBE_ROWS; a leaf's word goes to each of its rows */
  if(words > 1 && rows <= TUBE_ROWS*words)
  { long mid = i_lo + words/2;
    tube_cell(b_lo, b_hi, i_lo, mid);
    tube_cell(b_lo, b_hi, mid, i_hi);
    return;
  }
  if(rows > TUBE_ROWS)
  { long mid = b_lo + (rows - 1)/2;
    tube_cell(b_lo, mid, i_lo, i_hi);
    tube_cell(mid + 1, b_hi, i_lo, i_hi);
    return;
  }
  { long b;
    for(b = b_lo; b <= b_hi; b++) { tube_emit(b, i_lo, i_lo + 1); }
  }
  return;
}

/* the ranges of the row cut to its real region: those of the band's row
   (from the cells) intersected with those of row_region, into rlo, rhi;
   their number returned */
static long row_cut(long row, long nreg, long *reglo, long *reghi, long *rlo, long *rhi)
{
  long nr = band_nr[row], i = 0, j = 0, n = 0;
  while(i < nr && j < nreg)
  { long lo = (band_lo[row][i] > reglo[j]) ? band_lo[row][i] : reglo[j];
    long hi = (band_hi[row][i] < reghi[j]) ? band_hi[row][i] : reghi[j];
    if(lo < hi) { rlo[n] = lo; rhi[n] = hi; n++; }
    if(band_hi[row][i] < reghi[j]) { i++; } else { j++; }
  }
  return(n);
}

/* for the cost model, from bands of 64 rows of a few a spread over the
   box: the words per row the tube leaves, as they are and cut to the
   real region, the cells it visits per row and the fraction of the rows
   with a word; the words per row of the real region alone, in how many
   ranges, and the fraction of the rows that meet it */
void tube_sample(double *tube_words, double *tube_words_cut, double *tube_cells, double *tube_rows,
                 double *reg_words, double *reg_ranges, double *reg_rows)
{
  long k, rows = 0, left = 0, tleft = 0, w_low = (-height)>>LONG_SHIFT, w_high = (height>>LONG_SHIFT) + 1;
  double twords = 0.0, tcut = 0.0, cells = 0.0, rwords = 0.0, ranges = 0.0;
  th = (double)dbound; th2 = th*th;
  w_low_all = w_low; w_high_all = w_high;
  if(band_nr == NULL) { init_thread_sieve(); }
  for(k = 0; k < 24; k++)
  { long a = (k * 7919L) % (height + 1), b_lo, b_hi, b;
    if(degree == 5 && mpz_cmp_si(&coeffs[5], 1) == 0)
    { a = a % (long)floor(sqrt((double)height) + 1); a = a*a; }
    b_lo = ((k * 104729L) % (2*height + 1)) - height;
    if(a == 0 && b_lo < 1) { b_lo = 1; }
    b_hi = b_lo + 63;
    if(b_hi > height) { b_hi = height; }
    tube_a_init(a);
    tube_ncells = 0;
    band_b0 = b_lo;
    for(b = b_lo; b <= b_hi; b++) { band_nr[b - b_lo] = 0; }
    if(dbounded) { tube_cell(b_lo, b_hi, w_low, w_high); }
    for(b = b_lo; b <= b_hi; b++)
    { long reglo[ROW_RANGES], reghi[ROW_RANGES], rlo[ROW_RANGES + TUBE_RANGES], rhi[ROW_RANGES + TUBE_RANGES];
      long nreg = row_region(a, b, reglo, reghi), n, r;
      for(r = 0; r < nreg; r++) { rwords += (double)(reghi[r] - reglo[r]); }
      ranges += (double)nreg;
      if(nreg > 0) { left++; }
      n = row_cut(b - b_lo, nreg, reglo, reghi, rlo, rhi);
      for(r = 0; r < n; r++) { tcut += (double)(rhi[r] - rlo[r]); }
      for(r = 0; r < band_nr[b - b_lo]; r++) { twords += (double)(band_hi[b - b_lo][r] - band_lo[b - b_lo][r]); }
      if(band_nr[b - b_lo] > 0) { tleft++; }
    }
    cells += (double)tube_ncells;
    rows += b_hi - b_lo + 1;
  }
  *tube_words = (rows > 0) ? twords / (double)rows : 1.0;
  *tube_words_cut = (rows > 0) ? tcut / (double)rows : 1.0;
  *tube_cells = (rows > 0) ? cells / (double)rows : 1.0;
  *tube_rows = (rows > 0) ? (double)tleft / (double)rows : 1.0;
  *reg_words = (rows > 0) ? rwords / (double)rows : 1.0;
  *reg_ranges = (rows > 0) ? ranges / (double)rows : 1.0;
  *reg_rows = (rows > 0) ? (double)left / (double)rows : 1.0;
  return;
}

/* the word i of the row (a, b) whose table rows are set: the mask, the
   ends of the row and the excluded c0, then the table primes one by
   one, then the bits; returns 1 when one point is enough and one was
   found */
static _Thread_local long cur_chunk = 0;
static int tube_word(long a, long b, long i, bit_array mask, long c0)
{
  chunk_spec *ch;
  long k, n;
  bit_array nums = mask;
  unsigned short *r0;
  while(i >= chunks[cur_chunk].w_high) { cur_chunk++; }
  while(i < chunks[cur_chunk].w_low) { cur_chunk--; }
  ch = &chunks[cur_chunk];
  k = i - ch->w_low;
  r0 = &res0[k*npr];
  if(i == w_low_all) { nums &= begmask; }
  if(i == w_high_all - 1) { nums &= endmask; }
  if((c0>>LONG_SHIFT) == i) { nums &= ~(1UL<<(c0 & LONG_MASK)); }
  /* all the table primes, without a test in between: the loads then
     overlap, and no branch is mispredicted at an exit (a test halfway
     was measured and lost) */
  for(n = 0; n < sieve_primes2; n++)
  { long r = r0[n] + ch->res[n];
    r -= (r >= pr[n]) ? pr[n] : 0;
    nums &= rowptr[n][r];
  }
  if(nums) { return(check_bits(nums, a, b, i, r0, ch->res)); }
  return(0);
}

static void row_setup(long a, long b);
static int sift_range(long, long, long, long, long, bit_array, long);

/* the rows of a, in the tube or the region: the multiples of m among the
   b (see row_step), from b_first on; band by band, the words that can
   hold a point collected, then the rows sieved in order, each collected
   word on its own (the tube, its words cut to the real region when the
   model says so) or each range by the passes (the region) */
int sift_bands(long a, long b_first, long m)
{
  long b_lo;
  int tube = (sieve_mode == 1);
  tube_a_init(a);
  for(b_lo = b_first; b_lo <= height; b_lo += TUBE_BAND)
  { long b_hi = b_lo + TUBE_BAND - 1, b, k;
    if(b_hi > height) { b_hi = height; }
    band_b0 = b_lo;
    if(tube)
    { for(k = 0; k <= b_hi - b_lo; k++) { band_nr[k] = 0; }
      tube_cell(b_lo, b_hi, w_low_all, w_high_all);
    }
    for(b = b_lo; b <= b_hi; b += m)
    { long reglo[ROW_RANGES], reghi[ROW_RANGES], clo[ROW_RANGES + TUBE_RANGES], chi[ROW_RANGES + TUBE_RANGES];
      long *rlo, *rhi, nr, c0, i;
      bit_array mask;
      if(tube && band_nr[b - b_lo] == 0) { continue; }
      if(use2 && a != 0 && !alive2[a & 63][b & 63]) { continue; }
      if(tube && !tube_cut)
      { nr = band_nr[b - b_lo]; rlo = band_lo[b - b_lo]; rhi = band_hi[b - b_lo]; }
      else if(tube)
      { nr = row_region(a, b, reglo, reghi);
        nr = row_cut(b - b_lo, nr, reglo, reghi, clo, chi);
        rlo = clo; rhi = chi;
      }
      else
      { nr = row_region(a, b, reglo, reghi); rlo = reglo; rhi = reghi; }
      if(nr == 0) { continue; }
      row_setup(a, b);
      c0 = height + LONG_LENGTH;
      if(a != 0 && (b*b)%(4*a) == 0) { c0 = (b*b)/(4*a); }
      mask = ~zero;
      if(use2 && a != 0) { mask = mask2[a & 63][b & 63]; }
      else if(!((a|b)&1)) { mask = HALF_MASK; }
      for(k = 0; k < nr; k++)
      { long lo = rlo[k], hi = rhi[k];
        if(tube)
        { num_surv1 += hi - lo;
          for(i = lo; i < hi; i++)
          { if(tube_word(a, b, i, mask, c0) && one_point) { return(1); } }
        }
        else
        { long c;
          while(lo >= chunks[cur_chunk].w_high) { cur_chunk++; }
          while(lo < chunks[cur_chunk].w_low) { cur_chunk--; }
          for(c = cur_chunk; c < num_chunks && chunks[c].w_low < hi; c++)
          { long l = (lo > chunks[c].w_low) ? lo : chunks[c].w_low;
            long h = (hi < chunks[c].w_high) ? hi : chunks[c].w_high;
            if(sift_range(a, b, c, l, h, mask, c0) && one_point) { return(1); }
          }
        }
      }
    }
  }
  return(0);
}

/**************************************************************************
 * The sieving procedure itself                                           *
 **************************************************************************/

/*------------------------------------------------------------------------+
 | The following procedure is the heart of the matter.                    |
 | The overall speed of the program highly depends on the quality         |
 | of the code the compiler produces for the innermost loops in           |
 | sift_range.                                                            |
 +------------------------------------------------------------------------*/

/* The walk of the n-th prime over the words [off, off + len) of the
   chunk: OP(word of the bit array, table word) for every word.  The
   scalar walk (VW = 1) reads the table row from its word f1 up to its
   end, then in full periods of p words, then the tail; the vector walk
   reads VW words at a time from the offset f1 on, wrapping at the
   period.  For the whole chunk the walk is the precomputed one; for a
   range it starts off words further on, and its parts are computed
   from the residue tables, without a division. */
static inline void walk_parts(long n, long p, walk_spec *w, long off, long len,
                              int whole, long *f1, long *hd, long *np, long *tl)
{
  if(whole) { *f1 = w->first; *hd = w->head; *np = w->nper; *tl = w->tail; return; }
  if(VW == 1)
  { long f = w->first + res0[off*npr + n], h, rem, r;
    if(f > p) { f -= p; }
    h = p - f;
    if(h > len) { h = len; }
    rem = len - h; r = res0[rem*npr + n];
    *f1 = f; *hd = h; *tl = r;
    *np = (long)((double)(rem - r) * invp[n] + 0.5);   /* exact: p | rem - r */
  }
  else
  { long t = w->first + ((period[n] == p) ? res0[off*npr + n] : off % period[n]);
    if(t >= period[n]) { t -= period[n]; }
    *f1 = t; *hd = *np = *tl = 0;
  }
  return;
}
#if VW == 1
#define WALK(OP) \
  { bit_array *surv = surv0; \
    bit_array *siv1 = &rowptr[n][f1]; \
    long j; \
    for(j = hd; j; j--) { OP(*surv, *siv1); surv++; siv1++; } \
    /* now siv1 points at the end of the table row, whose word p repeats \
       its word 0, if there is anything left to do */ \
    for(j = np; j; j--) \
    { bit_array *siv0 = siv1; \
      siv1 -= p; \
      do { OP(*surv, *siv1); surv++; siv1++; } while(siv1 != siv0); \
    } \
    if((j = tl)) \
    { siv1 -= p; \
      for( ; j; j--) { OP(*surv, *siv1); surv++; siv1++; } \
    } \
  }
#define AND_OP(d, s) ((d) &= (s))
#define FILL_OP(d, s) ((d) = mask & (s))
#else
#define WALK(OP) \
  { vecu *surv = (vecu *)surv0; \
    bit_array *row = rowptr[n]; \
    long t = f1, Pd = period[n], j; \
    (void)p; (void)hd; (void)np; (void)tl; \
    for(j = nv; j; j--) \
    { vec s = *(vecu *)(row + t); \
      OP(*surv, s); surv++; \
      t += VW; t -= (t >= Pd) ? Pd : 0; \
    } \
  }
#define AND_OP(d, s) ((d) &= (s))
#define FILL_OP(d, s) ((d) = vmask & (s))
#endif

/* Sieve the words [lo, hi) of the k-th chunk of the row (a, b) -- the
   whole chunk in the box: the bit array set from the mask (all bits, or
   the odd c when a and b are even) and the table of the first prime in
   one pass, the ends of the row and the excluded c0 taken out, the other
   first-stage primes ANDed in, the surviving words passed to the second
   stage. */
static int sift_range(long a, long b, long k, long lo, long hi, bit_array mask, long c0)
{
  chunk_spec *ch = &chunks[k];
  long n, range = ch->w_high - ch->w_low, off = lo - ch->w_low, len = hi - lo;
  long f1, hd, np, tl, cw = (c0>>LONG_SHIFT) - ch->w_low;
  int whole;
  bit_array *surv0;
#if VW > 1
  vec vmask;
  long nv;
  for(n = 0; n < VW; n++) { vmask[n] = mask; }
  /* whole vectors: the range widened to them (the words beyond the chunk
     cleared before the scan) */
  { long o = off & ~(VW - 1);
    len += off - o; off = o;
    len = (len + VW - 1) & ~(VW - 1);
    if(off + len > ch->nvec*VW) { len = ch->nvec*VW - off; }
    nv = len/VW;
  }
#endif
  whole = (off == 0 && len >= range);
  surv0 = survivors + off;
  /* the fill, merged with the first prime's pass if there is one */
  if(sieve_primes1 == 0)
  { bit_array *surv = surv0;
    long i;
    for(i = len; i; i--) { *surv++ = mask; }
    n = 0;
  }
  else
  { long p = pr[0];
    walk_spec *w = &ch->walk[0];
    n = 0;
    walk_parts(0, p, w, off, len, whole, &f1, &hd, &np, &tl);
    WALK(FILL_OP);
    n = 1;
  }
  if(k == 0 && off == 0) { survivors[0] &= begmask; }
  if(k == num_chunks - 1 && off + len >= range) { survivors[range-1] &= endmask; }
  if(off <= cw && cw < off + len) { survivors[cw] &= ~(1UL<<(c0 & LONG_MASK)); }
#if (DEBUG >= 1)
  printf("\n sift_range(%ld, %ld)\n", lo, hi);
#endif
  /* now do the sieving (fast!) */
  for( ; n < sieve_primes1; n++)
  { long p = pr[n];
    walk_spec *w = &ch->walk[n];
    walk_parts(n, p, w, off, len, whole, &f1, &hd, &np, &tl);
    WALK(AND_OP);
  }
#if (DEBUG >= 3)
  { long i;
    for(i = off; i < off + len; i++) { printf(" %8.8lx",survivors[i]); }
    printf("\n");
  }
#endif
  /* Check the points that have survived the sieve if they really are points */
#if VW == 1
  { bit_array *surv = surv0;
    bit_array nums;
    long i;
    for(i = off; i < off + len; i++)
    { if((nums = *surv++))
      { if(check_point(nums, a, b, ch->w_low + i, i, ch->res) && one_point)
        { return(1); }
  } } }
#else
  { /* the padding beyond the chunk cleared, then a vector at a time */
    vecu *sv = (vecu *)surv0;
    long i, j;
    for(i = range; i < off + len; i++) { survivors[i] = 0; }
    for(j = 0; j < nv; j++)
    { vec v = sv[j];
      bit_array any = 0;
      for(i = 0; i < VW; i++) { any |= v[i]; }
      if(any)
      { for(i = 0; i < VW; i++)
        { bit_array nums = v[i];
          long pos = off + j*VW + i;
          if(nums && check_point(nums, a, b, ch->w_low + pos, pos, ch->res)
             && one_point)
          { return(1); }
        }
      }
    }
  }
#endif
  return(0);
}

/* the residues of b modulo the sieving primes and the table rows of
   (a, b): b steps by a constant from one row to the next, so they are
   kept and advanced, and computed by division only when a changes or
   the caller jumps back; a step of any size is reduced mod p by the
   reciprocal in double precision (exact: the error of d/p is far below
   1/p) */
static void row_setup(long a, long b)
{
  static _Thread_local long last_a = -1, last_b = 0, last_d = 0;
  static _Thread_local long bres[NUM_PRIMES], dres[NUM_PRIMES], drow[NUM_PRIMES];
  long n, d = b - last_b;
  if(a == last_a && d > 0 && d < 3)
  { /* a step below every prime: the residues step by d */
    for(n = 0; n < sieve_primes2; n++)
    { bres[n] += d;
      rowptr[n] += d*rowlen[n];
      if(bres[n] >= pr[n]) { bres[n] -= pr[n]; rowptr[n] -= blocklen[n]; }
    }
  }
  else if(a == last_a && d > 0)
  { if(d != last_d)
    { /* the step d and its residues, once per change of the step */
      for(n = 0; n < sieve_primes2; n++)
      { long p = pr[n], r = d;
        if(r >= p)
        { r = d - p * (long)((double)d * invp[n]);
          if(r >= p) { r -= p; }
          if(r < 0) { r += p; }
        }
        dres[n] = r;
        drow[n] = r*rowlen[n];
      }
      last_d = d;
    }
    for(n = 0; n < sieve_primes2; n++)
    { bres[n] += dres[n];
      rowptr[n] += drow[n];
      if(bres[n] >= pr[n]) { bres[n] -= pr[n]; rowptr[n] -= blocklen[n]; }
    }
  }
  else
  { for(n = 0; n < sieve_primes2; n++)
    { long p = pr[n], ap = a % p, bp = b % p;
      if(bp < 0) { bp += p; }
      bres[n] = bp;
      rowptr[n] = &sieve_tab[n][(ap*p + bp)*rowlen[n]];
    }
    last_d = 0;
  }
  last_a = a; last_b = b;
  return;
}

int sift(long a, long b)
/* the row (a, b) in the box: every chunk whole */
{
  long k;
  /* c0 is value of c that has to be excluded */
  long c0 = height + LONG_LENGTH;
  bit_array mask;
  if(a != 0 && (b*b)%(4*a) == 0) { c0 = (b*b)/(4*a); }
#if (DEBUG >= 1)
  printf("\n sift(a = %ld, b = %ld)\n", a, b);
#endif
  row_setup(a, b);
  /* the fill: the condition at 2 for the row, or just "not all even" */
  mask = ~zero;
  if(use2 && a != 0) { mask = mask2[a & 63][b & 63]; }
  else if(!((a|b)&1)) { mask = HALF_MASK; }
  /* Now the chunks of longwords (= bit_arrays) */
  for(k = 0; k < num_chunks; k++)
  { if(sift_range(a, b, k, chunks[k].w_low, chunks[k].w_high, mask, c0) && one_point)
    { return(1); }
  }
  return(0);
}
