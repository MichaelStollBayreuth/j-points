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
{ long first;  /* the table word the chunk begins with, in [1, p] */
  long head;   /* the words up to the end of the table row (or the chunk) */
  long nper;   /* the full periods of p words after the head */
  long tail;   /* the words after the last full period */
} walk_spec;

typedef struct
{ long w_low, w_high;  /* the words of the chunk: w_low <= i < w_high */
  walk_spec *walk;     /* one per first-stage prime */
  long *res;           /* per second-stage prime: w_low mod p */
} chunk_spec;

static chunk_spec *chunks;
static long num_chunks;

static long npr;    /* the number of sieving primes, all stages */
static unsigned short *res0;
/* res0[k*npr + n] = k mod pr[n] for 0 <= k < array_size: with the residue
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
static unsigned char *tab3[NUM_PRIMES];  /* per third-stage prime: the row */
static long mult3[NUM_PRIMES];           /* ... and the factor */
static unsigned long recip3[NUM_PRIMES];
/* recip3[m] = ceil(2^32 / p) for the m-th third-stage prime p: for
   0 <= x < 2^16, (x * recip3[m]) >> 32 is exactly x / p, since x / p is
   either an integer or at least 1/p away from one, and the error of the
   product is below 2^-16.  So x mod p costs two multiplications. */
#define MOD3(x, m) ((x) - pr[sieve_primes2 + (m)] * (long)(((x) * recip3[m]) >> 32))
static long row3_a = -1, row3_b = 0;     /* the row the data is for */

/* the sieving primes, the lengths of their table rows and of the block
   of rows of one residue class of a, and the table rows of the current
   (a, b): rowptr[n] for the n-th prime */
static long pr[NUM_PRIMES], rowlen[NUM_PRIMES], blocklen[NUM_PRIMES];
static bit_array *rowptr[NUM_PRIMES];
static long w_low_all, w_high_all;   /* the words of a whole row */


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
    rowlen[n] = pr[n] + 1;
    blocklen[n] = pr[n]*rowlen[n];
  }
  for(n = 0; n < MAX_PRIME_EVEN; n++) { all_ones[n] = 1; }
  for(n = 0; n < sieve_primes3 - sieve_primes2; n++)
  { unsigned long p = pr[sieve_primes2 + n];
    recip3[n] = ((1UL << 32) + p - 1) / p;
  }
  num_chunks = CEIL(w_high - w_low, array_size);
  survivors = (bit_array *)malloc(array_size*sizeof(bit_array));
  chunks = (chunk_spec *)malloc(num_chunks*sizeof(chunk_spec));
  walk = (walk_spec *)malloc((num_chunks*np1 + 1)*sizeof(walk_spec));
  res = (long *)malloc((num_chunks*npr + 1)*sizeof(long));
  res0 = (unsigned short *)malloc((array_size*npr + 1)*sizeof(unsigned short));
  if(survivors == NULL || chunks == NULL || walk == NULL || res == NULL
       || res0 == NULL)
  { error(7); }
  for(k = 0; k < num_chunks; k++)
  { chunk_spec *ch = &chunks[k];
    long range;
    ch->w_low = w_low + k*array_size;
    ch->w_high = ch->w_low + array_size;
    if(ch->w_high > w_high) { ch->w_high = w_high; }
    range = ch->w_high - ch->w_low;
    ch->walk = &walk[k*np1];
    ch->res = &res[k*npr];
    for(n = 0; n < np1; n++)
    { long p = prime[pnn[n]], start = ch->w_low % p;
      walk_spec *w = &ch->walk[n];
      if(start < 0) { start += p; }
      w->first = (start == 0) ? p : start;
      w->head = p - w->first;
      if(w->head > range) { w->head = range; }
      w->nper = (range - w->head)/p;
      w->tail = range - w->head - w->nper*p;
    }
    for(n = 0; n < npr; n++)
    { long p = pr[n], start = ch->w_low % p;
      if(start < 0) { start += p; }
      ch->res[n] = start;
    }
  }
  for(n = 0; n < npr; n++)
  { long p = pr[n];
    for(k = 0; k < array_size; k++) { res0[k*npr + n] = k % p; }
  }
  return;
}

/* the data of the third stage for the row (a, b), see above; the part
   that depends on a alone (its inverse mod p, or that p | a) is kept
   from row to row while a stays */
static long inva3[NUM_PRIMES];  /* 1/a mod p, or 0 when p | a */
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
 * The tube                                                               *
 **************************************************************************/

/* A plain run bounds the fourth coordinate d as well, and d is a root of
   Q(d) = k2 d^2 + k1 d + k0, so a point of height at most h needs a root
   of Q in [-h, h]: either Q(h) and Q(-h) differ in sign (or one is 0),
   or both roots lie inside, which needs |k1| <= 2 h |k2| and Q(h), Q(-h)
   on the side of k2.  This holds only in a thin tube around the plane
   section d = 0 of the Kummer surface (or near the origin): at h = 2000
   on a few of the 63 words of a row, in many rows on none.  (The weaker
   |k0| <= h |k1| + h^2 |k2| is far from thin when h^2 k2 is the large
   term, as on curves of degree 5.)

   So the rows are not sieved by passes over a bit array.  For each a,
   the plane of (b, c) is cut into cells of rows times words, and a cell
   is dropped when the bounds show that no d in [-h, h] exists anywhere
   in it: Q(h), Q(-h), k1, k2 are polynomials in (b, c), and with the
   Taylor expansion sum p_jk s^j t^k at the centre (bm, cm), |s| <= rb,
   |t| <= rc, the value keeps the sign of p_00 throughout when |p_00| >
   sum' |p_jk| rb^j rc^k; the rounding errors of the coefficients, which
   may cancel, are covered by TUBE_EPS times the absolute polynomial at
   (|bm| + rb, |cm| + rc), so that no word that can hold a point is
   lost.  A cell that survives is halved, in rows or in words, until it
   is at most TUBE_ROWS rows by one word; that word is then collected
   for each of the rows.  The rows of a band are then sieved in order:
   each collected word tested against all the table primes, one by one,
   and its surviving bits against the third stage and the exact check. */

#define TUBE_EPS 1.5e-5   /* 15 coefficients, each with a margin of 1e-6 */
#define TUBE_ROWS 64      /* rows of a leaf cell (16 to 128 measured much the same) */
#define TUBE_BAND 4096    /* rows analysed and then sieved together */
#define TUBE_RANGES 8     /* word ranges kept per row; more are merged */

/* for the current a, as polynomials in (b, c) with P[j][k] the coefficient
   of b^j c^k: k0 (degree 4), k1 (3), k2 (2), and the absolute values of
   the coefficients of |k0| + h |k1| + h^2 |k2|, the scale of the errors */
static double K0[25], K1[25], K2[25], KA[25];
static double th, th2;              /* h and h^2 */
static long tube_a;                 /* the current a */

/* the words collected for the rows of the band: per row up to
   TUBE_RANGES ranges [lo, hi) */
static long band_b0;                              /* the first row */
static long band_nr[TUBE_BAND];                   /* ranges per row */
static long band_lo[TUBE_BAND][TUBE_RANGES], band_hi[TUBE_BAND][TUBE_RANGES];
static int tube_counting;            /* count the words only */
static long tube_n, tube_ncells;

static void tube_a_init(long a)
{
  double A = (double)a, A2 = A*A, A3 = A2*A, A4 = A3*A;
  long i;
  for(i = 0; i < 25; i++) { K0[i] = 0.0; K1[i] = 0.0; K2[i] = 0.0; }
#define P(j, k) [5*(j) + (k)]
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
#undef P
  for(i = 0; i < 25; i++)
  { KA[i] = fabs(K0[i]) + th*fabs(K1[i]) + th2*fabs(K2[i]); }
  tube_a = a;
  return;
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

/* the shift of P[j][k] (b^j c^k, j + k <= deg, row length 5) to (bm, cm) */
static inline void shift2(double *P, long deg, double bm, double cm)
{
  long j, k;
  for(k = 0; k <= deg; k++) { shift1(P + k, deg - k, 5, bm); }
  for(j = 0; j <= deg; j++) { shift1(P + 5*j, deg - j, 1, cm); }
  return;
}

/* the sum over (j, k) != (0, 0), j + k <= deg, of |P[j][k]| rp[5j + k],
   rp holding the products rb^j rc^k */
static inline double tube_sum(double *P, long deg, double *rp)
{
  double sum = 0.0;
  long j, k;
  for(j = 0; j <= deg; j++)
    for(k = (j == 0) ? 1 : 0; k <= deg - j; k++)
    { sum += fabs(P[5*j + k]) * rp[5*j + k]; }
  return(sum);
}

/* the value of the polynomial with the coefficients PA at (x, y) */
static inline double tube_eval(double *PA, long deg, double x, double y)
{
  double sum = 0.0, xj = 1.0;
  long j, k;
  for(j = 0; j <= deg; j++)
  { double yk = 1.0;
    for(k = 0; k <= deg - j; k++) { sum += PA[5*j + k] * xj * yk; yk *= y; }
    xj *= x;
  }
  return(sum);
}

/* the word i can hold a point of the row b: count it, or collect it into
   the row's ranges (merged with the last one when they touch; beyond
   TUBE_RANGES ranges the row keeps their hull) */
static void tube_emit(long b, long i)
{
  long row = b - band_b0, nr;
  if(tube_counting) { tube_n++; return; }
  nr = band_nr[row];
  if(nr > 0 && band_hi[row][nr-1] >= i)
  { if(band_hi[row][nr-1] < i + 1) { band_hi[row][nr-1] = i + 1; }
    return;
  }
  if(nr == TUBE_RANGES)
  { band_hi[row][0] = i + 1; band_nr[row] = 1; return; }
  band_lo[row][nr] = i; band_hi[row][nr] = i + 1;
  band_nr[row] = nr + 1;
  return;
}

/* the cell of the rows b_lo..b_hi and the words [i_lo, i_hi) */
static void tube_cell(long b_lo, long b_hi, long i_lo, long i_hi)
{
  static double rp[25], last_rb = -1.0, last_rc = -1.0;
  tube_ncells++;
  long rows = b_hi - b_lo + 1, words = i_hi - i_lo;
  double p0[25], p1[25], p2[25], p[25], m[25];
  double bm = 0.5*(double)(b_lo + b_hi), rb = 0.5*(double)(b_hi - b_lo);
  double cm = 32.0*(double)(i_lo + i_hi) - 0.5;
  double rc = 32.0*(double)(i_hi - i_lo) - 0.5;
  double ma, m1, m2, lowp, lowm, low2, up2;
  long j, k;
  /* the products rb^j rc^k, the same for cells of one size */
  if(rb != last_rb || rc != last_rc)
  { double rbj;
    for(j = 0, rbj = 1.0; j < 5; j++, rbj *= rb)
    { double rck = rbj;
      for(k = 0; k < 5; k++, rck *= rc) { rp[5*j + k] = rck; }
    }
    last_rb = rb; last_rc = rc;
  }
  for(j = 0; j < 25; j++) { p0[j] = K0[j]; p1[j] = K1[j]; p2[j] = K2[j]; }
  shift2(p0, 4, bm, cm); shift2(p1, 3, bm, cm); shift2(p2, 2, bm, cm);
  for(j = 0; j < 25; j++)
  { p[j] = p0[j] + th*p1[j] + th2*p2[j];    /* Q(h) */
    m[j] = p0[j] - th*p1[j] + th2*p2[j];    /* Q(-h) */
  }
  /* the scales of the rounding errors: |k0| + h |k1| + h^2 |k2| over
     the cell bounds each of them */
  ma = tube_eval(KA, 4, fabs(bm) + rb, fabs(cm) + rc);
  m1 = ma / th; m2 = ma / th2;
  /* a root of Q in [-h, h] by a sign change between -h and h: out when
     Q(h) and Q(-h) keep one and the same sign throughout the cell */
  lowp = fabs(p[0]) - tube_sum(p, 4, rp) - TUBE_EPS*ma;
  lowm = fabs(m[0]) - tube_sum(m, 4, rp) - TUBE_EPS*ma;
  if(lowp > 0.0 && lowm > 0.0 && (p[0] > 0.0) == (m[0] > 0.0))
  { /* both roots inside: k2 must have the sign of Q(h), Q(-h), and
       |k1| <= 2 h |k2|; out when either fails throughout */
    low2 = fabs(p2[0]) - tube_sum(p2, 2, rp) - TUBE_EPS*m2;
    if(low2 > 0.0 && (p2[0] > 0.0) != (p[0] > 0.0)) { return; }
    up2 = fabs(p2[0]) + tube_sum(p2, 2, rp) + TUBE_EPS*m2;
    if(fabs(p1[0]) - tube_sum(p1, 3, rp) - TUBE_EPS*m1 > 2.0*th*up2) { return; }
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
    for(b = b_lo; b <= b_hi; b++) { tube_emit(b, i_lo); }
  }
  return;
}

/* the words per row the tube leaves and the cells per row its analysis
   visits, for the cost model, from bands of 64 rows of a few a spread
   over the box */
void tube_sample(double *words_per_row, double *cells_per_row)
{
  long k, rows = 0, w_low = (-height)>>LONG_SHIFT, w_high = (height>>LONG_SHIFT) + 1;
  double words = 0.0, cells = 0.0;
  th = (double)height; th2 = th*th;
  tube_counting = 1;
  for(k = 0; k < 24; k++)
  { long a = (k * 7919L) % (height + 1), b_lo, b_hi;
    if(degree == 5 && mpz_cmp_si(&coeffs[5], 1) == 0)
    { a = a % (long)floor(sqrt((double)height) + 1); a = a*a; }
    b_lo = ((k * 104729L) % (2*height + 1)) - height;
    if(a == 0 && b_lo < 1) { b_lo = 1; }
    b_hi = b_lo + 63;
    if(b_hi > height) { b_hi = height; }
    tube_a_init(a);
    tube_n = 0; tube_ncells = 0;
    band_b0 = b_lo;
    tube_cell(b_lo, b_hi, w_low, w_high);
    words += (double)tube_n;
    cells += (double)tube_ncells;
    rows += b_hi - b_lo + 1;
  }
  tube_counting = 0;
  *words_per_row = (rows > 0) ? words / (double)rows : 1.0;
  *cells_per_row = (rows > 0) ? cells / (double)rows : 1.0;
  return;
}

/* the word i of the row (a, b) whose table rows are set: the mask, the
   ends of the row and the excluded c0, then the table primes one by
   one, then the bits; returns 1 when one point is enough and one was
   found */
static int tube_word(long a, long b, long i, bit_array mask, long c0)
{
  static long cur_chunk = 0;
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

/* the rows of a, in the tube: the multiples of m among the b (see
   row_step), from b_first on; band by band, the words that can hold a
   point collected, then the rows sieved in order */
int sift_tube(long a, long b_first, long m)
{
  long b_lo;
  th = (double)height; th2 = th*th;
  tube_a_init(a);
  for(b_lo = b_first; b_lo <= height; b_lo += TUBE_BAND)
  { long b_hi = b_lo + TUBE_BAND - 1, b, k;
    if(b_hi > height) { b_hi = height; }
    band_b0 = b_lo;
    for(k = 0; k <= b_hi - b_lo; k++) { band_nr[k] = 0; }
    tube_cell(b_lo, b_hi, w_low_all, w_high_all);
    for(b = b_lo; b <= b_hi; b += m)
    { long nr = band_nr[b - b_lo], c0, i;
      bit_array mask;
      if(nr == 0) { continue; }
      row_setup(a, b);
      c0 = height + LONG_LENGTH;
      if(a != 0 && (b*b)%(4*a) == 0) { c0 = (b*b)/(4*a); }
      mask = ~zero;
      if(!((a|b)&1)) { mask = HALF_MASK; }
      for(k = 0; k < nr; k++)
      { num_surv1 += band_hi[b - b_lo][k] - band_lo[b - b_lo][k];
        for(i = band_lo[b - b_lo][k]; i < band_hi[b - b_lo][k]; i++)
        { if(tube_word(a, b, i, mask, c0) && one_point) { return(1); } }
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
 | of the code the compiler produces for the innermost loops in sift0.    |
 +------------------------------------------------------------------------*/

/* The walk of the n-th prime over the chunk: OP(word of the bit array,
   table word) for every word, the table row read from its word w->first
   up to its end, then in full periods of p words, then the tail. */
#define WALK(OP) \
  { bit_array *surv = survivors; \
    bit_array *siv1 = &rowptr[n][w->first]; \
    long j; \
    for(j = w->head; j; j--) { OP(*surv, *siv1); surv++; siv1++; } \
    /* now siv1 points at the end of the table row, whose word p repeats \
       its word 0, if there is anything left to do */ \
    for(j = w->nper; j; j--) \
    { bit_array *siv0 = siv1; \
      siv1 -= p; \
      do { OP(*surv, *siv1); surv++; siv1++; } while(siv1 != siv0); \
    } \
    if((j = w->tail)) \
    { siv1 -= p; \
      for( ; j; j--) { OP(*surv, *siv1); surv++; siv1++; } \
    } \
  }
#define AND_OP(d, s) ((d) &= (s))
#define FILL_OP(d, s) ((d) = mask & (s))

/* Sieve the k-th chunk of the row (a, b): the bit array set from the mask
   (all bits, or the odd c when a and b are even) and the table of the
   first prime in one pass, the ends of the row and the excluded c0 taken
   out, the other first-stage primes ANDed in, the surviving words passed
   to the second stage. */
static int sift0(long a, long b, long k, bit_array mask, long c0)
{
  chunk_spec *ch = &chunks[k];
  long n, range = ch->w_high - ch->w_low;
  /* the fill, merged with the first prime's pass if there is one */
  if(sieve_primes1 == 0)
  { bit_array *surv = survivors;
    long i;
    for(i = range; i; i--) { *surv++ = mask; }
    n = 0;
  }
  else
  { long p = pr[0];
    walk_spec *w = &ch->walk[0];
    n = 0;
    WALK(FILL_OP);
    n = 1;
  }
  if(k == 0) { survivors[0] &= begmask; }
  if(k == num_chunks - 1) { survivors[range-1] &= endmask; }
  if(ch->w_low <= (c0>>LONG_SHIFT) && (c0>>LONG_SHIFT) < ch->w_high)
  { survivors[(c0>>LONG_SHIFT) - ch->w_low] &= ~(1UL<<(c0 & LONG_MASK)); }
#if (DEBUG >= 3)
  { long i;
    for(i = 0; i < range; i++) printf(" %8.8lx",survivors[i]);
  }
#endif
#if (DEBUG >= 1)
  printf("\n sift0(%ld, %ld)\n", ch->w_low, ch->w_high);
#endif
  /* now do the sieving (fast!) */
  for( ; n < sieve_primes1; n++)
  { long p = pr[n];
    walk_spec *w = &ch->walk[n];
    WALK(AND_OP);
  }
#if (DEBUG >= 3)
  { long i;
    for(i = 0; i < range; i++) { printf(" %8.8lx",survivors[i]); }
    printf("\n");
  }
#endif
  /* Check the points that have survived the sieve if they really are points */
  { bit_array *surv0 = &survivors[0];
    bit_array nums;
    long i;
    for(i = 0; i < range; i++)
    { if((nums = *surv0++))
      { if(check_point(nums, a, b, ch->w_low + i, i, ch->res) && one_point)
        { return(1); }
  } } }
  return(0);
}

/* the residues of b modulo the sieving primes and the table rows of
   (a, b): b steps by a constant from one row to the next, so they are
   kept and advanced, and computed by division only when a changes or
   the caller jumps back; a step of any size is reduced mod p by the
   reciprocal in double precision (exact: the error of b/p is far below
   1/p) */
static void row_setup(long a, long b)
{
  static long last_a = -1, last_b = 0, last_d = 0;
  static long bres[NUM_PRIMES], dres[NUM_PRIMES], drow[NUM_PRIMES];
  static double invp[NUM_PRIMES];
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
      invp[n] = 1.0 / (double)p;
    }
    last_d = 0;
  }
  last_a = a; last_b = b;
  return;
}

int sift(long a, long b)
/* print points surviving sieve */
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
  mask = ~zero;
  if(!((a|b)&1)) { mask = HALF_MASK; }
  /* Now the chunks of longwords (= bit_arrays) */
  for(k = 0; k < num_chunks; k++)
  { if(sift0(a, b, k, mask, c0) && one_point) { return(1); } }
  return(0);
}
