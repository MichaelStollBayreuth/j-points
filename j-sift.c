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
 *    the surviving bits handed to the exact check; the bookkeeping   *
 *    of a row done without integer divisions                          *
 ***********************************************************************/

#include <stdlib.h>
#include <stdio.h>

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

static long np2;                  /* the number of second-stage primes */
static long p2[NUM_PRIMES];       /* the second-stage primes */
static unsigned short *res0;
/* res0[k*np2 + m] = k mod p2[m] for 0 <= k < array_size: with the residue
   of the chunk's first word, the residue of any word of the chunk costs an
   addition and a comparison instead of a division */

/* the sieving primes, the lengths of their table rows and of the block
   of rows of one residue class of a, and the table rows of the current
   (a, b): rowptr[n] for the n-th prime, row2 = &rowptr[sieve_primes1]
   for the second-stage primes */
static long pr[NUM_PRIMES], rowlen[NUM_PRIMES], blocklen[NUM_PRIMES];
static bit_array *rowptr[NUM_PRIMES];
static bit_array **row2;


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

  np2 = sieve_primes2 - sieve_primes1;
  for(n = 0; n < sieve_primes2; n++)
  { pr[n] = prime[pnn[n]];
    rowlen[n] = pr[n] + 1;
    blocklen[n] = pr[n]*rowlen[n];
  }
  row2 = &rowptr[np1];
  num_chunks = CEIL(w_high - w_low, array_size);
  survivors = (bit_array *)malloc(array_size*sizeof(bit_array));
  chunks = (chunk_spec *)malloc(num_chunks*sizeof(chunk_spec));
  walk = (walk_spec *)malloc((num_chunks*np1 + 1)*sizeof(walk_spec));
  res = (long *)malloc((num_chunks*np2 + 1)*sizeof(long));
  res0 = (unsigned short *)malloc((array_size*np2 + 1)*sizeof(unsigned short));
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
    ch->res = &res[k*np2];
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
    for(n = 0; n < np2; n++)
    { long p = prime[pnn[np1 + n]], start = ch->w_low % p;
      if(start < 0) { start += p; }
      ch->res[n] = start;
    }
  }
  for(n = 0; n < np2; n++)
  { long p = prime[pnn[np1 + n]];
    p2[n] = p;
    for(k = 0; k < array_size; k++) { res0[k*np2 + n] = k % p; }
  }
  return;
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
  num_surv1++;
  { unsigned short *r0 = &res0[k*np2];
    long m;
    for(m = 0; m < np2 && nums; m++)
    { long r = r0[m] + res[m];
      if(r >= p2[m]) { r -= p2[m]; }
      nums &= row2[m][r];
    }
  }
  if(nums)
  { long c;

    /* c will be the coordinate corresponding to the selected bit */
    for(c = i<<LONG_SHIFT; nums; nums >>= 1, c++)
    {/* test one bit */
      if(nums & 1)
      { num_surv2++;
        if(relprime3(a, b, c) && check_one_point(a, b, c) && one_point)
        { return(1); }
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

int sift(long a, long b)
/* print points surviving sieve */
{
  /* the residues of b modulo the sieving primes and the table rows of
     (a, b): b steps by a constant from one row to the next, so they are
     kept and advanced, and computed by division only when a or the step
     changes or the caller jumps back */
  static long last_a = -1, last_b = 0, last_d = 0;
  static long bres[NUM_PRIMES], dres[NUM_PRIMES], drow[NUM_PRIMES];
  long n, k, d = b - last_b;
  /* c0 is value of c that has to be excluded */
  long c0 = height + LONG_LENGTH;
  bit_array mask;
  if(a != 0 && (b*b)%(4*a) == 0) { c0 = (b*b)/(4*a); }
#if (DEBUG >= 1)
  printf("\n sift(a = %ld, b = %ld)\n", a, b);
#endif

  if(a == last_a && d > 0)
  { if(d != last_d)
    { /* the step d and its residues, at most once per a */
      for(n = 0; n < sieve_primes2; n++)
      { dres[n] = d % pr[n];
        drow[n] = dres[n]*rowlen[n];
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

  /* Now the chunks of longwords (= bit_arrays) */
  mask = ~zero;
  if(!((a|b)&1)) { mask = HALF_MASK; }
  for(k = 0; k < num_chunks; k++)
  { if(sift0(a, b, k, mask, c0) && one_point) { return(1); } }
  return(0);
}
