/***********************************************************************
 * j-points-2.1                                                        *
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
 *    the surviving bits handed to the exact check                     *
 ***********************************************************************/

#include <stdlib.h>
#include <stdio.h>

#include "j-points.h"


sieve_spec sieves1[NUM_PRIMES];
sieve_spec sieves2p[NUM_PRIMES];
sieve_spec sieves2n[NUM_PRIMES];


/**************************************************************************
 * Some helper functions                                                  *
 **************************************************************************/

// static inline long long m, long n)
// {
//   if(n < 0) n = -n;
//   /* m is always nonnegative here */
//   while(1)
//   { if(n == 0) return(m);
//     m %= n;
//     if(m == 0) return(n);
//     n %= m;
//   }
// }

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

static inline int check_point(bit_array nums, long a, long b, long i)
{ /* to be more precise, this checks a complete bit_array of possible
     survivors */
  sieve_spec *ssp = (i < 0) ? &sieves2n[0] : &sieves2p[0];

  num_surv1++;
  { long n;
    long ii = i;
    for(n = sieve_primes2 - sieve_primes1; n && nums; n--)
    { nums &= ssp->ptr[ii%(ssp->p)];
      ssp++;
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

int sift0(long a, long b, long w_low, long w_high)
{
   /* now do the sieving (fast!) */
  long n, i, range = w_high - w_low;
  for(n = 0; n < sieve_primes1; n++)
  {
    long p = prime[pnn[n]];
    bit_array *sieve_n = sieves1[n].ptr;
    long p_low = CEIL(w_low, p), p_high = FLOOR(w_high, p);
    bit_array *surv = survivors;

    if(p_high < p_low)
    {
      bit_array *siv1;
      long i;
      siv1 = &sieve_n[w_low - p * p_high];
      for(i = range; i ; i--) { *surv++ &= *siv1++; }
    }
    else
    {
      bit_array *siv1;
      long j = p * p_low - w_low;
      siv1 = &sieve_n[p-j];
      for( ; j; j--) { *surv++ &= *siv1++; }
      j = p_high - p_low;
      p_high *= p;
      for( ; j; j--)
      {
        bit_array *siv0;
        siv0 = siv1;
        siv1 -= p;
        do {*surv++ &= *siv1++;} while(siv1 != siv0);
      }
      siv1 -= p;
      for(j = w_high - p_high; j; j--) { *surv++ &= *siv1++; }
  } }
#if (DEBUG >= 3)
  for(i = 0; i < range; i++) { printf(" %8.8lx",survivors[i]); }
  printf("\n");
#endif
  /* Check the points that have survived the sieve if they really are points */
  { bit_array *surv0 = &survivors[0];
    bit_array nums;
    for(i = w_low; i < w_high; i++)
    { if((nums = *surv0++))
      { if(check_point(nums, a, b, i) && one_point) { return(1); }
  } } }
  return(0);
}

int sift(long a, long b)
/* print points surviving sieve */
{
  long range;
  long w_low, w_high;
  long w_low0, w_high0;
  /* c0 is value of c that has to be excluded */
  long c0 = height + LONG_LENGTH;
  bit_array mask;
  if(a != 0 && (b*b)%(4*a) == 0) { c0 = (b*b)/(4*a); }
#if (DEBUG >= 1)
  printf("\n sift(a = %ld, b = %ld)\n", a, b);
#endif

  { long n, m;
    for(n = 0; n < sieve_primes1; n++)
    { long pn = pnn[n], p = prime[pn], bp = b%p;
      if(bp < 0) { bp += p; }
      sieves1[n].ptr = &sieve[n][a%p][bp][0];
    }
    for(m = 0 ; n < sieve_primes2; n++, m++)
    { long pn = pnn[n], p = prime[pn], bp = b%p;
      if(bp < 0) { bp += p; }
      sieves2n[m].ptr = (sieves2p[m].ptr = &sieve[n][a%p][bp][0]) + p;
    }
  }
  /* Now the range of longwords (= bit_arrays) */
  mask = ~zero;
  if(!((a|b)&1)) { mask = HALF_MASK; }
  w_low = (-height)>>LONG_SHIFT;    /* FLOOR(-height, LONG_LENGTH); */
  w_high = (height>>LONG_SHIFT)+1;  /* CEIL(height+1, LONG_LENGTH); */
  for(w_low0 = w_low; w_low0 < w_high; w_low0 += array_size)
  { w_high0 = w_low0 + array_size;
    if(w_high0 > w_high) { w_high0 = w_high; }
    range = w_high0 - w_low0;
    /* initialise the bits */
    {
      bit_array *surv;
      long i;
      surv = &survivors[0];
      for(i = range; i; i--) { *surv++ = mask; }
    }
    if(w_low0 == w_low) { survivors[0] &= begmask; }
    if(w_high0 == w_high) { survivors[range-1] &= endmask; }
    if(w_low0 <= (c0>>LONG_SHIFT) && (c0>>LONG_SHIFT) < w_high0)
    { survivors[(c0>>LONG_SHIFT) - w_low0] &= ~(1UL<<(c0 & LONG_MASK)); }
#if (DEBUG >= 3)
    { long i;
      for(i = 0; i < range; i++) printf(" %8.8lx",survivors[i]);
    }
#endif
#if (DEBUG >= 1)
    printf("\n sift0(%ld, %ld)\n", w_low0, w_high0);
#endif
    if(sift0(a, b, w_low0, w_high0) && one_point) { return(1); }
  }
  return(0);
}

