/******* j-sift0-1.0.c ***************************************************/

#include <stdlib.h>
#include <stdio.h>

#include "j-points.h"


sieve_spec sieves1[NUM_PRIMES];
sieve_spec sieves2p[NUM_PRIMES];
sieve_spec sieves2n[NUM_PRIMES];


/**************************************************************************
 * compute the gcd of two integers                                        *
 **************************************************************************/

static inline long gcd(long m, long n)
{
  if(n < 0) n = -n;
  /* m is always nonnegative here */
  while(1)
  { if(n == 0) return(m);
    m %= n;
    if(m == 0) return(n);
    n %= m;
  }
}

static inline int check_point(bit_array nums, long a, long b, long i) 
{ /* to be more precise, this checks a complete bit_array of possible
     survivors */   
  sieve_spec *ssp = &sieves2p[0];
  if(i < 0) ssp = &sieves2n[0];
  num_surv1++;
  { long n;
    register long ii = i;
    for(n = sieve_primes2 - sieve_primes1; n && nums; n--)
    { nums &= ssp->ptr[ii%(ssp->p)];
      ssp++;
    }
  }
  if(nums)
  { long c;
    bit_array test; /* test has one bit set and is used to isolate
                       the corresponding bit in nums */
    /* c will be the coordinate corresponding to the selected bit */
    for(c = i<<LONG_SHIFT, test = 1 ; test; c++, test <<= 1)
      /* test one bit */
      if(nums & test)
      { num_surv2++;
        if(gcd(gcd(a, b), c) == 1 && check_one_point(a, b, c) && one_point)
          return(1);
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
 | It is advisable to force the compiler to put the variables             |
 | surv, siv0, siv1 into registers.                                       |
 | With gcc-2.7.2.1, I got a speedup of more than 30% ! -- MS             |
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
    register bit_array *surv;
    surv = survivors;
    
    if(p_high < p_low)
    { 
      register bit_array *siv1;
      register long i;
      siv1 = &sieve_n[w_low - p * p_high];
      for(i = range; i ; i--) *surv++ &= *siv1++;
    }
    else
    {
      register bit_array *siv1;
      { register long j;
        j = p * p_low - w_low;
        siv1 = &sieve_n[p-j];
        for( ; j; j--) *surv++ &= *siv1++;
      }
      {
        register long i = p_high - p_low;
        p_high *= p;
        for( ; i; i--)
        {
          register bit_array *siv0;
          siv0 = siv1;
          siv1 -= p;
          do {*surv++ &= *siv1++;} while(siv1 != siv0);
        }
      }
      siv1 -= p;
      { register long j;
        for(j = w_high - p_high; j; j--) *surv++ &= *siv1++;
  } } }
#if (DEBUG >= 3)
  for(i = 0; i < range; i++) printf(" %8.8lx",survivors[i]);
  printf("\n");
#endif
  /* Check the points that have survived the sieve if they really are points */
  { bit_array *surv0 = &survivors[0];
    bit_array nums;
    for(i = w_low; i < w_high; i++)
      if((nums = *surv0++))
        if(check_point(nums, a, b, i) && one_point)
          return(1);
  }
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
  if(a != 0 && (b*b)%(4*a) == 0) c0 = (b*b)/(4*a);
#if (DEBUG >= 1)
  printf("\n sift(a = %ld, b = %ld)\n", a, b);
#endif
  
  { long n, m;
    for(n = 0; n < sieve_primes1; n++)
    { long pn = pnn[n], p = prime[pn], bp = b%p;
      if(bp < 0) bp += p;
      sieves1[n].ptr = &sieve[n][a%p][bp][0];
    }
    for(m = 0 ; n < sieve_primes2; n++, m++)
    { long pn = pnn[n], p = prime[pn], bp = b%p;
      if(bp < 0) bp += p;
      sieves2n[m].ptr = (sieves2p[m].ptr = &sieve[n][a%p][bp][0]) + p;
    }
  }
  /* Now the range of longwords (= bit_arrays) */
  mask = ~zero;
  if(!((a|b)&1)) mask = HALF_MASK; 
  w_low = (-height)>>LONG_SHIFT;    /* FLOOR(-height, LONG_LENGTH); */
  w_high = (height>>LONG_SHIFT)+1;  /* CEIL(height+1, LONG_LENGTH); */
  for(w_low0 = w_low; w_low0 < w_high; w_low0 += array_size)
  { w_high0 = w_low0 + array_size;
    if(w_high0 > w_high) w_high0 = w_high;
    range = w_high0 - w_low0;
    /* initialise the bits */
    {
      register bit_array *surv;
      register long i;
      register long mask1;
      surv = &survivors[0];
      mask1 = mask;
      for(i = range; i; i--) *surv++ = mask1;
    }
    if(w_low0 == w_low) survivors[0] &= begmask;
    if(w_high0 == w_high) survivors[range-1] &= endmask;
    if(w_low0 <= (c0>>LONG_SHIFT) && (c0>>LONG_SHIFT) < w_high0)
      survivors[(c0>>LONG_SHIFT) - w_low0] &= ~(1UL<<(c0 & LONG_MASK));
#if (DEBUG >= 3)
    { long i;
      for(i = 0; i < range; i++) printf(" %8.8lx",survivors[i]);
    }
#endif
#if (DEBUG >= 1)
    printf("\n sift0(%ld, %ld)\n", w_low0, w_high0);
#endif
    if(sift0(a, b, w_low0, w_high0) && one_point) return(1);
  }
  return(0);
}

