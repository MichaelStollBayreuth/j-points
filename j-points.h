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
 * j-points.h                                                          *
 *  - the constants, types and declarations shared by j-points.c and   *
 *    j-sift.c                                                         *
 ***********************************************************************/

#ifndef DEBUG
#define DEBUG 0
#endif

/* MAX_PRIME_EVEN is the least power of 2 that is larger than MAX_PRIME */
/* NUM_PRIMES odd primes up to MAX_PRIME are in the table; the primes up
   to MAX_TABLE_PRIME can get a sieve table (p*p*(p+1) words) for the
   first two stages, the larger ones serve the third stage, per bit */
#if (DEBUG > 0)
# define NUM_PRIMES 5
# define MAX_PRIME 13
# define MAX_PRIME_EVEN 16
# define MAX_TABLE_PRIME 13
#else
# define NUM_PRIMES 53
# define MAX_PRIME 251
# define MAX_PRIME_EVEN 256
# define MAX_TABLE_PRIME 127
#endif /* DEBUG > 0 */

#define FLOOR(a,b) (((a) < 0) ? -(1 + (-(a)-1) / (b)) : (a) / (b))
#define CEIL(a,b) (((a) <= 0) ? -(-(a) / (b)) : 1 + ((a)-1) / (b))

#define LONG_LENGTH (8 * sizeof(unsigned long))
   /* number of bits in an unsigned long */
#define LONG_SHIFT ((LONG_LENGTH == 16) ? 4 : \
                    (LONG_LENGTH == 32) ? 5 : \
		    (LONG_LENGTH == 64) ? 6 : 0)
#define LONG_MASK (~(-1L<<LONG_SHIFT))
#define MAX_HEIGHT ((LONG_LENGTH == 16) ? 181UL : \
                    (LONG_LENGTH == 32) ? 46340UL : \
		    (LONG_LENGTH == 64) ? 3037000499UL : 0UL)

#define HALF_MASK ((bit_array)(~((unsigned long)(-1L) / 3)))

typedef unsigned long bit_array;
#define zero ((bit_array)0);

extern long prime[];
extern long pnn[];
extern int one_point;        /* A flag saying if one point is enough */
extern long sieve_primes1;
   /* The number of primes used for the first sieving stage */
extern long sieve_primes2;
   /* The number of primes used for the first two sieving stages */
extern long sieve_primes3;
   /* The number of primes used for all three sieving stages */
extern bit_array *sieve_tab[NUM_PRIMES];
   /* The sieve table of each sieving prime, see j-points.c */
extern long inverses[NUM_PRIMES][MAX_PRIME_EVEN];
extern unsigned char has_infinity[NUM_PRIMES];
extern unsigned char is_f_square[NUM_PRIMES][MAX_PRIME_EVEN];
extern unsigned char is_point_on_j[NUM_PRIMES][MAX_PRIME_EVEN][MAX_PRIME_EVEN];
   /* The tables of the third stage, see j-points.c */

extern double kd400, kd310, kd301, kd220, kd211, kd202, kd130, kd121, kd112,
              kd103, kd040, kd031, kd022, kd013, kd004;
extern double fd[7];
extern int tube_mode;    /* the tube of plain runs, see j-sift.c */
extern void tube_sample(double *, double *);
extern long degree;
extern MP_INT coeffs[7];
extern int sift_tube(long, long, long);

extern long num_surv1;   /* Used to count the survivors of the first stage */
extern long num_surv2;   /* Used to count the survivors of the second stage */
extern long num_surv3;   /* Used to count the survivors of the third stage */
extern bit_array *survivors; /* In this array the sieving takes place */
extern long height;          /* The height bound */
extern long array_size;      /* The size of the survivors array (in longs) */
extern bit_array begmask, endmask;
   /* Bit masks for the beginning and end of the sieving array */

extern int check_one_point(long, long, long);
extern void error(long);

extern void init_sift(void);
extern int sift(long, long);
