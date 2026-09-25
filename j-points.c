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
 * j-points.c                                                          *
 *  - the driver: reads the curve and the options, sets up the tables  *
 *    and the Kummer equation, finds the double points and the points  *
 *    with a = 0, runs the sieve over the other triples (a, b) and     *
 *    checks its survivors exactly                                     *
 ***********************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <gmp.h>

#include "j-points.h"

/**************************************************************************
 * define                                                                 *
 **************************************************************************/

#define DEFAULT_SIZE 10     /* Default value for the -s option */

#define J_POINTS_VERSION \
  "This is j-points-3.0 by Michael Stoll (2026-09-25).\n\n" \
  "Please acknowledge use of the program in published work.\n"

/**************************************************************************
 * global variables                                                       *
 **************************************************************************/

#if (DEBUG > 0)
long prime[NUM_PRIMES+1] =
{3,5,7,11,13};
#else
long prime[NUM_PRIMES+1] =
{3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101,
 103,107,109,113,127,131,137,139,149,151,157,163,167,173,179,181,191,193,197,
 199,211,223,227,229,233,239,241,251};
#endif
long pnn[NUM_PRIMES+1];  /* This array holds the numbers (= index in prime[])
                            of the sieving primes, stage by stage */

long num_primes = NUM_PRIMES; /* the primes considered, from the beginning
                                 of the table; -p sets it */

/* The cost model that chooses the sieving primes and their stages (see
   choose_primes): the costs in cycles of one core, measured on the curve
   of test1 at 2000 with the sieve of 3.0.  Only their ratios matter for
   the choice, and those are much the same on any current machine. */
#define COST_AND      0.9  /* one word of one first-stage pass ... */
#define COST_SIZE   270.0  /* ... times 1 + p/COST_SIZE: the tables of the
                              larger primes spill out of the cache (12%
                              more per word for a set of primes with mean
                              84 than for one with mean 51) */
#define COST_WORD    56.0  /* a word surviving the first stage: found and
                              handed to the second stage */
#define COST_TEST2    8.0  /* one second-stage test of a word */
#define COST_BIT     30.0  /* a bit surviving the second stage: found,
                              its coprimality tested */
#define COST_ROW3    40.0  /* per third-stage prime, the data of a row
                              that has such a bit (see row3_setup) */
#define COST_TEST3   15.0  /* one third-stage test of a bit */
#define COST_EXACT 1300.0  /* the exact check of a triple */
#define COST_ROW    300.0  /* a row of the box sieve besides its passes:
                              the fill, the scan, the set-up */
#define COST_CELL   320.0  /* a cell of the tube's analysis */
#define COST_TWORD   12.0  /* a word the tube leaves, besides its tests */
#define COST_TTEST    4.0  /* one test of such a word (no early exit) */
#define COST_TROW    45.0  /* a row in the tube: the set-up */
#define COST_TABLE   13.0  /* one word of a sieve table at set-up */
#define COST_INIT   110.0  /* one class (1, b, c) of the exact table of a
                              prime whose rate was sampled, at set-up */
#define SAMPLE_PAIRS 4096  /* the classes (1, b, c) sampled for the rate
                              of a prime, see init_fmodpsquare */

bit_array bits[LONG_LENGTH]; /* An array of bit masks */

typedef struct {long p; double r; long n; long np0;} entry;

long inverses[NUM_PRIMES][MAX_PRIME_EVEN];
/* inverses[pn][a] = b such that a*b = 1 mod p = prime[pn]; 0 < a, b < p. */
int squares[NUM_PRIMES][MAX_PRIME];
/* squares[pn][x] = 1 if x is a square mod prime[pn], 0 if not */
unsigned char has_infinity[NUM_PRIMES];
/* has_infinity[pn] = 1 if the curve has rational points at infinity,
   0 if not */
unsigned char is_f_square[NUM_PRIMES][MAX_PRIME_EVEN];
/* is_f_square[pn][x] = 1 if f(x) is a square mod prime[pn], 0 if not */
unsigned char is_point_on_j[NUM_PRIMES][MAX_PRIME_EVEN][MAX_PRIME_EVEN];
/* is_point_on_j[pn][x][y] = 1 if there are points on J mod prime[pn]
   with first three coordinates on the Kummer surface (1,x,y), 0 if not */

bit_array *sieve_tab[NUM_PRIMES];
/* The sieve table of the pn-th sieving prime p = prime[pnn[pn]], allocated
   by init_sieve: p*p rows of p+1 words, the row of (a1, b1) starting at
   index (a1*p + b1)*(p+1).  Bit k of word c1 of that row is 0 iff (a,b,c)
   is excluded mod p when a = a1, b = b1 and c = c1*LONG_LENGTH + k mod p.
   The row is periodic in c1 with period p; its word p repeats its word 0. */

MP_INT coeffs[7];  /* The coefficients of f */
MP_INT bc[7]; /* A helper array */
MP_INT fff, tmp, tmp2, tmp3, ddd;   /* Some multi-precision integer variables */
MP_INT k400, k310, k301, k220, k211, k202, k130, k121, k112, k103,
       k040, k031, k022, k013, k004;
        /* coefficients for Kummer equation */
MP_INT x12, x22, x32;
MP_INT kummer[3];

/* The tube (see j-sift.c): the coefficients of the Kummer equation and
   of f as doubles, for the bounds per row, and whether the tube is used
   (plain runs: the bound on the fourth coordinate cuts the box to a
   tube around the plane section d = 0 of the Kummer surface; -a has no
   such bound, and the unsieved run must not use it) */
double kd400, kd310, kd301, kd220, kd211, kd202, kd130, kd121, kd112,
       kd103, kd040, kd031, kd022, kd013, kd004;
double fd[7];
int tube_mode = 0;
double tube_words = 0.0;   /* words per row the tube leaves, sampled */
double tube_cells = 0.0;   /* cells per row its analysis visits, sampled */
MP_INT cpf1, cpf2, cpf3, cpfa, cpfb, cpfc;

long degree;
long coeffs_mod_p[NUM_PRIMES][8];
                         /* The coefficients of f reduced modulo the various
                            primes */
entry prec[NUM_PRIMES];  /* This array is used for sorting in order to
                            determine the `best' sieving primes. */

long height;          /* The height bound */
long sieve_primes1;   /* The number of primes used for the first sieving stage */
long sieve_primes2;   /* The number of primes used for the first two stages */
long sieve_primes3;   /* The number of primes used for all three stages */
int quiet;            /* A flag saying whether to suppress messages */
int one_point;        /* A flag saying if one point is enough */
int all_points;       /* Indicates that the `-a' option was given */
char *print_format;   /* The printf format for printing points */
long array_size;      /* The size of the survivors array (in longs) */

bit_array *survivors; /* In this array the sieving takes place */

long num_surv1 = 0;   /* Used to count the survivors of the first stage */
long num_surv2 = 0;   /* Used to count the survivors of the second stage */
long num_surv3 = 0;   /* Used to count the survivors of the third stage */

long total = 0;       /* Counts the points found */

long bound;

bit_array begmask, endmask;  /* Bit masks for the beginning and end of
                                the sieving array */

/**************************************************************************
 * prototypes                                                             *
 **************************************************************************/

void init_main(void);
void find_points(void);
void read_input(long, char *argv[]);
char *scan_mpz(char*, MP_INT*);
void init_inverses(void);
void init_squares(void);
void init_fmodpsquare(void);
void choose_primes(void);
void init_sieve(void);
void kummer_init(void);
void tube_init(void);
static inline int relprime(long, long);
int check_one_point_final(long, long, long, MP_INT *, MP_INT *);
int check_one_point(long, long, long);
void print_poly(MP_INT*, long);
void message(long, long);
void error(long);

/**************************************************************************
 * main                                                                   *
 **************************************************************************/

int main(int argc, char *argv[])
{
  long s;

  if(LONG_SHIFT == 0) error(1);
  init_main();
  /* read input */
  if(argc < 3) error(2);
  read_input(argc-1, &argv[0]);
  if(!quiet)
  { message(0, 0);
    message(5, degree);
    message(6, height);
    message(3, 0);
  }
  begmask = (~0UL)<<((-height) & LONG_MASK);
  endmask = (~0UL)>>((~height) & LONG_MASK);
  bound = (all_points) ? MAX_HEIGHT : height;
  s = 2*CEIL(height+1, LONG_LENGTH);
  array_size <<= 13 - LONG_SHIFT; /* from kbytes to longs */
  if(s < array_size) array_size = s;
  /* initialise data for equations */
  kummer_init();
  tube_init();
  /* find and count points */
  find_points();
  if(!quiet) { message(12, 0); message(2, total); }
  return(0);
}

/**************************************************************************
 * procedures                                                             *
 **************************************************************************/

/**************************************************************************
 * get at the input                                                       *
 **************************************************************************/

void read_input(long argc, char *argv[])
{
  { char *s = argv[1];
    degree = 0;
    while((degree <= 6) && (s = scan_mpz(s, &coeffs[degree]))) degree++;
    degree--;
    if(scan_mpz(s, &fff)) error(3);
    /* zero leading coefficients do not count */
    while(degree > 0 && mpz_sgn(&coeffs[degree]) == 0) degree--;
  }
  if(degree < 5) error(5);
  if(degree == 5) mpz_set_si(&coeffs[6], 0);
  if(sscanf(argv[2], " %ld", &height) != 1 || height < 1 || height > MAX_HEIGHT)
  { error(4); }
  /* Set global variables to their default values */
  sieve_primes1 = -1;  /* automatic determination of this number */
  sieve_primes2 = -1;  /* automatic determination of this number */
  sieve_primes3 = -1;  /* automatic determination of this number */
  quiet = 0;           /* don't be quiet */
  one_point = 0;       /* look for all points */
  all_points = 0;      /* take naive Kummer height */
  print_format = "(%ld, %ld, %ld, %ld)\n";
  array_size = DEFAULT_SIZE;
  /* recognise optional args */
  { long i = 3;
    while(i <= argc)
    {
      if(*(argv[i]) != '-') error(6);
      switch(argv[i][1])
      { case 'n': /* number of primes used for first stage of sieving */
          if(argc == i) error(6);
          i++;
          if(sscanf(argv[i], " %ld", &sieve_primes1) != 1) error(6);
          if(sieve_primes1 < 0) sieve_primes1 = 0;
          i++;
          break;
        case 'M': /* number of primes used for the first two stages */
          if(argc == i) error(6);
          i++;
          if(sscanf(argv[i], " %ld", &sieve_primes2) != 1) error(6);
          if(sieve_primes2 < 0) sieve_primes2 = 0;
          i++;
          break;
        case 'N': /* number of primes used for sieving altogether */
          if(argc == i) error(6);
          i++;
          if(sscanf(argv[i], " %ld", &sieve_primes3) != 1) error(6);
          if(sieve_primes3 < 0) sieve_primes3 = 0;
          i++;
          break;
        case 's': /* size of survivors array in kbytes */
          if(argc == i) error(6);
          i++;
          if(sscanf(argv[i], " %ld", &array_size) != 1) error(6);
          if(array_size <= 0) error(6);
          i++;
          break;
        case 'f': /* printing format */
          if(argc == i) error(6);
          i++;
          { long l = strlen(argv[i]);
            print_format = malloc((l+1)*sizeof(char));
            strcpy(print_format, argv[i]);
          }
          i++;
          break;
        case 'q': /* quiet */
          quiet = 1;
          i++;
          break;
        case 'a': /* all points */
          all_points = 1;
          i++;
          break;
        case '1': /* only one point */
          one_point = 1;
          i++;
          break;
        case 'p': /* set number of primes used */
          if(argc == i) error(6);
          i++;
          if(sscanf(argv[i], " %ld", &num_primes) != 1) error(6);
          if(num_primes <= 0) error(6);
          if(num_primes > NUM_PRIMES) { num_primes = NUM_PRIMES; }
          i++;
          break;
        default: error(6);
  } } }
  /* the numbers of sieving primes cannot exceed the number of primes
     considered, and a stage cannot have more primes than the stages
     after it, whatever the order of -n, -M, -N and -p */
  if(sieve_primes3 > num_primes) sieve_primes3 = num_primes;
  if(sieve_primes2 > num_primes) sieve_primes2 = num_primes;
  if(sieve_primes3 >= 0 && sieve_primes2 > sieve_primes3)
    sieve_primes2 = sieve_primes3;
  if(sieve_primes1 > num_primes) sieve_primes1 = num_primes;
  if(sieve_primes3 >= 0 && sieve_primes1 > sieve_primes3)
    sieve_primes1 = sieve_primes3;
  if(sieve_primes2 >= 0 && sieve_primes1 > sieve_primes2)
    sieve_primes1 = sieve_primes2;
}

/* Read in a long long long integer. Should really be in the library. */
char *scan_mpz(char *s, MP_INT *x)
{
  long neg = 0;
  if(s == NULL || *s == 0) return NULL;
  while(*s == ' ') s++;
  if(*s == 0) return NULL;
  if(*s == '-') {neg = 1; s++;}
  else if(*s == '+') s++;
  mpz_set_si(&tmp2, 0);
  while('0' <= *s && *s <= '9')
  { mpz_mul_ui(&tmp2, &tmp2, 10);
    mpz_add_ui(&tmp2, &tmp2, (long)(*s - '0'));
    s++; }
  if(neg) mpz_neg(&tmp2, &tmp2);
  mpz_set(x, &tmp2);
  return s;
}

/**************************************************************************
 * initialisations                                                        *
 **************************************************************************/

void init_main(void)
{
  bit_array bit;
  long n;
  /* initialise multi-precision integer variables */
  for(n = 0; n <= 6 ; n++)
  { mpz_init(&coeffs[n]); mpz_init(&bc[n]); }
  for(n = 0; n < 3; n++) mpz_init(&kummer[n]);
  mpz_init(&fff);
  mpz_init(&tmp); mpz_init(&tmp2); mpz_init(&tmp3); mpz_init(&ddd);
  mpz_init(&k400); mpz_init(&k310); mpz_init(&k301); mpz_init(&k220);
  mpz_init(&k211); mpz_init(&k202); mpz_init(&k130); mpz_init(&k121);
  mpz_init(&k112); mpz_init(&k103); mpz_init(&k040); mpz_init(&k031);
  mpz_init(&k022); mpz_init(&k013); mpz_init(&k004);
  mpz_init(&x12); mpz_init(&x22); mpz_init(&x32);
  mpz_init(&cpf1); mpz_init(&cpf2); mpz_init(&cpf3);
  mpz_init(&cpfa); mpz_init(&cpfb); mpz_init(&cpfc);

  /* intialise bits[] */
  bit = (bit_array)1;
  for(n = 0; n < LONG_LENGTH; n++) { bits[n] = bit; bit <<= 1; }
  /* initialise inverses[][] */
  init_inverses();
  /* initialise squares[][] */
  init_squares();
}

/**************************************************************************
 * initialise the tables                                                  *
 **************************************************************************/

long invert(long a, long p)
{ /* compute a^(p-2) mod p */
  long n, r = 1, b = a;
  for(n = p-2; n > 1; n >>= 1)
  { if(n&1) { r *= b; r %= p; }
    b *= b; b %= p;
  }
  return (r*b) % p;
}

void init_inverses(void)
{
  long a, pn, p;
  for(pn = 0; pn < NUM_PRIMES; pn++)
  { p = prime[pn];
    for(a = 1; a < p; a++) inverses[pn][a] = invert(a, p);
} }

void init_squares(void)
/* initialise squares[][] */
{
  long a, pn, p;
  for(pn = 0; pn < NUM_PRIMES; pn++)
  {
    p = prime[pn];
    for(a = 0; a < p; a++) squares[pn][a] = p;
    for(a = 0; a < p; a += 2) squares[pn][(a*a)%p] = a;
         /* if p==2 then this does not work! */
  }
  return;
}

/* This is a comparison function needed for sorting in order to determine
   the `best' primes for sieving. */
int compare_entries(const void *a, const void *b)
{
  double diff = (((entry *)a)->r - ((entry *)b)->r);
  return (diff > 0) ? 1 : (diff < 0) ? -1 : 0;
}

/* The rate of a prime p: the fraction of the residue triples (a, b, c)
   mod p that the sieve admits, from the number of admitted classes
   (1, b, c) and the a = 0 part, as init_fmodpsquare counts them. */
static double rate_of(long p, long np0, double np1)
{ return(((double)(p-1) * ((double)np0 + np1) + 1.0)
           / ((double)p * (double)p * (double)p)); }

/* Is (1, b, c) mod p = prime[pn] the image of a point of J(F_p)?  The
   remainder s1*x + s2 of f modulo x^2 - b*x + c gives f(x1)*f(x2) = s2^2
   + b*s1*s2 + c*s1^2 for the roots x1, x2, which must be a square t^2,
   and then one of (y1 +- y2)^2 = 2*s2 + b*s1 +- 2*t must be a square.
   The reductions mod p by the reciprocal, see MODP in the table; every
   operand stays below 2^24.  (The diagonal c = b^2/4 is always admitted.) */
#define MODP(x) ((x) - p * (long)(((x) * recip) >> 32))
static inline int point_on_j(long pn, long p, unsigned long recip,
                             long b, long c)
{
  long s1 = coeffs_mod_p[pn][degree], s2 = coeffs_mod_p[pn][degree-1], t, n;
  for(n = degree-2; n >= 0; n--)
  { t = s1;
    s1 = MODP(s1*b + s2);
    s2 = MODP(coeffs_mod_p[pn][n] + (p-c)*t);
  }
  t = MODP(s2*s2 + MODP(b*s1)*s2 + MODP(c*s1)*s1);
  t = squares[pn][t];
  if(t == p) { return(0); }
  { long u1 = MODP(2*s2 + b*s1 + 2*t), u2 = MODP(u1 + 4*(p - t));
    return(squares[pn][u1] != p || squares[pn][u2] != p);
  }
}

/* Fill is_point_on_j[pn] for the prime p = prime[pn], and return the
   number of admitted classes (1, b, c). */
static long init_table(long pn)
{
  long p = prime[pn], b, c, count = 0;
  unsigned long recip = ((1UL << 32) + p - 1) / p;
  for(b = 0; b < p; b++)
  { long bb = ( ((1 + p*(3-(p&2)))/4) * b*b ) % p;   /* b^2/4 mod p */
    for(c = 0; c < p; c++)
    { if(c == bb || point_on_j(pn, p, recip, b, c))
      { is_point_on_j[pn][b][c] = 1; count++; }
      else
      { is_point_on_j[pn][b][c] = 0; }
    }
  }
  return(count);
}

/* The rates of the primes come from a sample of SAMPLE_PAIRS classes
   (1, b, c) -- all of them for the primes up to 61 -- and the tables are
   filled by init_table for the primes chosen; so the 53 primes of the
   table cost little at set-up.  The sample: every stride-th pair in the
   order b*p + c, the stride below p, so that c runs through all residues
   evenly. */

void init_fmodpsquare(void)
/* initialise is_f_square[][], the rates of the primes, and the tables
   is_point_on_j[][][] of the primes chosen */
{
  long a, pn, p, n, s, np0;
  int sampled[NUM_PRIMES], tabled[NUM_PRIMES];

  for(pn = 0; pn < num_primes; pn++)
  {
    unsigned long recip;
    long i, pp, stride, total, count;
    p = prime[pn];
    pp = p*p;
    recip = ((1UL << 32) + p - 1) / p;
    /* compute coefficients mod p */
    for(n = 0; n <= 6; n++)
      coeffs_mod_p[pn][n] = mpz_fdiv_r_ui(&tmp, &coeffs[n], p);
    /* deal with (0,0,1) */
    has_infinity[pn] = (squares[pn][coeffs_mod_p[pn][6]] == p) ? 0 : 1;
    np0 = 1;
    /* deal with (0,1,a) */
    if((is_f_square[pn][0] = (squares[pn][coeffs_mod_p[pn][0]] == p) ? 0 : 1))
    { np0++; }
    for(a = 1 ; a < p; a++)
    {
      s = coeffs_mod_p[pn][degree];
      for(n = degree-1 ; n >= 0 ; n--)
      { s *= a;
        s += coeffs_mod_p[pn][n];
        s %= p;
      }
      if((is_f_square[pn][a] = (squares[pn][s] == p) ? 0 : 1)) np0++;
    }
    /* deal with (1,b,c): the whole table, or a sample */
    stride = (pp + SAMPLE_PAIRS - 1) / SAMPLE_PAIRS;
    sampled[pn] = (stride > 1);
    tabled[pn] = 0;
    if(stride == 1)
    { count = init_table(pn); total = pp; tabled[pn] = 1; }
    else
    { count = 0; total = 0;
      for(i = 0; i < pp; i += stride)
      { long b = i / p, c = i - b*p;
        long bb = ( ((1 + p*(3-(p&2)))/4) * b*b ) % p;
        total++;
        if(c == bb || point_on_j(pn, p, recip, b, c)) { count++; }
      }
    }
    /* Fill array with info for p */
    prec[pn].p = p;
    prec[pn].n = pn;
    prec[pn].r = rate_of(p, np0, (double)count * (double)pp / (double)total);
    prec[pn].np0 = np0;
  }
  /* a sampled rate of 0 is no proof that no class is admitted: those
     primes get the exact count (a rate of 0 means that there are no
     points at all, see find_points) */
  for(pn = 0; pn < num_primes; pn++)
  { if(sampled[pn] && prec[pn].r == 0.0)
    { prec[pn].r = rate_of(prime[pn], prec[pn].np0, (double)init_table(pn));
      tabled[pn] = 1;
    }
  }
  /* sort the array to get at the best primes */
  qsort(prec, num_primes, sizeof(entry), compare_entries);
  choose_primes();
  /* the tables of the primes chosen */
  for(n = 0; n < sieve_primes3; n++)
  { if(!tabled[pnn[n]]) { init_table(pnn[n]); tabled[pnn[n]] = 1; } }
  if(!quiet)
  { message(4, 0);
    if(sieve_primes3 > 0)
    { message(7, 0); message(8, 0); }
  }
  return;
}

/* The candidates of a stage, ranked by their value for it: the
   reduction of the survivors, -log r, per cycle the prime costs in that
   stage per row at this height (see choose_primes). */
typedef struct { long n; double v; } cand;   /* n: the index in prec */
static int compare_cands(const void *a, const void *b)
{ double diff = ((cand *)b)->v - ((cand *)a)->v;
  return (diff > 0) ? 1 : (diff < 0) ? -1 : 0;
}

/* Choose the primes of the three stages and their numbers, unless -n,
   -M, -N pinned the numbers, and fill pnn[] stage by stage.  The model
   (the COST_ constants) prices a run: per row of W words, n1 passes of W
   word-ANDs, each with the size penalty; the words surviving the first
   stage found and tested one prime at a time until they die; the bits
   surviving the second stage found and tested one prime at a time, and
   the data of the third stage for the rows that have such a bit; the
   exact check of the bits surviving all stages; and at set-up the sieve
   tables (stages 1 and 2, p <= MAX_TABLE_PRIME) and the exact tables
   is_point_on_j of the primes whose rate was sampled.  The rates are
   taken as the densities of the admitted c in a row, the primes as
   independent.  For each number n1 the first stage takes the n1 best
   table primes by their value for the first stage (-log r per cycle of
   a pass plus the amortised table); the second stage the best n2 - n1
   of the table primes left, by their value for the second stage (per
   cycle of the tests on the words the first stage leaves); the third
   stage the best n3 - n2 of all primes left, by their value for it.
   The numbers minimise the modelled cost. */
static double choose_primes_mode(int tube);

/* the two ways of a plain run priced, the cheaper taken (the tube
   itself, and the cost of a row of the box, are not among the terms
   that the choice of the primes sees, so they are added here) */
void choose_primes(void)
{
  double box, tube;
  long pin1 = sieve_primes1, pin2 = sieve_primes2, pin3 = sieve_primes3;
  long p1, p2, p3, pn[NUM_PRIMES], n;
  box = choose_primes_mode(0);
  if(all_points) { return; }
  /* the box's choice kept aside, the pins restored for the tube's */
  p1 = sieve_primes1; p2 = sieve_primes2; p3 = sieve_primes3;
  for(n = 0; n < p3; n++) { pn[n] = pnn[n]; }
  sieve_primes1 = pin1; sieve_primes2 = pin2; sieve_primes3 = pin3;
  tube = choose_primes_mode(1);
  if(tube < box) { tube_mode = 1; return; }
  sieve_primes1 = p1; sieve_primes2 = p2; sieve_primes3 = p3;
  for(n = 0; n < p3; n++) { pnn[n] = pn[n]; }
  return;
}

static double choose_primes_mode(int tube)
{
  long ne = 0, n, n1, n2, n3, k;
  cand list1[NUM_PRIMES], list2[NUM_PRIMES], list3[NUM_PRIMES];
  long best1[NUM_PRIMES], best2[NUM_PRIMES], best3[NUM_PRIMES];
  double W = (double)(2*(height>>LONG_SHIFT) + 2);  /* words per row */
  double rows = (double)(height + 1) * (double)(2*height + 1);
  double bits;                                      /* bits per row */
  double best = -1.0, fixed;
  long b1 = 0, b2 = 0, b3 = 0;
  long n1lo, n1hi, pin1 = sieve_primes1;
  double logr[NUM_PRIMES], tabcost[NUM_PRIMES], initcost[NUM_PRIMES];
  int used[NUM_PRIMES];

  /* in the tube there is no first stage: the words the tube leaves are
     tested one by one, like the survivors of a first stage but without
     an early exit; the cost of the analysis and of a row are fixed */
  if(tube)
  { W = (tube_words > 0.01) ? tube_words : 0.01;
    pin1 = 0;
    fixed = tube_cells * COST_CELL + COST_TROW + W * COST_TWORD;
  }
  else
  { fixed = COST_ROW; }
  bits = (double)LONG_LENGTH * W;
  for(n = 0; n < num_primes; n++)
  { double p = (double)prec[n].p;
    logr[n] = (prec[n].r > 0.0) ? -log(prec[n].r) : 1.0e9;
    tabcost[n] = COST_TABLE * p*p*(p + 1.0);
    initcost[n] = (p*p > (double)SAMPLE_PAIRS) ? COST_INIT * p*p : 0.0;
  }
  /* the first stage's ranking of the table primes */
  for(n = 0; n < num_primes; n++)
  { if(prec[n].p <= MAX_TABLE_PRIME)
    { list1[ne].n = n;
      list1[ne].v = logr[n] / (W * COST_AND * (1.0 + (double)prec[n].p / COST_SIZE)
                               + (tabcost[n] + initcost[n]) / rows);
      ne++;
    }
  }
  qsort(list1, ne, sizeof(cand), compare_cands);
  /* the range of n1: pinned, or 0..ne, never above pinned n2 or n3 */
  n1lo = 0; n1hi = ne;
  if(sieve_primes2 >= 0 && n1hi > sieve_primes2) { n1hi = sieve_primes2; }
  if(sieve_primes3 >= 0 && n1hi > sieve_primes3) { n1hi = sieve_primes3; }
  if(pin1 >= 0)
  { n1lo = n1hi = (pin1 < n1hi) ? pin1 : n1hi; }
  for(n1 = n1lo; n1 <= n1hi; n1++)
  { double cost1 = 0.0, setup1 = 0.0, rho1 = 1.0, s1;
    long n2lo, n2hi, ne2 = 0;
    for(n = 0; n < num_primes; n++) { used[n] = 0; }
    for(k = 0; k < n1; k++)
    { n = list1[k].n;
      used[n] = 1;
      cost1 += W * COST_AND * (1.0 + (double)prec[n].p / COST_SIZE);
      setup1 += tabcost[n] + initcost[n];
      rho1 *= prec[n].r;
    }
    s1 = 1.0 - pow(1.0 - rho1, (double)LONG_LENGTH);  /* words surviving */
    /* the second stage's ranking of the table primes left */
    for(k = 0; k < ne; k++)
    { n = list1[k].n;
      if(!used[n])
      { list2[ne2].n = n;
        list2[ne2].v = logr[n] / (W * s1 * COST_TEST2
                                  + (tabcost[n] + initcost[n]) / rows);
        ne2++;
      }
    }
    qsort(list2, ne2, sizeof(cand), compare_cands);
    n2lo = n1; n2hi = n1 + ne2;
    if(sieve_primes3 >= 0 && n2hi > sieve_primes3) { n2hi = sieve_primes3; }
    if(sieve_primes2 >= 0)
    { n2lo = n2hi = (sieve_primes2 < n2hi) ? sieve_primes2 : n2hi;
      if(n2lo < n1) { n2lo = n2hi = n1; }
    }
    for(n2 = n2lo; n2 <= n2hi; n2++)
    { double cost2, setup2 = setup1, rho2 = rho1, w = s1, rowbit;
      long n3lo, n3hi, ne3 = 0;
      /* the words the first stage leaves, tested until they die */
      cost2 = tube ? 0.0 : W * w * COST_WORD;
      for(k = 0; k < n2 - n1; k++)
      { n = list2[k].n;
        cost2 += tube ? W * COST_TTEST : W * w * COST_TEST2;
        setup2 += tabcost[n] + initcost[n];
        rho2 *= prec[n].r;
        w = 1.0 - pow(1.0 - rho2, (double)LONG_LENGTH);
      }
      rowbit = 1.0 - pow(1.0 - rho2, bits);  /* rows with a bit left */
      /* the third stage's ranking of all the primes left */
      for(n = 0; n < num_primes; n++)
      { int in2 = 0;
        for(k = 0; k < n2 - n1; k++) { if(list2[k].n == n) { in2 = 1; } }
        if(!used[n] && !in2)
        { list3[ne3].n = n;
          list3[ne3].v = logr[n] / (bits * rho2 * COST_TEST3
                                          + rowbit * COST_ROW3
                                          + initcost[n] / rows);
          ne3++;
        }
      }
      qsort(list3, ne3, sizeof(cand), compare_cands);
      n3lo = n2; n3hi = n2 + ne3;
      if(sieve_primes3 >= 0)
      { n3lo = n3hi = (sieve_primes3 < n3hi) ? sieve_primes3 : n3hi;
        if(n3lo < n2) { n3lo = n3hi = n2; }
      }
      { double cost3 = bits * rho2 * COST_BIT, setup3 = setup2, rho = rho2;
        for(n3 = n2; n3 <= n3hi; n3++)
        { if(n3 >= n3lo)
          { double total = rows * (fixed + cost1 + cost2 + cost3
                                   + rowbit * (double)(n3 - n2) * COST_ROW3
                                   + bits * rho * COST_EXACT)
                           + setup3;
            if(best < 0.0 || total < best)
            { best = total; b1 = n1; b2 = n2; b3 = n3;
              for(k = 0; k < n1; k++) { best1[k] = list1[k].n; }
              for(k = 0; k < n2 - n1; k++) { best2[k] = list2[k].n; }
              for(k = 0; k < n3 - n2; k++) { best3[k] = list3[k].n; }
            }
          }
          if(n3 < n3hi)
          { n = list3[n3 - n2].n;
            cost3 += bits * rho * COST_TEST3;
            setup3 += initcost[n];
            rho *= prec[n].r;
          }
        }
      }
    }
  }
  sieve_primes1 = b1; sieve_primes2 = b2; sieve_primes3 = b3;
  for(k = 0; k < b1; k++) { pnn[k] = prec[best1[k]].n; }
  for(k = 0; k < b2 - b1; k++) { pnn[b1 + k] = prec[best2[k]].n; }
  for(k = 0; k < b3 - b2; k++) { pnn[b2 + k] = prec[best3[k]].n; }
  return(best);
}

/* help is an array for temporarily storing the sieving information.
   Bit j of help[a][b][c0] says whether (a : b : c) are the first three
   coordinates of a point on K(F_p), where c = c0*LONG_LENGTH + j.
   c runs from 0 to k*p - 1, where k*p is the least multiple of p exceeding
   LONG_LENGTH. */
bit_array help[MAX_TABLE_PRIME][MAX_TABLE_PRIME][MAX_TABLE_PRIME / LONG_LENGTH + 2];

/* allocate and initalise the sieve tables */
void init_sieve(void)
{
  long a, b, c, i, pn, p, aa, ab, k, n, kp;

#if (DEBUG >= 2)
  printf("\n sieve:\n");
#endif
  for(pn = 0; pn < sieve_primes2; pn++)
  {
    n = pnn[pn];
    p = prime[n];
    k = (LONG_LENGTH + p - 1)/p;
    kp = k*p;
    sieve_tab[pn] = (bit_array *)malloc(p*p*(p+1)*sizeof(bit_array));
    if(sieve_tab[pn] == NULL) { error(7); }
    /* determine which triples (a,b,c) are excluded mod p */
    aa = p>>LONG_SHIFT;
    if(aa == 0) aa = 1;

#if (DEBUG >= 2)
    printf("  p = %ld, aa = %ld\n", p, aa);
#endif

    for(a = 0; a < p; a++)
      for(b = 0; b < p; b++)
        for(c = 0; c <= aa; c++) help[a][b][c] = zero;
    /* a = 0, b = 0 */
    for(c = 0; c <= aa; c++) help[0][0][c] = ~zero;
    /* a = 0, b != 0 */
    if(has_infinity[n])
      for(c = 0; c < p; c++)
        if(is_f_square[n][c])
          for(b = 1; b < p; b++)
            for(i = (b*c)%p; i < kp; i += p)
              help[0][b][i>>LONG_SHIFT] |= bits[i & LONG_MASK];
    /* a != 0 */
    for(b = 0; b < p; b++)
      for(c = 0; c < p; c++)
        if(is_point_on_j[n][b][c])
          for(a = 1; a < p; a++)
          { ab = (a*b)%p;
            for(i = (a*c)%p; i < kp; i += p)
              help[a][ab][i>>LONG_SHIFT] |= bits[i & LONG_MASK];
	  }

#if (DEBUG >= 3)
    printf(" help(%ld):\n", p);
    for(a = 0; a < p; a++)
    { printf(" a = %3ld:\n", a);
      for(b = 0; b < p; b++)
      { printf("  b = %3ld: ", b);
        for(c = 0; c <= aa; c++)
	  printf(" %8.8lx", help[a][b][c]);
        printf("\n");
    } }
#endif

    /* fill the bit pattern from help[][][] into the table of p.
	Word c0 of the row of (a, b) has the same semantics as help[a][b][c0],
	but here, c0 runs from 0 to p-1 and all bits are filled. */
    for(a = 0; a < p; a++)
      for(b = 0; b < p; b++)
      { bit_array *si = &sieve_tab[pn][(a*p + b)*(p+1)];
        bit_array *he = help[a][b];
        long p1 = (LONG_LENGTH/p + 1) * p;
        long diff_shift = p1 & LONG_MASK;
        long diff = LONG_LENGTH - diff_shift;
        bit_array diff_mask = ~(bit_array)(-1L<<diff);
        long c1;
        long wp = p1>>LONG_SHIFT;
        /* copy the first chunk from help[a][b][] into the row */
        for(c = 0; c < wp; c++) si[c] = he[c];
        /* now keep repeating the bit pattern, rotating it in help */
        for(c1 = c ; c < p; c++)
        {
	  he[c1] |= (he[(c1 == wp) ? 0 : c1 + 1] & diff_mask)<<diff_shift;
	  si[c] = he[c1];
	  if(c1 == wp) c1 = 0; else c1++;
	  he[c1] >>= diff;
        }
        si[p] = si[0];
      }

#if (DEBUG >= 3)
    printf(" sieve(%ld):\n", p);
    for(a = 0; a < p; a++)
    { printf(" a = %3ld:\n", a);
      for(b = 0; b < p; b++)
      { printf("  b = %3ld: ", b);
        for(c = 0; c < p; c++)
	  printf(" %8.8lx", sieve_tab[pn][(a*p + b)*(p+1) + c]);
        printf("\n");
    } }
#endif

  } /* end for pn */
  return;
}

/**************************************************************************
 * find points by looping over the first two coordinates and sieving      *
 * on the third                                                           *
 **************************************************************************/

int find_double_points(void)
{ /* Find `double points', i.e. points of the form [2*P - O] */
  /* Use a simple version of the sieve, since height will be no more
     than ~ 100 */
  long wheight = floor(sqrt(height));
  long a, b, n;
  struct stype {long p; long ap; unsigned char *ptr;} ssp0[NUM_PRIMES];
  struct stype *ssp;
  for(n = 0; n < sieve_primes3; n++)
    ssp0[n].p = prime[pnn[n]];
  if(degree == 5 && mpz_cmp_si(&coeffs[5], 1) == 0)
  { /* can use only squares for a */
    long aa;
    for(aa = 1; (a = aa*aa) <= wheight; aa++)
    { long bmax;
      for(ssp = &ssp0[0], n = 0; n < sieve_primes3; n++, ssp++)
      { long pn = pnn[n], ap0 = a%(ssp->p);
        if(ap0 == 0)
        { ssp->ap = 0;
          ssp->ptr = &has_infinity[pn];
        }
        else
        { ssp->ap = inverses[pn][ap0];
          ssp->ptr = &is_f_square[pn][0];
      } }
      bmax = height/(2*a);
      if(wheight < bmax) bmax = wheight;
      for(b = -bmax; b <= bmax; b++)
      { if(relprime(a, b)) /* a is positive */
        { for(ssp = &ssp0[0], n = sieve_primes3; n; n--, ssp++)
          { long p = ssp->p, bp = b%p;
            if(bp < 0) bp += p;
            if(!ssp->ptr[((ssp->ap)*bp)%p]) goto nextb1;
          }
          /* Check if we really have got a point */
          if(check_one_point(a*a, 2*a*b, b*b) && one_point) return(1);
        }
nextb1: ;
      }
    }
  }
  else
    for(a = 1; a <= wheight; a++)
    { long bmax;
      for(ssp = &ssp0[0], n = 0; n < sieve_primes3; n++, ssp++)
      { long pn = pnn[n], ap0 = a%(ssp->p);
        if(ap0 == 0)
        { if(!has_infinity[pn]) goto nexta2;
          ssp->ap = 0;
          ssp->ptr = &has_infinity[pn];
        }
        else
        { ssp->ap = inverses[pn][ap0];
          ssp->ptr = &is_f_square[pn][0];
      } }
      bmax = height/(2*a);
      if(wheight < bmax) bmax = wheight;
      for(b = -bmax; b <= bmax; b++)
      { if(relprime(a, b)) /* a is positive */
        { for(ssp = &ssp0[0], n = sieve_primes3; n; n--, ssp++)
          { long p = ssp->p, bp = b%p;
            if(bp < 0) bp += p;
            if(!ssp->ptr[((ssp->ap)*bp)%p]) goto nextb2;
          }
          /* Check if we really have got a point */
          if(check_one_point(a*a, 2*a*b, b*b) && one_point) return(1);
        }
nextb2: ;
      }
nexta2: ;
    }
  return(0);
}

/* The rows (a, b) with p | a and p not | b, for a sieving prime p modulo
   which the leading coefficient is not a square: no point of J reduces to
   (0 : b : c) mod p with b nonzero then (the curve has no point at
   infinity over F_p), and the sieve would empty every such row.  So for
   a given a the loop over b takes only the multiples of the product of
   these primes, which this function returns -- or height + 1 when the
   product exceeds the height bound, so that b = 0 is the only row (and
   a = 0, whose b starts at 1, has none). */
static long row_step(long a)
{
  long m = 1, n;
  for(n = 0; n < sieve_primes3; n++)
  { long pn = pnn[n], p = prime[pn];
    if(!has_infinity[pn] && a % p == 0)
    { m *= p;
      if(m > height) { return(height + 1); }
    }
  }
  return(m);
}

void find_points(void)
{
  long a, b;

  /* initialise is_f_square[][] */
  init_fmodpsquare();
  /* allocate and initalise the sieve tables */
  init_sieve();
  if(sieve_primes3 > 0 && prec[0].r == 0.0)
  { if(!quiet) message(1,0); return; }
  if(sieve_primes3 == 0) { tube_mode = 0; }
  /* the bit array and the tables of the chunks */
  init_sift();
  /* deal with (0, 0, c, d) */
  if(degree == 6 && mpz_perfect_square_p(&coeffs[6]))
  { mpz_mul(&fff, &coeffs[5], &coeffs[5]);
    mpz_mul_ui(&tmp, &coeffs[6], 4);
    mpz_mul(&tmp2, &tmp, &coeffs[4]);
    mpz_sub(&fff, &fff, &tmp2);
    if(check_one_point_final(0, 0, 1, &fff, &tmp) && one_point) return;
  }
  /* deal with `double points' */
  if(find_double_points() && one_point) return;
  if(degree == 5 && mpz_cmp_si(&coeffs[5], 1) == 0)
  { /* can take only squares for the first coordinate */
    long aa;
    for(a = 0; (aa = a*a) <= height; a++)
    { long m = row_step(aa), b0 = (a == 0) ? m : -(height/m)*m;
      if(tube_mode)
      { if(sift_tube(aa, b0, m) && one_point) return; }
      else
      { for(b = b0; b <= height; b += m)
        {
#ifdef VERBOSE
          printf(" a = %ld, b = %ld\n", aa, b);
#endif
          if(sift(aa, b) && one_point) return;
        }
      }
    }
  }
  else
  { for(a = 0; a <= height; a++)
    { long m = row_step(a), b0 = (a == 0) ? m : -(height/m)*m;
      if(tube_mode)
      { if(sift_tube(a, b0, m) && one_point) return; }
      else
      { for(b = b0; b <= height; b += m)
        {
#ifdef VERBOSE
          printf(" a = %ld, b = %ld\n", a, b);
#endif
          if(sift(a, b) && one_point) return;
        }
      }
    }
  }
  return;
}

/**************************************************************************
 * Compute coefficients of the Kummer surface equation.                   *
 * k<i1><i2><i3> is the coefficient of x1^i1*x2^i2*x3^i3.                 *
 **************************************************************************/

void kummer_init()
{ /* -4*x1^4*f0*f2 + x1^4*f1^2 */
  mpz_mul(&k400, &coeffs[1], &coeffs[1]); mpz_mul(&tmp, &coeffs[0], &coeffs[2]);
  mpz_mul_ui(&tmp, &tmp, 4); mpz_sub(&k400, &k400, &tmp);
  /* - 4*x1^3*x2*f0*f3 */
  mpz_mul(&k310, &coeffs[0], &coeffs[3]); mpz_mul_si(&k310, &k310, -4);
  /* - 2*x1^3*x3*f1*f3 */
  mpz_mul(&k301, &coeffs[1], &coeffs[3]); mpz_mul_si(&k301, &k301, -2);
  /* - 4*x1^2*x2^2*f0*f4 */
  mpz_mul(&k220, &coeffs[0], &coeffs[4]); mpz_mul_si(&k220, &k220, -4);
  /* + 4*x1^2*x2*x3*f0*f5 - 4*x1^2*x2*x3*f1*f4 */
  mpz_mul(&k211, &coeffs[0], &coeffs[5]); mpz_mul(&tmp, &coeffs[1], &coeffs[4]);
  mpz_sub(&k211, &k211, &tmp); mpz_mul_ui(&k211, &k211, 4);
  /* - 4*x1^2*x3^2*f0*f6 + 2*x1^2*x3^2*f1*f5 - 4*x1^2*x3^2*f2*f4
     + x1^2*x3^2*f3^2 */
  mpz_mul(&k202, &coeffs[3], &coeffs[3]); mpz_mul(&tmp, &coeffs[2], &coeffs[4]);
  mpz_mul_ui(&tmp, &tmp, 4), mpz_sub(&k202, &k202, &tmp);
  mpz_mul(&tmp, &coeffs[1], &coeffs[5]); mpz_mul_ui(&tmp, &tmp, 2);
  mpz_add(&k202, &k202, &tmp); mpz_mul(&tmp, &coeffs[0], &coeffs[6]);
  mpz_mul_ui(&tmp, &tmp, 4); mpz_sub(&k202, &k202, &tmp);
  /* - 4*x1*x2^3*f0*f5 */
  mpz_mul(&k130, &coeffs[0], &coeffs[5]); mpz_mul_si(&k130, &k130, -4);
  /* + 8*x1*x2^2*x3*f0*f6 - 4*x1*x2^2*x3*f1*f5 */
  mpz_mul(&k121, &coeffs[0], &coeffs[6]); mpz_mul_ui(&k121, &k121, 8);
  mpz_mul(&tmp, &coeffs[1], &coeffs[5]); mpz_mul_ui(&tmp, &tmp, 4);
  mpz_sub(&k121, &k121, &tmp);
  /* + 4*x1*x2*x3^2*f1*f6 - 4*x1*x2*x3^2*f2*f5 */
  mpz_mul(&k112, &coeffs[1], &coeffs[6]); mpz_mul(&tmp, &coeffs[2], &coeffs[5]);
  mpz_sub(&k112, &k112, &tmp); mpz_mul_ui(&k112, &k112, 4);
  /* - 2*x1*x3^3*f3*f5 */
  mpz_mul(&k103, &coeffs[3], &coeffs[5]); mpz_mul_si(&k103, &k103, -2);
  /* - 4*x2^4*f0*f6 */
  mpz_mul(&k040, &coeffs[0], &coeffs[6]); mpz_mul_si(&k040, &k040, -4);
  /* - 4*x2^3*x3*f1*f6 */
  mpz_mul(&k031, &coeffs[1], &coeffs[6]); mpz_mul_si(&k031, &k031, -4);
  /* - 4*x2^2*x3^2*f2*f6 */
  mpz_mul(&k022, &coeffs[2], &coeffs[6]); mpz_mul_si(&k022, &k022, -4);
  /* - 4*x2*x3^3*f3*f6 */
  mpz_mul(&k013, &coeffs[3], &coeffs[6]); mpz_mul_si(&k013, &k013, -4);
  /* - 4*x3^4*f4*f6 + x3^4*f5^2 */
  mpz_mul(&k004, &coeffs[5], &coeffs[5]); mpz_mul(&tmp, &coeffs[4], &coeffs[6]);
  mpz_mul_ui(&tmp, &tmp, 4); mpz_sub(&k004, &k004, &tmp);
}

/* the coefficients as doubles for the tube, and whether it is used; the
   words per row it leaves are sampled over a few hundred rows spread
   over the box (deterministically), for the cost model */
void tube_init(void)
{
  long n;
  kd400 = mpz_get_d(&k400); kd310 = mpz_get_d(&k310);
  kd301 = mpz_get_d(&k301); kd220 = mpz_get_d(&k220);
  kd211 = mpz_get_d(&k211); kd202 = mpz_get_d(&k202);
  kd130 = mpz_get_d(&k130); kd121 = mpz_get_d(&k121);
  kd112 = mpz_get_d(&k112); kd103 = mpz_get_d(&k103);
  kd040 = mpz_get_d(&k040); kd031 = mpz_get_d(&k031);
  kd022 = mpz_get_d(&k022); kd013 = mpz_get_d(&k013);
  kd004 = mpz_get_d(&k004);
  for(n = 0; n <= 6; n++) { fd[n] = mpz_get_d(&coeffs[n]); }
  /* a plain run may use the tube; choose_primes decides */
  tube_mode = 0;
  if(!all_points) { tube_sample(&tube_words, &tube_cells); }
  return;
}

/* Compute the coefficients of the quadratic equation for the fourth
 * Kummer coordinate, given the first three coordinates a, b, c.
 * We require that max{|a|^2, |b|^2, |c|^2} < 2^(LONG_LENGTH-1),
 * so that products of two coordinates fit into a (signed) long. */
static inline void kummer_eqn(long a, long b, long c)
{ long aa = a*a, bb = b*b, cc = c*c, ab = a*b, ac = a*c, bc = b*c;
  /* Constant term in kummer[0]: */
  /* k400*a^4 + k310*a^3*b + k301*a^3*c + k220*a^2*b^2 + k211*a^2*b*c + k202*a^2*c^2
   *  + k130*a*b^3 + k121*a*b^2*c + k112*a*b*c^2 + k103*a*c^3
   *  + k040*b^4 + k031*b^3*c + k022*b^2*c^2 + k013*b*c^3 + k004*c^4
   * = (k400*a^2 + k310*a*b + k301*a*c + k220*b^2 + k211*b*c + k202*c^2)*a^2
   *    + (k130*a*b + k121*a*c + k040*b^2 + k031*b*c + k022*c^2)*b^2
   *    + (k112*a*b + k103*a*c + k013*b*c + k004*c^2)*c^2 */
  mpz_mul_si(&kummer[0], &k400, aa);
  mpz_mul_si(&tmp, &k310, ab); mpz_add(&kummer[0], &kummer[0], &tmp);
  mpz_mul_si(&tmp, &k301, ac); mpz_add(&kummer[0], &kummer[0], &tmp);
  mpz_mul_si(&tmp, &k220, bb); mpz_add(&kummer[0], &kummer[0], &tmp);
  mpz_mul_si(&tmp, &k211, bc); mpz_add(&kummer[0], &kummer[0], &tmp);
  mpz_mul_si(&tmp, &k202, cc); mpz_add(&kummer[0], &kummer[0], &tmp);
  mpz_mul_si(&kummer[0], &kummer[0], aa);

  mpz_mul_si(&kummer[1], &k130, ab);
  mpz_mul_si(&tmp, &k121, ac); mpz_add(&kummer[1], &kummer[1], &tmp);
  mpz_mul_si(&tmp, &k040, bb); mpz_add(&kummer[1], &kummer[1], &tmp);
  mpz_mul_si(&tmp, &k031, bc); mpz_add(&kummer[1], &kummer[1], &tmp);
  mpz_mul_si(&tmp, &k022, cc); mpz_add(&kummer[1], &kummer[1], &tmp);
  mpz_mul_si(&kummer[1], &kummer[1], bb);

  mpz_mul_si(&kummer[2], &k112, ab);
  mpz_mul_si(&tmp, &k103, ac); mpz_add(&kummer[2], &kummer[2], &tmp);
  mpz_mul_si(&tmp, &k013, bc); mpz_add(&kummer[2], &kummer[2], &tmp);
  mpz_mul_si(&tmp, &k004, cc); mpz_add(&kummer[2], &kummer[2], &tmp);
  mpz_mul_si(&kummer[2], &kummer[2], cc);

  mpz_add(&kummer[0], &kummer[0], &kummer[1]);
  mpz_add(&kummer[0], &kummer[0], &kummer[2]);

  /* Coefficient of linear term: */
  /* -4*a^3*f0 - 2*a^2*b*f1 - 4*a^2*c*f2 - 2*a*b*c*f3 - 4*a*c^2*f4 - 2*b*c^2*f5 - 4*c^3*f6
   *  = (f0*a^2 + f2*a*c)*(-4*a) + (f1*a^2 + f3*a*c + f5*c^2)*(-2*b) + (f4*ac + f6*c^2)*(-4*c) */
  mpz_mul_si(&kummer[1], &coeffs[0], aa);
  mpz_mul_si(&tmp, &coeffs[2], ac); mpz_add(&kummer[1], &kummer[1], &tmp);
  mpz_mul_si(&kummer[1], &kummer[1], -4*a);

  mpz_mul_si(&kummer[2], &coeffs[1], aa);
  mpz_mul_si(&tmp, &coeffs[3], ac); mpz_add(&kummer[2], &kummer[2], &tmp);
  mpz_mul_si(&tmp, &coeffs[5], cc); mpz_add(&kummer[2], &kummer[2], &tmp);
  mpz_mul_si(&kummer[2], &kummer[2], -2*b);
  mpz_add(&kummer[1], &kummer[1], &kummer[2]);

  mpz_mul_si(&kummer[2], &coeffs[4], ac);
  mpz_mul_si(&tmp, &coeffs[6], cc); mpz_add(&kummer[2], &kummer[2], &tmp);
  mpz_mul_si(&kummer[2], &kummer[2], -4*c);
  mpz_add(&kummer[1], &kummer[1], &kummer[2]);

  /* Leading coefficient: */
  /* b^2-4*a*c */
  mpz_set_si(&tmp, ac); mpz_mul_ui(&tmp, &tmp, 4);
  mpz_set_si(&kummer[2], bb);
  mpz_sub(&kummer[2], &kummer[2], &tmp);
}

/**************************************************************************
 * Check if two integers are coprime.                                     *
 **************************************************************************/

static inline int relprime(long n, long m)
{
  /* n is always positive here */
  if(m == 0) { return(n == 1); }
  if(m < 0) { m = -m; }
  if(!(m & 1)) /* m is even */
  { if(!(n & 1)) { return(0); } /* n is also even */
    m >>= 1; while(!(m & 1)) { m >>= 1; } /* n odd: replace m by odd part */
  }
  while(!(n & 1)) { n >>= 1; } /* replace n by odd part */
  /* successively subtract the smaller from the larger
   * and replace the result by its odd part,
   * until both are equal (to their gcd) */
  while(n != m)
  { if(n > m)
    { n -= m; n >>= 1; while(!(n & 1)) { n >>= 1; } }
    else
    { m -= n; m >>= 1; while(!(m & 1)) { m >>= 1; } }
  }
  return(m == 1);
}


/**************************************************************************
 * check a `survivor' of the sieve if it really gives a point             *
 **************************************************************************/

int check_lifts(long a, long b, long c, long d)
{ /* test whether (a : b : c : d) in K(Q) comes from a point on J(Q) */
  long k;

  if(a == 0)
    if(b == 0)
      if(c == 0)
        /* point is the origin */
        return(1);
      else
        /* (0 : 0 : c : d) : must have f6 a square */
        return(mpz_perfect_square_p(&coeffs[6]));
    else /* b /= 0 */
    { /* (0 : b : c : d): f6 must be zero or a square, and
         f(c/b) must be a square */
      if(!(mpz_cmp_si(&coeffs[6], 0) == 0 || mpz_perfect_square_p(&coeffs[6])))
        return(0);
      /* compute entries bc[k] = coeffs[k] * b^(degree-k), k < degree */
      mpz_set_si(&tmp, 1);
      for(k = 5; k >= 0; k--)
      { mpz_mul_ui(&tmp, &tmp, (unsigned long)b);
        mpz_mul(&bc[k], &coeffs[k], &tmp);
      }
      /* now compute b^6*f(c/b) = f6*c^6 + f5*c^5*b + ... + f0*b^6 */
      mpz_set(&fff, &coeffs[6]);
      for(k = 5; k >= 0; k--)
      {
        mpz_mul_si(&fff, &fff, c);
        mpz_add(&fff, &fff, &bc[k]);
      }
      return(mpz_cmp_si(&fff, 0) == 0 || mpz_perfect_square_p(&fff));
    }
  else /* a /= 0 */
  { /* (a : b : c : d): must have A^2 == 0 && (B^2 == 0 or a sqaure)
       or A^2 a non-zero square, where
       A^2 = k1^3*k4 + f2*k1^4 + f3*k1^3*k2 + f4*k1^2*k2^2
             + f5*k1*k2*(k2^2-k1*k3) + f6*(k2^2-k1*k3)^2     and
       B^2 = k1^2*k3*k4 + f0*k1^4 + f4*k1^2*k3^2 + f5*k1*k2*k3^2
              + f6*k2^2*k3^2 */
    /* compute A^2 */
    mpz_set_si(&x12, a*a); mpz_set_si(&x22, b*b); mpz_set_si(&x32, c*c);

    mpz_mul_si(&tmp, &x12, a*d);

    mpz_mul(&tmp2, &x12, &x12); mpz_mul(&tmp2, &tmp2, &coeffs[2]);
    mpz_add(&tmp, &tmp, &tmp2);

    mpz_mul_si(&tmp2, &x12, a*b); mpz_mul(&tmp2, &tmp2, &coeffs[3]);
    mpz_add(&tmp, &tmp, &tmp2);

    mpz_mul(&tmp2, &x12, &x22); mpz_mul(&tmp2, &tmp2, &coeffs[4]);
    mpz_add(&tmp, &tmp, &tmp2);

    /* b^2 - a*c in tmp3: it can exceed a long, unlike the products */
    mpz_set_si(&tmp3, b); mpz_mul_si(&tmp3, &tmp3, b);
    mpz_set_si(&tmp2, a); mpz_mul_si(&tmp2, &tmp2, c);
    mpz_sub(&tmp3, &tmp3, &tmp2);
    mpz_mul(&tmp2, &coeffs[5], &tmp3);
    mpz_mul_si(&tmp2, &tmp2, a*b); mpz_add(&tmp, &tmp, &tmp2);

    mpz_mul(&tmp2, &coeffs[6], &tmp3); mpz_mul(&tmp2, &tmp2, &tmp3);
    mpz_add(&tmp, &tmp, &tmp2);
#ifdef VERBOSE
    printf("  A^2 = %s\n", mpz_get_str((char *) 0, 10, &tmp));
#endif
    /* test if A^2 is a square */
    if(mpz_cmp_si(&tmp, 0) != 0)
      return(mpz_perfect_square_p(&tmp));
    /* A == 0: Compute B^2 */
    mpz_mul_si(&tmp, &x12, c*d);

    mpz_mul(&tmp2, &x12, &x12); mpz_mul(&tmp2, &tmp2, &coeffs[0]);
    mpz_add(&tmp, &tmp, &tmp2);

    mpz_mul(&tmp2, &x12, &x32); mpz_mul(&tmp2, &tmp2, &coeffs[4]);
    mpz_add(&tmp, &tmp, &tmp2);

    mpz_mul(&tmp2, &coeffs[5], &x32); mpz_mul_si(&tmp2, &tmp2, a*b);
    mpz_add(&tmp, &tmp, &tmp2);

    mpz_mul(&tmp2, &x22, &x32); mpz_mul(&tmp2, &tmp2, &coeffs[6]);
    mpz_add(&tmp, &tmp, &tmp2);
#ifdef VERBOSE
    printf("  B^2 = %s\n", mpz_get_str((char *) 0, 10, &tmp));
#endif
    return(mpz_cmp_si(&tmp, 0) == 0 || mpz_perfect_square_p(&tmp));
  }
}

int check_lifts_mpz(MP_INT *a, MP_INT *b, MP_INT *c, MP_INT *d)
{ /* test whether (a : b : c : d) in K(Q) comes from a point on J(Q) */
  long k;

  if(mpz_cmp_si(a, 0) == 0)
    if(mpz_cmp_si(b, 0) == 0)
      if(mpz_cmp_si(c, 0) == 0)
        /* point is the origin */
        return(1);
      else
        /* (0 : 0 : c : d) : must have f6 a square */
        return(mpz_perfect_square_p(&coeffs[6]));
    else /* b /= 0 */
    { /* (0 : b : c : d): f6 must be zero or a square, and
         f(c/b) must be a square */
      if(!(mpz_cmp_si(&coeffs[6], 0) == 0 || mpz_perfect_square_p(&coeffs[6])))
        return(0);
      /* compute entries bc[k] = coeffs[k] * b^(degree-k), k < degree */
      mpz_set_si(&tmp, 1);
      for(k = 5; k >= 0; k--)
      { mpz_mul(&tmp, &tmp, b);
        mpz_mul(&bc[k], &coeffs[k], &tmp);
      }
      /* now compute b^6*f(c/b) = f6*c^6 + f5*c^5*b + ... + f0*b^6 */
      mpz_set(&fff, &coeffs[6]);
      for(k = 5; k >= 0; k--)
      {
        mpz_mul(&fff, &fff, c);
        mpz_add(&fff, &fff, &bc[k]);
      }
      return(mpz_cmp_si(&fff, 0) == 0 || mpz_perfect_square_p(&fff));
    }
  else /* a /= 0 */
  { /* (a : b : c : d): must have A^2 == 0 && (B^2 == 0 or a sqaure)
       or A^2 a non-zero square, where
       A^2 = k1^3*k4 + f2*k1^4 + f3*k1^3*k2 + f4*k1^2*k2^2
             + f5*k1*k2*(k2^2-k1*k3) + f6*(k2^2-k1*k3)^2     and
       B^2 = k1^2*k3*k4 + f0*k1^4 + f4*k1^2*k3^2 + f5*k1*k2*k3^2
              + f6*k2^2*k3^2 */
    /* compute A^2 */
    mpz_mul(&x12, a, a); mpz_mul(&x22, b, b); mpz_mul(&x32, c, c);

    mpz_mul(&tmp, &x12, a); mpz_mul(&tmp, &tmp, d);

    mpz_mul(&tmp2, &x12, &x12); mpz_mul(&tmp2, &tmp2, &coeffs[2]);
    mpz_add(&tmp, &tmp, &tmp2);

    mpz_mul(&tmp2, &x12, a); mpz_mul(&tmp2, &tmp2, b); mpz_mul(&tmp2, &tmp2, &coeffs[3]);
    mpz_add(&tmp, &tmp, &tmp2);

    mpz_mul(&tmp2, &x12, &x22); mpz_mul(&tmp2, &tmp2, &coeffs[4]);
    mpz_add(&tmp, &tmp, &tmp2);

    mpz_mul(&tmp3, a, c); mpz_sub(&tmp3, &x22, &tmp3); /* k = b*b - a*c; */
    mpz_mul(&tmp2, &coeffs[5], &tmp3);
    mpz_mul(&tmp2, &tmp2, a); mpz_mul(&tmp2, &tmp2, b); mpz_add(&tmp, &tmp, &tmp2);

    mpz_mul(&tmp2, &coeffs[6], &tmp3); mpz_mul(&tmp2, &tmp2, &tmp3);
    mpz_add(&tmp, &tmp, &tmp2);
#ifdef VERBOSE
    printf("  A^2 = %s\n", mpz_get_str((char *) 0, 10, &tmp));
#endif
    /* test if A^2 is a square */
    if(mpz_cmp_si(&tmp, 0) != 0)
      return(mpz_perfect_square_p(&tmp));
    /* A == 0: Compute B^2 */
    mpz_mul(&tmp, &x12, c); mpz_mul(&tmp, &tmp, d);

    mpz_mul(&tmp2, &x12, &x12); mpz_mul(&tmp2, &tmp2, &coeffs[0]);
    mpz_add(&tmp, &tmp, &tmp2);

    mpz_mul(&tmp2, &x12, &x32); mpz_mul(&tmp2, &tmp2, &coeffs[4]);
    mpz_add(&tmp, &tmp, &tmp2);

    mpz_mul(&tmp2, &coeffs[5], &x32); mpz_mul(&tmp2, &tmp2, a); mpz_mul(&tmp2, &tmp2, b);
    mpz_add(&tmp, &tmp, &tmp2);

    mpz_mul(&tmp2, &x22, &x32); mpz_mul(&tmp2, &tmp2, &coeffs[6]);
    mpz_add(&tmp, &tmp, &tmp2);
#ifdef VERBOSE
    printf("  B^2 = %s\n", mpz_get_str((char *) 0, 10, &tmp));
#endif
    return(mpz_cmp_si(&tmp, 0) == 0 || mpz_perfect_square_p(&tmp));
  }
}

void printf_mpz(const char *format, MP_INT *a, MP_INT *b, MP_INT *c, MP_INT *d)
{ /* print the four coordinates with the user's format, in which each %ld
     becomes the %Zd of gmp_printf (the two are of the same length) */
  long l = strlen(format), i, j;
  char *fmt = malloc(l + 1);
  for(i = 0, j = 0; i < l; i++, j++)
  { fmt[j] = format[i];
    if(format[i] == '%' && format[i+1] == 'l' && format[i+2] == 'd')
    { fmt[++j] = 'Z'; fmt[++j] = 'd'; i += 2; }
  }
  fmt[j] = 0;
  gmp_printf(fmt, a, b, c, d);
  free(fmt);
}

int check_one_point_final(long a, long b, long c, MP_INT *d1, MP_INT *d2)
{ /* Given a, b, c and d = d1/d2, check if this gives a point satisfying
     the height condition. If so, print and count it. */
  long m = (a > labs(b)) ? ((a > labs(c)) ? a : labs(c))
                         : ((labs(b) > labs(c)) ? labs(b) : labs(c));
  long h = bound/m;
  mpz_gcd(&cpf3, d1, d2);
  mpz_divexact(&cpf1, d1, &cpf3);
  mpz_divexact(&cpf2, d2, &cpf3);
  if(mpz_cmp_si(&cpf1, bound) <= 0 && mpz_cmp_si(&cpf1, -bound) >= 0
      && mpz_cmp_si(&cpf2, h) <= 0 && mpz_cmp_si(&cpf2, -h) >= 0)
  { long g = mpz_get_si(&cpf2), d = mpz_get_si(&cpf1);
    if(g < 0) { g = -g; d = -d; }
    a *= g; b *= g; c *= g;
#ifdef VERBOSE
    printf("  coordinates = (%ld : %ld : %ld : %ld): ", a, b, c, d);
#endif
    if(check_lifts(a, b, c, d))
    {
#ifdef VERBOSE
      printf("lifts.\n");
#endif
      printf(print_format, a, b, c, d);
      total++;
      return(1);
    }
    else
    {
#ifdef VERBOSE
      printf("does not lift.\n");
#endif
      return(0);
    }
  }
  else if(all_points) /* no bound for the height and machine-size is not sufficient */
  {
    /* scale by denominator &cpf2; fourth coordinate is numerator &cpf1 */
#ifdef VERBOSE
    printf("  coordinates exceed machine size ");
#endif
    mpz_mul_si(&cpfa, &cpf2, a);
    mpz_mul_si(&cpfb, &cpf2, b);
    mpz_mul_si(&cpfc, &cpf2, c);
    if(check_lifts_mpz(&cpfa, &cpfb, &cpfc, &cpf1))
    {
#ifdef VERBOSE
      printf("lifts.\n");
#endif
      printf_mpz(print_format, &cpfa, &cpfb, &cpfc, &cpf1);
      total++;
      return(1);
    }
    else
    {
#ifdef VERBOSE
      printf("does not lift.\n");
#endif
      return(0);
    }
  }
  else
  {
#ifdef VERBOSE
    printf("  height is too big.\n");
#endif
    return(0);
  }
}

int check_one_point(long a, long b, long c)
{ /* Check if there is some d such that (a : b : c : d) is a point on
     K  that lifts to J(Q). If so, print the point(s) that satisfy the
     height bound and return 1, if we are only looking for one point,
     else return 0. */
#ifdef VERBOSE
  printf("\n check_one_point(%ld, %ld, %ld):\n", a, b, c);
#endif
  kummer_eqn(a, b, c);
  mpz_mul(&tmp, &kummer[1], &kummer[1]);
  mpz_mul(&tmp2, &kummer[0], &kummer[2]);
  mpz_mul_ui(&tmp2, &tmp2, 4);
  mpz_sub(&tmp, &tmp, &tmp2); /* the discriminant */
#ifdef VERBOSE
  printf("  Kummer eqn: ");
  print_poly(kummer, 2);
  printf("\n  Discriminant = %s\n", mpz_get_str((char *) 0, 10, &tmp));
#endif
  if(mpz_cmp_si(&tmp, 0) == 0)
  { /* one point -> -k[1]/(2*k[2])*/
#ifdef VERBOSE
    printf("  Disc = 0 ==> one point\n");
#endif
    if(mpz_cmp_si(&kummer[2], 0) == 0 && mpz_cmp_si(&kummer[1], 0) == 0)
    {
#ifdef VERBOSE
      printf("  Kummer eqn. is constant ==> no point.\n");
#endif
      return(0);
    }
    mpz_mul_ui(&ddd, &kummer[2], 2);
    mpz_neg(&tmp, &kummer[1]);
    return(check_one_point_final(a, b, c, &tmp, &ddd));
  }
  else if(mpz_perfect_square_p(&tmp))
  { /* two points */
#ifdef VERBOSE
    printf("  Disc is a square ==> two points\n");
#endif
    if(mpz_cmp_si(&kummer[2], 0) == 0)
    { /* one point, after all */
#ifdef VERBOSE
      printf("  Leading coeff = 0 ==> one (finite) point\n");
#endif
      mpz_neg(&tmp, &kummer[0]);
      return(check_one_point_final(a, b, c, &tmp, &kummer[1]));
    }
    else /* kummer[2] != 0 */
    { /* two points, really */
      mpz_sqrt(&ddd, &tmp);
      mpz_sub(&tmp, &ddd, &kummer[1]);
      mpz_mul_ui(&kummer[2], &kummer[2], 2);
#ifdef VERBOSE
      printf("  sqrt(disc): %s\n", mpz_get_str((char *) 0, 10, &ddd));
      printf("  %s/%s\n", mpz_get_str((char *) 0, 10, &tmp),
                          mpz_get_str((char *) 0, 10, &kummer[2]));
#endif
      if(check_one_point_final(a, b, c, &tmp, &kummer[2]) && one_point)
        return(1);
      mpz_add(&tmp, &ddd, &kummer[1]);
      mpz_neg(&tmp, &tmp);
#ifdef VERBOSE
      printf("  sqrt(disc): %s\n", mpz_get_str((char *) 0, 10, &ddd));
      printf("  %s/%s\n", mpz_get_str((char *) 0, 10, &tmp),
                          mpz_get_str((char *) 0, 10, &kummer[2]));
#endif
      return(check_one_point_final(a, b, c, &tmp, &kummer[2]));
    } /* if kummer[2] == 0 */
  } /* if disc is a square */
  else
  {
#ifdef VERBOSE
    printf("  Disc is a non-square ==> no points\n\n");
#endif
    /* no point */
    return(0);
  }
}


/**************************************************************************
 * output routines                                                        *
 **************************************************************************/

void print_poly(MP_INT *coeffs, long degree)
{
  int flag = 0;
  int i;
  for(i = degree; i >= 0; i--)
  { mpz_set(&tmp3, &coeffs[i]);
    if(mpz_cmp_si(&tmp3, 0) != 0)
    { if(mpz_cmp_si(&tmp3, 0) > 0)
      { printf(flag ? " + " : ""); }
      else
      { printf(flag ? " - " : "- ");
        mpz_neg(&tmp3, &tmp3);
      }
      flag = 1;
      switch(i)
      { case 0: printf("%s", mpz_get_str((char *) 0, 10, &tmp3)); break;
        case 1: if(mpz_cmp_si(&tmp3, 1) == 0)
                  printf("x");
                else
                  printf("%s x", mpz_get_str((char *) 0, 10, &tmp3));
                break;
        default: if(mpz_cmp_si(&tmp3, 1) == 0)
                   printf("x^%d", i);
                 else
                   printf("%s x^%d", mpz_get_str((char *) 0, 10, &tmp3), i);
                 break;
  } } }
  printf("\n");
  fflush(stdout);
}

void message(long n, long total)
{
  switch(n)
  { case 0: printf("\n%s\n", J_POINTS_VERSION); break;
    case 1: printf("\nprob = 0, hence no solutions.\n"); break;
    case 2: printf("\nFound %ld rational points on K lifting to J.\n", total);
            break;
    case 4: if(tube_mode)
            { printf("Sieving in the tube of the bound on the fourth coordinate:\n");
              printf("about %.2f words per row, %.2f cells of the analysis per row.\n",
                     tube_words, tube_cells);
            }
            printf("%ld primes used for the first stage of sieving,\n",
                   sieve_primes1);
            printf("%ld primes used for the first two stages together,\n",
                   sieve_primes2);
            printf("%ld primes used for all three stages together.\n",
                   sieve_primes3);
            break;
    case 5: printf("\ny^2 = "); print_poly(coeffs, total); printf("\n"); break;
    case 6: printf("max. Height = %ld\n", total); break;
    case 7: { long i;
              printf("Sieving primes:\n First stage: ");
              for(i = 0; i < sieve_primes1; i++)
              { printf("%ld", prime[pnn[i]]);
                if(i < sieve_primes1 - 1) printf(", ");
              }
              printf("\n Second stage: ");
              for( ; i < sieve_primes2; i++)
              { printf("%ld", prime[pnn[i]]);
                if(i < sieve_primes2 - 1) printf(", ");
              }
              printf("\n Third stage: ");
              for( ; i < sieve_primes3; i++)
              { printf("%ld", prime[pnn[i]]);
                if(i < sieve_primes3 - 1) printf(", ");
              }
              printf("\n");
              break;
            }
    case 8:
      { /* the survival rates: the best prime, the last of each stage, the
           worst prime considered */
        long i;
        printf("Probabilities: Min(%ld) = %f", prec[0].p, prec[0].r);
        for(i = 0; i < num_primes; i++)
        { if(i + 1 == sieve_primes1 || i + 1 == sieve_primes2
             || i + 1 == sieve_primes3)
          { long j;
            for(j = 0; prec[j].n != pnn[i]; j++) ;
            printf(", Cut%d(%ld) = %f",
                   (i + 1 == sieve_primes1) ? 1 :
                   (i + 1 == sieve_primes2) ? 2 : 3,
                   prec[j].p, prec[j].r);
          }
        }
        printf(", Max(%ld) = %f\n\n",
               prec[num_primes-1].p, prec[num_primes-1].r);
        break;
      }
    case 12: printf("\n%ld candidates survived the first stage,\n", num_surv1);
             printf("%ld candidates survived the second stage,\n", num_surv2);
             printf("%ld candidates survived the third stage.\n", num_surv3);
             break;
  }
  fflush(stdout);
}

void error(long errno)
{
  switch(errno)
  { case 1:
      printf("\nUnusual size of `unsigned long' type: %d\n\n",
             (int)LONG_LENGTH);
      break;
    case 3: printf("\nToo many coefficients.\n\n"); break;
    case 4: printf("\nIncorrect height argument.\n");
            printf("  Height must be in [1, %ld].\n\n", MAX_HEIGHT); break;
    case 5: printf("\nThe polynomial must have degree at least 5.\n\n"); break;
    case 7: printf("\nNot enough memory.\n\n"); break;
    case 6: printf("\nWrong syntax for optional arguments:\n\n");
    case 2:
      printf("\n");
      printf("Usage: j-points 'a_0 a_1 ... a_d' max_height\n");
      printf("                [-n num_primes1] [-M num_primes2] [-N num_primes3]\n");
      printf("                [-p num_primes] [-s size] [-f format] [-1] [-q] [-a]\n");
      break;
  }
  fflush(stdout);
  exit(errno);
}
