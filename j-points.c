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

#define _GNU_SOURCE   /* asprintf */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <pthread.h>
#include <gmp.h>

#include "j-points.h"

/**************************************************************************
 * define                                                                 *
 **************************************************************************/

#define DEFAULT_SIZE 10     /* Default value for the -s option */

#define J_POINTS_VERSION \
  "This is j-points-3.0 by Michael Stoll (2026-09-26).\n\n" \
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

/* The cost model that chooses the sieving primes, their stages and the
   mode of the sieve (see choose_primes): the costs in cycles of one core,
   measured one operation at a time on the curve of test1 at 2000 with
   the sieve of 3.0, and four of them -- the pass, its size penalty, the
   second-stage test, the tube's test -- then adjusted end to end by
   "make tune" on a laptop and a desktop of 2023 for each width VW of the
   passes: a pass costs less per word the wider its vectors, and the
   balance between the stages shifts with it (the change log of 3.0).
   Only their ratios matter for the choice.  These are the defaults:
   "make tune" measures the machine's own and writes them to tuning.mk,
   whose -DCOST_X=value replaces a default here, and -c X=value replaces
   it for one run. */
#if VW == 4
#define COST_AND_VW     0.45
#define COST_SIZE_VW  450.0
#define COST_TEST2_VW   4.0
#else                        /* the scalar build, and SSE2 (two words per
                                AND buy less than the sieve's other work
                                gives back: measured the same set) */
#define COST_AND_VW     0.63
#define COST_SIZE_VW  756.0
#define COST_TEST2_VW   5.6
#endif
#ifndef COST_AND
#define COST_AND   COST_AND_VW   /* one word of one first-stage pass ... */
#endif
#ifndef COST_SIZE
#define COST_SIZE  COST_SIZE_VW  /* ... times 1 + p/COST_SIZE: the tables
                                    of the larger primes spill out of the
                                    cache */
#endif
#ifndef COST_WORD
#define COST_WORD    56.0  /* a word surviving the first stage: found and
                              handed to the second stage */
#endif
#ifndef COST_TEST2
#define COST_TEST2 COST_TEST2_VW /* one second-stage test of a word */
#endif
#ifndef COST_BIT
#define COST_BIT     30.0  /* a bit surviving the second stage: found,
                              its coprimality tested */
#endif
#ifndef COST_ROW3
#define COST_ROW3    40.0  /* per third-stage prime, the data of a row
                              that has such a bit (see row3_setup) */
#endif
#ifndef COST_TEST3
#define COST_TEST3   15.0  /* one third-stage test of a bit */
#endif
#ifndef COST_EXACT
#define COST_EXACT 1300.0  /* the exact check of a triple */
#endif
#ifndef COST_ROW
#define COST_ROW    300.0  /* a row of the box sieve besides its passes:
                              the fill, the scan, the set-up */
#endif
#ifndef COST_CELL
#define COST_CELL   320.0  /* a cell of the tube's analysis */
#endif
#ifndef COST_TWORD
#define COST_TWORD   12.0  /* a word the tube leaves, besides its tests */
#endif
#ifndef COST_TTEST
#define COST_TTEST    2.0  /* one test of such a word (no early exit) */
#endif
#ifndef COST_TROW
#define COST_TROW    45.0  /* a row with words in the tube: the set-up */
#endif
#ifndef COST_RROW
#define COST_RROW    40.0  /* a row's real region (row_region) ... */
#endif
#ifndef COST_REXCL
#define COST_REXCL   15.0  /* ... plus this per interval of c it excludes
                              (two per interval where f is negative) */
#endif
#ifndef COST_RCUT
#define COST_RCUT    30.0  /* the cut of a row's tube words to it (row_cut) */
#endif
#ifndef COST_RANGE
#define COST_RANGE   20.0  /* a range of words of the region, besides its passes */
#endif
#ifndef COST_PASS
#define COST_PASS     5.0  /* one pass over such a range, besides its words */
#endif
#ifndef COST_TABLE
#define COST_TABLE   13.0  /* one word of a sieve table at set-up */
#endif
#ifndef COST_INIT
#define COST_INIT   110.0  /* one class (1, b, c) of the exact table of a
                              prime whose rate was sampled, at set-up */
#endif
#define SAMPLE_PAIRS 4096  /* the classes (1, b, c) sampled for the rate
                              of a prime, see init_fmodpsquare */

/* The constants in force: the defaults, or what -c set (set = 1 then) */
static double cost_and = COST_AND, cost_size = COST_SIZE, cost_word = COST_WORD,
  cost_test2 = COST_TEST2, cost_bit = COST_BIT, cost_row3 = COST_ROW3,
  cost_test3 = COST_TEST3, cost_exact = COST_EXACT, cost_row = COST_ROW,
  cost_cell = COST_CELL, cost_tword = COST_TWORD, cost_ttest = COST_TTEST,
  cost_trow = COST_TROW, cost_rrow = COST_RROW, cost_rexcl = COST_REXCL,
  cost_rcut = COST_RCUT, cost_range = COST_RANGE, cost_pass = COST_PASS,
  cost_table = COST_TABLE, cost_init = COST_INIT;
static struct { const char *name; double *value; int set; } cost_names[] =
{ {"AND", &cost_and, 0}, {"SIZE", &cost_size, 0}, {"WORD", &cost_word, 0},
  {"TEST2", &cost_test2, 0}, {"BIT", &cost_bit, 0}, {"ROW3", &cost_row3, 0},
  {"TEST3", &cost_test3, 0}, {"EXACT", &cost_exact, 0}, {"ROW", &cost_row, 0},
  {"CELL", &cost_cell, 0}, {"TWORD", &cost_tword, 0}, {"TTEST", &cost_ttest, 0},
  {"TROW", &cost_trow, 0}, {"RROW", &cost_rrow, 0}, {"REXCL", &cost_rexcl, 0},
  {"RCUT", &cost_rcut, 0}, {"RANGE", &cost_range, 0}, {"PASS", &cost_pass, 0},
  {"TABLE", &cost_table, 0}, {"INIT", &cost_init, 0} };
#define NUM_COSTS (sizeof(cost_names)/sizeof(cost_names[0]))
static int costs_set = 0;   /* whether -c set any of them */

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
   by init_sieve: p*p rows of sieve_rowlen(p) words, the row of (a1, b1)
   starting at index (a1*p + b1)*sieve_rowlen(p).  Bit k of word c1 of
   that row is 0 iff (a,b,c) is excluded mod p when a = a1, b = b1 and
   c = c1*LONG_LENGTH + k mod p.  The row is periodic in c1 with period p,
   and continues the pattern beyond the period the walk uses. */

MP_INT coeffs[7];  /* The coefficients of f */
/* the multi-precision temporaries of the lifting test, one set per
   thread (init_thread_mpz) */
_Thread_local MP_INT bc[7]; /* A helper array */
_Thread_local MP_INT fff, tmp, tmp2, tmp3, ddd;   /* Some multi-precision integer variables */
MP_INT k400, k310, k301, k220, k211, k202, k130, k121, k112, k103,
       k040, k031, k022, k013, k004;
        /* coefficients for Kummer equation */
_Thread_local MP_INT x12, x22, x32;
_Thread_local MP_INT kummer[3];

/* The analysis of the plane of (b, c) (see j-sift.c): the coefficients
   of the Kummer equation and of f as doubles, for the bounds per cell,
   and how the rows are sieved: the box (every word of every row), the
   tube (the words the analysis leaves, one by one: a plain run's thin
   tube) or the region (the passes of the box over the ranges of words
   the analysis leaves: the real region of -a); the unsieved run must
   not use the analysis */
double kd400, kd310, kd301, kd220, kd211, kd202, kd130, kd121, kd112,
       kd103, kd040, kd031, kd022, kd013, kd004;
double fd[7];
int sieve_mode = 0;
int tube_cut = 0;          /* the tube's words cut to the real region */
/* sampled: the words per row the tube leaves, as they are and cut to
   the real region, the cells it visits per row and the fraction of the
   rows with a word; the words per row of the real region, in how many
   ranges, and the fraction of the rows that meet it */
double tube_words = 0.0, tube_words_cut = 0.0, tube_cells = 0.0, tube_rows = 1.0;
double reg_words = 0.0, reg_ranges = 0.0, reg_rows = 1.0;
/* The real roots of f, for condition (3) of the analysis: the open
   intervals on which f is negative, each given grown by the isolating
   intervals of the roots at its ends (f may be negative there) and
   shrunk by them (f is negative there for sure); see real_roots_init */
long num_neg = 0;
double neg_grown_lo[8], neg_grown_hi[8], neg_shrunk_lo[8], neg_shrunk_hi[8];

/* The condition at 2 (see twoadic_init): per class of (a, b) mod 64 the
   word of the admitted c mod 64, the fill of the bit array of a row, and
   whether the class admits any c at all; the density of the admitted
   triples and the fraction of empty classes for the cost model */
bit_array mask2[64][64];
unsigned char alive2[64][64];
int use2 = 0;
double density2 = 1.0, empty2 = 0.0;
_Thread_local MP_INT cpf1, cpf2, cpf3, cpfa, cpfb, cpfc;

long degree;
/* The curve is searched in a form that makes the search cheaper when
   there is one, see normalise_curve: reversed (x -> 1/x) and negated
   (x -> -x); the points are printed for the curve given. */
int reversed = 0, negated = 0;
long coeffs_mod_p[NUM_PRIMES][8];
                         /* The coefficients of f reduced modulo the various
                            primes */
entry prec[NUM_PRIMES];  /* This array is used for sorting in order to
                            determine the `best' sieving primes. */

long height;          /* The height bound */
long dbound = -1;     /* The bound on the fourth coordinate: the height
                         bound, or the value of -w */
int dbounded = 1;     /* whether the fourth coordinate is bounded at all
                         (not with -a, unless -w gave a bound) */
long sieve_primes1;   /* The number of primes used for the first sieving stage */
long sieve_primes2;   /* The number of primes used for the first two stages */
long sieve_primes3;   /* The number of primes used for all three stages */
int quiet;            /* A flag saying whether to suppress messages */
int one_point;        /* A flag saying if one point is enough */
int all_points;       /* Indicates that the `-a' option was given */
char *print_format;   /* The printf format for printing points */
long array_size;      /* The size of the survivors array (in longs) */

_Thread_local bit_array *survivors; /* In this array the sieving takes place */

/* the counters, per thread, and their sums over the threads */
_Thread_local long num_surv1 = 0;   /* Used to count the survivors of the first stage */
_Thread_local long num_surv2 = 0;   /* Used to count the survivors of the second stage */
_Thread_local long num_surv3 = 0;   /* Used to count the survivors of the third stage */
_Thread_local long total = 0;       /* Counts the points found */
long tot_surv1 = 0, tot_surv2 = 0, tot_surv3 = 0, tot_points = 0;

/* The threads (-t): the values of a are the units of work, handed out in
   order from a counter; each thread sieves its a with its own state and
   collects the points of that a in a buffer, and the buffers are printed
   in the order of the a, so that the output is that of one thread. */
long num_threads = 1;
static pthread_mutex_t out_lock = PTHREAD_MUTEX_INITIALIZER;
static long next_work = 0;           /* the next unit to hand out */
static long next_print = 0;          /* the next unit to print */
typedef struct outnode { long unit; char *text; struct outnode *next; } outnode;
static outnode *out_list = NULL;     /* finished units not yet printed, by unit */
static _Thread_local char *obuf = NULL;       /* the buffer of the unit in work */
static _Thread_local size_t olen = 0, ocap = 0;

bit_array begmask, endmask;  /* Bit masks for the beginning and end of
                                the sieving array */

/**************************************************************************
 * prototypes                                                             *
 **************************************************************************/

void init_main(void);
void init_thread_mpz(void);
void normalise_curve(void);
void find_points(void);
void read_input(long, char *argv[]);
char *scan_mpz(char*, MP_INT*);
static int squarefree(void);
static char *unescape(const char *);
void init_inverses(void);
void init_squares(void);
void init_fmodpsquare(void);
void choose_primes(void);
void init_sieve(void);
void kummer_init(void);
void tube_init(void);
static void real_roots_init(void);
void twoadic_init(void);
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
  /* "j-points -c list" prints the constants of the cost model in force
     in this build, which tune.sh measures against */
  if(argc == 3 && strcmp(argv[1], "-c") == 0 && strcmp(argv[2], "list") == 0)
  { size_t k;
    for(k = 0; k < NUM_COSTS; k++)
    { printf("%s %g\n", cost_names[k].name, *cost_names[k].value); }
    return(0);
  }
  /* read input */
  if(argc < 3) error(2);
  read_input(argc-1, &argv[0]);
  if(!quiet)
  { message(0, 0);
    message(5, degree);
    message(6, height);
    message(3, 0);
  }
  normalise_curve();
  begmask = (~0UL)<<((-height) & LONG_MASK);
  endmask = (~0UL)>>((~height) & LONG_MASK);
  s = 2*CEIL(height+1, LONG_LENGTH);
  array_size <<= 13 - LONG_SHIFT; /* from kbytes to longs */
  if(s < array_size) array_size = s;
  /* initialise data for equations */
  kummer_init();
  tube_init();
  /* the condition at 2, when the run is long enough to pay for its table */
  if(height >= 500) { twoadic_init(); }
  /* find and count points */
  find_points();
  if(!quiet) { message(12, 0); message(2, 0); }
  return(0);
}

/**************************************************************************
 * procedures                                                             *
 **************************************************************************/

/**************************************************************************
 * get at the input                                                       *
 **************************************************************************/

/* The argument of -c: NAME=value, several separated by commas, the names
   those of the COST_ constants without the prefix, the values positive;
   0 if it is not of that form */
static int set_costs(char *s)
{
  while(*s)
  { char *eq = strchr(s, '='), *end;
    size_t k;
    double v;
    if(eq == NULL) { return(0); }
    for(k = 0; k < NUM_COSTS; k++)
    { if(strlen(cost_names[k].name) == (size_t)(eq - s)
         && strncmp(cost_names[k].name, s, eq - s) == 0) { break; } }
    if(k == NUM_COSTS) { return(0); }
    v = strtod(eq + 1, &end);
    if(end == eq + 1 || !(v > 0.0) || (*end != '\0' && *end != ','))
    { return(0); }
    *cost_names[k].value = v;
    cost_names[k].set = 1;
    costs_set = 1;
    s = (*end == ',') ? end + 1 : end;
  }
  return(1);
}

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
  if(!squarefree()) error(9);
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
        case 'f': /* printing format, its backslash escapes interpreted */
          if(argc == i) error(6);
          i++;
          print_format = unescape(argv[i]);
          i++;
          break;
        case 'w': /* the bound on the fourth coordinate */
          if(argc == i) error(6);
          i++;
          if(sscanf(argv[i], " %ld", &dbound) != 1) error(6);
          if(dbound < 1) error(6);
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
        case 't': /* the number of threads */
          if(argc == i) error(6);
          i++;
          if(sscanf(argv[i], " %ld", &num_threads) != 1) error(6);
          if(num_threads < 1) error(6);
          i++;
          break;
        case 'c': /* the constants of the cost model: NAME=value, several
                     separated by commas */
          if(argc == i) error(6);
          i++;
          if(!set_costs(argv[i])) error(10);
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
  /* one point is found by one thread: the first in the order of the
     search */
  if(one_point) { num_threads = 1; }
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
  /* the bound on the fourth coordinate: -w's, else the height bound, and
     none with -a unless -w gave one */
  dbounded = (dbound > 0) || !all_points;
  if(dbound <= 0) { dbound = height; }
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

/* The search is cheapest for a monic quintic, whose points have a square
   first coordinate: the units of work are the squares up to h then,
   sqrt(h) instead of h.  A sextic with f0 = 0 and f1 = 1 becomes one
   under x -> 1/x, y -> y/x^3, which takes y^2 = f(x) to y^2 = x^6 f(1/x)
   (the coefficients in the opposite order), and a quintic with leading
   coefficient -1 under x -> -x (the odd coefficients negated); a sextic
   with f0 = 0 and f1 = -1 needs both.  The Kummer coordinates of a point
   go to (c : b : a : d) under the first map (the fourth coordinate
   x4 = (F0(x, u) - 2 y v) / (x - u)^2 in the Cassels-Flynn normalisation
   is unchanged, since F0 reverses with f) and to (a : -b : c : d) under
   the second, so both maps preserve the height, and the points of the
   curve given are those of the curve searched, transformed back when
   they are printed.  The messages name the curve searched. */
void normalise_curve(void)
{
  long i;
  if(degree == 6 && mpz_sgn(&coeffs[0]) == 0 && mpz_cmpabs_ui(&coeffs[1], 1) == 0)
  { for(i = 0; i < 3; i++) { mpz_swap(&coeffs[i], &coeffs[6-i]); }
    degree = 5;
    reversed = 1;
  }
  if(degree == 5 && mpz_cmp_si(&coeffs[5], -1) == 0)
  { for(i = 1; i <= 5; i += 2) { mpz_neg(&coeffs[i], &coeffs[i]); }
    negated = 1;
  }
  if((reversed || negated) && !quiet)
  { printf("The search runs on y^2 = "); print_poly(coeffs, degree);
    printf("(the curve under %s%s%s), a monic quintic, and the points are\n"
           "transformed back.\n\n",
           reversed ? "x -> 1/x" : "", (reversed && negated) ? " and " : "",
           negated ? "x -> -x" : "");
  }
  return;
}

/* Is f squarefree?  It must be, for a curve of genus 2: the gcd of f and
   f' over Q, by Euclid's algorithm on polynomials with rational
   coefficients (of degree at most 6), must be a constant. */
static int squarefree(void)
{
  mpq_t u[7], v[7], q, t, *pu = u, *pv = v, *ps;
  long du = degree, dv = degree - 1, i;
  int result;
  for(i = 0; i <= 6; i++) { mpq_init(u[i]); mpq_init(v[i]); }
  mpq_init(q); mpq_init(t);
  for(i = 0; i <= degree; i++) { mpq_set_z(u[i], &coeffs[i]); }
  for(i = 1; i <= degree; i++)
  { mpq_set_z(v[i-1], &coeffs[i]); mpq_set_ui(t, (unsigned long)i, 1);
    mpq_mul(v[i-1], v[i-1], t);
  }
  /* u of degree du, v of degree dv <= du: replace u by its remainder
     modulo v, then swap, until v is 0; u is then the gcd */
  while(dv >= 0)
  { while(du >= dv)
    { mpq_div(q, pu[du], pv[dv]);
      for(i = 0; i <= dv; i++)
      { mpq_mul(t, q, pv[i]); mpq_sub(pu[i + du - dv], pu[i + du - dv], t); }
      /* the leading term of u is gone (exactly) */
      du--;
      while(du >= 0 && mpq_sgn(pu[du]) == 0) { du--; }
    }
    ps = pu; pu = pv; pv = ps;
    i = du; du = dv; dv = i;
  }
  result = (du == 0);
  for(i = 0; i <= 6; i++) { mpq_clear(u[i]); mpq_clear(v[i]); }
  mpq_clear(q); mpq_clear(t);
  return(result);
}

/* the format of -f with its backslash escapes \n, \t and \\ interpreted;
   any other backslash is kept */
static char *unescape(const char *s)
{
  char *r = malloc(strlen(s) + 1), *d = r;
  if(r == NULL) { error(7); }
  for( ; *s; s++)
  { if(*s == '\\' && (s[1] == 'n' || s[1] == 't' || s[1] == '\\'))
    { *d++ = (s[1] == 'n') ? '\n' : (s[1] == 't') ? '\t' : '\\'; s++; }
    else { *d++ = *s; }
  }
  *d = 0;
  return(r);
}

/**************************************************************************
 * initialisations                                                        *
 **************************************************************************/

/* the multi-precision temporaries of the calling thread */
void init_thread_mpz(void)
{
  long n;
  for(n = 0; n <= 6 ; n++) { mpz_init(&bc[n]); }
  for(n = 0; n < 3; n++) mpz_init(&kummer[n]);
  mpz_init(&fff);
  mpz_init(&tmp); mpz_init(&tmp2); mpz_init(&tmp3); mpz_init(&ddd);
  mpz_init(&x12); mpz_init(&x22); mpz_init(&x32);
  mpz_init(&cpf1); mpz_init(&cpf2); mpz_init(&cpf3);
  mpz_init(&cpfa); mpz_init(&cpfb); mpz_init(&cpfc);
  return;
}

void init_main(void)
{
  bit_array bit;
  long n;
  /* initialise multi-precision integer variables */
  for(n = 0; n <= 6 ; n++) { mpz_init(&coeffs[n]); }
  mpz_init(&k400); mpz_init(&k310); mpz_init(&k301); mpz_init(&k220);
  mpz_init(&k211); mpz_init(&k202); mpz_init(&k130); mpz_init(&k121);
  mpz_init(&k112); mpz_init(&k103); mpz_init(&k040); mpz_init(&k031);
  mpz_init(&k022); mpz_init(&k013); mpz_init(&k004);
  init_thread_mpz();

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
static double choose_primes_mode(int mode, int cut);

/* the ways of a run priced, the cheapest taken with its primes: the box,
   the region, and for a run with a bound on the fourth coordinate the
   tube, with or without the cut of its words to the real region (the
   analysis, the cost of a row and of a range are not among the terms
   that the choice of the primes sees, so they are added there) */
void choose_primes(void)
{
  long pin1 = sieve_primes1, pin2 = sieve_primes2, pin3 = sieve_primes3;
  long b1 = 0, b2 = 0, b3 = 0, pn[NUM_PRIMES], n, k;
  int modes[4] = {0, 2, 1, 1}, cuts[4] = {0, 0, 0, 1}, nmodes = dbounded ? 4 : 2;
  double best = -1.0;
  for(k = 0; k < nmodes; k++)
  { double cost;
    sieve_primes1 = pin1; sieve_primes2 = pin2; sieve_primes3 = pin3;
    cost = choose_primes_mode(modes[k], cuts[k]);
    if(best < 0.0 || cost < best)
    { best = cost; sieve_mode = modes[k]; tube_cut = cuts[k];
      b1 = sieve_primes1; b2 = sieve_primes2; b3 = sieve_primes3;
      for(n = 0; n < b3; n++) { pn[n] = pnn[n]; }
    }
  }
  sieve_primes1 = b1; sieve_primes2 = b2; sieve_primes3 = b3;
  for(n = 0; n < b3; n++) { pnn[n] = pn[n]; }
  return;
}

static double choose_primes_mode(int mode, int cut)
{
  int tube = (mode == 1);
  long ne = 0, n, n1, n2, n3, k;
  cand list1[NUM_PRIMES], list2[NUM_PRIMES], list3[NUM_PRIMES];
  long best1[NUM_PRIMES], best2[NUM_PRIMES], best3[NUM_PRIMES];
  double W = (double)(2*(height>>LONG_SHIFT) + 2);  /* words per row */
  double rows = (double)(2*height + 1)   /* the rows: only the square a of a monic quintic */
                * ((degree == 5 && mpz_cmp_si(&coeffs[5], 1) == 0)
                   ? floor(sqrt((double)height)) + 1.0 : (double)(height + 1));
  double bits;                                      /* bits per row */
  double best = -1.0, fixed, passx = 0.0;   /* per row; passx per pass */
  double rrow = cost_rrow + cost_rexcl * 2.0 * (double)num_neg;   /* a row's real region */
  long b1 = 0, b2 = 0, b3 = 0;
  long n1lo, n1hi, pin1 = sieve_primes1;
  double logr[NUM_PRIMES], tabcost[NUM_PRIMES], initcost[NUM_PRIMES];
  int used[NUM_PRIMES];

  /* in the tube there is no first stage: the words the tube leaves are
     tested one by one, like the survivors of a first stage but without
     an early exit; the cost of the analysis and of a row are fixed.  In
     the region the passes run over the words left, range by range, with
     the row's cost for the rows that meet it. */
  if(tube)
  { W = cut ? tube_words_cut : tube_words;
    if(W < 0.01) { W = 0.01; }
    pin1 = 0;
    fixed = tube_cells * cost_cell + tube_rows * (cost_trow + (cut ? rrow + cost_rcut : 0.0))
            + W * cost_tword;
  }
  else if(mode == 2)
  { W = (reg_words > 0.01) ? reg_words : 0.01;
    fixed = rrow + reg_rows * cost_row + reg_ranges * cost_range;
    passx = reg_ranges * cost_pass;
  }
  else
  { fixed = cost_row; }
  bits = (double)LONG_LENGTH * W;
  /* the condition at 2: the rows its empty classes remove, the density
     the survivors start from */
  if(use2) { rows *= 1.0 - empty2; }
  for(n = 0; n < num_primes; n++)
  { double p = (double)prec[n].p;
    logr[n] = (prec[n].r > 0.0) ? -log(prec[n].r) : 1.0e9;
    tabcost[n] = cost_table * p*p*(p + 1.0);
    initcost[n] = (p*p > (double)SAMPLE_PAIRS) ? cost_init * p*p : 0.0;
  }
  /* the first stage's ranking of the table primes */
  for(n = 0; n < num_primes; n++)
  { if(prec[n].p <= MAX_TABLE_PRIME)
    { list1[ne].n = n;
      list1[ne].v = logr[n] / (W * cost_and * (1.0 + (double)prec[n].p / cost_size)
                               + passx + (tabcost[n] + initcost[n]) / rows);
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
  { double cost1 = 0.0, setup1 = 0.0, rho1 = use2 ? density2 : 1.0, s1;
    long n2lo, n2hi, ne2 = 0;
    for(n = 0; n < num_primes; n++) { used[n] = 0; }
    for(k = 0; k < n1; k++)
    { n = list1[k].n;
      used[n] = 1;
      cost1 += W * cost_and * (1.0 + (double)prec[n].p / cost_size) + passx;
      setup1 += tabcost[n] + initcost[n];
      rho1 *= prec[n].r;
    }
    s1 = 1.0 - pow(1.0 - rho1, (double)LONG_LENGTH);  /* words surviving */
    /* the second stage's ranking of the table primes left */
    for(k = 0; k < ne; k++)
    { n = list1[k].n;
      if(!used[n])
      { list2[ne2].n = n;
        list2[ne2].v = logr[n] / (W * s1 * cost_test2
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
      cost2 = tube ? 0.0 : W * w * cost_word;
      for(k = 0; k < n2 - n1; k++)
      { n = list2[k].n;
        cost2 += tube ? W * cost_ttest : W * w * cost_test2;
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
          list3[ne3].v = logr[n] / (bits * rho2 * cost_test3
                                          + rowbit * cost_row3
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
      { double cost3 = bits * rho2 * cost_bit, setup3 = setup2, rho = rho2;
        for(n3 = n2; n3 <= n3hi; n3++)
        { if(n3 >= n3lo)
          { double total = rows * (fixed + cost1 + cost2 + cost3
                                   + rowbit * (double)(n3 - n2) * cost_row3
                                   + bits * rho * cost_exact)
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
            cost3 += bits * rho * cost_test3;
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
    sieve_tab[pn] = (bit_array *)malloc(p*p*sieve_rowlen(p)*sizeof(bit_array));
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
      { bit_array *si = &sieve_tab[pn][(a*p + b)*sieve_rowlen(p)];
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
        /* the continuation of the pattern beyond one period */
        for( ; c < sieve_rowlen(p); c++) si[c] = si[c - p];
      }

#if (DEBUG >= 3)
    printf(" sieve(%ld):\n", p);
    for(a = 0; a < p; a++)
    { printf(" a = %3ld:\n", a);
      for(b = 0; b < p; b++)
      { printf("  b = %3ld: ", b);
        for(c = 0; c < p; c++)
	  printf(" %8.8lx", sieve_tab[pn][(a*p + b)*sieve_rowlen(p) + c]);
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

static long unit_a(long);
static void run_units(long);
static void emit(char *);

void find_points(void)
{

  /* initialise is_f_square[][] */
  init_fmodpsquare();
  /* allocate and initalise the sieve tables */
  init_sieve();
  if(sieve_primes3 > 0 && prec[0].r == 0.0)
  { if(!quiet) message(1,0); return; }
  if(sieve_primes3 == 0) { sieve_mode = 0; use2 = 0; }
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
  /* the main search, the units of work being the values of a */
  { long units = (degree == 5 && mpz_cmp_si(&coeffs[5], 1) == 0)
                 ? (long)floor(sqrt((double)height)) + 1 : height + 1;
    while(units > 0 && unit_a(units - 1) > height) { units--; }
    run_units(units);
  }
  return;
}

/* the first coordinate of the unit of work k: k itself, or k^2 when f is
   monic of degree 5 */
static long unit_a(long k)
{ return((degree == 5 && mpz_cmp_si(&coeffs[5], 1) == 0) ? k*k : k); }

/* the rows of one a: those the primes without a point at infinity leave
   (row_step) and the condition at 2 leaves, in the tube or the box;
   returns 1 when one point is enough and one was found */
static int run_a(long a)
{
  long m = row_step(a), b0 = (a == 0) ? m : -(height/m)*m, b;
  if(sieve_mode) { return(sift_bands(a, b0, m)); }
  for(b = b0; b <= height; b += m)
  {
#ifdef VERBOSE
    printf(" a = %ld, b = %ld\n", a, b);
#endif
    if(use2 && a != 0 && !alive2[a & 63][b & 63]) { continue; }
    if(sift(a, b) && one_point) { return(1); }
  }
  return(0);
}

/* a unit of work done by a thread: its points, if any, handed to the
   printer, which prints the finished units in order */
static void finish_unit(long unit)
{
  outnode *node = (outnode *)malloc(sizeof(outnode)), **pp;
  if(node == NULL) { error(7); }
  node->unit = unit;
  node->text = obuf; obuf = NULL; olen = ocap = 0;
  pthread_mutex_lock(&out_lock);
  for(pp = &out_list; *pp != NULL && (*pp)->unit < unit; pp = &(*pp)->next) ;
  node->next = *pp; *pp = node;
  while(out_list != NULL && out_list->unit == next_print)
  { outnode *done = out_list;
    if(done->text != NULL) { fputs(done->text, stdout); free(done->text); }
    out_list = done->next;
    free(done);
    next_print++;
  }
  pthread_mutex_unlock(&out_lock);
  return;
}

/* a thread: units from the counter until they run out */
static long work_units;
static void *worker(void *arg)
{
  long k;
  if(arg != NULL)   /* a thread of its own, not the main one */
  { init_thread_mpz(); init_thread_sieve(); }
  while((k = __atomic_fetch_add(&next_work, 1, __ATOMIC_RELAXED)) < work_units)
  { run_a(unit_a(k));
    finish_unit(k);
  }
  pthread_mutex_lock(&out_lock);
  tot_surv1 += num_surv1; tot_surv2 += num_surv2; tot_surv3 += num_surv3;
  tot_points += total;
  pthread_mutex_unlock(&out_lock);
  return(NULL);
}

static void run_units(long units)
{
  work_units = units;
  if(num_threads == 1)
  { long k;
    for(k = 0; k < units; k++)
    { if(run_a(unit_a(k)) && one_point) { break; } }
    tot_surv1 = num_surv1; tot_surv2 = num_surv2; tot_surv3 = num_surv3;
    tot_points = total;
    return;
  }
  { pthread_t *th = (pthread_t *)malloc(num_threads*sizeof(pthread_t));
    long t;
    if(th == NULL) { error(7); }
    for(t = 1; t < num_threads; t++)
    { if(pthread_create(&th[t], NULL, worker, (void *)th) != 0) { error(8); } }
    worker(NULL);
    for(t = 1; t < num_threads; t++) { pthread_join(th[t], NULL); }
    free(th);
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

/* the coefficients as doubles for the analysis of the plane of (b, c);
   what it leaves per row is sampled over a few hundred rows spread over
   the box (deterministically), for the cost model, which then chooses
   the way to sieve */
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
  sieve_mode = 0;
  real_roots_init();
  tube_sample(&tube_words, &tube_words_cut, &tube_cells, &tube_rows, &reg_words, &reg_ranges, &reg_rows);
  return;
}

/* The real roots of f by Sturm's theorem in exact arithmetic: with the
   chain p0 = f, p1 = f', p_i = -(p_{i-2} mod p_{i-1}) over Q, the
   number of roots in (x, y] is the number of sign changes of the chain
   at x less that at y (zeros skipped; the roots are simple, f being
   squarefree).  The roots are isolated by bisection of [-B, B], B the
   Cauchy bound, and each isolating interval is bisected on to a relative
   width of 2^-40; between consecutive roots f has one sign, alternating
   from the sign at +infinity.  The negative intervals go to the arrays
   above, their ends rounded outwards (grown) or inwards (shrunk). */
static mpq_t sturm[7][7];      /* the chain, sturm[i][k] the coefficient of x^k of p_i */
static long sturm_deg[7], sturm_len;
static mpq_t root_lo[6], root_hi[6];   /* the isolating intervals, in order */
static long num_roots;

static long sturm_changes(mpq_t x)
{
  mpq_t v, t;
  long i, k, changes = 0;
  int last = 0;
  mpq_init(v); mpq_init(t);
  for(i = 0; i < sturm_len; i++)
  { int s;
    mpq_set(v, sturm[i][sturm_deg[i]]);
    for(k = sturm_deg[i] - 1; k >= 0; k--)
    { mpq_mul(v, v, x); mpq_add(v, v, sturm[i][k]); }
    s = mpq_sgn(v);
    if(s != 0)
    { if(last != 0 && s != last) { changes++; }
      last = s;
    }
  }
  mpq_clear(v); mpq_clear(t);
  return(changes);
}

/* the roots in (lo, hi], nlo - nhi of them, isolated and refined */
static void sturm_isolate(mpq_t lo, mpq_t hi, long nlo, long nhi)
{
  mpq_t mid;
  long nmid;
  if(nlo == nhi) { return; }
  mpq_init(mid);
  if(nlo - nhi == 1)
  { /* one root: bisect to a relative width of 2^-40 */
    mpq_t w;
    mpq_init(w);
    for(;;)
    { double d = mpq_get_d(hi) - mpq_get_d(lo), s = fabs(mpq_get_d(lo)) + fabs(mpq_get_d(hi));
      if(d <= 0.5*ldexp(1.0, -40)*(s + 1.0)) { break; }
      mpq_add(mid, lo, hi); mpq_div_2exp(mid, mid, 1);
      nmid = sturm_changes(mid);
      if(nlo - nmid == 1) { mpq_set(hi, mid); nhi = nmid; }
      else { mpq_set(lo, mid); nlo = nmid; }
    }
    mpq_clear(w);
    if(num_roots < 6)
    { mpq_set(root_lo[num_roots], lo); mpq_set(root_hi[num_roots], hi); num_roots++; }
    mpq_clear(mid);
    return;
  }
  mpq_add(mid, lo, hi); mpq_div_2exp(mid, mid, 1);
  nmid = sturm_changes(mid);
  { mpq_t l2, h2;
    mpq_init(l2); mpq_init(h2);
    mpq_set(l2, lo); mpq_set(h2, mid);
    sturm_isolate(l2, h2, nlo, nmid);
    mpq_set(l2, mid); mpq_set(h2, hi);
    sturm_isolate(l2, h2, nmid, nhi);
    mpq_clear(l2); mpq_clear(h2);
  }
  mpq_clear(mid);
  return;
}

static void real_roots_init(void)
{
  mpq_t q, lo, hi;
  long i, k, sgn_inf, nlo, nhi;
  double B = 0.0;
  for(i = 0; i < 7; i++) { for(k = 0; k < 7; k++) { mpq_init(sturm[i][k]); } }
  for(i = 0; i < 6; i++) { mpq_init(root_lo[i]); mpq_init(root_hi[i]); }
  mpq_init(q); mpq_init(lo); mpq_init(hi);
  /* p0 = f, p1 = f' */
  for(k = 0; k <= degree; k++) { mpq_set_z(sturm[0][k], &coeffs[k]); }
  sturm_deg[0] = degree;
  for(k = 1; k <= degree; k++)
  { mpq_set_z(sturm[1][k-1], &coeffs[k]); mpq_set_ui(q, (unsigned long)k, 1);
    mpq_mul(sturm[1][k-1], sturm[1][k-1], q);
  }
  sturm_deg[1] = degree - 1;
  sturm_len = 2;
  /* p_i = -(p_{i-2} mod p_{i-1}) */
  while(sturm_deg[sturm_len-1] > 0)
  { long i0 = sturm_len - 2, i1 = sturm_len - 1, i2 = sturm_len, d;
    for(k = 0; k < 7; k++) { mpq_set(sturm[i2][k], sturm[i0][k]); }
    d = sturm_deg[i0];
    while(d >= sturm_deg[i1])
    { mpq_div(q, sturm[i2][d], sturm[i1][sturm_deg[i1]]);
      for(k = 0; k <= sturm_deg[i1]; k++)
      { mpq_t t;
        mpq_init(t);
        mpq_mul(t, q, sturm[i1][k]);
        mpq_sub(sturm[i2][k + d - sturm_deg[i1]], sturm[i2][k + d - sturm_deg[i1]], t);
        mpq_clear(t);
      }
      d--;
      while(d >= 0 && mpq_sgn(sturm[i2][d]) == 0) { d--; }
    }
    if(d < 0) { break; }   /* the remainder is 0 */
    for(k = 0; k <= d; k++) { mpq_neg(sturm[i2][k], sturm[i2][k]); }
    sturm_deg[i2] = d;
    sturm_len++;
  }
  /* the Cauchy bound: every root is below 1 + max |f_k / f_d| */
  for(k = 0; k < degree; k++)
  { double r = fabs(mpz_get_d(&coeffs[k]) / mpz_get_d(&coeffs[degree]));
    if(r > B) { B = r; }
  }
  B = ceil(B) + 2.0;
  mpq_set_d(lo, -B); mpq_set_d(hi, B);
  nlo = sturm_changes(lo); nhi = sturm_changes(hi);
  num_roots = 0;
  sturm_isolate(lo, hi, nlo, nhi);
  /* the sign of f beyond the last root, and on each interval down from
     there, alternating; a negative interval is recorded */
  sgn_inf = mpz_sgn(&coeffs[degree]);
  num_neg = 0;
  for(i = num_roots; i >= 0; i--)
  { long s = ((num_roots - i) % 2 == 0) ? sgn_inf : -sgn_inf;   /* on (root i-1, root i) */
    if(s < 0)
    { double glo, ghi, slo, shi;
      if(i == 0) { glo = slo = -HUGE_VAL; }
      else
      { glo = nextafter(mpq_get_d(root_lo[i-1]), -HUGE_VAL);
        shi = 0.0; slo = nextafter(mpq_get_d(root_hi[i-1]), HUGE_VAL);
      }
      if(i == num_roots) { ghi = shi = HUGE_VAL; }
      else
      { ghi = nextafter(mpq_get_d(root_hi[i]), HUGE_VAL);
        shi = nextafter(mpq_get_d(root_lo[i]), -HUGE_VAL);
      }
      neg_grown_lo[num_neg] = glo; neg_grown_hi[num_neg] = ghi;
      neg_shrunk_lo[num_neg] = slo; neg_shrunk_hi[num_neg] = shi;
      num_neg++;
    }
  }
  mpq_clear(q); mpq_clear(lo); mpq_clear(hi);
  for(i = 0; i < 7; i++) { for(k = 0; k < 7; k++) { mpq_clear(sturm[i][k]); } }
  for(i = 0; i < 6; i++) { mpq_clear(root_lo[i]); mpq_clear(root_hi[i]); }
  return;
}

/* The condition at 2.  A point (a' : b' : c' : d') of K with coprime
   integer coordinates that lifts to J satisfies the Kummer equation, and
   A^2 = a'^3 d' + A0(a', b', c') is a square or 0, with A0 = f2 a^4 + f3
   a^3 b + f4 a^2 b^2 + f5 a b (b^2 - a c) + f6 (b^2 - a c)^2 (the
   identity needs a' != 0).  The sieve's triple (a, b, c) is coprime, and
   the point's coordinates are (v a, v b, v c, u) for d = u/v.  So, mod
   64, (a, b, c) must admit either some d with k2 d^2 + k1 d + k0 = 0 and
   a^3 d + A0 a square residue (v odd), or, for v = 2^e, some odd t with
   k2 t^2 + 2^e k1 t + 4^e k0 = 0 and, when e = 1, 8 (a^3 t + 2 A0) a
   square residue (for e >= 2 the square condition is empty).  Modulo 64
   the admitted c form one word, which is the fill of the bit array for
   the row (a, b); rows whose class admits no c are skipped.  The table:
   sol[k2][k1][k0] holds the d with k2 d^2 + k1 d + k0 = 0 mod 64, built
   by running over (k2, k1, d); the cases v = 2^e use it with k1, k0
   scaled.  Rows with a = 0 keep the plain condition (not all even). */
/* the table for the modulus m (16 or 64; a word holds the 64 residues
   of c, the pattern repeated for m = 16); returns the number of admitted
   triples not all even, and counts the empty classes of (a, b) */
static long twoadic_table(long m, bit_array masks[64][64],
                          unsigned char alive[64][64], long *empty)
{
  bit_array *sol = (bit_array *)malloc(m*m*m*sizeof(bit_array));
  long fm[7], km[5][5][5], k2, k1, d, i, j, k, a2, b2, c2, mm = m - 1;
  long count = 0;
  int sq[64];
  bit_array odd = 0, all = (m == 64) ? ~(bit_array)0 : ((1UL << m) - 1);
  if(sol == NULL) { error(7); }
  for(d = 1; d < m; d += 2) { odd |= (1UL << d); }
  for(i = 0; i < m; i++) { sq[i] = 0; }
  for(i = 0; i < m; i++) { sq[(i*i) & mm] = 1; }
  for(i = 0; i < m*m*m; i++) { sol[i] = 0; }
  for(k2 = 0; k2 < m; k2++)
    for(k1 = 0; k1 < m; k1++)
      for(d = 0; d < m; d++)
      { long v = (k2*d*d + k1*d) & mm;
        sol[(k2*m + k1)*m + ((m - v) & mm)] |= (1UL << d);
      }
  for(i = 0; i <= 6; i++) { fm[i] = mpz_fdiv_ui(&coeffs[i], m); }
  for(i = 0; i < 5; i++) for(j = 0; j < 5; j++) for(k = 0; k < 5; k++) { km[i][j][k] = 0; }
  km[4][0][0] = mpz_fdiv_ui(&k400, m); km[3][1][0] = mpz_fdiv_ui(&k310, m);
  km[3][0][1] = mpz_fdiv_ui(&k301, m); km[2][2][0] = mpz_fdiv_ui(&k220, m);
  km[2][1][1] = mpz_fdiv_ui(&k211, m); km[2][0][2] = mpz_fdiv_ui(&k202, m);
  km[1][3][0] = mpz_fdiv_ui(&k130, m); km[1][2][1] = mpz_fdiv_ui(&k121, m);
  km[1][1][2] = mpz_fdiv_ui(&k112, m); km[1][0][3] = mpz_fdiv_ui(&k103, m);
  km[0][4][0] = mpz_fdiv_ui(&k040, m); km[0][3][1] = mpz_fdiv_ui(&k031, m);
  km[0][2][2] = mpz_fdiv_ui(&k022, m); km[0][1][3] = mpz_fdiv_ui(&k013, m);
  km[0][0][4] = mpz_fdiv_ui(&k004, m);
  *empty = 0;
  for(a2 = 0; a2 < m; a2++)
  { long ap[5], bp[5], cp[5];
    ap[0] = 1; for(i = 1; i < 5; i++) { ap[i] = (ap[i-1]*a2) & mm; }
    for(b2 = 0; b2 < m; b2++)
    { bit_array mask = 0;
      bp[0] = 1; for(i = 1; i < 5; i++) { bp[i] = (bp[i-1]*b2) & mm; }
      for(c2 = 0; c2 < m; c2++)
      { long k0 = 0, kk, A0, a3, e;
        int ok = 0;
        bit_array S;
        if(((a2 | b2 | c2) & 1) == 0) { continue; }
        cp[0] = 1; for(i = 1; i < 5; i++) { cp[i] = (cp[i-1]*c2) & mm; }
        for(i = 0; i < 5; i++)
          for(j = 0; i + j < 5; j++)
          { k = 4 - i - j;
            k0 = (k0 + km[i][j][k]*ap[i]*bp[j]*cp[k]) & mm;
          }
        k1 = (m - ((4*ap[3]*fm[0] + 2*ap[2]*b2*fm[1] + 4*ap[2]*c2*fm[2]
                    + 2*a2*b2*c2*fm[3] + 4*a2*cp[2]*fm[4] + 2*b2*cp[2]*fm[5]
                    + 4*cp[3]*fm[6]) & mm)) & mm;
        k2 = (bp[2] + m*4 - 4*a2*c2) & mm;
        kk = (bp[2] + m*4 - a2*c2) & mm;
        A0 = (fm[2]*ap[4] + fm[3]*ap[3]*b2 + fm[4]*ap[2]*bp[2] + fm[5]*a2*b2*kk
              + fm[6]*kk*kk) & mm;
        a3 = ap[3];
        /* v odd: some d */
        S = sol[(k2*m + k1)*m + k0];
        for(d = 0; d < m && !ok; d++)
        { if((S >> d) & 1) { ok = sq[(a3*d + A0) & mm]; } }
        /* v = 2: some odd t */
        if(!ok)
        { S = sol[(k2*m + ((2*k1) & mm))*m + ((4*k0) & mm)] & odd;
          for(d = 1; d < m && !ok; d += 2)
          { if((S >> d) & 1) { ok = sq[(8*(a3*d + 2*A0)) & mm]; } }
        }
        /* v = 4, ..., m: some odd t, no square condition */
        for(e = 2; e <= 6 && !ok; e++)
        { S = sol[(k2*m + ((k1 << e) & mm))*m + ((k0 << (2*e)) & mm)] & odd;
          ok = (S != 0);
        }
        if(ok) { mask |= (1UL << c2); count++; }
      }
      /* the pattern over the 64 residues of a word */
      if(m < 64) { for(i = m; i < 64; i += m) { mask |= mask << m; } }
      masks[a2][b2] = mask & ((m < 64) ? ~(bit_array)0 : all);
      alive[a2][b2] = (mask != 0);
      if(mask == 0) { (*empty)++; }
    }
  }
  free(sol);
  return(count);
}

void twoadic_init(void)
{
  static bit_array masks16[64][64];
  static unsigned char alive16[64][64];
  long count, empty;
  /* the modulus 16 first, cheap: when it admits every class (f a square
     modulo 4, as on the record curve), the modulus 64 does too */
  count = twoadic_table(16, masks16, alive16, &empty);
  if(count == 16*16*16 - 8*8*8) { use2 = 0; return; }
  count = twoadic_table(64, mask2, alive2, &empty);
  density2 = (double)count / (double)(64*64*64 - 32*32*32);
  empty2 = (double)empty / (double)(64*64);
  use2 = 1;
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
  /* the point of the curve given, its first nonzero coordinate positive */
  if(reversed) { MP_INT *x = a; a = c; c = x; }
  if(negated) { mpz_neg(b, b); }
  if((reversed || negated)
     && (mpz_sgn(a) < 0 || (mpz_sgn(a) == 0 && (mpz_sgn(b) < 0 || (mpz_sgn(b) == 0 && mpz_sgn(c) < 0)))))
  { mpz_neg(a, a); mpz_neg(b, b); mpz_neg(c, c); mpz_neg(d, d); }
  if(num_threads == 1) { gmp_printf(fmt, a, b, c, d); }
  else
  { char *s;
    if(gmp_asprintf(&s, fmt, a, b, c, d) < 0) { error(7); }
    emit(s);
    free(s);
  }
  free(fmt);
}

/* a point found: printed, or, with threads, appended to the buffer of
   the unit of work in hand (see find_points) */
static void emit(char *s)
{
  size_t l = strlen(s);
  if(num_threads == 1) { fputs(s, stdout); return; }
  if(olen + l + 1 > ocap)
  { ocap = 2*(olen + l + 1) + 256;
    obuf = (char *)realloc(obuf, ocap);
    if(obuf == NULL) { error(7); }
  }
  memcpy(obuf + olen, s, l + 1);
  olen += l;
  return;
}

static void emit_point(long a, long b, long c, long d)
{
  char *s;
  /* the point of the curve given, its first nonzero coordinate positive */
  if(reversed) { long x = a; a = c; c = x; }
  if(negated) { b = -b; }
  if(a < 0 || (a == 0 && (b < 0 || (b == 0 && c < 0))))
  { a = -a; b = -b; c = -c; d = -d; }
  if(num_threads == 1) { printf(print_format, a, b, c, d); return; }
  if(asprintf(&s, print_format, a, b, c, d) < 0) { error(7); }
  emit(s);
  free(s);
  return;
}

int check_one_point_final(long a, long b, long c, MP_INT *d1, MP_INT *d2)
{ /* Given a, b, c and d = d1/d2, check if this gives a point satisfying
     the height condition. If so, print and count it. */
  long m = (a > labs(b)) ? ((a > labs(c)) ? a : labs(c))
                         : ((labs(b) > labs(c)) ? labs(b) : labs(c));
  /* the coordinates of the point are (g a, g b, g c, d1) with g = d2 once
     d is in lowest terms: g m is bounded by the height bound (not with
     -a), d1 by the bound on the fourth coordinate (with -a only when -w
     gave one).  A point is printed as machine words when its coordinates
     are at most MAX_HEIGHT, so that the products of two of them fit a
     long in the lifting test; otherwise through gmp. */
  long h = (all_points ? MAX_HEIGHT : height)/m;
  long dmax = (dbounded && dbound < MAX_HEIGHT) ? dbound : MAX_HEIGHT;
  int dok;
  mpz_gcd(&cpf3, d1, d2);
  mpz_divexact(&cpf1, d1, &cpf3);
  mpz_divexact(&cpf2, d2, &cpf3);
  dok = !dbounded
        || (mpz_cmp_si(&cpf1, dbound) <= 0 && mpz_cmp_si(&cpf1, -dbound) >= 0);
  if(mpz_cmp_si(&cpf1, dmax) <= 0 && mpz_cmp_si(&cpf1, -dmax) >= 0
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
      emit_point(a, b, c, d);
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
  else if(dok && (all_points
                  || (mpz_cmp_si(&cpf2, h) <= 0 && mpz_cmp_si(&cpf2, -h) >= 0)))
  { /* the coordinates exceed a machine word (with -a, or with a bound on
       the fourth coordinate above MAX_HEIGHT) */
    /* scale by denominator &cpf2; fourth coordinate is numerator &cpf1;
       the sign of the scaling made positive, as for machine words */
#ifdef VERBOSE
    printf("  coordinates exceed machine size ");
#endif
    if(mpz_sgn(&cpf2) < 0) { mpz_neg(&cpf2, &cpf2); mpz_neg(&cpf1, &cpf1); }
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
    case 2: printf("\nFound %ld rational points on K lifting to J.\n", tot_points);
            break;
    case 4: if(sieve_mode == 1)
            { printf("Sieving in the tube of the bound on the fourth coordinate%s:\n",
                     tube_cut ? ", cut to the real region" : "");
              printf("about %.2f words per row, %.2f cells of the analysis per row.\n",
                     tube_cut ? tube_words_cut : tube_words, tube_cells);
            }
            if(sieve_mode == 2)
            { printf("Sieving the real region by passes over its ranges of words:\n");
              printf("about %.2f words per row in %.2f ranges, %.3f of the rows.\n",
                     reg_words, reg_ranges, reg_rows);
            }
            if(use2)
            { printf("The condition at 2 admits %.3f of the classes mod 64 and empties %.3f of the rows.\n",
                     density2, empty2);
            }
            printf("%ld primes used for the first stage of sieving,\n",
                   sieve_primes1);
            printf("%ld primes used for the first two stages together,\n",
                   sieve_primes2);
            printf("%ld primes used for all three stages together.\n",
                   sieve_primes3);
            if(costs_set)
            { size_t k;
              printf("Constants of the cost model set by -c:");
              for(k = 0; k < NUM_COSTS; k++)
              { if(cost_names[k].set)
                { printf(" %s=%g", cost_names[k].name, *cost_names[k].value); }
              }
              printf("\n");
            }
            break;
    case 5: printf("\ny^2 = "); print_poly(coeffs, total); printf("\n"); break;
    case 6: printf("max. Height = %ld\n", total);
            if(dbounded && dbound != height)
            { printf("Bound on the fourth coordinate = %ld\n", dbound); }
            break;
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
    case 12: printf("\n%ld candidates survived the first stage,\n", tot_surv1);
             printf("%ld candidates survived the second stage,\n", tot_surv2);
             printf("%ld candidates survived the third stage.\n", tot_surv3);
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
    case 8: printf("\nCould not start a thread.\n\n"); break;
    case 9: printf("\nThe polynomial must be squarefree.\n\n"); break;
    case 10: { size_t k;
               printf("\nThe argument of -c must be NAME=value, or several of them separated\n");
               printf("by commas, with a positive value and one of the names\n");
               for(k = 0; k < NUM_COSTS; k++)
               { printf(" %s", cost_names[k].name); }
               printf(".\n\n");
               break;
             }
    case 6: printf("\nWrong syntax for optional arguments:\n\n");
    case 2:
      printf("\n");
      printf("Usage: j-points 'a_0 a_1 ... a_d' max_height\n");
      printf("                [-n num_primes1] [-M num_primes2] [-N num_primes3]\n");
      printf("                [-p num_primes] [-s size] [-t threads] [-f format]\n");
      printf("                [-w bound4] [-c constants] [-1] [-q] [-a]\n");
      break;
  }
  fflush(stdout);
  exit(errno);
}
