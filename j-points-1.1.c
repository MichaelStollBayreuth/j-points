/******* j-points-1.1.c **************************************************
 *                                                                       * 
 * Change Log:                                                           *
 *                                                                       *
 * + version 1.1, 2006-08-09, Michael Stoll:                             *
 *   Removed a call 'mpz_mul_ui(&fff, &fff, (unsigned long)b)'           *
 *   in check_lifts that let the program miss points of the form         *
 *   (0 : b : c : d) when b was not a square.                            *
 *   (Was computing b^7*f(c/b) instead of b^6*f(c/b).)                   *
 *                                                                       *
 *************************************************************************/

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
  "This is j-points-1.1 by Michael Stoll (2006-08-10).\n\n" \
  "Please acknowledge use of the program in published work.\n"


#define mpz_mul_si(ROP, OP1, OP2) \
  { long mpz_mul_si_temp = (OP2); \
    if(mpz_mul_si_temp >= 0) \
      mpz_mul_ui((ROP), (OP1), (unsigned long)(mpz_mul_si_temp)); \
    else \
    { mpz_mul_ui((ROP), (OP1), (unsigned long)(-mpz_mul_si_temp)); \
      mpz_neg((ROP), (ROP)); \
  } }

/**************************************************************************
 * global variables                                                       *
 **************************************************************************/

#if (DEBUG > 0)
long prime[NUM_PRIMES+1] =  
{3,5,7,11,13};
#else
long prime[NUM_PRIMES+1] =  
{3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61};
#endif
long pnn[NUM_PRIMES+1];  /* This array holds the numbers (= index in prime[])
                            of the sieving primes */

double ratio1_def = 100;
double ratio2_def = 2;

double ratio1 = 0.0;  /* This will be set if the -r option is given */
double ratio2 = 0.0;  /* This will be set if the -R option is given */

bit_array bits[LONG_LENGTH]; /* An array of bit masks */

typedef struct {long p; double r; long n;} entry;

long inverses[NUM_PRIMES][MAX_PRIME_EVEN];
/* inverses[pn][a] = b such that a*b = 1 mod p = prime[pn]; 0 < a, b < p. */
char squares[NUM_PRIMES][MAX_PRIME];
/* squares[pn][x] = 1 if x is a square mod prime[pn], 0 if not */
char has_infinity[NUM_PRIMES];
/* has_infinity[pn] = 1 if the curve has rational points at infinity,
   0 if not */
char is_f_square[NUM_PRIMES][MAX_PRIME_EVEN];
/* is_f_square[pn][x] = 1 if f(x) is a square mod prime[pn], 0 if not */
char is_point_on_j[NUM_PRIMES][MAX_PRIME_EVEN][MAX_PRIME_EVEN];
/* is_point_on_j[pn][x][y] = 1 if there are points on J mod prime[pn]
   with first three coordinates on the Kummer surface (1,x,y), 0 if not */

bit_array sieve[NUM_PRIMES][MAX_PRIME_EVEN][MAX_PRIME_EVEN][MAX_PRIME_EVEN];
/* bit k of sieve[pn][a1][b1][c1] = 0
   iff (a,b,c) is excluded mod p = prime[pnn[pn]],
   when  c = c1*LONG_LENGTH + k mod p, b = b1 mod p and a = a1 mod p. */

MP_INT coeffs[7];  /* The coefficients of f */
MP_INT bc[7]; /* A helper array */
MP_INT fff, tmp, tmp2, tmp3, ddd;   /* Some multi-precision integer variables */
MP_INT k400, k310, k301, k220, k211, k202, k130, k121, k112, k103,
       k040, k031, k022, k013, k004;
        /* coefficients for Kummer equation */
MP_INT x12, x22, x32;
MP_INT kummer[3];
MP_INT cpf1, cpf2, cpf3;

long degree;
long coeffs_mod_p[NUM_PRIMES][8];
                         /* The coefficients of f reduced modulo the various
                            primes */
entry prec[NUM_PRIMES];  /* This array is used for sorting in order to
                            determine the `best' sieving primes. */

long height;          /* The height bound */
long w_height;        /* The height bound divided by the word length in bits */
long sieve_primes1;   /* The number of primes used for the first sieving stage */
long sieve_primes2;   /* The number of primes used for both sieving stages */
int quiet;            /* A flag saying whether to suppress messages */
int one_point;        /* A flag saying if one point is enough */
int all_points;       /* Indicates that the `-a' option was given */
char *print_format;   /* The printf format for printing points */
long array_size;      /* The size of the survivors array (in longs) */

bit_array *survivors; /* In this array the sieving takes place */

long num_surv1 = 0;   /* Used to count the survivors of the first stage */
long num_surv2 = 0;   /* Used to count the survivors of the second stage */

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
void init_sieve(void);
void kummer_init(void);
inline long gcd(long, long);
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
  }
  if(degree < 5) error(5);
  if(degree == 5) mpz_set_si(&coeffs[6], 0);
  if(sscanf(argv[2], " %ld", &height) != 1 || height < 1 || height > MAX_HEIGHT)
  { error(4); }
  w_height = ((height - 1) / LONG_LENGTH) + 1;
  /* Set global variables to their default values */
  sieve_primes1 = -1;  /* automatic determination of this number */
  sieve_primes2 = -1;  /* automatic determination of this number */
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
          else
          { if(sieve_primes2 >= 0 && sieve_primes1 > sieve_primes2)
              sieve_primes1 = sieve_primes2;
            else { if(sieve_primes1 > NUM_PRIMES) sieve_primes1 = NUM_PRIMES; }
          }
          i++;
          break;
        case 'N': /* number of primes used for sieving altogether */
          if(argc == i) error(6);
          i++;
          if(sscanf(argv[i], " %ld", &sieve_primes2) != 1) error(6);
          if(sieve_primes2 < 0) sieve_primes2 = 0;
          else
          { if(sieve_primes1 >= 0 && sieve_primes1 > sieve_primes2)
              sieve_primes2 = sieve_primes1;
            else { if(sieve_primes2 > NUM_PRIMES) sieve_primes2 = NUM_PRIMES; }
          }
          i++;
          break;
        case 'r': /* speed ratio 1 */
          if(argc == i) error(6);
          i++;
          if(sscanf(argv[i], " %lf", &ratio1) != 1) error(6);
          if(ratio1 <= 0.0) error(6);
          i++;
          break;
        case 'R': /* speed ratio 2 */
          if(argc == i) error(6);
          i++;
          if(sscanf(argv[i], " %lf", &ratio2) != 1) error(6);
          if(ratio2 <= 0.0) error(6);
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
        default: error(6);
  } } }
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

void init_fmodpsquare(void)
/* initialise is_f_square[][] and is_point_on_j[][][] */
{
  long a, b, c, pn, p, n, s, np, s1, s2, t;

  for(pn = 0; pn < NUM_PRIMES; pn++)
  {
    p = prime[pn];
    np = 1;
    /* compute coefficients mod p */
    for(n = 0; n <= 6; n++)
      coeffs_mod_p[pn][n] = mpz_fdiv_r_ui(&tmp, &coeffs[n], p);
    /* deal with (0,0,1) */
    has_infinity[pn] = (squares[pn][coeffs_mod_p[pn][6]] == p) ? 0 : 1;
    /* deal with (0,1,a) */
    if((is_f_square[pn][0] = (squares[pn][coeffs_mod_p[pn][0]] == p) ? 0 : 1))
    { np++; }
    for(a = 1 ; a < p; a++)
    {
      s = coeffs_mod_p[pn][degree];
      for(n = degree-1 ; n >= 0 ; n--)
      { s *= a;
        s += coeffs_mod_p[pn][n];
        s %= p;
      }
      if((is_f_square[pn][a] = (squares[pn][s] == p) ? 0 : 1)) np++;
    }
    /* deal with (1,b,c) */
    for(b = 0; b < p; b++)
    { long bb = ( ((1 + p*(3-(p&2)))/4) * b*b ) % p;
      for(c = 0; c < p; c++)
      { if(c == bb)
        { is_point_on_j[pn][b][c] = 1;
          np++;
          continue;
        }
        s1 = coeffs_mod_p[pn][degree];
        s2 = coeffs_mod_p[pn][degree-1];
        for(n = degree-2; n >= 0; n--)
        { t = s1;
          s1 *= b;
          s1 += s2;
          s1 %= p;
          s2 = (coeffs_mod_p[pn][n] + (p-c)*t)%p;
        }
        t = (s2*s2 + b*s1*s2 + c*s1*s1)%p;
        t = squares[pn][t];
        if(t == p)
        { is_point_on_j[pn][b][c] = 0; }
        else
        { long u1 = (2*s2 + b*s1 + 2*t)%p,
               u2 = (u1 + 4*(p - t))%p;
          if(squares[pn][u1] == p && squares[pn][u2] == p)
          { is_point_on_j[pn][b][c] = 0; }
          else
          { is_point_on_j[pn][b][c] = 1;
            np++;
    } } } }
    np *= p-1;
    np++; /* count origin */
    /* Fill array with info for p */
    prec[pn].p = p;
    prec[pn].n = pn;
    prec[pn].r = (double)np/(double)(p*p*p);
  }
  /* sort the array to get at the best primes */
  qsort(prec, NUM_PRIMES, sizeof(entry), compare_entries);
  /* Determine the number of sieving primes:
     First stage: take smallest n such that prod(j<=n) prob(j) * ratio1 < 1.
     Second stage: take smallest N such that (1-prob(N+1)) * ratio2 < 1. */
  if(sieve_primes1 < 0)
  { long n = 0, m = (sieve_primes2 >= 0) ? sieve_primes2-1 : NUM_PRIMES-1;
    double prod;
    if(ratio1 == 0.0) ratio1 = ratio1_def;
    prod = ratio1;
    for(n = 0; n < m; n++)
    { prod *= prec[n].r;
      if(prod < 1.0) break;
    }
    sieve_primes1 = n + 1;
  }
  if(sieve_primes2 < 0)
  { long n;
    if(ratio2 == 0.0) ratio2 = ratio2_def;
    for(n = sieve_primes1; n < NUM_PRIMES; n++)
      if(ratio2*(1.0 - prec[n].r) < 1.0) break;
    sieve_primes2 = n;
  }
  for(n = 0; n < sieve_primes1; n++)
  { pnn[n] = prec[n].n;
    sieves1[n].p = prime[pnn[n]];
  }
  for( ; n < sieve_primes2; n++)
  { pnn[n] = prec[n].n;
    sieves2p[n-sieve_primes1].p = sieves2n[n-sieve_primes1].p = prime[pnn[n]];
  }
  if(!quiet)
  { message(4, 0);
    if(sieve_primes2 > 0)
    { message(7, 0); message(8, 0); }
  }
  return;
}

/* help is an array for temporarily storing the sieving information.
   Bit j of help[a][b][c0] says whether (a : b : c) are the first three
   coordinates of a point on K(F_p), where c = c0*LONG_LENGTH + j. 
   c runs from 0 to k*p - 1, where k*p is the least multiple of p exceeding 
   LONG_LENGTH. */
bit_array help[MAX_PRIME][MAX_PRIME_EVEN][MAX_PRIME / LONG_LENGTH + 2];

/* initalise sieve[][][][] */
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

    /* fill the bit pattern from help[][][] into sieve[pn][][][].
	sieve[pn][a][b][c0] has the same semantics as help[a][b][c0],
	but here, c0 runs from 0 to p-1 and all bits are filled. */
    for(a = 0; a < p; a++)
      for(b = 0; b < p; b++)
      { bit_array *si = sieve[pn][a][b];
        bit_array *he = help[a][b];
        long p1 = (LONG_LENGTH/p + 1) * p;
        long diff_shift = p1 & LONG_MASK;
        long diff = LONG_LENGTH - diff_shift;
        bit_array diff_mask = ~(bit_array)(-1L<<diff);
        long c1;
        long wp = p1>>LONG_SHIFT;
        /* copy the first chunk from help[a][b][] into sieve[pn][a][b][] */
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
	  printf(" %8.8lx", sieve[pn][a][b][c]);
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
  struct stype {long p; long ap; char *ptr;} ssp0[NUM_PRIMES];
  struct stype *ssp;
  for(n = 0; n < sieve_primes2; n++)
    ssp0[n].p = prime[pnn[n]];
  if(degree == 5 && mpz_cmp_si(&coeffs[5], 1) == 0)
  { /* can use only squares for a */
    long aa;
    for(aa = 1; (a = aa*aa) <= wheight; aa++)
    { long bmax;
      for(ssp = &ssp0[0], n = 0; n < sieve_primes2; n++, ssp++)
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
      { if(gcd(a, b) == 1)
        { for(ssp = &ssp0[0], n = sieve_primes2; n; n--, ssp++)
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
      for(ssp = &ssp0[0], n = 0; n < sieve_primes2; n++, ssp++)
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
      { if(gcd(a, b) == 1)
        { for(ssp = &ssp0[0], n = sieve_primes2; n; n--, ssp++)
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

void find_points(void)
{
  long a, b;
  
  /* initialise is_f_square[][] */
  init_fmodpsquare();
  /* initalise sieve[][][][] */
  init_sieve();
  if(sieve_primes2 > 0 && prec[0].r == 0.0)
  { if(!quiet) message(1,0); return; }
  survivors = (bit_array *)malloc(array_size*sizeof(bit_array));
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
      for(b = (a == 0) ? 1 : -height; b <= height; b++)
      { 
#ifdef VERBOSE
        printf(" a = %ld, b = %ld\n", aa, b);
#endif
        if(sift(aa, b) && one_point) return;
      }
  }
  else
  { for(a = 0; a <= height; a++)
      for(b = (a == 0) ? 1 : -height; b <= height; b++)
      { 
#ifdef VERBOSE
        printf(" a = %ld, b = %ld\n", a, b);
#endif
        if(sift(a, b) && one_point) return;
      }
  }
  return;
}

/**************************************************************************
 * some Kummer surface formula stuff                                      *
 **************************************************************************/

void kummer_init()
{ /* -4*x1^4*f0*f2 + x1^4*f1^2 */
  mpz_mul(&k400, &coeffs[1], &coeffs[1]); mpz_mul(&tmp, &coeffs[0], &coeffs[2]);
  mpz_mul_ui(&tmp, &tmp, 4); mpz_sub(&k400, &k400, &tmp);
  /* - 4*x1^3*x2*f0*f3 */
  mpz_mul(&k310, &coeffs[0], &coeffs[3]); mpz_mul_ui(&k310, &k310, 4);
  mpz_neg(&k310, &k310);
  /* - 2*x1^3*x3*f1*f3 */
  mpz_mul(&k301, &coeffs[1], &coeffs[3]); mpz_mul_ui(&k301, &k301, 2);
  mpz_neg(&k301, &k301);
  /* - 4*x1^2*x2^2*f0*f4 */
  mpz_mul(&k220, &coeffs[0], &coeffs[4]); mpz_mul_ui(&k220, &k220, 4);
  mpz_neg(&k220, &k220);
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
  mpz_mul(&k130, &coeffs[0], &coeffs[5]); mpz_mul_ui(&k130, &k130, 4);
  mpz_neg(&k130, &k130);
  /* + 8*x1*x2^2*x3*f0*f6 - 4*x1*x2^2*x3*f1*f5 */
  mpz_mul(&k121, &coeffs[0], &coeffs[6]); mpz_mul_ui(&k121, &k121, 8);
  mpz_mul(&tmp, &coeffs[1], &coeffs[5]); mpz_mul_ui(&tmp, &tmp, 4);
  mpz_sub(&k121, &k121, &tmp);
  /* + 4*x1*x2*x3^2*f1*f6 - 4*x1*x2*x3^2*f2*f5 */
  mpz_mul(&k112, &coeffs[1], &coeffs[6]); mpz_mul(&tmp, &coeffs[2], &coeffs[5]);
  mpz_sub(&k112, &k112, &tmp); mpz_mul_ui(&k112, &k112, 4);
  /* - 2*x1*x3^3*f3*f5 */
  mpz_mul(&k103, &coeffs[3], &coeffs[5]); mpz_mul_ui(&k103, &k103, 2);
  mpz_neg(&k103, &k103);
  /* - 4*x2^4*f0*f6 */
  mpz_mul(&k040, &coeffs[0], &coeffs[6]); mpz_mul_ui(&k040, &k040, 4);
  mpz_neg(&k040, &k040);
  /* - 4*x2^3*x3*f1*f6 */
  mpz_mul(&k031, &coeffs[1], &coeffs[6]); mpz_mul_ui(&k031, &k031, 4);
  mpz_neg(&k031, &k031);
  /* - 4*x2^2*x3^2*f2*f6 */
  mpz_mul(&k022, &coeffs[2], &coeffs[6]); mpz_mul_ui(&k022, &k022, 4);
  mpz_neg(&k022, &k022); 
  /* - 4*x2*x3^3*f3*f6 */
  mpz_mul(&k013, &coeffs[3], &coeffs[6]); mpz_mul_ui(&k013, &k013, 4);
  mpz_neg(&k013, &k013);
  /* - 4*x3^4*f4*f6 + x3^4*f5^2 */
  mpz_mul(&k004, &coeffs[5], &coeffs[5]); mpz_mul(&tmp, &coeffs[4], &coeffs[6]);
  mpz_mul_ui(&tmp, &tmp, 4); mpz_sub(&k004, &k004, &tmp);
}

inline void kummer_eqn(long a, long b, long c)
{ mpz_set_si(&x12, a*a); mpz_set_si(&x22, b*b); mpz_set_si(&x32, c*c);
  /* -4*x1^4*f0*f2 + x1^4*f1^2 - 4*x1^3*x2*f0*f3 - 2*x1^3*x3*f1*f3
     - 4*x1^2*x2^2*f0*f4 + 4*x1^2*x2*x3*f0*f5 - 4*x1^2*x2*x3*f1*f4
     - 4*x1^2*x3^2*f0*f6 + 2*x1^2*x3^2*f1*f5 - 4*x1^2*x3^2*f2*f4 
     + x1^2*x3^2*f3^2 - 4*x1*x2^3*f0*f5 + 8*x1*x2^2*x3*f0*f6 
     - 4*x1*x2^2*x3*f1*f5 + 4*x1*x2*x3^2*f1*f6 - 4*x1*x2*x3^2*f2*f5 
     - 2*x1*x3^3*f3*f5 - 4*x2^4*f0*f6 - 4*x2^3*x3*f1*f6 - 4*x2^2*x3^2*f2*f6 
     - 4*x2*x3^3*f3*f6 - 4*x3^4*f4*f6 + x3^4*f5^2 */
  mpz_mul(&tmp, &k400, &x12); mpz_mul(&kummer[0], &tmp, &x12);
  
  mpz_mul(&tmp, &k310, &x12); mpz_mul_si(&tmp, &tmp, a*b);
  mpz_add(&kummer[0], &kummer[0], &tmp);
  
  mpz_mul(&tmp, &k301, &x12); mpz_mul_si(&tmp, &tmp, a*c);
  mpz_add(&kummer[0], &kummer[0], &tmp);
  
  mpz_mul(&tmp, &k220, &x12); mpz_mul(&tmp, &tmp, &x22);
  mpz_add(&kummer[0], &kummer[0], &tmp);
  
  mpz_mul(&tmp, &k211, &x12); mpz_mul_si(&tmp, &tmp, b*c);
  mpz_add(&kummer[0], &kummer[0], &tmp);
  
  mpz_mul(&tmp, &k202, &x12); mpz_mul(&tmp, &tmp, &x32);
  mpz_add(&kummer[0], &kummer[0], &tmp);
  
  mpz_mul(&tmp, &k130, &x22); mpz_mul_si(&tmp, &tmp, a*b);
  mpz_add(&kummer[0], &kummer[0], &tmp);
  
  mpz_mul(&tmp, &k121, &x22); mpz_mul_si(&tmp, &tmp, a*c);
  mpz_add(&kummer[0], &kummer[0], &tmp);
  
  mpz_mul(&tmp, &k112, &x32); mpz_mul_si(&tmp, &tmp, a*b);
  mpz_add(&kummer[0], &kummer[0], &tmp);
  
  mpz_mul(&tmp, &k103, &x32); mpz_mul_si(&tmp, &tmp, a*c);
  mpz_add(&kummer[0], &kummer[0], &tmp);
  
  mpz_mul(&tmp, &k040, &x22); mpz_mul(&tmp, &tmp, &x22);
  mpz_add(&kummer[0], &kummer[0], &tmp);
  
  mpz_mul(&tmp, &k031, &x22); mpz_mul_si(&tmp, &tmp, b*c);
  mpz_add(&kummer[0], &kummer[0], &tmp);
  
  mpz_mul(&tmp, &k022, &x22); mpz_mul(&tmp, &tmp, &x32);
  mpz_add(&kummer[0], &kummer[0], &tmp);
  
  mpz_mul(&tmp, &k013, &x32); mpz_mul_si(&tmp, &tmp, b*c);
  mpz_add(&kummer[0], &kummer[0], &tmp);
  
  mpz_mul(&tmp, &k004, &x32); mpz_mul(&tmp, &tmp, &x32);
  mpz_add(&kummer[0], &kummer[0], &tmp);
  /* -4*x1^3*f0 - 2*x1^2*x2*f1 - 4*x1^2*x3*f2 - 2*x1*x2*x3*f3 
     - 4*x1*x3^2*f4 - 2*x2*x3^2*f5 - 4*x3^3*f6 */
  mpz_mul_si(&tmp, &coeffs[0], a); mpz_mul(&kummer[1], &tmp, &x12);
  mpz_mul(&tmp, &coeffs[2], &x12); mpz_mul_si(&tmp, &tmp, c);
  mpz_add(&kummer[1], &kummer[1], &tmp); mpz_mul_si(&tmp, &coeffs[4], a);
  mpz_mul(&tmp, &tmp, &x32); mpz_add(&kummer[1], &kummer[1], &tmp);
  mpz_mul_si(&tmp, &coeffs[6], c); mpz_mul(&tmp, &tmp, &x32);
  mpz_add(&kummer[1], &kummer[1], &tmp); mpz_mul_ui(&kummer[1], &kummer[1], 2);
  mpz_mul(&tmp, &coeffs[1], &x12); mpz_mul_si(&tmp, &tmp, b);
  mpz_add(&kummer[1], &kummer[1], &tmp); mpz_mul_si(&tmp, &coeffs[3], a);
  mpz_mul_si(&tmp, &tmp, b*c);
  mpz_add(&kummer[1], &kummer[1], &tmp); mpz_mul_si(&tmp, &coeffs[5], b);
  mpz_mul(&tmp, &tmp, &x32); mpz_add(&kummer[1], &kummer[1], &tmp);
  mpz_mul_ui(&kummer[1], &kummer[1], 2); mpz_neg(&kummer[1], &kummer[1]);
  /* x2^2-4*x1*x3 */
  mpz_set_si(&tmp, a*c); mpz_mul_ui(&tmp, &tmp, 4);
  mpz_sub(&kummer[2], &x22, &tmp);
}

/**************************************************************************
 * compute the gcd of two integers                                        *
 **************************************************************************/

inline long gcd(long m, long n)
{
  if(m < 0) m = -m;
  if(n < 0) n = -n;
  while(1)
  { if(n == 0) return(m);
    m %= n;
    if(m == 0) return(n);
    n %= m;
  }
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
      { /* Why isn't there a `mpz_mul_si' ? */
        if(c > 0)
          mpz_mul_ui(&fff, &fff, (unsigned long)c);
        else
        { mpz_mul_ui(&fff, &fff, (unsigned long)(-c));
          mpz_neg(&fff, &fff);
        }
        mpz_add(&fff, &fff, &bc[k]);
      }
      /* mpz_mul_ui(&fff, &fff, (unsigned long)b); removed for version 1.1 */
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
    
    k = b*b - a*c; mpz_mul_si(&tmp2, &coeffs[5], k);
    mpz_mul_si(&tmp2, &tmp2, a*b); mpz_add(&tmp, &tmp, &tmp2);
    
    mpz_mul_si(&tmp2, &coeffs[6], k); mpz_mul_si(&tmp2, &tmp2, k);
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

int check_one_point_final(long a, long b, long c, MP_INT *d1, MP_INT *d2)
{ /* Given a, b, c and d = d1/d2, check if this gives a point satisfying
     the height condition. If so, print and count it. */
  long m = (a > abs(b)) ? ((a > abs(c)) ? a :abs(c))
                        : ((abs(b) > abs(c)) ? abs(b) : abs(c));
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
    case 4: printf("%ld primes used for first stage of sieving,\n", 
                   sieve_primes1);
            printf("%ld primes used for both stages of sieving together.\n",
                   sieve_primes2);
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
              printf("\n");
              break;
            }
    case 8:
      printf("Probabilities: Min(%ld) = %f, Cut1(%ld) = %f, ",
              prec[0].p, prec[0].r, 
              prec[sieve_primes1-1].p, prec[sieve_primes1-1].r);
      printf("Cut2(%ld) = %f, Max(%ld) = %f\n\n",
              prec[sieve_primes2-1].p, prec[sieve_primes2-1].r,
              prec[NUM_PRIMES-1].p, prec[NUM_PRIMES-1].r);  break;
    case 9: printf("Using speed ratios %f and %f\n", ratio1, ratio2); break;
    case 12: printf("\n%ld candidates survived the first stage,\n", num_surv1);
             printf("%ld candidates survived the second stage.\n", num_surv2);
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
    case 6: printf("\nWrong syntax for optional arguments:\n\n");
    case 2:
      printf("\n");
      printf("Usage: j-points 'a_0 a_1 ... a_n' max_height\n");
      printf("                [-n num_primes1] [-N num_primes2]\n");
      printf("                [-r ratio1] [-R ratio2] [-f format] [-s size]\n");
      break;
  }
  fflush(stdout);
  exit(errno);
}
