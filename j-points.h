/******* j-points.h ******************************************************/

#ifndef DEBUG
#define DEBUG 0
#endif

/* MAX_PRIME_EVEN is the least power of 2 that is larger than MAX_PRIME */
#if (DEBUG > 0)
# define NUM_PRIMES 5
# define MAX_PRIME 13
# define MAX_PRIME_EVEN 16
#else
# define NUM_PRIMES 17
# define MAX_PRIME 61
# define MAX_PRIME_EVEN 64
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

typedef struct {long p; bit_array *ptr;} sieve_spec;

extern long prime[];
extern long pnn[];
extern int one_point;        /* A flag saying if one point is enough */
extern long sieve_primes1;
   /* The number of primes used for the first sieving stage */
extern long sieve_primes2;
   /* The number of primes used for both sieving stages */
extern bit_array 
       sieve[NUM_PRIMES][MAX_PRIME_EVEN][MAX_PRIME_EVEN][MAX_PRIME_EVEN];

extern long num_surv1;   /* Used to count the survivors of the first stage */
extern long num_surv2;   /* Used to count the survivors of the second stage */
extern bit_array *survivors; /* In this array the sieving takes place */
extern long height;          /* The height bound */
extern long array_size;      /* The size of the survivors array (in longs) */
extern bit_array begmask, endmask;
   /* Bit masks for the beginning and end of the sieving array */

extern sieve_spec sieves1[NUM_PRIMES];
extern sieve_spec sieves2p[NUM_PRIMES];
extern sieve_spec sieves2n[NUM_PRIMES];

extern int check_one_point(long, long, long);

extern int sift(long, long);
