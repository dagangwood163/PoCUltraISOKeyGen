/*
   touch lippar.h
   compile lip.c to produce lip.o
   compile and run timer.c (using lip.o), produces non-empty lippar.h
   recompile lip.c to produce new lip.o
*/


#include <stddef.h>
#include <stdio.h>
#include <math.h>
#include <malloc.h>
#include <signal.h>
#ifdef WIN32
#include <time.h>
#else
#include <sys/time.h>
#include <sys/resource.h>
#endif
#include "lip.h"


/* do NUMBER timings of integers of size FIRST + i * INCR bits,
for 0 <= i < NUMBER */

long FIRST = 59;
long INCR = 120;
long NUMBER = 5;

#define ITER	10000		/* # of iterations for timings */
long SEED=11L;			/* for randomness */

/* do not worry about anything below this point */


#ifndef KAR_DEPTH
#define KAR_DEPTH 20
#endif
long K_M_C = 20;
long K_S_C = 20;

/*Function prototypes for static internal functions.*/


static void zhalt(
	char *c
	);

#define zaddmulone(ama, amb) \
{ \
	register long lami; \
	register long lams = 0; \
	register verylong lama = (ama); \
	register verylong lamb = (amb); \
 \
	lams = 0; \
	for (lami = (*lamb++); lami > 0; lami--) \
	{ \
		lams += (*lama + *lamb++); \
		*lama++ = lams & RADIXM; \
		lams >>= NBITS; \
	} \
	*lama += lams; \
}


/* for long division */
static double log10rad = -1.0;
static double epsilon;
static double fradix = (double)RADIX;
static double fudge = -1.0;
static double fudge1 = -1.0;
#ifdef ALPHA
static double fudge2 = -1.0;
#endif
#ifdef ALPHA50
static double alpha50fudge = -1.0;
static double alpha50fradix = (double) ALPHA50RADIX;
#endif
/* for random generator */
static verylong zseed = 0;
static verylong zranp = 0;
static verylong zprroot = 0;
/* for convenience */
static long oner[] = {1, 1, 1};
static long glosho[] = {1, 1, 0};
static verylong one = &oner[1];
/* for karatsuba */
static verylong kar_mem[5*KAR_DEPTH];
static long kar_mem_initialized = 0;


static void
zhalt(
        char *c
	)
{
#ifdef NOHALT
	fprintf(stderr,"error:\n   %s\ncontinue...\n",c);
#else
	fprintf(stderr,"fatal error:\n   %s\nexit...\n",c);
	(void)exit((int)0);
#endif
}



/* SINGLE_MUL=0, PLAIN=0, KARAT=1 */
void
zmul3(
        verylong a,
        verylong b,
        verylong *c
        )
{       /* output not input */
        register long aneg;
        register long bneg;
        verylong olda;
        verylong oldb;

        if (ALLOCATE && (!a || !b))
        {
                zzero(c);
                return;
        }
        if (a == b)
        {
                zsq(a, c);
                return;
        }
        if (!kar_mem_initialized)
        {
                kar_mem_initialized = 1;
                for (aneg = (5 * KAR_DEPTH) - 1; aneg >= 0; aneg--)
                        kar_mem[aneg] = (verylong) 0;
        }
        olda = a;
        oldb = b;
        if (aneg = (*a < 0))
                a[0] = -a[0];
        if (bneg = (*b < 0))
                b[0] = -b[0];
        if (*a > *b)
                kar_mul3(a, b, c, (long) 0);
        else
                kar_mul3(b, a, c, (long) 0);
        if (aneg != bneg && ((*c)[1] || (*c)[0] != 1))
                (*c)[0] = -(*c)[0];
        if (aneg)
                olda[0] = -olda[0];
        if (bneg)
                oldb[0] = -oldb[0];
}


void
zsq3(
        verylong a,
        verylong *c
        )
{       /* output is not input */
        register long aneg;

        if (ALLOCATE && !a)
        {
                zzero(c);
                return;
        }
        if (!kar_mem_initialized)
        {
                kar_mem_initialized = 1;
                for (aneg = (5 * KAR_DEPTH) - 1; aneg >= 0; aneg--)
                        kar_mem[aneg] = (verylong) 0;
        }
        if (aneg = (*a < 0))
                a[0] = -a[0];
        kar_sq3(a, c, (long) 0);
        if (aneg)
                a[0] = -a[0];
}


#define zaddmulp3(  a,  b,  d,  t) { \
	register unsigned long aa= (*a) + (*t); \
        register unsigned long b1 = b & RADIXROOTM; \
        register unsigned long d1 = d & RADIXROOTM; \
        register unsigned long bd,b1d1,m; \
        register unsigned long ld = (d>>NBITSH); \
        register unsigned long lb = (b>>NBITSH); \
 \
        bd=lb*ld; \
        b1d1=b1*d1; \
        m=(lb+b1)*(ld+d1) - bd - b1d1; \
        aa += ( b1d1+ ((m&RADIXROOTM)<<NBITSH)); \
        (*t) = (aa >> NBITS) + bd + (m>>NBITSH); \
        (*a) = aa & RADIXM; \
}


#define zaddmul3(ams, ama, amb) \
{ \
        register long lami; \
        register long lams = (ams); \
        register verylong lama = (ama); \
        register verylong lamb = (amb); \
        long lamcarry = 0; \
 \
        for (lami = (*lamb++); lami > 0; lami--) \
        { \
                zaddmulp3(lama, *lamb, lams, &lamcarry); \
        /* Be careful, the last lama is unnormalized */ \
                lama++; \
                lamb++; \
        } \
        *lama += lamcarry; \
}


kar_mul3(
	verylong a,
	verylong b,
	verylong *c,
	long shi
	)
{
	register long al;
	register long hal;
	register long i;
	register long restoreb0 = b[0];
	register verylong pc;
	register long bbig = 1;
	verylong *a0;
	verylong *a1;
	verylong *a2;
	verylong *a3;
	verylong *a4;

	zsetlength(c, (hal = (al = a[0]) + (i = b[0])), "in kar_mul, third argument");
	if ((shi >= (5 * KAR_DEPTH)) || (al < K_M_C) || (i < K_M_C))
	{
		pc = &(*c)[1];
		for (i = hal; i > 0; i--)
			*pc++ = 0;
		pc = &(*c)[1];
		if (al <= *b)
			for (i = al; i; i--)
			{
				zaddmul3(*(++a), pc++, b);
			}
		else
			for (i = *b; i; i--)
			{
				zaddmul3(*(++b), pc++, a);
			}
		while ((hal > 1) && (!((*c)[hal])))
			hal--;
		(*c)[0] = hal;
		return;
	}
	a0 = &(kar_mem[shi]);
	a1 = &(kar_mem[shi + 1]);
	a2 = &(kar_mem[shi + 2]);
	a3 = &(kar_mem[shi + 3]);
	a4 = &(kar_mem[shi + 4]);
	hal = ((al + 1) >> 1);
	zsetlength(a0, al, "in kar_mul, locals\n");
	zsetlength(a1, al, "");
	zsetlength(a2, al, "");
	zsetlength(a3, al + hal, "");
	zsetlength(a4, al + 2, "");
	i = hal;
	while ((i > 1) && (!(a[i])))
		i--;
	a[0] = i;
	if (hal >= b[0])
		bbig = 0;
	else
	{
		i = hal;
		while ((i > 1) && (!(b[i])))
			i--;
		b[0] = i;
	}
	for (i = hal + 1; i <= al; i++)
		(*a1)[i - hal] = a[i];
	(*a1)[0] = al - hal;
	if (bbig)
	{
		for (i = hal + 1; i <= restoreb0; i++)
			(*a3)[i - hal] = b[i];
		(*a3)[0] = restoreb0 - hal;
	}
	kar_mul3(a, b, a4, shi + 5);
	zadd(a, (*a1), a0);
	a[0] = al;
	if (bbig)
	{
		kar_mul3((*a1), (*a3), c, shi + 5);
		zadd(b, (*a3), a2);
		b[0] = restoreb0;
		kar_mul3((*a0), (*a2), a3, shi + 5);
	}
	else
		kar_mul3((*a0), b, a3, shi + 5);
	zsubpos((*a3), (*a4), a3);
	if (bbig)
		zsubpos((*a3), *c, a3);
	zlshift((*a3), hal * NBITS, a3);
	hal <<= 1;
	if (bbig)
	{
		for (i = (*c)[0]; i; i--)
			(*c)[i + hal] = (*c)[i];
		for (i = hal; i > (*a4)[0]; i--)
			(*c)[i] = 0;
		for (; i; i--)
			(*c)[i] = (*a4)[i];
		(*c)[0] += hal;
	}
	else
	{
		for (i = (*a4)[0]; i >= 0; i--)
			(*c)[i] = (*a4)[i];
	}
	zadd(*c, (*a3), c);
}


#define zaddmulsq3(sql, sqa, sqb) \
{ \
        register long lsqi = (sql); \
        register long lsqs = *(sqb); \
        register verylong lsqa = (sqa); \
        register verylong lsqb = (sqb); \
        long lsqcarry = 0; \
 \
        lsqb++; \
        for (; lsqi > 0; lsqi--) \
        { \
                zaddmulp3(lsqa, *lsqb, lsqs, &lsqcarry); \
                lsqa++; \
                lsqb++; \
        } \
        *lsqa += lsqcarry; \
/* Be careful, the last lama is unnormalized */ \
}


#define zaddmulpsq3(_a, _b, _t) \
{ \
        register long lb = (_b); \
        register long b1 = (_b) & RADIXROOTM; \
        register long aa = *(_a) + b1 * b1; \
 \
        b1 = (b1 * (lb >>= NBITSH) << 1) + (aa >> NBITSH); \
        aa = (aa & RADIXROOTM) + ((b1 & RADIXROOTM) << NBITSH); \
        *(_t) = lb * lb + (b1 >> NBITSH) + (aa >> NBITS); \
        *(_a) = (aa & RADIXM); \
}



kar_sq3(
	verylong a,
	verylong *c,
	long shi
	)
{
	register long al;
	register long hal;
	register long i;
	register verylong pc;
	verylong *a0;
	verylong *a1;
	verylong *a2;

	zsetlength(c, (i = ((al = a[0]) << 1)), "in kar_sq, second argument");
	if ((shi >= (3 * KAR_DEPTH)) || (al < K_S_C))
	{
		register unsigned long uncar;
		long carry = 0;
		pc = &(*c)[1];
		for (; i > 0; i--)
			*pc++ = 0;
		for (hal = 1; hal <= al; hal++)
		{
			i += 2;
			{
				zaddmulsq3(al - hal, &((*c)[i]), &(a[hal]));
			}
			uncar = ((*c)[i - 1] << 1) + carry;
			(*c)[i - 1] = uncar & RADIXM;
			uncar = ((*c)[i] << 1) + (uncar >> NBITS);
			{
				zaddmulpsq3(&(*c)[i - 1], a[hal], &carry);
			}
			uncar += carry;
			carry = uncar >> NBITS;
			(*c)[i] = uncar & RADIXM;
		}
		while ((i > 1) && (!((*c)[i])))
			i--;
		(*c)[0] = i;
		return;
	}
	a0 = &(kar_mem[shi]);
	a1 = &(kar_mem[shi + 1]);
	a2 = &(kar_mem[shi + 2]);
	hal = ((al + 1) >> 1);
	zsetlength(a0, al + hal + 2, "in kar_sq, locals\n");
	zsetlength(a1, al + 2, "");
	zsetlength(a2, al, "");
	i = hal;
	while ((i > 1) && (!(a[i])))
		i--;
	a[0] = i;
	for (i = hal + 1; i <= al; i++)
		(*a0)[i - hal] = a[i];
	(*a0)[0] = al - hal;
	kar_sq3(a, a1, shi + 3);
	zadd(a, (*a0), a2);
	kar_sq3((*a0), c, shi + 3);
	a[0] = al;
	kar_sq3((*a2), a0, shi + 3);
	zsubpos((*a0), (*a1), a0);
	zsubpos((*a0), *c, a0);
	zlshift((*a0), hal * NBITS, a0);
	hal <<= 1;
	for (i = (*c)[0]; i; i--)
		(*c)[i + hal] = (*c)[i];
	for (i = hal; i > (*a1)[0]; i--)
		(*c)[i] = 0;
	for (; i; i--)
		(*c)[i] = (*a1)[i];
	(*c)[0] += hal;
	zadd(*c, (*a0), c);
}


/* SINGLE_MUL=0, PLAIN=0, KARAT=0 */

void
zmul0(
        verylong a,
        verylong b,
        verylong *c
        )
{       /* output not input */
        register long aneg;
        register long bneg;
        verylong olda;
        verylong oldb;

/*
printf("zmul0\n"); fflush(stdout);
*/
        if (ALLOCATE && (!a || !b))
        {
                zzero(c);
                return;
        }
        if (a == b)
        {
                zsq(a, c);
                return;
        }
        if (!kar_mem_initialized)
        {
                kar_mem_initialized = 1;
                for (aneg = (5 * KAR_DEPTH) - 1; aneg >= 0; aneg--)
                        kar_mem[aneg] = (verylong) 0;
        }
        olda = a;
        oldb = b;
        if (aneg = (*a < 0))
                a[0] = -a[0];
        if (bneg = (*b < 0))
                b[0] = -b[0];
        if (*a > *b)
                kar_mul0(a, b, c, (long) 0);
        else
                kar_mul0(b, a, c, (long) 0);
        if (aneg != bneg && ((*c)[1] || (*c)[0] != 1))
                (*c)[0] = -(*c)[0];
        if (aneg)
                olda[0] = -olda[0];
        if (bneg)
                oldb[0] = -oldb[0];
}


void
zsq0(
        verylong a,
        verylong *c
        )
{       /* output is not input */
        register long aneg;

        if (ALLOCATE && !a)
        {
                zzero(c);
                return;
        }
        if (!kar_mem_initialized)
        {
                kar_mem_initialized = 1;
                for (aneg = (5 * KAR_DEPTH) - 1; aneg >= 0; aneg--)
                        kar_mem[aneg] = (verylong) 0;
        }
        if (aneg = (*a < 0))
                a[0] = -a[0];
        kar_sq0(a, c, (long) 0);
        if (aneg)
                a[0] = -a[0];
}

static double fradix_inv0 = 1.0 / RADIX; /* Global constant */
#define zaddmulp0(_a, _b, _d, _t) \
{ \
        register at = *(_a) + *(_t); \
        register long aa = (at + (_b) * (_d)) & RADIXM; \
 \
        *(_t) = (long) (0.25 + fradix_inv0 * (((double) (at - aa)) \
                       + ((double) (_b)) * ((double) (_d)))); \
        *(_a) = aa; \
}


#define zaddmul0(ams, ama, amb) \
{ \
        register long lami; \
        register long lams = (ams); \
        register verylong lama = (ama); \
        register verylong lamb = (amb); \
        long lamcarry = 0; \
 \
        for (lami = (*lamb++); lami > 0; lami--) \
        { \
                zaddmulp0(lama, *lamb, lams, &lamcarry); \
        /* Be careful, the last lama is unnormalized */ \
                lama++; \
                lamb++; \
        } \
        *lama += lamcarry; \
}


kar_mul0(
	verylong a,
	verylong b,
	verylong *c,
	long shi
	)
{
	register long al;
	register long hal;
	register long i;
	register long restoreb0 = b[0];
	register verylong pc;
	register long bbig = 1;
	verylong *a0;
	verylong *a1;
	verylong *a2;
	verylong *a3;
	verylong *a4;

/*
printf("karmul0\n"); fflush(stdout);
*/
	zsetlength(c, (hal = (al = a[0]) + (i = b[0])), "in kar_mul, third argument");
	if ((shi >= (5 * KAR_DEPTH)) || (al < K_M_C) || (i < K_M_C))
	{
		pc = &(*c)[1];
		for (i = hal; i > 0; i--)
			*pc++ = 0;
		pc = &(*c)[1];
		if (al <= *b)
			for (i = al; i; i--)
			{
				zaddmul0(*(++a), pc++, b);
			}
		else
			for (i = *b; i; i--)
			{
				zaddmul0(*(++b), pc++, a);
			}
		while ((hal > 1) && (!((*c)[hal])))
			hal--;
		(*c)[0] = hal;
		return;
	}
	a0 = &(kar_mem[shi]);
	a1 = &(kar_mem[shi + 1]);
	a2 = &(kar_mem[shi + 2]);
	a3 = &(kar_mem[shi + 3]);
	a4 = &(kar_mem[shi + 4]);
	hal = ((al + 1) >> 1);
	zsetlength(a0, al, "in kar_mul, locals\n");
	zsetlength(a1, al, "");
	zsetlength(a2, al, "");
	zsetlength(a3, al + hal, "");
	zsetlength(a4, al + 2, "");
	i = hal;
	while ((i > 1) && (!(a[i])))
		i--;
	a[0] = i;
	if (hal >= b[0])
		bbig = 0;
	else
	{
		i = hal;
		while ((i > 1) && (!(b[i])))
			i--;
		b[0] = i;
	}
	for (i = hal + 1; i <= al; i++)
		(*a1)[i - hal] = a[i];
	(*a1)[0] = al - hal;
	if (bbig)
	{
		for (i = hal + 1; i <= restoreb0; i++)
			(*a3)[i - hal] = b[i];
		(*a3)[0] = restoreb0 - hal;
	}
	kar_mul0(a, b, a4, shi + 5);
	zadd(a, (*a1), a0);
	a[0] = al;
	if (bbig)
	{
		kar_mul0((*a1), (*a3), c, shi + 5);
		zadd(b, (*a3), a2);
		b[0] = restoreb0;
		kar_mul0((*a0), (*a2), a3, shi + 5);
	}
	else
		kar_mul0((*a0), b, a3, shi + 5);
	zsubpos((*a3), (*a4), a3);
	if (bbig)
		zsubpos((*a3), *c, a3);
	zlshift((*a3), hal * NBITS, a3);
	hal <<= 1;
	if (bbig)
	{
		for (i = (*c)[0]; i; i--)
			(*c)[i + hal] = (*c)[i];
		for (i = hal; i > (*a4)[0]; i--)
			(*c)[i] = 0;
		for (; i; i--)
			(*c)[i] = (*a4)[i];
		(*c)[0] += hal;
	}
	else
	{
		for (i = (*a4)[0]; i >= 0; i--)
			(*c)[i] = (*a4)[i];
	}
	zadd(*c, (*a3), c);
}


#define zaddmulsq0(sql, sqa, sqb) \
{ \
        register long lsqi = (sql); \
        register long lsqs = *(sqb); \
        register verylong lsqa = (sqa); \
        register verylong lsqb = (sqb); \
        long lsqcarry = 0; \
 \
        lsqb++; \
        for (; lsqi > 0; lsqi--) \
        { \
                zaddmulp0(lsqa, *lsqb, lsqs, &lsqcarry); \
                lsqa++; \
                lsqb++; \
        } \
        *lsqa += lsqcarry; \
/* Be careful, the last lama is unnormalized */ \
}


#define zaddmulpsq0(_a, _b, _t) \
{ \
        register at = *(_a); \
        register long aa = (at + (_b) * (_b)) & RADIXM; \
 \
        *(_t) = (long) (0.25 + fradix_inv0 * (((double) (at - aa)) \
			      + ((double) (_b)) * ((double) (_b)))); \
        *(_a) = aa; \
}


kar_sq0(
	verylong a,
	verylong *c,
	long shi
	)
{
	register long al;
	register long hal;
	register long i;
	register verylong pc;
	verylong *a0;
	verylong *a1;
	verylong *a2;

	zsetlength(c, (i = ((al = a[0]) << 1)), "in kar_sq, second argument");
	if ((shi >= (3 * KAR_DEPTH)) || (al < K_S_C))
	{
		register unsigned long uncar;
		long carry = 0;
		pc = &(*c)[1];
		for (; i > 0; i--)
			*pc++ = 0;
		for (hal = 1; hal <= al; hal++)
		{
			i += 2;
			{
				zaddmulsq0(al - hal, &((*c)[i]), &(a[hal]));
			}
			uncar = ((*c)[i - 1] << 1) + carry;
			(*c)[i - 1] = uncar & RADIXM;
			uncar = ((*c)[i] << 1) + (uncar >> NBITS);
			{
				zaddmulpsq0(&(*c)[i - 1], a[hal], &carry);
			}
			uncar += carry;
			carry = uncar >> NBITS;
			(*c)[i] = uncar & RADIXM;
		}
		while ((i > 1) && (!((*c)[i])))
			i--;
		(*c)[0] = i;
		return;
	}
	a0 = &(kar_mem[shi]);
	a1 = &(kar_mem[shi + 1]);
	a2 = &(kar_mem[shi + 2]);
	hal = ((al + 1) >> 1);
	zsetlength(a0, al + hal + 2, "in kar_sq, locals\n");
	zsetlength(a1, al + 2, "");
	zsetlength(a2, al, "");
	i = hal;
	while ((i > 1) && (!(a[i])))
		i--;
	a[0] = i;
	for (i = hal + 1; i <= al; i++)
		(*a0)[i - hal] = a[i];
	(*a0)[0] = al - hal;
	kar_sq0(a, a1, shi + 3);
	zadd(a, (*a0), a2);
	kar_sq0((*a0), c, shi + 3);
	a[0] = al;
	kar_sq0((*a2), a0, shi + 3);
	zsubpos((*a0), (*a1), a0);
	zsubpos((*a0), *c, a0);
	zlshift((*a0), hal * NBITS, a0);
	hal <<= 1;
	for (i = (*c)[0]; i; i--)
		(*c)[i + hal] = (*c)[i];
	for (i = hal; i > (*a1)[0]; i--)
		(*c)[i] = 0;
	for (; i; i--)
		(*c)[i] = (*a1)[i];
	(*c)[0] += hal;
	zadd(*c, (*a0), c);
}
/* SINGLE_MUL = 1, PLAIN = 0, KARAT = 0 */
/* warning: must first convert to 26 bit radix */


void
zmul1(
        verylong a,
        verylong b,
        verylong *c
        )
{       /* output not input */
        register long aneg;
        register long bneg;
        verylong olda;
        verylong oldb;

        if (ALLOCATE && (!a || !b))
        {
                zzero(c);
                return;
        }
        if (a == b)
        {
                zsq(a, c);
                return;
        }
        if (!kar_mem_initialized)
        {
                kar_mem_initialized = 1;
                for (aneg = (5 * KAR_DEPTH) - 1; aneg >= 0; aneg--)
                        kar_mem[aneg] = (verylong) 0;
        }
        olda = a;
        oldb = b;
        if (aneg = (*a < 0))
                a[0] = -a[0];
        if (bneg = (*b < 0))
                b[0] = -b[0];
        if (*a > *b)
                kar_mul1(a, b, c, (long) 0);
        else
                kar_mul1(b, a, c, (long) 0);
        if (aneg != bneg && ((*c)[1] || (*c)[0] != 1))
                (*c)[0] = -(*c)[0];
        if (aneg)
                olda[0] = -olda[0];
        if (bneg)
                oldb[0] = -oldb[0];
}


void
zsq1(
        verylong a,
        verylong *c
        )
{       /* output is not input */
        register long aneg;

        if (ALLOCATE && !a)
        {
                zzero(c);
                return;
        }
        if (!kar_mem_initialized)
        {
                kar_mem_initialized = 1;
                for (aneg = (5 * KAR_DEPTH) - 1; aneg >= 0; aneg--)
                        kar_mem[aneg] = (verylong) 0;
        }
        if (aneg = (*a < 0))
                a[0] = -a[0];
        kar_sq1(a, c, (long) 0);
        if (aneg)
                a[0] = -a[0];
}

#define LO1(x) (((unsigned long *) &x)[1])
#define HI1(x) (((unsigned long *) &x)[0])
static void zsubpos1();

typedef union { double d; unsigned long rep[2]; } d_or_rep;

#define FetchHiLo1(hi,lo,x) \
{ \
   d_or_rep _l_xx; \
   _l_xx.d = x; \
   hi = _l_xx.rep[0]; \
   lo = _l_xx.rep[1]; \
}

void zaddmul1(
        long ams,
        long *ama,
        long *amb
        )
{
        register long carry = 0;
        register long i = (*amb++);
        register double dams = (double) ams;
        register double xx;
        register double yy;
        register unsigned long lo_wd, lo;
        register unsigned long hi_wd, hi;

        xx  =  ((double) (*amb))*dams + 4503599627370496.0;
        for (; i > 1; i--)
        {
                yy =  ((double) (*(++amb)))*dams +4503599627370496.0;
                FetchHiLo1(hi_wd, lo_wd, xx);
                lo = lo_wd & 0x3FFFFFF;
                hi = ((hi_wd<<6)|(lo_wd>>26)) & 0x3FFFFFF;
                lo = lo + (*ama) + carry;
                *ama = lo & 0x3FFFFFF;
                carry = hi + (lo >> 26);
                ama++;
                xx = yy;
        }

        FetchHiLo1(hi_wd, lo_wd, xx);
        lo = lo_wd & 0x3FFFFFF;
        hi = ((hi_wd<<6)|(lo_wd>>26)) & 0x3FFFFFF;
        lo = lo + (*ama) + carry;
        *ama = lo & 0x3FFFFFF;
        carry = hi + (lo >> 26);
        ama++;
        *ama += carry;
}


static void
zaddls1(
        verylong a,
        verylong b,
        verylong c
        )
{
/* low level c = a+b, b shorter than a, output can be input */
        register long sa = *a;
        register long sb = *b;
        register long carry = 0;
        register long i;
        register verylong pc;

 /* we know that sa and sb are both >0 or both <0 */
        pc = &c[0];
        if (sa < 0)
        {
                sa = -sa;
                sb = -sb;
        }
        for (i = 1; i <= sb; i++)
        {
                if ((*(++pc) = (*(++a)) + (*(++b)) + carry) < (1L<<26))
                        carry = 0;
                else
                {
                        *pc -= (1L<<26);
                        carry = 1;
                }
        }
        for (; i <= sa; i++)
        {
                if ((*(++pc) = (*(++a)) + carry) < (1L<<26))
                        carry = 0;
                else
                {
                        *pc -= (1L<<26);
                        carry = 1;
                }
        }
        if (carry)
        {
                c[0] = sa + 1;
                *(++pc) = 1;
        }
        else
                c[0] = sa;
}

void
z2mul1(
        verylong n,
        verylong *rres
        )
{
        register long sn;
        register long i;
        register long carry = 0;
        verylong res = *rres;

        if (ALLOCATE && !n)
        {
                zzero(rres);
                return;
        }
        if ((!n[1]) && (n[0] == 1))
        {
                zzero(rres);
                return;
        }
        if ((sn = n[0]) < 0)
                sn = -sn;
        zsetlength(&res, sn + 1, "in z2mul, second argument");
        if (n == *rres) n = res;
        *rres = res;
        for (i = 1; i <= sn; i++)
        {
                if ((res[i] = (n[i] << 1) + carry) >= (1L<<26))
                {
                        res[i] -= (1L<<26);
                        carry = 1;
                }
                else
                        carry = 0;
        }
        if (carry)
                res[++sn] = 1;
        if (n[0] < 0)
                res[0] = -sn;
        else
                res[0] = sn;
}

void
zlshift1(
	verylong n,
	long k,
	verylong *rres
	)
{
	register long big;
	register long small;
	register long sn;
	register long i;
	register long cosmall;
	verylong res = *rres;

	if (ALLOCATE && !n)
	{
		zzero(rres);
		return;
	}
	if ((!n[1]) && (n[0] == 1))
	{
		zzero(rres);
		return;
	}
	if (!k)
	{
		if (n != *rres)
			zcopy(n, rres);
		return;
	}
	if (k < 0)
	{
		zrshift(n, -k, rres);
		return;
	}
	if (k == 1)
	{
		z2mul1(n, rres);
		return;
	}
	if ((sn = n[0]) < 0)
		sn = -sn;
	i = sn + (big = k / 26);
	if (small = k - big * 26)
	{
		zsetlength(&res, i + 1, "in zlshift, third argument");
		if (n == *rres) n = res;
		*rres = res;
		res[i + 1] = n[sn] >> (cosmall = 26 - small);
		for (i = sn; i > 1; i--)
			res[i + big] = ((n[i] << small) & ((1L<<26)-1)) + (n[i - 1] >> cosmall);
		res[big + 1] = (n[1] << small) & ((1L<<26)-1);
		for (i = big; i; i--)
			res[i] = 0;
		if (res[sn + big + 1])
			big++;
	}
	else
	{
		zsetlength(&res, i, "in zlshift, third argument");
		if (n == *rres) n = res;
		*rres = res;
		for (i = sn; i; i--)
			res[i + big] = n[i];
		for (i = big; i; i--)
			res[i] = 0;
	}
	if (n[0] > 0)
		res[0] = n[0] + big;
	else
		res[0] = n[0] - big;
}


void
zadd1(
        verylong a,
        verylong b,
        verylong *cc
        )
{
        register long sa;
        register long sb;
        register long anegative;
        verylong c = *cc;

        if (ALLOCATE && !a)
        {
                if (b)
                        zcopy(b, cc);
                else
                        zzero(cc);
                return;
        }
        if (ALLOCATE && !b)
        {
                zcopy(a, cc);
                return;
        }
        if ((anegative = ((sa = a[0]) < 0)) == ((sb = b[0]) < 0))
        {
        /* signs a and b are the same */
                register long i;
                if (anegative)
                {
                        sa = -sa;
                        sb = -sb;
                }
                zsetlength(&c, (sa > sb ? sa : sb) + 1, "in zadd, third argument");
                if (a == *cc) a = c;
                if (b == *cc) b = c;
                *cc = c;
                if (sa == sb)
                {
                        register verylong pc;
                        pc = &c[0];
                        i = 0;
                        for (; sa; sa--)
                        {
                                if ((*(++pc) = (*(++a)) + (*(++b)) + i) < (1L<<26))i = 0;
                                else
                                {
                                        *pc -= (1L<<26);
                                        i = 1;
                                }
                        }
                        if (i)
                        {
                                c[0] = sb + 1;
                                *(++pc) = 1;
                        }
                        else
                                c[0] = sb;
                }
                else if (sa > sb)
                        zaddls1(a, b, c);
                else
                        zaddls1(b, a, c);
                if (anegative)
                        c[0] = -c[0];
        /* if anegative, then c cannot be zero */
        }
        else
        {
        /* signs a and b are different */
                verylong old;
                verylong oldc;

                oldc = c;
                if (anegative)
                {
                        a[0] = -a[0];
                        old = a;
                }
                else
                {
                        b[0] = -b[0];
                        old = b;
                }
        /* the one that's negative cannot be zero */
                if (!(sa = zcompare(a, b)))
                {
                        zzero(&c);
                        *cc = c;
                        if (anegative)
                        {
                                if (old != oldc)
                                        a[0] = -a[0];
                        }
                        else if (old != oldc)
                                b[0] = -b[0];
                }
                else if (sa > 0)
                {
                        zsubpos1(a, b, &c);
                        *cc = c;
                        if (anegative)
                        {
                                c[0] = -c[0];
                                if (old != oldc)
                                        a[0] = -a[0];
                        }
                        else if (old != oldc)
                                b[0] = -b[0];
                }
                else
                {
                        zsubpos1(b, a, &c);
                        *cc = c;
                        if (!anegative)
                        {
                                c[0] = -c[0];
                                if (old != oldc)
                                        b[0] = -b[0];
                        }
                        else if (old != oldc)
                                a[0] = -a[0];
                }
        }
}


void
zsub1(
	verylong a,
	verylong b,
	verylong *cc
	)
{
	register long sa;
	register long sb;
	register long anegative;
	verylong c = *cc;

	if (ALLOCATE && !b)
	{
		if (a)
			zcopy(a, cc);
		else
			zzero(cc);
		return;
	}
	if (ALLOCATE && !a)
	{
		zcopy(b, cc);
		znegate(cc);
		return;
	}
	if ((anegative = ((sa = a[0]) < 0)) == ((sb = b[0]) < 0))
	{
	/* signs a and b are the same */
		register long carry = 0;
		register long i;
		register verylong pc;
		register long agrb;

		if (!(agrb = zcompare(a, b)))
		{
			zzero(cc);
			return;
		}
		if ((agrb > 0 && anegative) || (agrb < 0 && !anegative))
		{
			pc = a;
			a = b;
			b = pc;
			sa = *a;
			sb = *b;
		}
		if (anegative)
		{
			sa = -sa;
			sb = -sb;
		}
		zsetlength(&c, sa, "in zsub, third argument");
                if (b == *cc) b = c;
		*cc = c;
		pc = &c[0];
		for (i = 1; i <= sb; i++)
		{
			if ((*(++pc) = (*(++a)) - (*(++b)) - carry) >= 0)
				carry = 0;
			else
			{
				*pc += (1L<<26);
				carry = 1;
			};
		}
		for (; i <= sa; i++)
		{
			if ((*(++pc) = (*(++a)) - carry) >= 0)
				carry = 0;
			else
			{
				*pc += (1L<<26);
				carry = 1;
			};
		}
		i = sa;
		while ((i > 1) && (!(*pc)))
		{
			i--;
			pc--;
		}
		if (agrb > 0)
			c[0] = i;
		else
			c[0] = -i;
	}
	else
	{
	/* signs a and b are different */
		verylong old;
		verylong oldc;

		oldc = c;
		if (anegative)
		{
			a[0] = -a[0];
			old = a;
		}
		else
		{
			b[0] = -b[0];
			old = b;
		}
	/* the one that's negative cannot be zero */
		zadd1(a, b, &c);
		*cc = c;
		if (anegative)
		{
			c[0] = -c[0];
			if (old != oldc)
				a[0] = -a[0];
		}
		else if (old != oldc)
			b[0] = -b[0];
	}
}



static void
zsubpos1(
        verylong a,
        verylong b,
        verylong *cc
        )
{
        register long sa = a[0];
        register long sb = b[0];
        register long carry = 0;
        register long i;
        register verylong pc;
        verylong c = *cc;

        if (ALLOCATE && !b)
        {
                if (a)
                        zcopy(a, cc);
                else
                        zzero(cc);
                return;
        }
        if (ALLOCATE && !a)
        {
                zzero(cc);
                return;
        }
        zsetlength(&c, sa, "in zsubpos, third argument");
        /* if *cc == a, then nothing will happen */
        /* if *cc == b, then b might point to freed space, so */
        if (b == *cc) b = c;
        *cc = c;
        pc = &c[0];
        for (i = 1; i <= sb; i++)
        {
                if ((*(++pc) = (*(++a)) - (*(++b)) - carry) >= 0)
                        carry = 0;
                else
                {
                        *pc += (1L<<26);
                        carry = 1;
                };
        }
        for (; i <= sa; i++)
        {
                if ((*(++pc) = (*(++a)) - carry) >= 0)
                        carry = 0;
                else
                {
                        *pc += (1L<<26);
                        carry = 1;
                };
        }
        i = sa;
        while ((i > 1) && (!(*pc)))
        {
                i--;
                pc--;
        }
        c[0] = i;
}

kar_mul1(
	verylong a,
	verylong b,
	verylong *c,
	long shi
	)
{
	register long al;
	register long hal;
	register long i;
	register long restoreb0 = b[0];
	register verylong pc;
	register long bbig = 1;
	verylong *a0;
	verylong *a1;
	verylong *a2;
	verylong *a3;
	verylong *a4;

	zsetlength(c, (hal = (al = a[0]) + (i = b[0])), "in kar_mul, third argument");
	if ((shi >= (5 * KAR_DEPTH)) || (al < K_M_C) || (i < K_M_C))
	{
		pc = &(*c)[1];
		for (i = hal; i > 0; i--)
			*pc++ = 0;
		pc = &(*c)[1];
		if (al <= *b)
			for (i = al; i; i--)
			{
				zaddmul1(*(++a), pc++, b);
			}
		else
			for (i = *b; i; i--)
			{
				zaddmul1(*(++b), pc++, a);
			}
		while ((hal > 1) && (!((*c)[hal])))
			hal--;
		(*c)[0] = hal;
		return;
	}
	a0 = &(kar_mem[shi]);
	a1 = &(kar_mem[shi + 1]);
	a2 = &(kar_mem[shi + 2]);
	a3 = &(kar_mem[shi + 3]);
	a4 = &(kar_mem[shi + 4]);
	hal = ((al + 1) >> 1);
	zsetlength(a0, al, "in kar_mul, locals\n");
	zsetlength(a1, al, "");
	zsetlength(a2, al, "");
	zsetlength(a3, al + hal, "");
	zsetlength(a4, al + 2, "");
	i = hal;
	while ((i > 1) && (!(a[i])))
		i--;
	a[0] = i;
	if (hal >= b[0])
		bbig = 0;
	else
	{
		i = hal;
		while ((i > 1) && (!(b[i])))
			i--;
		b[0] = i;
	}
	for (i = hal + 1; i <= al; i++)
		(*a1)[i - hal] = a[i];
	(*a1)[0] = al - hal;
	if (bbig)
	{
		for (i = hal + 1; i <= restoreb0; i++)
			(*a3)[i - hal] = b[i];
		(*a3)[0] = restoreb0 - hal;
	}
	kar_mul1(a, b, a4, shi + 5);
	zadd1(a, (*a1), a0);
	a[0] = al;
	if (bbig)
	{
		kar_mul1((*a1), (*a3), c, shi + 5);
		zadd1(b, (*a3), a2);
		b[0] = restoreb0;
		kar_mul1((*a0), (*a2), a3, shi + 5);
	}
	else
		kar_mul1((*a0), b, a3, shi + 5);
	zsubpos1((*a3), (*a4), a3);
	if (bbig)
		zsubpos1((*a3), *c, a3);
	zlshift1((*a3), hal * 26, a3);
	hal <<= 1;
	if (bbig)
	{
		for (i = (*c)[0]; i; i--)
			(*c)[i + hal] = (*c)[i];
		for (i = hal; i > (*a4)[0]; i--)
			(*c)[i] = 0;
		for (; i; i--)
			(*c)[i] = (*a4)[i];
		(*c)[0] += hal;
	}
	else
	{
		for (i = (*a4)[0]; i >= 0; i--)
			(*c)[i] = (*a4)[i];
	}
	zadd1(*c, (*a3), c);
}


#define ExtractHiLo1(Hi,Lo,x) \
{ \
double t=x+(4503599627370496.0); \
unsigned long *it = (unsigned long *)&t; \
Lo = it[1]; \
Hi = ((it[0]<<6)|(it[1]>>26)); \
Lo &= 0x3FFFFFF; \
Hi &= 0x3FFFFFF; \
}

#define zaddmulpsq1(a,b,t) \
{ \
   double __lx = ((double) (*a)) + ((double) b)*((double) b); \
   register long  __lhi = 0, __llo = 0;\
   ExtractHiLo1(__lhi,__llo,__lx);\
   (*a) = __llo;\
   (*t) = __lhi;\
}


void zaddmulsq1(
        long ams,
        long *ama,
        long *amb
        )
{
        register long carry = 0;
        register long i = ams;
        register double dams = (double) (*amb);
        double xx;
        register double yy;
        register unsigned long lo;
        register unsigned long hi;

        xx =  ((double) (*(++amb)))*dams + 4503599627370496.0;
        for (; i > 1; i--)
        {
                yy =  ((double) (*(++amb)))*dams + 4503599627370496.0;
                lo = LO1(xx) & 0x3FFFFFF;
                hi = ((HI1(xx)<<6)|(LO1(xx)>>26)) & 0x3FFFFFF;
                lo = lo + (*ama) + carry;
                *ama = lo & 0x3FFFFFF;
                carry = hi + (lo >> 26);
                ama++;
                xx = yy;
        }
        if (i==1)
        {
                lo = LO1(xx) & 0x3FFFFFF;
                hi = ((HI1(xx)<<6)|(LO1(xx)>>26)) & 0x3FFFFFF;
                lo = lo + (*ama) + carry;
                *ama = lo & 0x3FFFFFF;
                carry = hi + (lo >> 26);
                ama++;
        }
        *ama += carry;
}


kar_sq1(
	verylong a,
	verylong *c,
	long shi
	)
{
	register long al;
	register long hal;
	register long i;
	register verylong pc;
	verylong *a0;
	verylong *a1;
	verylong *a2;

	zsetlength(c, (i = ((al = a[0]) << 1)), "in kar_sq, second argument");
	if ((shi >= (3 * KAR_DEPTH)) || (al < K_S_C))
	{
		register unsigned long uncar;
		long carry = 0;
		pc = &(*c)[1];
		for (; i > 0; i--)
			*pc++ = 0;
		for (hal = 1; hal <= al; hal++)
		{
			i += 2;
			{
				zaddmulsq1(al - hal, &((*c)[i]), &(a[hal]));
			}
			uncar = ((*c)[i - 1] << 1) + carry;
			(*c)[i - 1] = uncar & ((1L<<26)-1);
			uncar = ((*c)[i] << 1) + (uncar >> 26);
			{
				zaddmulpsq1(&(*c)[i - 1], a[hal], &carry);
			}
			uncar += carry;
			carry = uncar >> 26;
			(*c)[i] = uncar & ((1L<<26)-1);
		}
		while ((i > 1) && (!((*c)[i])))
			i--;
		(*c)[0] = i;
		return;
	}
	a0 = &(kar_mem[shi]);
	a1 = &(kar_mem[shi + 1]);
	a2 = &(kar_mem[shi + 2]);
	hal = ((al + 1) >> 1);
	zsetlength(a0, al + hal + 2, "in kar_sq, locals\n");
	zsetlength(a1, al + 2, "");
	zsetlength(a2, al, "");
	i = hal;
	while ((i > 1) && (!(a[i])))
		i--;
	a[0] = i;
	for (i = hal + 1; i <= al; i++)
		(*a0)[i - hal] = a[i];
	(*a0)[0] = al - hal;
	kar_sq1(a, a1, shi + 3);
	zadd(a, (*a0), a2);
	kar_sq1((*a0), c, shi + 3);
	a[0] = al;
	kar_sq1((*a2), a0, shi + 3);
	zsubpos1((*a0), (*a1), a0);
	zsubpos1((*a0), *c, a0);
	zlshift1((*a0), hal * 26, a0);
	hal <<= 1;
	for (i = (*c)[0]; i; i--)
		(*c)[i + hal] = (*c)[i];
	for (i = hal; i > (*a1)[0]; i--)
		(*c)[i] = 0;
	for (; i; i--)
		(*c)[i] = (*a1)[i];
	(*c)[0] += hal;
	zadd(*c, (*a0), c);
}


/* SINGLE_MUL = 0, PLAIN = 1, KARAT = 0 */

void
zmul2(
        verylong a,
        verylong b,
        verylong *c
        )
{       /* output not input */
        register long aneg;
        register long bneg;
        verylong olda;
        verylong oldb;

        if (ALLOCATE && (!a || !b))
        {
                zzero(c);
                return;
        }
        if (a == b)
        {
                zsq(a, c);
                return;
        }
        if (!kar_mem_initialized)
        {
                kar_mem_initialized = 1;
                for (aneg = (5 * KAR_DEPTH) - 1; aneg >= 0; aneg--)
                        kar_mem[aneg] = (verylong) 0;
        }
        olda = a;
        oldb = b;
        if (aneg = (*a < 0))
                a[0] = -a[0];
        if (bneg = (*b < 0))
                b[0] = -b[0];
        if (*a > *b)
                kar_mul2(a, b, c, (long) 0);
        else
                kar_mul2(b, a, c, (long) 0);
        if (aneg != bneg && ((*c)[1] || (*c)[0] != 1))
                (*c)[0] = -(*c)[0];
        if (aneg)
                olda[0] = -olda[0];
        if (bneg)
                oldb[0] = -oldb[0];
}

void
zsq2(
        verylong a,
        verylong *c
        )
{       /* output is not input */
        register long aneg;

        if (ALLOCATE && !a)
        {
                zzero(c);
                return;
        }
        if (!kar_mem_initialized)
        {
                kar_mem_initialized = 1;
                for (aneg = (5 * KAR_DEPTH) - 1; aneg >= 0; aneg--)
                        kar_mem[aneg] = (verylong) 0;
        }
        if (aneg = (*a < 0))
                a[0] = -a[0];
        kar_sq2(a, c, (long) 0);
        if (aneg)
                a[0] = -a[0];
}


#  define zaddmulp2(_a, _b, _d, _t) \
{ \
        register long lb = (_b); \
        register long ld = (_d); \
        register long b1 = (_b) & RADIXROOTM; \
        register long d1 = (_d) & RADIXROOTM; \
        register long aa = *(_a) + b1 * d1; \
 \
        b1 = b1 * (ld >>= NBITSH) + d1 * (lb >>= NBITSH) + (aa >> NBITSH); \
        aa = (aa & RADIXROOTM) + ((b1 & RADIXROOTM) << NBITSH) + *(_t); \
        *(_t) = ld * lb + (b1 >> NBITSH) + (aa >> NBITS); \
        *(_a) = (aa & RADIXM); \
}


#define zaddmul2(ams, ama, amb) \
{ \
        register long lami; \
        register long lams = (ams); \
        register verylong lama = (ama); \
        register verylong lamb = (amb); \
        long lamcarry = 0; \
 \
        for (lami = (*lamb++); lami > 0; lami--) \
        { \
                zaddmulp2(lama, *lamb, lams, &lamcarry); \
        /* Be careful, the last lama is unnormalized */ \
                lama++; \
                lamb++; \
        } \
        *lama += lamcarry; \
}

kar_mul2(
	verylong a,
	verylong b,
	verylong *c,
	long shi
	)
{
	register long al;
	register long hal;
	register long i;
	register long restoreb0 = b[0];
	register verylong pc;
	register long bbig = 1;
	verylong *a0;
	verylong *a1;
	verylong *a2;
	verylong *a3;
	verylong *a4;

	zsetlength(c, (hal = (al = a[0]) + (i = b[0])), "in kar_mul, third argument");
	if ((shi >= (5 * KAR_DEPTH)) || (al < K_M_C) || (i < K_M_C))
	{
		pc = &(*c)[1];
		for (i = hal; i > 0; i--)
			*pc++ = 0;
		pc = &(*c)[1];
		if (al <= *b)
			for (i = al; i; i--)
			{
				zaddmul2(*(++a), pc++, b);
			}
		else
			for (i = *b; i; i--)
			{
				zaddmul2(*(++b), pc++, a);
			}
		while ((hal > 1) && (!((*c)[hal])))
			hal--;
		(*c)[0] = hal;
		return;
	}
	a0 = &(kar_mem[shi]);
	a1 = &(kar_mem[shi + 1]);
	a2 = &(kar_mem[shi + 2]);
	a3 = &(kar_mem[shi + 3]);
	a4 = &(kar_mem[shi + 4]);
	hal = ((al + 1) >> 1);
	zsetlength(a0, al, "in kar_mul, locals\n");
	zsetlength(a1, al, "");
	zsetlength(a2, al, "");
	zsetlength(a3, al + hal, "");
	zsetlength(a4, al + 2, "");
	i = hal;
	while ((i > 1) && (!(a[i])))
		i--;
	a[0] = i;
	if (hal >= b[0])
		bbig = 0;
	else
	{
		i = hal;
		while ((i > 1) && (!(b[i])))
			i--;
		b[0] = i;
	}
	for (i = hal + 1; i <= al; i++)
		(*a1)[i - hal] = a[i];
	(*a1)[0] = al - hal;
	if (bbig)
	{
		for (i = hal + 1; i <= restoreb0; i++)
			(*a3)[i - hal] = b[i];
		(*a3)[0] = restoreb0 - hal;
	}
	kar_mul2(a, b, a4, shi + 5);
	zadd(a, (*a1), a0);
	a[0] = al;
	if (bbig)
	{
		kar_mul2((*a1), (*a3), c, shi + 5);
		zadd(b, (*a3), a2);
		b[0] = restoreb0;
		kar_mul2((*a0), (*a2), a3, shi + 5);
	}
	else
		kar_mul2((*a0), b, a3, shi + 5);
	zsubpos((*a3), (*a4), a3);
	if (bbig)
		zsubpos((*a3), *c, a3);
	zlshift((*a3), hal * NBITS, a3);
	hal <<= 1;
	if (bbig)
	{
		for (i = (*c)[0]; i; i--)
			(*c)[i + hal] = (*c)[i];
		for (i = hal; i > (*a4)[0]; i--)
			(*c)[i] = 0;
		for (; i; i--)
			(*c)[i] = (*a4)[i];
		(*c)[0] += hal;
	}
	else
	{
		for (i = (*a4)[0]; i >= 0; i--)
			(*c)[i] = (*a4)[i];
	}
	zadd(*c, (*a3), c);
}


#define zaddmulsq2(sql, sqa, sqb) \
{ \
        register long lsqi = (sql); \
        register long lsqs = *(sqb); \
        register verylong lsqa = (sqa); \
        register verylong lsqb = (sqb); \
        long lsqcarry = 0; \
 \
        lsqb++; \
        for (; lsqi > 0; lsqi--) \
        { \
                zaddmulp2(lsqa, *lsqb, lsqs, &lsqcarry); \
                lsqa++; \
                lsqb++; \
        } \
        *lsqa += lsqcarry; \
/* Be careful, the last lama is unnormalized */ \
}


#define zaddmulpsq2(_a, _b, _t) \
{ \
        register long lb = (_b); \
        register long b1 = (_b) & RADIXROOTM; \
        register long aa = *(_a) + b1 * b1; \
 \
        b1 = (b1 * (lb >>= NBITSH) << 1) + (aa >> NBITSH); \
        aa = (aa & RADIXROOTM) + ((b1 & RADIXROOTM) << NBITSH); \
        *(_t) = lb * lb + (b1 >> NBITSH) + (aa >> NBITS); \
        *(_a) = (aa & RADIXM); \
}



kar_sq2(
	verylong a,
	verylong *c,
	long shi
	)
{
	register long al;
	register long hal;
	register long i;
	register verylong pc;
	verylong *a0;
	verylong *a1;
	verylong *a2;

	zsetlength(c, (i = ((al = a[0]) << 1)), "in kar_sq, second argument");
	if ((shi >= (3 * KAR_DEPTH)) || (al < K_S_C))
	{
		register unsigned long uncar;
		long carry = 0;
		pc = &(*c)[1];
		for (; i > 0; i--)
			*pc++ = 0;
		for (hal = 1; hal <= al; hal++)
		{
			i += 2;
			{
				zaddmulsq2(al - hal, &((*c)[i]), &(a[hal]));
			}
			uncar = ((*c)[i - 1] << 1) + carry;
			(*c)[i - 1] = uncar & RADIXM;
			uncar = ((*c)[i] << 1) + (uncar >> NBITS);
			{
				zaddmulpsq2(&(*c)[i - 1], a[hal], &carry);
			}
			uncar += carry;
			carry = uncar >> NBITS;
			(*c)[i] = uncar & RADIXM;
		}
		while ((i > 1) && (!((*c)[i])))
			i--;
		(*c)[0] = i;
		return;
	}
	a0 = &(kar_mem[shi]);
	a1 = &(kar_mem[shi + 1]);
	a2 = &(kar_mem[shi + 2]);
	hal = ((al + 1) >> 1);
	zsetlength(a0, al + hal + 2, "in kar_sq, locals\n");
	zsetlength(a1, al + 2, "");
	zsetlength(a2, al, "");
	i = hal;
	while ((i > 1) && (!(a[i])))
		i--;
	a[0] = i;
	for (i = hal + 1; i <= al; i++)
		(*a0)[i - hal] = a[i];
	(*a0)[0] = al - hal;
	kar_sq2(a, a1, shi + 3);
	zadd(a, (*a0), a2);
	kar_sq2((*a0), c, shi + 3);
	a[0] = al;
	kar_sq2((*a2), a0, shi + 3);
	zsubpos((*a0), (*a1), a0);
	zsubpos((*a0), *c, a0);
	zlshift((*a0), hal * NBITS, a0);
	hal <<= 1;
	for (i = (*c)[0]; i; i--)
		(*c)[i + hal] = (*c)[i];
	for (i = hal; i > (*a1)[0]; i--)
		(*c)[i] = 0;
	for (; i; i--)
		(*c)[i] = (*a1)[i];
	(*c)[0] += hal;
	zadd(*c, (*a0), c);
}


#define lower_radix1(in,out) { \
register long i=1,j=0,rb=0,nrb=0,bits; \
register verylong pin=in; \
long size; \
 \
zsetlength( (&out), (long) (30.0/26.0 * in[0] + 1),""); \
if ((size= *pin)==2) { \
        pin++; \
        out[1] = *pin & ((1L<<26)-1); \
        out[2] = *pin >> 26; \
        pin++; \
        out[2] |= (*pin & ((1L << (2*26-30))-1L)) << (30 - \
         26); \
        bits = *pin >> (2*26-30); \
        if (bits) { \
                out[3] = bits; \
                out[0]=3; \
                } \
        else \
                out[0]=2; \
        } \
else { \
        for (j=0; i <= size ; j++) { \
                if (nrb >= 26) { \
                        bits = rb & ((1L<<26)-1); \
                        rb >>= 26; \
                        nrb -= 26; \
                        } \
                else { \
                        bits = rb | ((((1L << (26 - nrb)) -1L) & \
                         pin[i]) << nrb); \
                        rb = pin[i] >> (26 - nrb); \
                        nrb += (30-26); \
                        i++; \
                        } \
                out[j+1]=bits; \
                } \
        if (rb) { \
                j++; \
                if (rb < (1L<<26)) \
                        out[j]= rb; \
                else { \
                        out[j]= rb & (( 1L << 26) -1); \
                        j++; \
                        out[j] = rb>>26; \
                        } \
                } \
        out[0]=j; \
        } \
}


void
zmul5(
        verylong a,
        verylong b,
        verylong *c
        )
{       /* output not input */
        register long aneg;
        register long bneg;
        verylong olda;
        verylong oldb;

        if (ALLOCATE && (!a || !b))
        {
                zzero(c);
                return;
        }
        if (a == b)
        {
                zsq(a, c);
                return;
        }
        if (!kar_mem_initialized)
        {
                kar_mem_initialized = 1;
                for (aneg = (5 * KAR_DEPTH) - 1; aneg >= 0; aneg--)
                        kar_mem[aneg] = (verylong) 0;
        }
        olda = a;
        oldb = b;
        if (aneg = (*a < 0))
                a[0] = -a[0];
        if (bneg = (*b < 0))
                b[0] = -b[0];
        if (*a > *b)
                kar_mul5(a, b, c, (long) 0);
        else
                kar_mul5(b, a, c, (long) 0);
        if (aneg != bneg && ((*c)[1] || (*c)[0] != 1))
                (*c)[0] = -(*c)[0];
        if (aneg)
                olda[0] = -olda[0];
        if (bneg)
                oldb[0] = -oldb[0];
}


void
zsq5(
        verylong a,
        verylong *c
        )
{       /* output is not input */
        register long aneg;

        if (ALLOCATE && !a)
        {
                zzero(c);
                return;
        }
        if (!kar_mem_initialized)
        {
                kar_mem_initialized = 1;
                for (aneg = (5 * KAR_DEPTH) - 1; aneg >= 0; aneg--)
                        kar_mem[aneg] = (verylong) 0;
        }
        if (aneg = (*a < 0))
                a[0] = -a[0];
        kar_sq5(a, c, (long) 0);
        if (aneg)
                a[0] = -a[0];
}

#define LO5(x) (((unsigned long *) &x)[0])
#define HI5(x) (((unsigned long *) &x)[1])

#define FetchHiLo5(hi,lo,x) \
{ \
   d_or_rep _l_xx; \
   _l_xx.d = x; \
   hi = _l_xx.rep[1]; \
   lo = _l_xx.rep[0]; \
}

void zaddmul5(
        long ams,
        long *ama,
        long *amb
        )
{
        register long carry = 0;
        register long i = (*amb++);
        register double dams = (double) ams;
        register double xx;
        register double yy;
        register unsigned long lo_wd, lo;
        register unsigned long hi_wd, hi;

        xx  =  ((double) (*amb))*dams + 4503599627370496.0;
        for (; i > 1; i--)
        {
                yy =  ((double) (*(++amb)))*dams +4503599627370496.0;
                FetchHiLo5(hi_wd, lo_wd, xx);
                lo = lo_wd & 0x3FFFFFF;
                hi = ((hi_wd<<6)|(lo_wd>>26)) & 0x3FFFFFF;
                lo = lo + (*ama) + carry;
                *ama = lo & 0x3FFFFFF;
                carry = hi + (lo >> 26);
                ama++;
                xx = yy;
        }

        FetchHiLo5(hi_wd, lo_wd, xx);
        lo = lo_wd & 0x3FFFFFF;
        hi = ((hi_wd<<6)|(lo_wd>>26)) & 0x3FFFFFF;
        lo = lo + (*ama) + carry;
        *ama = lo & 0x3FFFFFF;
        carry = hi + (lo >> 26);
        ama++;
        *ama += carry;
}


kar_mul5(
	verylong a,
	verylong b,
	verylong *c,
	long shi
	)
{
	register long al;
	register long hal;
	register long i;
	register long restoreb0 = b[0];
	register verylong pc;
	register long bbig = 1;
	verylong *a0;
	verylong *a1;
	verylong *a2;
	verylong *a3;
	verylong *a4;

	zsetlength(c, (hal = (al = a[0]) + (i = b[0])), "in kar_mul, third argument");
	if ((shi >= (5 * KAR_DEPTH)) || (al < K_M_C) || (i < K_M_C))
	{
		pc = &(*c)[1];
		for (i = hal; i > 0; i--)
			*pc++ = 0;
		pc = &(*c)[1];
		if (al <= *b)
			for (i = al; i; i--)
			{
				zaddmul5(*(++a), pc++, b);
			}
		else
			for (i = *b; i; i--)
			{
				zaddmul5(*(++b), pc++, a);
			}
		while ((hal > 1) && (!((*c)[hal])))
			hal--;
		(*c)[0] = hal;
		return;
	}
	a0 = &(kar_mem[shi]);
	a1 = &(kar_mem[shi + 1]);
	a2 = &(kar_mem[shi + 2]);
	a3 = &(kar_mem[shi + 3]);
	a4 = &(kar_mem[shi + 4]);
	hal = ((al + 1) >> 1);
	zsetlength(a0, al, "in kar_mul, locals\n");
	zsetlength(a1, al, "");
	zsetlength(a2, al, "");
	zsetlength(a3, al + hal, "");
	zsetlength(a4, al + 2, "");
	i = hal;
	while ((i > 1) && (!(a[i])))
		i--;
	a[0] = i;
	if (hal >= b[0])
		bbig = 0;
	else
	{
		i = hal;
		while ((i > 1) && (!(b[i])))
			i--;
		b[0] = i;
	}
	for (i = hal + 1; i <= al; i++)
		(*a1)[i - hal] = a[i];
	(*a1)[0] = al - hal;
	if (bbig)
	{
		for (i = hal + 1; i <= restoreb0; i++)
			(*a3)[i - hal] = b[i];
		(*a3)[0] = restoreb0 - hal;
	}
	kar_mul5(a, b, a4, shi + 5);
	zadd(a, (*a1), a0);
	a[0] = al;
	if (bbig)
	{
		kar_mul5((*a1), (*a3), c, shi + 5);
		zadd(b, (*a3), a2);
		b[0] = restoreb0;
		kar_mul5((*a0), (*a2), a3, shi + 5);
	}
	else
		kar_mul5((*a0), b, a3, shi + 5);
	zsubpos((*a3), (*a4), a3);
	if (bbig)
		zsubpos((*a3), *c, a3);
	zlshift((*a3), hal * NBITS, a3);
	hal <<= 1;
	if (bbig)
	{
		for (i = (*c)[0]; i; i--)
			(*c)[i + hal] = (*c)[i];
		for (i = hal; i > (*a4)[0]; i--)
			(*c)[i] = 0;
		for (; i; i--)
			(*c)[i] = (*a4)[i];
		(*c)[0] += hal;
	}
	else
	{
		for (i = (*a4)[0]; i >= 0; i--)
			(*c)[i] = (*a4)[i];
	}
	zadd(*c, (*a3), c);
}


#define ExtractHiLo5(Hi,Lo,x) \
{ \
double t=x+(4503599627370496.0); \
unsigned long *it = (unsigned long *)&t; \
Lo = it[0]; \
Hi = ((it[1]<<6)|(it[0]>>26)); \
Lo &= 0x3FFFFFF; \
Hi &= 0x3FFFFFF; \
}


#define zaddmulpsq5(a,b,t) \
{ \
   double __lx = ((double) (*a)) + ((double) b)*((double) b); \
   register long  __lhi = 0, __llo = 0;\
   ExtractHiLo5(__lhi,__llo,__lx);\
   (*a) = __llo;\
   (*t) = __lhi;\
}


void zaddmulsq5(
        long ams,
        long *ama,
        long *amb
        )
{
        register long carry = 0;
        register long i = ams;
        register double dams = (double) (*amb);
        double xx;
        register double yy;
        register unsigned long lo;
        register unsigned long hi;

        xx =  ((double) (*(++amb)))*dams + 4503599627370496.0;
        for (; i > 1; i--)
        {
                yy =  ((double) (*(++amb)))*dams + 4503599627370496.0;
                lo = LO5(xx) & 0x3FFFFFF;
                hi = ((HI5(xx)<<6)|(LO5(xx)>>26)) & 0x3FFFFFF;
                lo = lo + (*ama) + carry;
                *ama = lo & 0x3FFFFFF;
                carry = hi + (lo >> 26);
                ama++;
                xx = yy;
        }
        if (i==1)
        {
                lo = LO5(xx) & 0x3FFFFFF;
                hi = ((HI5(xx)<<6)|(LO5(xx)>>26)) & 0x3FFFFFF;
                lo = lo + (*ama) + carry;
                *ama = lo & 0x3FFFFFF;
                carry = hi + (lo >> 26);
                ama++;
        }
        *ama += carry;
}


kar_sq5(
	verylong a,
	verylong *c,
	long shi
	)
{
	register long al;
	register long hal;
	register long i;
	register verylong pc;
	verylong *a0;
	verylong *a1;
	verylong *a2;

	zsetlength(c, (i = ((al = a[0]) << 1)), "in kar_sq, second argument");
	if ((shi >= (3 * KAR_DEPTH)) || (al < K_S_C))
	{
		register unsigned long uncar;
		long carry = 0;
		pc = &(*c)[1];
		for (; i > 0; i--)
			*pc++ = 0;
		for (hal = 1; hal <= al; hal++)
		{
			i += 2;
			{
				zaddmulsq5(al - hal, &((*c)[i]), &(a[hal]));
			}
			uncar = ((*c)[i - 1] << 1) + carry;
			(*c)[i - 1] = uncar & RADIXM;
			uncar = ((*c)[i] << 1) + (uncar >> NBITS);
			{
				zaddmulpsq5(&(*c)[i - 1], a[hal], &carry);
			}
			uncar += carry;
			carry = uncar >> NBITS;
			(*c)[i] = uncar & RADIXM;
		}
		while ((i > 1) && (!((*c)[i])))
			i--;
		(*c)[0] = i;
		return;
	}
	a0 = &(kar_mem[shi]);
	a1 = &(kar_mem[shi + 1]);
	a2 = &(kar_mem[shi + 2]);
	hal = ((al + 1) >> 1);
	zsetlength(a0, al + hal + 2, "in kar_sq, locals\n");
	zsetlength(a1, al + 2, "");
	zsetlength(a2, al, "");
	i = hal;
	while ((i > 1) && (!(a[i])))
		i--;
	a[0] = i;
	for (i = hal + 1; i <= al; i++)
		(*a0)[i - hal] = a[i];
	(*a0)[0] = al - hal;
	kar_sq5(a, a1, shi + 3);
	zadd(a, (*a0), a2);
	kar_sq5((*a0), c, shi + 3);
	a[0] = al;
	kar_sq5((*a2), a0, shi + 3);
	zsubpos((*a0), (*a1), a0);
	zsubpos((*a0), *c, a0);
	zlshift((*a0), hal * NBITS, a0);
	hal <<= 1;
	for (i = (*c)[0]; i; i--)
		(*c)[i + hal] = (*c)[i];
	for (i = hal; i > (*a1)[0]; i--)
		(*c)[i] = 0;
	for (; i; i--)
		(*c)[i] = (*a1)[i];
	(*c)[0] += hal;
	zadd(*c, (*a0), c);
}


static void
zsubmul0 (
        long r,
        verylong a,
        verylong b
        )
{
        register long rd = RADIX - r;
        register long i;
        long carry = RADIX;

        for (i = (*b++); i > 0; i--)
        {
                zaddmulp0(a, *b, rd, &carry);
                a++;
                carry += RADIXM - (*b++);
        }
        *a += carry - RADIX;    /* unnormalized */
}


void
zstart0()
{
        double local_one = (double) 1;
        double local_half = 1 / (double) 2;
        epsilon = local_one;
        fudge = local_one + epsilon;

 /* the following loop is sick, but we have to: */
 /* comparing one to one + epsilon does not have */
 /* the desired effect on some machines, because */
 /* one + epsilon is Not cast to double (as it should) */
        while (local_one != fudge)
        {
                epsilon = epsilon * local_half;
                fudge = local_one + epsilon;
        }
        epsilon += epsilon;
#ifndef ALPHA_OR_ALPHA50
        if ((epsilon * RADIX) > local_one)
        {
                fprintf(stderr,
                        "zstart failure: recompile with smaller NBITS\n")
;
                exit(0);
        }
#endif
        epsilon *= 3;
        fudge = fradix + epsilon * fradix;
#ifdef ALPHA
        fudge2 = fradix - epsilon * fradix;
#endif
#ifdef ALPHA50
        alpha50fudge = alpha50fradix + epsilon * alpha50fradix;
#endif
}

void
zdiv0(
	verylong in_a,
	verylong in_b,
	verylong *qqq,
	verylong *rrr
	)

{
	register long i;
	register long qq;
	long sa;
	long sb;
	long sq;
	verylong p;
	verylong pc;
	long sign;
	static verylong a = 0;
	static verylong b = 0;
	static verylong c = 0;
	static verylong d = 0;
	double btopinv;
	double aux;
	verylong q = *qqq;
	verylong r = *rrr;

#ifndef START
	if (fudge < 0)
		zstart0();
#endif

	if (ALLOCATE && !in_a)
	{
		zzero(qqq);
		zzero(rrr);
		return;
	}
	if ((!in_b) || (((sb = in_b[0]) == 1) && (!in_b[1])))
	{
		zhalt("division by zero in zdiv");
		return;
	}
	zcopy(in_a,&a);

	zcopy(in_b,&b);

	sign = (*a < 0 ? 2 : 0) | (*b < 0 ? 1 : 0);
	if (*a < 0)
		(*a) = (-(*a));
	sa = (*a);
	if (*b < 0)
		(*b) = (-(*b));
	sb = (*b);
	zsetlength(&c, (i = (sa > sb ? sa : sb) + 1), "in zdiv, locals\n");
	zsetlength(&d, i, "");
	zsetlength(&q, i, "in zdiv, third argument");
	*qqq = q;
	zsetlength(&r, sb + 1, "in zdiv, fourth argument");
	*rrr = r;
	p = &b[sb];
	if ((sb == 1) && (*p < RADIX))
		zintoz(zsdiv(a, *p, &q), &r);
	else
	{
		sq = sa - sb;	/* size of quotient */
		btopinv = (double) (*p) * fradix;
		if (sb > 1)
			btopinv += (*(--p));
		btopinv *= fradix;
		if (sb > 2)
			btopinv += (*(p - 1));
		btopinv = fudge / btopinv;
		p = &a[1];
		pc = &c[0];
		for (i = sa; i > 0; i--)
			*pc++ = *p++;
		p = pc;
		sa = 1 - sb;
		for (i = (-sq); i > 0; i--)
			*p++ = 0;
		*pc = 0;
		d[1] = 0;
		p = &d[sq + 1];
		for (i = sq; i >= 0; i--)
		{

			aux = fradix * (fradix * (*pc) + (*(pc - 1))) + 1.0;
			if (i > sa)
				aux += (*(pc - 2));
			qq = (long) (btopinv * aux + 0.5);
		/* dirty, but safe. On most machines (long)(btopinv*aux) */
		/* is correct, or one too big; on some however it becomes */
		/* too small. Could change zstart, but +0.5 and a while */
		/* instead of one single if is safer */
			if (qq > 0)
			{
				if (qq >= RADIX)
					qq = RADIX-1;
				zsubmul0(qq, &c[i], &b[0]);
				while (*pc < 0)
				{
					qq--;
					{
						zaddmulone(&c[i], &b[0]);
					}
				}
			}
			pc--;
			*p-- = qq;
		}
		sb--;
		while ((sb > 0) && (!(c[sb])))
			sb--;
		sb++;
		r[0] = sb;
		p = &r[1];
		pc = &c[0];
		for (i = sb; i > 0; i--)
			*p++ = *pc++;
		if (sq < 0)
		{
			q[0] = 1;
			q[1] = d[1];
		}
		else
		{
			sq++;
			while ((sq > 1) && (!(d[sq])))
				sq--;
			q[0] = sq;
			p = &q[1];
			pc = &d[1];
			for (i = sq; i > 0; i--)
				*p++ = *pc++;
		}
	}	/* non-trivial case */

	if (sign)
	{
		if (sign <= 2)
		{
			if (!(r[1]) && (r[0] == 1))
				q[0] = -(q[0]);
			else
			{
				zadd(q, one, &q);
				q[0] = -q[0];
				if (sign == 1)
					zsub(r, b, &r);
				else
					zsub(b, r, &r);
			}
		}
		else
			znegate(&r);
	}
	*qqq = q;
	*rrr = r;

}


#define zaddmulone1(ama, amb) \
{ \
        register long lami; \
        register long lams = 0; \
        register verylong lama = (ama); \
        register verylong lamb = (amb); \
 \
        lams = 0; \
        for (lami = (*lamb++); lami > 0; lami--) \
        { \
                lams += (*lama + *lamb++); \
                *lama++ = lams & ((1L<<26)-1); \
                lams >>= 26; \
        } \
        *lama += lams; \
}


static
void zsubmul1(
        long ams,
        long *ama,
        long *amb
        )
{
        register long carry = (1L<<26);
        register long i = (*amb++);
        register double dams = (double) ((1L<<26)-ams);
        double xx;
        register double yy;
        register unsigned long lo;
        register unsigned long hi;

        xx =  ((double) (*amb))*dams + 4503599627370496.0;
        for (; i > 1; i--)
        {
                yy =  ((double) (*(amb+1)))*dams + 4503599627370496.0;
                lo = LO1(xx) & 0x3FFFFFF;
                hi = ((HI1(xx)<<6)|(LO1(xx)>>26)) & 0x3FFFFFF;
                lo = lo + (*ama) + carry;
                *ama = lo & 0x3FFFFFF;
                carry = hi + (lo >> 26);
                ama++;
                carry += ((1L<<26)-1) - (*(amb++));
                xx = yy;
        }
        lo = LO1(xx) & 0x3FFFFFF;
        hi = ((HI1(xx)<<6)|(LO1(xx)>>26)) & 0x3FFFFFF;
        lo = lo + (*ama) + carry;
        *ama = lo & 0x3FFFFFF;
        carry = hi + (lo >> 26);
        ama++;
        carry += ((1L<<26)-1) - (*amb);
        *ama += carry - (1L<<26);
}


void
zstart1()
{
        double local_one = (double) 1;
        double local_half = 1 / (double) 2;
        epsilon = local_one;
        fudge1 = local_one + epsilon;

 /* the following loop is sick, but we have to: */
 /* comparing one to one + epsilon does not have */
 /* the desired effect on some machines, because */
 /* one + epsilon is Not cast to double (as it should) */
        while (local_one != fudge1)
        {
                epsilon = epsilon * local_half;
                fudge1 = local_one + epsilon;
        }
        epsilon += epsilon;
        if ((epsilon * (1L<<26)) > local_one)
        {
                fprintf(stderr,
                        "zstart failure: recompile with smaller NBITS\n");
                exit(0);
        }
        epsilon *= 3;
        fudge1 = 67108864.0 + epsilon * 67108864.0;
}


void
zdiv1(
	verylong in_a,
	verylong in_b,
	verylong *qqq,
	verylong *rrr
	)

{
	register long i;
	register long qq;
	long sa;
	long sb;
	long sq;
	verylong p;
	verylong pc;
	long sign;
	static verylong a = 0;
	static verylong b = 0;
	static verylong c = 0;
	static verylong d = 0;
	double btopinv;
	double aux;
	verylong q = *qqq;
	verylong r = *rrr;

#ifndef START
	if (fudge1 < 0)
		zstart1();
#endif

	if (ALLOCATE && !in_a)
	{
		zzero(qqq);
		zzero(rrr);
		return;
	}
	if ((!in_b) || (((sb = in_b[0]) == 1) && (!in_b[1])))
	{
		zhalt("division by zero in zdiv");
		return;
	}
	zcopy(in_a,&a);

	zcopy(in_b,&b);

	sign = (*a < 0 ? 2 : 0) | (*b < 0 ? 1 : 0);
	if (*a < 0)
		(*a) = (-(*a));
	sa = (*a);
	if (*b < 0)
		(*b) = (-(*b));
	sb = (*b);
	zsetlength(&c, (i = (sa > sb ? sa : sb) + 1), "in zdiv, locals\n");
	zsetlength(&d, i, "");
	zsetlength(&q, i, "in zdiv, third argument");
	*qqq = q;
	zsetlength(&r, sb + 1, "in zdiv, fourth argument");
	*rrr = r;
	p = &b[sb];
	if ((sb == 1) && (*p < (1L<<26)))
		zintoz(zsdiv(a, *p, &q), &r);
	else
	{
		sq = sa - sb;	/* size of quotient */
		btopinv = (double) (*p) * 67108864.0;
		if (sb > 1)
			btopinv += (*(--p));
		btopinv *= 67108864.0;
		if (sb > 2)
			btopinv += (*(p - 1));
		btopinv = fudge1 / btopinv;
		p = &a[1];
		pc = &c[0];
		for (i = sa; i > 0; i--)
			*pc++ = *p++;
		p = pc;
		sa = 1 - sb;
		for (i = (-sq); i > 0; i--)
			*p++ = 0;
		*pc = 0;
		d[1] = 0;
		p = &d[sq + 1];
		for (i = sq; i >= 0; i--)
		{

			aux = 67108864.0 * (67108864.0 * (*pc) + (*(pc - 1))) + 1.0;
			if (i > sa)
				aux += (*(pc - 2));
			qq = (long) (btopinv * aux + 0.5);
		/* dirty, but safe. On most machines (long)(btopinv*aux) */
		/* is correct, or one too big; on some however it becomes */
		/* too small. Could change zstart, but +0.5 and a while */
		/* instead of one single if is safer */
			if (qq > 0)
			{
				if (qq >= (1L<<26))
					qq = (1L<<26)-1;
				zsubmul1(qq, &c[i], &b[0]);
				while (*pc < 0)
				{
					qq--;
					{
						zaddmulone1(&c[i], &b[0]);
					}
				}
			}
			pc--;
			*p-- = qq;
		}
		sb--;
		while ((sb > 0) && (!(c[sb])))
			sb--;
		sb++;
		r[0] = sb;
		p = &r[1];
		pc = &c[0];
		for (i = sb; i > 0; i--)
			*p++ = *pc++;
		if (sq < 0)
		{
			q[0] = 1;
			q[1] = d[1];
		}
		else
		{
			sq++;
			while ((sq > 1) && (!(d[sq])))
				sq--;
			q[0] = sq;
			p = &q[1];
			pc = &d[1];
			for (i = sq; i > 0; i--)
				*p++ = *pc++;
		}
	}	/* non-trivial case */

	if (sign)
	{
		if (sign <= 2)
		{
			if (!(r[1]) && (r[0] == 1))
				q[0] = -(q[0]);
			else
			{
				zadd1(q, one, &q);
				q[0] = -q[0];
				if (sign == 1)
					zsub1(r, b, &r);
				else
					zsub1(b, r, &r);
			}
		}
		else
			znegate(&r);
	}
	*qqq = q;
	*rrr = r;

}


static void
zsubmul2(
        long r,
        verylong a,
        verylong b
        )
{
        register long rd = RADIX - r;
        register long i;
        long carry = RADIX;

        for (i = (*b++); i > 0; i--)
        {
                zaddmulp2(a, *b, rd, &carry);
                a++;
                carry += RADIXM - (*b++);
        }
        *a += carry - RADIX;    /* unnormalized */
}

#define correct2( q, x1, x2, x3, y1, y2, btopinv) { \
        register long ix1=x1,ix2=x2,ix3=x3,iy1=y1,iy2=y2; \
        long lprodlow=0,lprodhigh=0; \
\
        zaddmulp2(&lprodlow, q, iy2, &lprodhigh); \
        if ((ix3 -= lprodlow) < 0) { \
                ix3+= RADIX; \
                ix2--; \
        } \
        lprodlow=0; \
        zaddmulp2(&lprodlow, q, iy1, &lprodhigh); \
        if ((ix2 -= lprodlow) < 0) { \
                ix2 += RADIX; \
                --ix1; \
        } \
        if (ix1 < lprodhigh) q--; \
        else if ((ix1 -= lprodhigh)) { \
                q += (long) ((fradix * (fradix * ix1 + ix2))*btopinv ); \
        } \
        else { \
                q += (long) ((fradix * ix2 + ix3)*btopinv ); \
        } \
}

#ifndef ALPHA50
void
zdiv2(
	verylong in_a,
	verylong in_b,
	verylong *qqq,
	verylong *rrr
	)

{
	register long i;
	register long qq;
	long sa;
	long sb;
	long sq;
	verylong p;
	verylong pc;
	long sign;
	static verylong a = 0;
	static verylong b = 0;
	static verylong c = 0;
	static verylong d = 0;
	double btopinv;
	double aux;
	verylong q = *qqq;
	verylong r = *rrr;
#ifdef ALPHA
        double btopinv2;
#endif

#ifndef START
	if (fudge < 0)
		zstart0();
#endif

	if (ALLOCATE && !in_a)
	{
		zzero(qqq);
		zzero(rrr);
		return;
	}
	if ((!in_b) || (((sb = in_b[0]) == 1) && (!in_b[1])))
	{
		zhalt("division by zero in zdiv");
		return;
	}
	zcopy(in_a,&a);

	zcopy(in_b,&b);

	sign = (*a < 0 ? 2 : 0) | (*b < 0 ? 1 : 0);
	if (*a < 0)
		(*a) = (-(*a));
	sa = (*a);
	if (*b < 0)
		(*b) = (-(*b));
	sb = (*b);
	zsetlength(&c, (i = (sa > sb ? sa : sb) + 1), "in zdiv, locals\n");
	zsetlength(&d, i, "");
	zsetlength(&q, i, "in zdiv, third argument");
	*qqq = q;
	zsetlength(&r, sb + 1, "in zdiv, fourth argument");
	*rrr = r;
	p = &b[sb];
	if ((sb == 1) && (*p < RADIX))
		zintoz(zsdiv(a, *p, &q), &r);
	else
	{
		sq = sa - sb;	/* size of quotient */
		btopinv = (double) (*p) * fradix;
		if (sb > 1)
			btopinv += (*(--p));
		btopinv *= fradix;
#ifdef ALPHA
                btopinv2=btopinv;
                btopinv2 = fudge2/btopinv2;
#else
                if (sb > 2)
                        btopinv += (*(p - 1));
#endif
		btopinv = fudge / btopinv;
		p = &a[1];
		pc = &c[0];
		for (i = sa; i > 0; i--)
			*pc++ = *p++;
		p = pc;
		sa = 1 - sb;
		for (i = (-sq); i > 0; i--)
			*p++ = 0;
		*pc = 0;
		d[1] = 0;
		p = &d[sq + 1];
		for (i = sq; i >= 0; i--)
		{

			aux = fradix * (fradix * (*pc) + (*(pc - 1))) + 1.0;
#ifndef ALPHA
                        if (i > sa)
                                aux += (*(pc - 2));
#endif
#ifdef ALPHA
                        qq = (long) (btopinv2 * aux + 0.5);
#else
                        qq = (long) (btopinv * aux + 0.5);
                /* dirty, but safe. On most machines (long)(btopinv*aux) */
                /* is correct, or one too big; on some however it becomes */
                /* too small. Could change zstart, but +0.5 and a while */
                /* instead of one single if is safer */
#endif
			if (qq > 0)
			{
#ifdef ALPHA
                                if (qq > (1L<<48)) {
                                        correct2(qq,*pc,*(pc-1),(i>sa) ? *(pc-2):0,b[sb],b[sb-1],btopinv);
                                        if (qq >= RADIX)
                                                qq = RADIX-1;
                                }
#else
                                if (qq >= RADIX)
                                        qq = RADIX-1;
#endif
				zsubmul2(qq, &c[i], &b[0]);
				while (*pc < 0)
				{
					qq--;
					{
						zaddmulone(&c[i], &b[0]);
					}
				}
			}
			pc--;
			*p-- = qq;
		}
		sb--;
		while ((sb > 0) && (!(c[sb])))
			sb--;
		sb++;
		r[0] = sb;
		p = &r[1];
		pc = &c[0];
		for (i = sb; i > 0; i--)
			*p++ = *pc++;
		if (sq < 0)
		{
			q[0] = 1;
			q[1] = d[1];
		}
		else
		{
			sq++;
			while ((sq > 1) && (!(d[sq])))
				sq--;
			q[0] = sq;
			p = &q[1];
			pc = &d[1];
			for (i = sq; i > 0; i--)
				*p++ = *pc++;
		}
	}	/* non-trivial case */

	if (sign)
	{
		if (sign <= 2)
		{
			if (!(r[1]) && (r[0] == 1))
				q[0] = -(q[0]);
			else
			{
				zadd(q, one, &q);
				q[0] = -q[0];
				if (sign == 1)
					zsub(r, b, &r);
				else
					zsub(b, r, &r);
			}
		}
		else
			znegate(&r);
	}
	*qqq = q;
	*rrr = r;

}

#else

#define lower_radix(in,out) { \
register long i=1,j=0,rb=0,nrb=0,bits; \
register verylong pin=in; \
long size; \
 \
if ((size= *pin)==2) { \
        pin++; \
        out[1] = *pin & ALPHA50RADIXM; \
        out[2] = *pin >> ALPHA50NBITS; \
        pin++; \
        out[2] |= (*pin & ((1L << (2*ALPHA50NBITS-NBITS))-1L)) << (NBITS - \
         ALPHA50NBITS); \
        bits = *pin >> (2*ALPHA50NBITS-NBITS); \
        if (bits) { \
                out[3] = bits; \
                out[0]=3; \
                } \
        else \
                out[0]=2; \
        } \
else { \
        for (j=0; i <= size ; j++) { \
                if (nrb >= ALPHA50NBITS) { \
                        bits = rb & ALPHA50RADIXM; \
                        rb >>= ALPHA50NBITS; \
                        nrb -= ALPHA50NBITS; \
                        } \
                else { \
                        bits = rb | ((((1L << (ALPHA50NBITS - nrb)) -1L) & \
                         pin[i]) << nrb); \
                        rb = pin[i] >> (ALPHA50NBITS - nrb); \
                        nrb += (NBITS-ALPHA50NBITS); \
                        i++; \
                        } \
                out[j+1]=bits; \
                } \
        if (rb) { \
                j++; \
                if (rb < (1L<<ALPHA50NBITS)) \
                        out[j]= rb; \
                else { \
                        out[j]= rb & (( 1L << ALPHA50NBITS) -1); \
                        j++; \
                        out[j] = rb>>ALPHA50NBITS; \
                        } \
                } \
        out[0]=j; \
        } \
}

/* convert from ALPHA50NBITS to NBITS radix */
/* this can probably be made faster */

#define raise_radix(size,in,out) { \
register long i,j,nrb,bits; \
 \
nrb=ALPHA50NBITS; \
for (i=0,j=1; i<size; j++) { \
        bits = (in[i++] >> (ALPHA50NBITS-nrb)); \
        if (i <size ) { \
                if (nrb + ALPHA50NBITS >= NBITS) { \
                        bits = bits | ((in[i] & ( (1L<< (NBITS-nrb)) - 1)) << \
                         nrb); \
                        nrb -= (NBITS-ALPHA50NBITS); \
                        if (!nrb) { i++; nrb=ALPHA50NBITS; } \
                        } \
                else { \
                        bits = bits | (in[i++] << nrb); \
                        if (i< size) \
                                bits = bits | ((in[i] & ((1L <<( (NBITS- \
                                 ALPHA50NBITS)-nrb)) -1)) << (ALPHA50NBITS \
                                 + nrb)); \
                        nrb = 2*ALPHA50NBITS-NBITS+nrb; \
                        if (!nrb) { i++; nrb=ALPHA50NBITS; } \
                        } \
                } \
        out[j]=bits; \
        } \
if (j>2 && !out[j-1]) j--; \
out[0]=j-1; \
}


#define raise_radix2(size,in,out) { \
register long i,j,nrb,bits; \
 \
nrb=ALPHA50NBITS; \
for (i=1,j=1; i<=size; j++) { \
        bits = (in[i++] >> (ALPHA50NBITS-nrb)); \
        if (i <=size ) { \
                if (nrb + ALPHA50NBITS >= NBITS) { \
                        bits = bits | ((in[i] & ( (1L<< (NBITS-nrb)) - 1)) << nrb); \
                        nrb -= (NBITS-ALPHA50NBITS); \
                        if (!nrb) { i++; nrb=ALPHA50NBITS; } \
                        } \
                else { \
                        bits = bits | (in[i++] << nrb); \
                        if (i<=size) \
                                bits = bits | ((in[i] & ((1L <<( (NBITS- \
                                 ALPHA50NBITS)-nrb)) -1)) << (ALPHA50NBITS+nrb)) ; \
                        nrb = 2*ALPHA50NBITS-NBITS+nrb; \
                        if (!nrb) { i++; nrb=ALPHA50NBITS; } \
                        } \
                } \
        out[j]=bits; \
        } \
if (j>2 && !out[j-1]) j--; \
out[0]= (j > 1) ? j-1:1 ; \
}

#define zalphaaddmulone(ama, amb) \
{ \
        register long lami; \
        register long lams = 0; \
        register verylong lama = (ama); \
        register verylong lamb = (amb); \
 \
        lams = 0; \
        for (lami = (*lamb++); lami > 0; lami--) \
        { \
                lams += (*lama + *lamb++); \
                *lama++ = lams & ALPHA50RADIXM; \
                lams >>= ALPHA50NBITS; \
        } \
        *lama += lams; \
}


#define zaldivaddmulp2(_a, _b, _d, _t) \
{ \
        register long lb = (_b); \
        register long ld = (_d); \
        register long b1 = (_b) & ALPHA50RADIXROOTM; \
        register long d1 = (_d) & ALPHA50RADIXROOTM; \
        register long aa = *(_a) + b1 * d1; \
 \
        b1 = b1 * (ld >>= ALPHA50NBITSH) + d1 * (lb >>= ALPHA50NBITSH) + (aa >> ALPHA50NBITSH); \
        aa = (aa & ALPHA50RADIXM) + ((b1 & RADIXROOTM) << ALPHA50NBITSH) + *(_t); \
        *(_t) = ld * lb + (b1 >> ALPHA50NBITSH) + (aa >> ALPHA50NBITS); \
        *(_a) = (aa & ALPHA50RADIXM); \
}


static void
zalphasubmul2(
          long r,
          verylong a,
          verylong b
          )
{
          register long rd = ALPHA50RADIX - r;
          register long i;
          long carry = ALPHA50RADIX;

          for (i = (*b++); i > 0; i--)
          {
                  zaldivaddmulp2(a, *b, rd, &carry);
                  a++;
                  carry += ALPHA50RADIXM - (*b++);
          }
          *a += carry - ALPHA50RADIX;    /* unnormalized */
}

void
zdiv2(
        verylong in_a,
        verylong in_b,
        verylong *qqq,
        verylong *rrr
        )

{
        register long i;
        register long qq;
        long sa;
        long sb;
        long sq;
        verylong p;
        verylong pc;
        long sign;
        static verylong a = 0;
        static verylong b = 0;
        static verylong c = 0;
        static verylong d = 0;
        double btopinv;
        double aux;
        verylong q = *qqq;
        verylong r = *rrr;

#ifndef START
        if (fudge < 0)
                zstart0();
#endif
        if (ALLOCATE && !in_a)
        {
                zzero(qqq);
                zzero(rrr);
                return;
        }
        if ((!in_b) || (((sb = in_b[0]) == 1) && (!in_b[1])))
        {
                zhalt("division by zero in zdiv");
                return;
        }

        sign = (*in_a < 0 ? 2 : 0) | (*in_b < 0 ? 1 : 0);
        if (*in_a < 0)
                (*in_a) = (-(*in_a));
        if (*in_b < 0)
                (*in_b) = (-(*in_b));

        if ((*in_b == 1) && (in_b[1] < ALPHA50RADIX)) {
                zintoz(zsdiv(in_a, in_b[1], qqq), rrr);
                goto done;
                }

        zsetlength(&a,(long) (NBITS/(double)ALPHA50NBITS*(double)in_a[0]+1.0), "");
        lower_radix(in_a,a);    /* test this for large integers */

        /* convert b */
        zsetlength(&b,(long) (NBITS/(double)ALPHA50NBITS*(double)in_b[0]+1.0), "");
        lower_radix(in_b,b);

        sa = (*a);
        sb = (*b);
        zsetlength(&c, (i = (sa > sb ? sa : sb) + 1), "in zdiv, locals\n");
        zsetlength(&d, i, "");
        zsetlength(&q, i, "in zdiv, third argument");
        *qqq = q;
        zsetlength(&r, sb + 1, "in zdiv, fourth argument");
        *rrr = r;
        p = &b[sb];

        sq = sa - sb;
        btopinv = (double) (*p) * alpha50fradix;
        if (sb > 1)
                btopinv += (*(--p));
        btopinv *= (double) alpha50fradix;
        if (sb > 2)
                btopinv += (*(p - 1));
        btopinv = alpha50fudge / btopinv;
        p = &a[1];
        pc = &c[0];
        for (i = sa; i > 0; i--)
                *pc++ = *p++;
        p = pc;
        sa = 1 - sb;
        for (i = (-sq); i > 0; i--)
                *p++ = 0;
        *pc = 0;
        d[1] = 0;
        p = &d[sq + 1];
        for (i = sq; i >= 0; i--)
        {
                aux = alpha50fradix * (alpha50fradix * (*pc) + (*(pc - 1))) + 1.0;
                if (i > sa)
                        aux += (*(pc - 2));
                qq = (long) (btopinv * aux + 0.5);
        /* dirty, but safe. On most machines (long)(btopinv*aux) */
        /* is correct, or one too big; on some however it becomes */
        /* too small. Could change zstart, but +0.5 and a while */
        /* instead of one single if is safer */
                if (qq >= ALPHA50RADIX)
                        qq = ALPHA50RADIX-1;
                if (qq > 0)
                {
                        zalphasubmul2(qq, &c[i], &b[0]);
                        while (*pc < 0)
                        {
                                qq--;
                                {
                                        zalphaaddmulone(&c[i], &b[0]);
                                }
                        }
                }
                pc--;
                *p-- = qq;
        }
        sb--;
        while ((sb > 0) && (!(c[sb])))
                sb--;
        sb++;
        raise_radix(sb,c,r);
        if (sq < 0)
        {
                q[0] = 1;
                q[1] = d[1];
        }
        else
        {
                sq++;
                while ((sq > 1) && (!(d[sq])))
                        sq--;
                raise_radix2(sq,d,q);
        }

        if (sign)
        {
                if (sign <= 2)
                {
                        if (!(r[1]) && (r[0] == 1))
                                q[0] = -(q[0]);
                        else
                        {
                                zadd(q, one, &q);
                                q[0] = -q[0];
                                if (sign == 1)
                                        zsub(r, b, &r);
                                else
                                        zsub(b, r, &r);
                        }
                }
                else
                        znegate(&r);
        }
        *qqq = q;
        *rrr = r;
done:   ;

}
#endif


static void
zsubmul3(
        long r,
        verylong a,
        verylong b
        )
{
        register long rd = RADIX - r;
        register long i;
        long carry = RADIX;

        for (i = (*b++); i > 0; i--)
        {
                zaddmulp3(a, *b, rd, &carry);
                a++;
                carry += RADIXM - (*b++);
        }
        *a += carry - RADIX;    /* unnormalized */
}

#define correct3( q, x1, x2, x3, y1, y2, btopinv) { \
        register long ix1=x1,ix2=x2,ix3=x3,iy1=y1,iy2=y2; \
        long lprodlow=0,lprodhigh=0; \
\
        zaddmulp3(&lprodlow, q, iy2, &lprodhigh); \
        if ((ix3 -= lprodlow) < 0) { \
                ix3+= RADIX; \
                ix2--; \
        } \
        lprodlow=0; \
        zaddmulp3(&lprodlow, q, iy1, &lprodhigh); \
        if ((ix2 -= lprodlow) < 0) { \
                ix2 += RADIX; \
                --ix1; \
        } \
        if (ix1 < lprodhigh) q--; \
        else if ((ix1 -= lprodhigh)) { \
                q += (long) ((fradix * (fradix * ix1 + ix2))*btopinv ); \
        } \
        else { \
                q += (long) ((fradix * ix2 + ix3)*btopinv ); \
        } \
}


#ifndef ALPHA50
void
zdiv3(
	verylong in_a,
	verylong in_b,
	verylong *qqq,
	verylong *rrr
	)

{
	register long i;
	register long qq;
	long sa;
	long sb;
	long sq;
	verylong p;
	verylong pc;
	long sign;
	static verylong a = 0;
	static verylong b = 0;
	static verylong c = 0;
	static verylong d = 0;
	double btopinv;
	double aux;
	verylong q = *qqq;
	verylong r = *rrr;
#ifdef ALPHA
        double btopinv2;
#endif

#ifndef START
	if (fudge < 0)
		zstart0();
#endif

	if (ALLOCATE && !in_a)
	{
		zzero(qqq);
		zzero(rrr);
		return;
	}
	if ((!in_b) || (((sb = in_b[0]) == 1) && (!in_b[1])))
	{
		zhalt("division by zero in zdiv");
		return;
	}
	zcopy(in_a,&a);

	zcopy(in_b,&b);

	sign = (*a < 0 ? 2 : 0) | (*b < 0 ? 1 : 0);
	if (*a < 0)
		(*a) = (-(*a));
	sa = (*a);
	if (*b < 0)
		(*b) = (-(*b));
	sb = (*b);
	zsetlength(&c, (i = (sa > sb ? sa : sb) + 1), "in zdiv, locals\n");
	zsetlength(&d, i, "");
	zsetlength(&q, i, "in zdiv, third argument");
	*qqq = q;
	zsetlength(&r, sb + 1, "in zdiv, fourth argument");
	*rrr = r;
	p = &b[sb];
	if ((sb == 1) && (*p < RADIX))
		zintoz(zsdiv(a, *p, &q), &r);
	else
	{
		sq = sa - sb;	/* size of quotient */
		btopinv = (double) (*p) * fradix;
		if (sb > 1)
			btopinv += (*(--p));
		btopinv *= fradix;
#ifdef ALPHA
                btopinv2=btopinv;
                btopinv2 = fudge2/btopinv2;
#else
		if (sb > 2)
			btopinv += (*(p - 1));
#endif
		btopinv = fudge / btopinv;
		p = &a[1];
		pc = &c[0];
		for (i = sa; i > 0; i--)
			*pc++ = *p++;
		p = pc;
		sa = 1 - sb;
		for (i = (-sq); i > 0; i--)
			*p++ = 0;
		*pc = 0;
		d[1] = 0;
		p = &d[sq + 1];
		for (i = sq; i >= 0; i--)
		{

			aux = fradix * (fradix * (*pc) + (*(pc - 1))) + 1.0;
#ifndef ALPHA
			if (i > sa)
				aux += (*(pc - 2));
#endif
#ifdef ALPHA
			qq = (long) (btopinv2 * aux + 0.5);
#else
			qq = (long) (btopinv * aux + 0.5);
		/* dirty, but safe. On most machines (long)(btopinv*aux) */
		/* is correct, or one too big; on some however it becomes */
		/* too small. Could change zstart, but +0.5 and a while */
		/* instead of one single if is safer */
#endif
			if (qq > 0)
			{
#ifdef ALPHA
                                if (qq > (1L<<48)) {
                                        correct3(qq,*pc,*(pc-1),(i>sa) ? *(pc-2):0,b[sb],b[sb-1],btopinv);
                                        if (qq >= RADIX)
                                                qq = RADIX-1;
                                }
#else
				if (qq >= RADIX)
					qq = RADIX-1;
#endif
				zsubmul3(qq, &c[i], &b[0]);
				while (*pc < 0)
				{
					qq--;
					{
						zaddmulone(&c[i], &b[0]);
					}
				}
			}
			pc--;
			*p-- = qq;
		}
		sb--;
		while ((sb > 0) && (!(c[sb])))
			sb--;
		sb++;
		r[0] = sb;
		p = &r[1];
		pc = &c[0];
		for (i = sb; i > 0; i--)
			*p++ = *pc++;
		if (sq < 0)
		{
			q[0] = 1;
			q[1] = d[1];
		}
		else
		{
			sq++;
			while ((sq > 1) && (!(d[sq])))
				sq--;
			q[0] = sq;
			p = &q[1];
			pc = &d[1];
			for (i = sq; i > 0; i--)
				*p++ = *pc++;
		}
	}	/* non-trivial case */

	if (sign)
	{
		if (sign <= 2)
		{
			if (!(r[1]) && (r[0] == 1))
				q[0] = -(q[0]);
			else
			{
				zadd(q, one, &q);
				q[0] = -q[0];
				if (sign == 1)
					zsub(r, b, &r);
				else
					zsub(b, r, &r);
			}
		}
		else
			znegate(&r);
	}
	*qqq = q;
	*rrr = r;

}

#else

#define zaldivaddmulp3(  a,  b,  d,  t) { \
        register long b1 = b & ALPHA50RADIXROOTM; \
        register long d1 = d & ALPHA50RADIXROOTM; \
        register long bd,b1d1,m,aa= (*a) + (*t); \
        register long ld = (d>>ALPHA50NBITSH); \
        register long lb = (b>>ALPHA50NBITSH); \
 \
        bd=lb*ld; \
        b1d1=b1*d1; \
        m=(lb+b1)*(ld+d1) - bd - b1d1; \
        aa += ( b1d1+ ((m&ALPHA50RADIXROOTM)<<ALPHA50NBITSH)); \
        (*t) = (aa >> ALPHA50NBITS) + bd + (m>>ALPHA50NBITSH); \
        (*a) = aa & ALPHA50RADIXM; \
}

static void
zalphasubmul3(
          long r,
          verylong a,
          verylong b
          )
{
          register long rd = ALPHA50RADIX - r;
          register long i;
          long carry = ALPHA50RADIX;

          for (i = (*b++); i > 0; i--)
          {
                  zaldivaddmulp3(a, *b, rd, &carry);
                  a++;
                  carry += ALPHA50RADIXM - (*b++);
          }
          *a += carry - ALPHA50RADIX;    /* unnormalized */
}

void
zdiv3(
        verylong in_a,
        verylong in_b,
        verylong *qqq,
        verylong *rrr
        )

{
        register long i;
        register long qq;
        long sa;
        long sb;
        long sq;
        verylong p;
        verylong pc;
        long sign;
        static verylong a = 0;
        static verylong b = 0;
        static verylong c = 0;
        static verylong d = 0;
        double btopinv;
        double aux;
        verylong q = *qqq;
        verylong r = *rrr;

/*printf("zdiv: "); zwrite(in_a); fprintf(stderr," div "); zwriteln(in_b);*/
#ifndef START
        if (fudge < 0)
                zstart();
#endif
        if (ALLOCATE && !in_a)
        {
                zzero(qqq);
                zzero(rrr);
                return;
        }
        if ((!in_b) || (((sb = in_b[0]) == 1) && (!in_b[1])))
        {
                zhalt("division by zero in zdiv");
                return;
        }

        sign = (*in_a < 0 ? 2 : 0) | (*in_b < 0 ? 1 : 0);
        if (*in_a < 0)
                (*in_a) = (-(*in_a));
        if (*in_b < 0)
                (*in_b) = (-(*in_b));

        if ((*in_b == 1) && (in_b[1] < ALPHA50RADIX)) {
                zintoz(zsdiv(in_a, in_b[1], qqq), rrr);
                goto done;
                }

        zsetlength(&a,(long) (NBITS/(double)ALPHA50NBITS*(double)in_a[0]+1.0),
                "");
        lower_radix(in_a,a);    /* test this for large integers */

        /* convert b */
        zsetlength(&b,(long) (NBITS/(double)ALPHA50NBITS*(double)in_b[0]+1.0),
                "");
        lower_radix(in_b,b);

        sa = (*a);
        sb = (*b);
        zsetlength(&c, (i = (sa > sb ? sa : sb) + 1), "in zdiv, locals\n");
        zsetlength(&d, i, "");
        zsetlength(&q, i, "in zdiv, third argument");
        *qqq = q;
        zsetlength(&r, sb + 1, "in zdiv, fourth argument");
        *rrr = r;
        p = &b[sb];

        sq = sa - sb;
        btopinv = (double) (*p) * alpha50fradix;
        if (sb > 1)
                btopinv += (*(--p));
        btopinv *= (double) alpha50fradix;
        if (sb > 2)
                btopinv += (*(p - 1));
        btopinv = alpha50fudge / btopinv;
        p = &a[1];
        pc = &c[0];
        for (i = sa; i > 0; i--)
                *pc++ = *p++;
        p = pc;
        sa = 1 - sb;
        for (i = (-sq); i > 0; i--)
                *p++ = 0;
        *pc = 0;
        d[1] = 0;
        p = &d[sq + 1];
        for (i = sq; i >= 0; i--)
        {
                aux = alpha50fradix * (alpha50fradix * (*pc) +
                    (*(pc - 1))) + 1.0;
                if (i > sa)
                        aux += (*(pc - 2));
                qq = (long) (btopinv * aux + 0.5);
        /* dirty, but safe. On most machines (long)(btopinv*aux) */
        /* is correct, or one too big; on some however it becomes */
        /* too small. Could change zstart, but +0.5 and a while */
        /* instead of one single if is safer */
                if (qq >= ALPHA50RADIX)
                        qq = ALPHA50RADIX-1;
                if (qq > 0)
                {
                        zalphasubmul3(qq, &c[i], &b[0]);
                        while (*pc < 0)
                        {
                                qq--;
                                {
                                        zalphaaddmulone(&c[i], &b[0]);
                                }
                        }
                }
                pc--;
                *p-- = qq;
        }
        sb--;
        while ((sb > 0) && (!(c[sb])))
                sb--;
        sb++;
        raise_radix(sb,c,r);

        if (sq < 0)
        {
                q[0] = 1;
                q[1] = d[1];
        }
        else
        {
                sq++;
                while ((sq > 1) && (!(d[sq])))
                        sq--;
                raise_radix2(sq,d,q);
        }

        if (sign)
        {
                if (sign <= 2)
                {
                        if (!(r[1]) && (r[0] == 1))
                                q[0] = -(q[0]);
                        else
                        {
                                zadd(q, one, &q);
                                q[0] = -q[0];
                                if (sign == 1)
                                        zsub(r, b, &r);
                                else
                                        zsub(b, r, &r);
                        }
                }
                else
                        znegate(&r);
        }
        *qqq = q;
        *rrr = r;
done:   ;

}
#endif


static
void zsubmul5(
        long ams,
        long *ama,
        long *amb
        )
{
        register long carry = (1L<<26);
        register long i = (*amb++);
        register double dams = (double) ((1L<<26)-ams);
        double xx;
        register double yy;
        register unsigned long lo;
        register unsigned long hi;

        xx =  ((double) (*amb))*dams + 4503599627370496.0;
        for (; i > 1; i--)
        {
                yy =  ((double) (*(amb+1)))*dams + 4503599627370496.0;
                lo = LO5(xx) & 0x3FFFFFF;
                hi = ((HI5(xx)<<6)|(LO5(xx)>>26)) & 0x3FFFFFF;
                lo = lo + (*ama) + carry;
                *ama = lo & 0x3FFFFFF;
                carry = hi + (lo >> 26);
                ama++;
                carry += ((1L<<26)-1) - (*(amb++));
                xx = yy;
        }
        lo = LO5(xx) & 0x3FFFFFF;
        hi = ((HI5(xx)<<6)|(LO5(xx)>>26)) & 0x3FFFFFF;
        lo = lo + (*ama) + carry;
        *ama = lo & 0x3FFFFFF;
        carry = hi + (lo >> 26);
        ama++;
        carry += ((1L<<26)-1) - (*amb);
        *ama += carry - (1L<<26);
}

void
zdiv5(
	verylong in_a,
	verylong in_b,
	verylong *qqq,
	verylong *rrr
	)

{
	register long i;
	register long qq;
	long sa;
	long sb;
	long sq;
	verylong p;
	verylong pc;
	long sign;
	static verylong a = 0;
	static verylong b = 0;
	static verylong c = 0;
	static verylong d = 0;
	double btopinv;
	double aux;
	verylong q = *qqq;
	verylong r = *rrr;

#ifndef START
	if (fudge1 < 0)
		zstart1();
#endif

	if (ALLOCATE && !in_a)
	{
		zzero(qqq);
		zzero(rrr);
		return;
	}
	if ((!in_b) || (((sb = in_b[0]) == 1) && (!in_b[1])))
	{
		zhalt("division by zero in zdiv");
		return;
	}
	zcopy(in_a,&a);

	zcopy(in_b,&b);

	sign = (*a < 0 ? 2 : 0) | (*b < 0 ? 1 : 0);
	if (*a < 0)
		(*a) = (-(*a));
	sa = (*a);
	if (*b < 0)
		(*b) = (-(*b));
	sb = (*b);
	zsetlength(&c, (i = (sa > sb ? sa : sb) + 1), "in zdiv, locals\n");
	zsetlength(&d, i, "");
	zsetlength(&q, i, "in zdiv, third argument");
	*qqq = q;
	zsetlength(&r, sb + 1, "in zdiv, fourth argument");
	*rrr = r;
	p = &b[sb];
	if ((sb == 1) && (*p < (1L<<26)))
		zintoz(zsdiv(a, *p, &q), &r);
	else
	{
		sq = sa - sb;	/* size of quotient */
		btopinv = (double) (*p) * 67108864.0;
		if (sb > 1)
			btopinv += (*(--p));
		btopinv *= 67108864.0;
		if (sb > 2)
			btopinv += (*(p - 1));
		btopinv = fudge1 / btopinv;
		p = &a[1];
		pc = &c[0];
		for (i = sa; i > 0; i--)
			*pc++ = *p++;
		p = pc;
		sa = 1 - sb;
		for (i = (-sq); i > 0; i--)
			*p++ = 0;
		*pc = 0;
		d[1] = 0;
		p = &d[sq + 1];
		for (i = sq; i >= 0; i--)
		{

			aux = 67108864.0 * (67108864.0 * (*pc) + (*(pc - 1))) + 1.0;
			if (i > sa)
				aux += (*(pc - 2));
			qq = (long) (btopinv * aux + 0.5);
		/* dirty, but safe. On most machines (long)(btopinv*aux) */
		/* is correct, or one too big; on some however it becomes */
		/* too small. Could change zstart, but +0.5 and a while */
		/* instead of one single if is safer */
			if (qq > 0)
			{
				if (qq >= (1L<<26))
					qq = (1L<<26)-1;
				zsubmul5(qq, &c[i], &b[0]);
				while (*pc < 0)
				{
					qq--;
					{
						zaddmulone1(&c[i], &b[0]);
					}
				}
			}
			pc--;
			*p-- = qq;
		}
		sb--;
		while ((sb > 0) && (!(c[sb])))
			sb--;
		sb++;
		r[0] = sb;
		p = &r[1];
		pc = &c[0];
		for (i = sb; i > 0; i--)
			*p++ = *pc++;
		if (sq < 0)
		{
			q[0] = 1;
			q[1] = d[1];
		}
		else
		{
			sq++;
			while ((sq > 1) && (!(d[sq])))
				sq--;
			q[0] = sq;
			p = &q[1];
			pc = &d[1];
			for (i = sq; i > 0; i--)
				*p++ = *pc++;
		}
	}	/* non-trivial case */

	if (sign)
	{
		if (sign <= 2)
		{
			if (!(r[1]) && (r[0] == 1))
				q[0] = -(q[0]);
			else
			{
				zadd1(q, one, &q);
				q[0] = -q[0];
				if (sign == 1)
					zsub1(r, b, &r);
				else
					zsub1(b, r, &r);
			}
		}
		else
			znegate(&r);
	}
	*qqq = q;
	*rrr = r;

}




verylong aa[8][4],a=0,a1=0;
verylong bb[8][4],b=0,b1=0;
verylong c=0;
verylong d=0;
verylong e=0, e1=0,e2=0,e3=0;
verylong f=0, f1=0,f2=0,f3=0;
verylong g=0;
long g0=1,g1=1,g2=1, g3=1;
double tt, tcnt, besttime,tot;
double ntime,ktime;
long bestmul= -1;
long iter;
FILE *fpar;

main () {
	register long i,j,k;
	FILE *fp;

	fp = fopen("TESTRANGE","r");
	if (fp != NULL) {
		long www,xxx,yyy,zzz;
		if (fscanf(fp,"%ld %ld %ld",&www,&xxx,&yyy,&zzz)==4) {
			NUMBER = www;
			FIRST = xxx;
			INCR =yyy;
			NUMBER = zzz;
		}
		fclose(fp);
	}

	fp = fopen("SEED","r");
	if (fp != NULL) {
		long xxx;
		fscanf(fp,"%ld",&xxx);
		SEED = xxx;
		fclose(fp);
	}

	fpar = fopen("lippar.h","w");

	test_multiplies();

	if (NUMBER > 10)
		NUMBER = 10;

	zrstarts(SEED);
	/* get random aa and bb */
	for (j=0; j<NUMBER; j++) {
		aa[j][0]=0; aa[j][1]=0; aa[j][2]=0; aa[j][3]=0;
		zrandoml(FIRST + j*INCR, &aa[j][0],zrandomb);
		zrandoml(FIRST + j*INCR, &aa[j][1],zrandomb);
		zrandoml(FIRST + j*INCR, &aa[j][2],zrandomb);
		zrandoml(FIRST + j*INCR, &aa[j][3],zrandomb);
		bb[j][0]=0; bb[j][1]=0; bb[j][2]=0; bb[j][3]=0;
		zrandoml(FIRST + j*INCR, &bb[j][0],zrandomb);
		zrandoml(FIRST + j*INCR, &bb[j][1],zrandomb);
		zrandoml(FIRST + j*INCR, &bb[j][2],zrandomb);
		zrandoml(FIRST + j*INCR, &bb[j][3],zrandomb);
		if (zcompare(aa[j][0],bb[j][0]) ==1)
			zswap(&aa[j][0],&bb[j][0]);
		if (zcompare(aa[j][1],bb[j][1]) ==1)
			zswap(&aa[j][1],&bb[j][1]);
		if (zcompare(aa[j][2],bb[j][2]) ==1)
			zswap(&aa[j][2],&bb[j][2]);
		if (zcompare(aa[j][3],bb[j][3]) ==1)
			zswap(&aa[j][3],&bb[j][3]);
	}
	fprintf(stderr,"doing %ld multiplications, divisions, squares\n\n",ITER);

	mult_time();
	squ_time();

	iter = ITER/4;

	if (iter < 10 )
		iter = 10;

	K_M_C = 30;
	K_S_C = 30;

	mult_crov();
	kar_crov();

	div_time();
	exit (0);
}





test_multiplies() {

#ifndef ALPHA_OR_ALPHA50
	if (NBITS != 30) {
		fprintf(stderr,"recompile, NBITS is %d, should be 30\n", NBITS);
		exit(1);
	}
#endif

	zsread("1152921504606846975",&a);
	zsread("1237940039285380274899124223",&b);
	zsread("1427247692705959879820345929011193356876775425",&d);

#ifndef ALPHA_OR_ALPHA50
	zmul0(a,b,&c);
	if (zcompare(c,d)) {
		fprintf(stderr,"error with SINGLE_MUL = 0, PLAIN = 0, KARAT = 0\n");
		g0=0;
	}

	zzero(&c);
	/* convert a and b to 26 bits e and f */
	zsetlength(&e,(long) (30.0/(double)26.0*(double)30.0+1.0), "");
	zsetlength(&e1,(long) (30.0/(double)26.0*(double)30.0+1.0), "");
	zsetlength(&e2,(long) (30.0/(double)26.0*(double)30.0+1.0), "");
	zsetlength(&e3,(long) (30.0/(double)26.0*(double)30.0+1.0), "");
	lower_radix1(a,e);
	zsetlength(&f,(long) (30.0/(double)26.0*(double)30.0+1.0), "");
	zsetlength(&f1,(long) (30.0/(double)26.0*(double)30.0+1.0), "");
	zsetlength(&f2,(long) (30.0/(double)26.0*(double)30.0+1.0), "");
	zsetlength(&f3,(long) (30.0/(double)26.0*(double)30.0+1.0), "");
	lower_radix1(b,f);
	zmul1(e,f,&c);
	/* convert d to 26 bits g */
	zsetlength(&g,(long) (30.0/(double)26.0*(double)d[0]+1.0), "");
	lower_radix1(d,g);
	if (zcompare(c,g)) {
		fprintf(stderr,"error with SINGLE_MUL = 1, PLAIN = 0, KARAT = 0\n");
		g1=0;
		fprintf(stderr,"trying DOUBLES_LOW_HIGH\n");
		zmul5(e,f,&c);
		if (!zcompare(c,g)) {
			g1=2;
		}
		else fprintf(stderr,"unsuccessful\n");
	}
#endif

	zzero(&c);
	zmul2(a,b,&c);
	if (zcompare(c,d)) {
		fprintf(stderr,"error with SINGLE_MUL = 0, PLAIN = 1, KARAT = 0\n");
		g2=0;
	}

	zzero(&c);
	zmul3(a,b,&c);
	if (zcompare(c,d)) {
		fprintf(stderr,"error with SINGLE_MUL = 0, PLAIN = 0, KARAT = 1\n");
		g3=0;
	}

	zsread("387294387387924879382974832874324",&a);
	zsread("8934732984782343834738247832473824732847834324",&b);
	zsread("3460371937815963399338956444782987599400115707528053812889435363520183465496976",&d);

#ifndef ALPHA_OR_ALPHA50
	if (g0) {
		zmul0(a,b,&c);
		if (zcompare(c,d)) {
			fprintf(stderr,"error with SINGLE_MUL = 0, PLAIN = 0, KARAT = 0\n");
			g0=0;
		}
	}

	if (g1==1) {
		zzero(&c);
		/* convert a and b to 26 bits e and f */
		zsetlength(&e,(long) (30.0/(double)26.0*(double)a[0]+1.0), "");
		lower_radix1(a,e);
		zsetlength(&f,(long) (30.0/(double)26.0*(double)b[0]+1.0), "");
		lower_radix1(b,f);
		zmul1(e,f,&c);
		/* convert d to 26 bits g */
		zsetlength(&g,(long) (30.0/(double)26.0*(double)d[0]+1.0), "");
		lower_radix1(d,g);
		if (zcompare(c,g)) {
			fprintf(stderr,"error with SINGLE_MUL = 1, PLAIN = 0, KARAT = 0\n");
			g1=0;
		}
	}
	else if (g1==2) {
		zzero(&c);
		/* convert a and b to 26 bits e and f */
		zsetlength(&e,(long) (30.0/(double)26.0*(double)a[0]+1.0), "");
		lower_radix1(a,e);
		zsetlength(&f,(long) (30.0/(double)26.0*(double)b[0]+1.0), "");
		lower_radix1(b,f);
		zmul5(e,f,&c);
		/* convert d to 26 bits g */
		zsetlength(&g,(long) (30.0/(double)26.0*(double)d[0]+1.0), "");
		lower_radix1(d,g);
		if (zcompare(c,g)) {
			fprintf(stderr,"error with SINGLE_MUL = 1, PLAIN = 0, KARAT = 0\n");
			g1=0;
		}
		else {
			fprintf(stderr,"successful\n");
			fprintf(stderr,"must define DOUBLES_LOW_HIGH for SINGLE_MUL = 1, PLAIN = 0, KARAT = 0\n");
		}
	}
#endif

	if (g2) {
		zzero(&c);
		zmul2(a,b,&c);
		if (zcompare(c,d)) {
			fprintf(stderr,"error with SINGLE_MUL = 0, PLAIN = 1, KARAT = 0\n");
			g2=0;
		}
	}

	if (g3) {
		zzero(&c);
		zmul3(a,b,&c);
		if (zcompare(c,d)) {
			fprintf(stderr,"error with SINGLE_MUL = 0, PLAIN = 0, KARAT = 1\n");
			g3=0;
		}
	}

}


mult_time() {

	register int i,j,k;
	double time;

	fprintf(stderr,"multiplication timings:\n\n");

#ifndef ALPHA_OR_ALPHA50
	if (g0) {
		fprintf(stderr,"timing SINGLE_MUL = 0, PLAIN = 0, KARAT = 0\n\n");
		fprintf(stderr,"\t");
		for (i=0; i<NUMBER; i++)
			fprintf(stderr,"\t%ldbit",FIRST + i*INCR);
		fprintf(stderr,"\n");
		fprintf(stderr,"\t");
		for (i=0; i<NUMBER; i++)
			fprintf(stderr,"\t-------");
		fprintf(stderr,"\n");
		tot = 0.0;
		for (j=0; j<NUMBER; j++) {
			i=FIRST + j*INCR;
			if (i<1000)
				fprintf(stderr,"%ldbit:\t",i);
			else
				fprintf(stderr,"%ldbit:",i);
			for (k=0; k<j ; k++)
				fprintf(stderr,"\t");
			for (;  k<NUMBER; k++) {
				if (bb[k][0][0] > K_M_C)
					K_M_C=bb[k][0][0] + 1;
				tt = gettime();
				for (i=0; i<ITER/4; i++) {
					zmul0(aa[j][0],bb[k][0],&c);
					zmul0(aa[j][1],bb[k][1],&c);
					zmul0(aa[j][2],bb[k][2],&c);
					zmul0(aa[j][3],bb[k][3],&c);
				}
				tcnt= gettime();
				time = tcnt-tt;
				fprintf(stderr,"\t%6.3lf",time);
				tot += time/(double) (i*(FIRST + k*INCR));
			}
			fprintf(stderr,"\n");
		}
		bestmul = 0;
		besttime = tot;
		fprintf(stderr,"\n");
	}

	if (g1) {
		fprintf(stderr,"timing SINGLE_MUL = 1, PLAIN = 0, KARAT = 0\n\n");
                fprintf(stderr,"\t");
                for (i=0; i<NUMBER; i++)
                        fprintf(stderr,"\t%ldbit",FIRST + i*INCR);
                fprintf(stderr,"\n");
                fprintf(stderr,"\t");
                for (i=0; i<NUMBER; i++)
                        fprintf(stderr,"\t-------");
                fprintf(stderr,"\n");
		tot = 0.0;
		for (j=0; j<NUMBER; j++) {
                        i=FIRST + j*INCR;
                        if (i<1000)
                                fprintf(stderr,"%ldbit:\t",i);
                        else
                                fprintf(stderr,"%ldbit:",i);
			zcopy(aa[j][0],&a);
			lower_radix1(a,e);
			zcopy(aa[j][1],&a);
			lower_radix1(a,e1);
			zcopy(aa[j][2],&a);
			lower_radix1(a,e2);
			zcopy(aa[j][3],&a);
			lower_radix1(a,e3);
			for (k=0; k<j ; k++)
				fprintf(stderr,"\t");
			for (;  k<NUMBER; k++) {
				lower_radix1((bb[k][0]),f);
				lower_radix1((bb[k][1]),f1);
				lower_radix1((bb[k][2]),f2);
				lower_radix1((bb[k][3]),f3);
				if (f[0] > K_M_C)
					K_M_C = f[0] + 1;
				tt = gettime();
				if (g1==1)
					for (i=0; i<ITER/4; i++) {
						zmul1(e,f,&c);
						zmul1(e3,f3,&c);
						zmul1(e2,f2,&c);
						zmul1(e1,f1,&c);
					}
				else
					for (i=0; i<ITER/4; i++) {
						zmul5(e,f,&c);
						zmul5(e3,f3,&c);
						zmul5(e2,f2,&c);
						zmul5(e1,f1,&c);
					}
				tcnt= gettime();
				time = tcnt-tt;
				fprintf(stderr,"\t%6.3lf",time);
				tot += time/(double) (i*(FIRST + k*INCR));
			}
			fprintf(stderr,"\n");
		}
		if (bestmul < 0 || tot < besttime) {
			bestmul = 1;
			besttime = tot;
		}
		fprintf(stderr,"\n");
	}
#endif

	if (g2) {
		fprintf(stderr,"timing SINGLE_MUL = 0, PLAIN = 1, KARAT = 0\n\n");
                fprintf(stderr,"\t");
                for (i=0; i<NUMBER; i++)
                        fprintf(stderr,"\t%ldbit",FIRST + i*INCR);
                fprintf(stderr,"\n");
                fprintf(stderr,"\t");
                for (i=0; i<NUMBER; i++)
                        fprintf(stderr,"\t-------");
                fprintf(stderr,"\n");
		tot = 0.0;
		for (j=0; j<NUMBER; j++) {
			i=FIRST + j*INCR;
                        if (i<1000)
                                fprintf(stderr,"%ldbit:\t",i);
                        else
                                fprintf(stderr,"%ldbit:",i);
			for (k=0; k<j ; k++)
				fprintf(stderr,"\t");
			for (;  k<NUMBER; k++) {
				tt = gettime();
				for (i=0; i<ITER/4; i++) {
					zmul2(aa[j][0],bb[k][0],&c);
					zmul2(aa[j][1],bb[k][1],&c);
					zmul2(aa[j][2],bb[k][2],&c);
					zmul2(aa[j][3],bb[k][3],&c);
				}
				tcnt= gettime();
				time = tcnt-tt;
				fprintf(stderr,"\t%6.3lf",time);
				tot += time/(double) (i*(FIRST + k*INCR));
			}
			fprintf(stderr,"\n");
		}
		if (bestmul < 0 || tot < besttime) {
			bestmul = 2;
			besttime = tot;
		}
		fprintf(stderr,"\n");
	}

	if (g3) {
		fprintf(stderr,"timing SINGLE_MUL = 0, PLAIN = 0, KARAT = 1\n\n");
		fprintf(stderr,"\t");
                for (i=0; i<NUMBER; i++)
                        fprintf(stderr,"\t%ldbit",FIRST + i*INCR);
                fprintf(stderr,"\n");
                fprintf(stderr,"\t");
                for (i=0; i<NUMBER; i++)
                        fprintf(stderr,"\t-------");
                fprintf(stderr,"\n");
		tot = 0.0;
		for (j=0; j<NUMBER; j++) {
			i=FIRST + j*INCR;
                        if (i<1000)
                                fprintf(stderr,"%ldbit:\t",i);
                        else
                                fprintf(stderr,"%ldbit:",i);
			for (k=0; k<j ; k++)
				fprintf(stderr,"\t");
			for (;  k<NUMBER; k++) {
				tt = gettime();
				for (i=0; i<ITER/4; i++) {
					zmul3(aa[j][0],bb[k][0],&c);
					zmul3(aa[j][1],bb[k][1],&c);
					zmul3(aa[j][2],bb[k][2],&c);
					zmul3(aa[j][3],bb[k][3],&c);
				}
				tcnt= gettime();
				time = tcnt-tt;
				fprintf(stderr,"\t%6.3lf",time);
				tot += time/(double) (i*(FIRST + k*INCR));
			}
			fprintf(stderr,"\n");
		}
		if (bestmul < 0 || tot < besttime) {
			bestmul = 3;
			besttime = tot;
		}
		fprintf(stderr,"\n");
	}

	switch (bestmul) {
		case -1: exit(1);
		case 0: break;
		case 1: fprintf(fpar,"#define SINGLE_MUL\t1\n");
			if (g1 == 2)
				fprintf(fpar,"#define DOUBLES_LOW_HIGH\t1\n");
			break;
		case 2: fprintf(fpar,"#define PLAIN\t1\n"); break;
		case 3: fprintf(fpar,"#define KARAT\t1\n"); break;
	}
}

squ_time() {

	register int i,j,k;
	double time;

	fprintf(stderr,"squaring timings:\n\n");
	K_S_C = K_M_C;

#ifndef ALPHA_OR_ALPHA50
        if (g0) {
                fprintf(stderr,"timing SINGLE_MUL = 0, PLAIN = 0, KARAT = 0\n") ;
                fprintf(stderr,"\t");
                for (i=0; i<NUMBER; i++)
                        fprintf(stderr,"\t%ldbit",FIRST + i*INCR);
                fprintf(stderr,"\n");
                fprintf(stderr,"\t");
                for (i=0; i<NUMBER; i++)
                        fprintf(stderr,"\t-------");
                fprintf(stderr,"\n");
		fprintf(stderr,"\t");
		for (j=0; j<NUMBER; j++) {
			tt = gettime();
			for (i=0; i<ITER/4; i++) {
				zsq0(bb[j][0],&c);
				zsq0(bb[j][1],&c);
				zsq0(bb[j][2],&c);
				zsq0(bb[j][3],&c);
			}
			tcnt= gettime();
			time = tcnt-tt;
			fprintf(stderr,"\t%6.3lf",time);
		}
		fprintf(stderr,"\n");
        }

        if (g1) {
                fprintf(stderr,"timing SINGLE_MUL = 1, PLAIN = 0, KARAT = 0\n") ;
                fprintf(stderr,"\t");
                for (i=0; i<NUMBER; i++)
                        fprintf(stderr,"\t%ldbit",FIRST + i*INCR);
                fprintf(stderr,"\n");
                fprintf(stderr,"\t");
                for (i=0; i<NUMBER; i++)
                        fprintf(stderr,"\t-------");
                fprintf(stderr,"\n\t");
		for (j=0; j<NUMBER; j++) {
			zcopy(bb[j][0],&b);
			lower_radix1(b,e);
			zcopy(bb[j][1],&b);
			lower_radix1(b,e1);
			zcopy(bb[j][2],&b);
			lower_radix1(b,e2);
			zcopy(bb[j][3],&b);
			lower_radix1(b,e3);
			tt = gettime();
			if (g1==1)
				for (i=0; i<ITER/4; i++) {
					zsq1(e,&c);
					zsq1(e1,&c);
					zsq1(e2,&c);
					zsq1(e3,&c);
				}
			else
				for (i=0; i<ITER/4; i++) {
					zsq5(e,&c);
					zsq5(e1,&c);
					zsq5(e2,&c);
					zsq5(e3,&c);
				}
			tcnt= gettime();
			time = tcnt-tt;
			fprintf(stderr,"\t%6.3lf",time);
		}
		fprintf(stderr,"\n");
        }
#endif

        if (g2) {
                fprintf(stderr,"timing SINGLE_MUL = 0, PLAIN = 1, KARAT = 0\n") ;
                fprintf(stderr,"\t");
                for (i=0; i<NUMBER; i++)
                        fprintf(stderr,"\t%ldbit",FIRST + i*INCR);
                fprintf(stderr,"\n");
                fprintf(stderr,"\t");
                for (i=0; i<NUMBER; i++)
                        fprintf(stderr,"\t-------");
                fprintf(stderr,"\n");
		fprintf(stderr,"\t");
		for (j=0; j<NUMBER; j++) {
			tt = gettime();
			for (i=0; i<ITER/4; i++) {
				zsq2(bb[j][0],&c);
				zsq2(bb[j][1],&c);
				zsq2(bb[j][2],&c);
				zsq2(bb[j][3],&c);
			}
			tcnt= gettime();
			time = tcnt-tt;
			fprintf(stderr,"\t%6.3lf",time);
		}
		fprintf(stderr,"\n");
        }

        if (g3) {
                fprintf(stderr,"timing SINGLE_MUL = 0, PLAIN = 0, KARAT = 1\n") ;
                fprintf(stderr,"\t");
                for (i=0; i<NUMBER; i++)
                        fprintf(stderr,"\t%ldbit",FIRST + i*INCR);
                fprintf(stderr,"\n");
                fprintf(stderr,"\t");
                for (i=0; i<NUMBER; i++)
                        fprintf(stderr,"\t-------");
                fprintf(stderr,"\n");
		fprintf(stderr,"\t");
		for (j=0; j<NUMBER; j++) {
			tt = gettime();
			for (i=0; i<ITER/4; i++) {
				zsq3(bb[j][0],&c);
				zsq3(bb[j][1],&c);
				zsq3(bb[j][2],&c);
				zsq3(bb[j][3],&c);
			}
			tcnt= gettime();
			time = tcnt-tt;
			fprintf(stderr,"\t%6.3lf",time);
		}
		fprintf(stderr,"\n");
        }

}


mult_crov() {
	register int i,j,k;
	double time;

	switch (bestmul) {
		case -1: exit(1);
		case 0:
			fprintf(stderr,"determining KAR_MUL_CROV for  SINGLE_MUL = 0, PLAIN = 0, KARAT = 0\n");
			while (a[0] < K_M_C) {
				zsmul(a,RADIXM    ,&a);
			}
			/* now a[0] == K_M_C */
			zsadd(a,1,&b);

			tt = gettime();
			for (i=0; i< iter; i++)
				zmul0(a,b,&c);
			tcnt= gettime();
			ktime = tcnt - tt;

			K_M_C++;
			tt = gettime();
			for (i=0; i< iter; i++)
				zmul0(a,b,&c);
			tcnt= gettime();
			ntime = tcnt - tt;
			if (ktime < ntime)
				while (ktime < ntime) {
					zsdiv(a,RADIXM    ,&a);
					zsadd(a,1,&b);
					K_M_C=a[0];
					tt = gettime();
					for (i=0; i< iter; i++)
						zmul0(a,b,&c);
					tcnt= gettime();
					ktime = tcnt - tt;

					K_M_C++;
					tt = gettime();
					for (i=0; i< iter; i++)
						zmul0(a,b,&c);
					tcnt= gettime();
					ntime = tcnt - tt;
				}
			else
				while (ntime < ktime) {
					zsmul(a,  RADIXM  ,&a);
					zsadd(a,1,&b);
					K_M_C=a[0];
					tt = gettime();
					for (i=0; i< iter; i++)
						zmul0(a,b,&c);
					tcnt= gettime();
					ktime = tcnt - tt;

					K_M_C++;
					tt = gettime();
					for (i=0; i< iter; i++)
						zmul0(a,b,&c);
					tcnt= gettime();
					ntime = tcnt - tt;
				}

			break;
		case 1:
			fprintf(stderr,"determining KAR_MUL_CROV for  SINGLE_MUL = 1, PLAIN = 0, KARAT = 0\n");
			lower_radix1(a,e);
			while (e[0] < K_M_C) {
				zsmul(a,(1L<<26)-1,&a);
				lower_radix1(a,e);
			}
			/* now e[0] == K_M_C */
			zcopy(e,&f);

			tt = gettime();
			if (g1==1)
				for (i=0; i< iter; i++)
					zmul1(e,f,&c);
			else
				for (i=0; i< iter; i++)
					zmul5(e,f,&c);
			tcnt= gettime();
			ktime = tcnt - tt;
			K_M_C++;
			tt = gettime();
			if (g1==1)
				for (i=0; i< iter; i++)
					zmul1(e,f,&c);
			else
				for (i=0; i< iter; i++)
					zmul5(e,f,&c);
			tcnt= gettime();
			ntime = tcnt - tt;

			if (ktime < ntime)
				while (ktime < ntime) {
					zsdiv(a,(1L<<26)-1,&a);
					lower_radix1(a,e);
					zcopy(e,&f);
					K_M_C=e[0];
					tt = gettime();
					if (g1==1)
						for (i=0; i< iter; i++)
							zmul1(e,f,&c);
					else
						for (i=0; i< iter; i++)
							zmul5(e,f,&c);
					tcnt= gettime();
					ktime = tcnt - tt;

					K_M_C++;
					tt = gettime();
					if (g1==1)
						for (i=0; i< iter; i++)
							zmul1(e,f,&c);
					else
						for (i=0; i< iter; i++)
							zmul5(e,f,&c);
					tcnt= gettime();
					ntime = tcnt - tt;
				}
			else
				while (ntime < ktime) {
					zsmul(a,(1L<<26)-1,&a);
					lower_radix1(a,e);
					zcopy(e,&f);
					K_M_C=e[0];
					tt = gettime();
					if (g1==1)
						for (i=0; i< iter; i++)
							zmul1(e,f,&c);
					else
						for (i=0; i< iter; i++)
							zmul5(e,f,&c);
					tcnt= gettime();
					ktime = tcnt - tt;

					K_M_C++;
					tt = gettime();
					if (g1==1)
						for (i=0; i< iter; i++)
							zmul1(e,f,&c);
					else
						for (i=0; i< iter; i++)
							zmul5(e,f,&c);
					tcnt= gettime();
					ntime = tcnt - tt;
				}

			break;
		case 2:
			fprintf(stderr,"determining KAR_MUL_CROV for  SINGLE_MUL = 0, PLAIN = 1, KARAT = 0\n");
			while (a[0] < K_M_C) {
				zsmul(a,  RADIXM  ,&a);
			}
			/* now a[0] == K_M_C */
			zsadd(a,1,&b);

			tt = gettime();
			for (i=0; i< iter; i++)
				zmul2(a,b,&c);
			tcnt= gettime();
			ktime = tcnt - tt;

			K_M_C++;
			tt = gettime();
			for (i=0; i< iter; i++)
				zmul2(a,b,&c);
			tcnt= gettime();
			ntime = tcnt - tt;

			if (ktime < ntime)
				while (ktime < ntime) {
					zsdiv(a,  RADIXM  ,&a);
					zsadd(a,1,&b);
					K_M_C=a[0];
					tt = gettime();
					for (i=0; i< iter; i++)
						zmul2(a,b,&c);
					tcnt= gettime();
					ktime = tcnt - tt;

					K_M_C++;
					tt = gettime();
					for (i=0; i< iter; i++)
						zmul2(a,b,&c);
					tcnt= gettime();
					ntime = tcnt - tt;
				}
			else
				while (ntime < ktime) {
					zsmul(a,  RADIXM  ,&a);
					zsadd(a,1,&b);
					K_M_C=a[0];
					tt = gettime();
					for (i=0; i< iter; i++)
						zmul2(a,b,&c);
					tcnt= gettime();
					ktime = tcnt - tt;

					K_M_C++;
					tt = gettime();
					for (i=0; i< iter; i++)
						zmul2(a,b,&c);
					tcnt= gettime();
					ntime = tcnt - tt;
				}

			break;
		case 3:
			fprintf(stderr,"determining KAR_MUL_CROV for  SINGLE_MUL = 0, PLAIN = 0, KARAT = 1\n");
			while (a[0] < K_M_C) {
				zsmul(a,  RADIXM  ,&a);
			}
			/* now a[0] == K_M_C */
			zsadd(a,1,&b);

			tt = gettime();
			for (i=0; i< iter; i++)
				zmul3(a,b,&c);
			tcnt= gettime();
			ktime = tcnt - tt;

			K_M_C++;
			tt = gettime();
			for (i=0; i< iter; i++)
				zmul3(a,b,&c);
			tcnt= gettime();
			ntime = tcnt - tt;

			if (ktime < ntime)
				while (ktime < ntime) {
					zsdiv(a,  RADIXM  ,&a);
					zsadd(a,1,&b);
					K_M_C=a[0];
					tt = gettime();
					for (i=0; i< iter; i++)
						zmul3(a,b,&c);
					tcnt= gettime();
					ktime = tcnt - tt;

					K_M_C++;
					tt = gettime();
					for (i=0; i< iter; i++)
						zmul3(a,b,&c);
					tcnt= gettime();
					ntime = tcnt - tt;
				}
			else
				while (ntime < ktime) {
					zsmul(a,  RADIXM  ,&a);
					zsadd(a,1,&b);
					K_M_C=a[0];
					tt = gettime();
					for (i=0; i< iter; i++)
						zmul3(a,b,&c);
					tcnt= gettime();
					ktime = tcnt - tt;

					K_M_C++;
					tt = gettime();
					for (i=0; i< iter; i++)
						zmul3(a,b,&c);
					tcnt= gettime();
					ntime = tcnt - tt;
				}

			break;
	}
	fprintf(stderr,"use KAR_MUL_CROV = %ld\n",K_M_C);
	fprintf(fpar,"#define KAR_MUL_CROV\t%d\n",K_M_C);

}

kar_crov() {
	double time;
	register int i,j,k;

	/* determine K_S_C */

	switch (bestmul) {
		case 0:
			fprintf(stderr,"determining KAR_SQU_CROV for  SINGLE_MUL = 0, PLAIN = 0, KARAT = 0\n");
			K_S_C = a[0];

			tt = gettime();
			for (i=0; i< iter; i++)
				zsq0(a,&c);
			tcnt= gettime();
			ktime = tcnt - tt;

			K_S_C++;
			tt = gettime();
			for (i=0; i< iter; i++)
				zsq0(a,&c);
			tcnt= gettime();
			ntime = tcnt - tt;

			if (ktime < ntime)
				while (ktime < ntime) {
					zsdiv(a,  RADIXM  ,&a);
					K_S_C=a[0];
					tt = gettime();
					for (i=0; i< iter; i++)
						zsq0(a,&c);
					tcnt= gettime();
					ktime = tcnt - tt;

					K_S_C++;
					tt = gettime();
					for (i=0; i< iter; i++)
						zsq0(a,&c);
					tcnt= gettime();
					ntime = tcnt - tt;
				}
			else
				while (ntime < ktime) {
					zsmul(a,  RADIXM  ,&a);
					K_S_C=a[0];
					tt = gettime();
					for (i=0; i< iter; i++)
						zsq0(a,&c);
					tcnt= gettime();
					ktime = tcnt - tt;

					K_S_C++;
					tt = gettime();
					for (i=0; i< iter; i++)
						zsq0(a,&c);
					tcnt= gettime();
					ntime = tcnt - tt;
				}

			break;
		case 1:
			fprintf(stderr,"determining KAR_SQU_CROV for  SINGLE_MUL = 1, PLAIN = 0, KARAT = 0\n");
			lower_radix1(a,e);
			K_S_C = e[0];

			tt = gettime();
			if (g1==1)
				for (i=0; i< iter; i++)
					zsq1(e,&c);
			else
				for (i=0; i< iter; i++)
					zsq5(e,&c);
			tcnt= gettime();
			ktime = tcnt - tt;
			K_S_C++;
			tt = gettime();
			if (g1==1)
				for (i=0; i< iter; i++)
					zsq1(e,&c);
			else
				for (i=0; i< iter; i++)
					zsq5(e,&c);
			tcnt= gettime();
			ntime = tcnt - tt;
			if (ktime < ntime)
				while (ktime < ntime) {
					zsdiv(a,(1L<<26)-1,&a);
					lower_radix1(a,e);
					K_S_C=e[0];
					tt = gettime();
					if (g1==1)
						for (i=0; i< iter; i++)
							zsq1(e,&c);
					else
						for (i=0; i< iter; i++)
							zsq5(e,&c);
					tcnt= gettime();
					ktime = tcnt - tt;

					K_S_C++;
					tt = gettime();
					if (g1==1)
						for (i=0; i< iter; i++)
							zsq1(e,&c);
					else
						for (i=0; i< iter; i++)
							zsq5(e,&c);
					tcnt= gettime();
					ntime = tcnt - tt;
				}
			else
				while (ntime < ktime) {
					zsmul(a,(1L<<26)-1,&a);
					lower_radix1(a,e);
					K_S_C=e[0];
					tt = gettime();
					if (g1==1)
						for (i=0; i< iter; i++)
							zsq1(e,&c);
					else
						for (i=0; i< iter; i++)
							zsq5(e,&c);
					tcnt= gettime();
					ktime = tcnt - tt;

					K_S_C++;
					tt = gettime();
					if (g1==1)
						for (i=0; i< iter; i++)
							zsq1(e,&c);
					else
						for (i=0; i< iter; i++)
							zsq5(e,&c);
					tcnt= gettime();
					ntime = tcnt - tt;
				}

			break;
		case 2:
			fprintf(stderr,"determining KAR_SQU_CROV for  SINGLE_MUL = 0, PLAIN = 1, KARAT = 0\n");
			K_S_C = a[0];

			tt = gettime();
			for (i=0; i< iter; i++)
				zsq2(a,&c);
			tcnt= gettime();
			ktime = tcnt - tt;

			K_S_C++;
			tt = gettime();
			for (i=0; i< iter; i++)
				zsq2(a,&c);
			tcnt= gettime();
			ntime = tcnt - tt;
			if (ktime < ntime)
				while (ktime < ntime) {
					zsdiv(a,  RADIXM  ,&a);
					K_S_C=a[0];
					tt = gettime();
					for (i=0; i< iter; i++)
						zsq2(a,&c);
					tcnt= gettime();
					ktime = tcnt - tt;

					K_S_C++;
					tt = gettime();
					for (i=0; i< iter; i++)
						zsq2(a,&c);
					tcnt= gettime();
					ntime = tcnt - tt;
				}
			else
				while (ntime < ktime) {
					zsmul(a,  RADIXM  ,&a);
					K_S_C=a[0];
					tt = gettime();
					for (i=0; i< iter; i++)
						zsq2(a,&c);
					tcnt= gettime();
					ktime = tcnt - tt;

					K_S_C++;
					tt = gettime();
					for (i=0; i< iter; i++)
						zsq2(a,&c);
					tcnt= gettime();
					ntime = tcnt - tt;
				}

			break;
		case 3:
			fprintf(stderr,"determining KAR_SQU_CROV for  SINGLE_MUL = 0, PLAIN = 0, KARAT = 1\n");
			K_S_C= a[0];

			tt = gettime();
			for (i=0; i< iter; i++)
				zsq3(a,&c);
			tcnt= gettime();
			ktime = tcnt - tt;

			K_S_C++;
			tt = gettime();
			for (i=0; i< iter; i++)
				zsq3(a,&c);
			tcnt= gettime();
			ntime = tcnt - tt;

			if (ktime < ntime)
				while (ktime < ntime) {
					zsdiv(a,  RADIXM  ,&a);
					K_S_C=a[0];
					tt = gettime();
					for (i=0; i< iter; i++)
						zsq3(a,&c);
					tcnt= gettime();
					ktime = tcnt - tt;

					K_S_C++;
					tt = gettime();
					for (i=0; i< iter; i++)
						zsq3(a,&c);
					tcnt= gettime();
					ntime = tcnt - tt;
				}
			else
				while (ntime < ktime) {
					zsmul(a,  RADIXM  ,&a);
					K_S_C=a[0];
					tt = gettime();
					for (i=0; i< iter; i++)
						zsq3(a,&c);
					tcnt= gettime();
					ktime = tcnt - tt;

					K_S_C++;
					tt = gettime();
					for (i=0; i< iter; i++)
						zsq3(a,&c);
					tcnt= gettime();
					ntime = tcnt - tt;
				}

			break;
	}
	fprintf(stderr,"use KAR_SQU_CROV = %ld\n",K_S_C);
	fprintf(fpar,"#define KAR_SQU_CROV\t%d\n",K_S_C);

}


div_time() {
	double time;
	register int i,j,k;

	fprintf(stderr,"\ndivision timings:\n\n");

#ifndef ALPHA_OR_ALPHA50
	if (g0) {
                fprintf(stderr,"timing SINGLE_MUL = 0, PLAIN = 0, KARAT = 0\n\n");
		fprintf(stderr,"\t");
                for (i=0; i<NUMBER; i++)
                        fprintf(stderr,"\t%ldbit",FIRST + i*INCR);
                fprintf(stderr,"\n");
                fprintf(stderr,"\t");
                for (i=0; i<NUMBER; i++)
                        fprintf(stderr,"\t-------");
                fprintf(stderr,"\n");
                tot = 0.0;
                for (j=0; j<NUMBER; j++) {
                        i=FIRST + j*INCR;
                        if (i<1000)
                                fprintf(stderr,"%ldbit:\t",i);
                        else
                                fprintf(stderr,"%ldbit:",i);
                        for (k=0; k<j ; k++)
                                fprintf(stderr,"\t");
                        for (;  k<NUMBER; k++) {
                                tt = gettime();
                                for (i=0; i<ITER/4; i++) {
                                        zdiv0(bb[k][0],aa[j][0],&f,&g);
                                        zdiv0(bb[k][1],aa[j][1],&f,&g);
                                        zdiv0(bb[k][2],aa[j][2],&f,&g);
                                        zdiv0(bb[k][3],aa[j][3],&f,&g);
				}
                                tcnt= gettime();
                                time = tcnt-tt;
                                fprintf(stderr,"\t%6.3lf",time);
                        }
                        fprintf(stderr,"\n");
                }
                fprintf(stderr,"\n");
        }

        if (g1) {
                fprintf(stderr,"timing SINGLE_MUL = 1, PLAIN = 0, KARAT = 0\n\n");
		fprintf(stderr,"\t");
                for (i=0; i<NUMBER; i++)
                        fprintf(stderr,"\t%ldbit",FIRST + i*INCR);
                fprintf(stderr,"\n");
                fprintf(stderr,"\t");
                for (i=0; i<NUMBER; i++)
                        fprintf(stderr,"\t-------");
                fprintf(stderr,"\n");
                for (j=0; j<NUMBER; j++) {
			i=FIRST + j*INCR;
                        if (i<1000)
                                fprintf(stderr,"%ldbit:\t",i);
                        else
                                fprintf(stderr,"%ldbit:",i);
                        zcopy(aa[j][0],&a);
                        lower_radix1(a,e);
                        zcopy(aa[j][1],&a);
                        lower_radix1(a,e1);
                        zcopy(aa[j][2],&a);
                        lower_radix1(a,e2);
                        zcopy(aa[j][3],&a);
                        lower_radix1(a,e3);
                        for (k=0; k<j ; k++)
                                fprintf(stderr,"\t");
                        for (;  k<NUMBER; k++) {
                                lower_radix1((bb[k][0]),f);
                                lower_radix1((bb[k][1]),f1);
                                lower_radix1((bb[k][2]),f2);
                                lower_radix1((bb[k][3]),f3);
                                tt = gettime();
				if (g1==1) {
					for (i=0; i<ITER/4; i++) {
						zdiv1(f,e,&d,&g);
						zdiv1(f3,e3,&d,&g);
						zdiv1(f1,e1,&d,&g);
						zdiv1(f2,e2,&d,&g);
					}
				}
				else {
					for (i=0; i<ITER/4; i++) {
						zdiv5(f3,e3,&d,&g);
						zdiv5(f1,e1,&d,&g);
						zdiv5(f2,e2,&d,&g);
						zdiv5(f,e,&d,&g);
					}
				}
                                tcnt= gettime();
                                time = tcnt-tt;
                                fprintf(stderr,"\t%6.3lf",time);
                        }
                        fprintf(stderr,"\n");
                }
                fprintf(stderr,"\n");
        }
#endif

#ifdef ALPHA
	fprintf(stderr,"using #define ALPHA\n\n");
#endif
#ifdef ALPHA50
	fprintf(stderr,"using #define ALPHA50\n\n");
#endif

        if (g2) {
                fprintf(stderr,"timing SINGLE_MUL = 0, PLAIN = 1, KARAT = 0\n\n");
                fprintf(stderr,"\t");
                for (i=0; i<NUMBER; i++)
                        fprintf(stderr,"\t%ldbit",FIRST + i*INCR);
                fprintf(stderr,"\n");
                fprintf(stderr,"\t");
                for (i=0; i<NUMBER; i++)
                        fprintf(stderr,"\t-------");
                fprintf(stderr,"\n");
                for (j=0; j<NUMBER; j++) {
			i=FIRST + j*INCR;
                        if (i<1000)
                                fprintf(stderr,"%ldbit:\t",i);
                        else
                                fprintf(stderr,"%ldbit:",i);
                        for (k=0; k<j ; k++)
                                fprintf(stderr,"\t");
                        for (;  k<NUMBER; k++) {
                                tt = gettime();
                                for (i=0; i<ITER/4; i++) {
                                        zdiv2(bb[k][0],aa[j][0],&f,&g);
                                        zdiv2(bb[k][1],aa[j][1],&f,&g);
                                        zdiv2(bb[k][2],aa[j][2],&f,&g);
                                        zdiv2(bb[k][3],aa[j][3],&f,&g);
				}
                                tcnt= gettime();
                                time = tcnt-tt;
                                fprintf(stderr,"\t%6.3lf",time);
                        }
                        fprintf(stderr,"\n");
                }
                fprintf(stderr,"\n");
        }

        if (g3) {
                fprintf(stderr,"timing SINGLE_MUL = 0, PLAIN = 0, KARAT = 1\n\n");
                fprintf(stderr,"\t");
                for (i=0; i<NUMBER; i++)
                        fprintf(stderr,"\t%ldbit",FIRST + i*INCR);
                fprintf(stderr,"\n");
                fprintf(stderr,"\t");
                for (i=0; i<NUMBER; i++)
                        fprintf(stderr,"\t-------");
                fprintf(stderr,"\n");
                for (j=0; j<NUMBER; j++) {
                        i=FIRST + j*INCR;
                        if (i<1000)
                                fprintf(stderr,"%ldbit:\t",i);
                        else
                                fprintf(stderr,"%ldbit:",i);
                        for (k=0; k<j ; k++)
                                fprintf(stderr,"\t");
                        for (;  k<NUMBER; k++) {
                                tt = gettime();
                                for (i=0; i<ITER/4; i++) {
                                        zdiv3(bb[k][0],aa[j][0],&f,&g);
                                        zdiv3(bb[k][1],aa[j][1],&f,&g);
                                        zdiv3(bb[k][2],aa[j][2],&f,&g);
                                        zdiv3(bb[k][3],aa[j][3],&f,&g);
				}
                                tcnt= gettime();
                                time = tcnt-tt;
                                fprintf(stderr,"\t%6.3lf",time);
                        }
                        fprintf(stderr,"\n");
                }
                fprintf(stderr,"\n");
        }
}


