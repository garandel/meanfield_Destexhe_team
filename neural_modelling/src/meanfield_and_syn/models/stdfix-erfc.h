/*
 * Copyright (c) 2017-2019 The University of Manchester
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

//! \file
//! \brief complementary error function in fix number arithmetics

#include <debug.h>
//#include <math.h>
//#include "../../../src/common/maths-util.h"
#include <stdfix-exp.h>
//#include <polynomial.h>
//#include <stdfix-full-iso.h>

/* origin: FreeBSD /usr/src/lib/msun/src/s_erf.c */
/*
 * ====================================================
 * Copyright (C) 1993 by Sun Microsystems, Inc. All rights reserved.
 *
 * Developed at SunPro, a Sun Microsystems, Inc. business.
 * Permission to use, copy, modify, and distribute this
 * software is freely granted, provided that this notice
 * is preserved.
 * ====================================================
 */
/* double erf(double x)
 * double erfc(double x)
 *                           x
 *                    2      |\
 *     erf(x)  =  ---------  | exp(-t*t)dt
 *                 sqrt(pi) \|
 *                           0
 *
 *     erfc(x) =  1-erf(x)
 *  Note that
 *              erf(-x) = -erf(x)
 *              erfc(-x) = 2 - erfc(x)
 *
 * Method:
 *      1. For |x| in [0, 0.84375]
 *          erf(x)  = x + x*R(x^2)
 *          erfc(x) = 1 - erf(x)           if x in [-.84375,0.25]
 *                  = 0.5 + ((0.5-x)-x*R)  if x in [0.25,0.84375]
 *         where R = P/Q where P is an odd poly of degree 8 and
 *         Q is an odd poly of degree 10.
 *                                               -57.90
 *                      | R - (erf(x)-x)/x | <= 2
 *
 *
 *         Remark. The formula is derived by noting
 *          erf(x) = (2/sqrt(pi))*(x - x^3/3 + x^5/10 - x^7/42 + ....)
 *         and that
 *          2/sqrt(pi) = 1.128379167095512573896158903121545171688
 *         is close to one. The interval is chosen because the fix
 *         point of erf(x) is near 0.6174 (i.e., erf(x)=x when x is
 *         near 0.6174), and by some experiment, 0.84375 is chosen to
 *         guarantee the error is less than one ulp for erf.
 *
 *      2. For |x| in [0.84375,1.25], let s = |x| - 1, and
 *         c = 0.84506291151 rounded to single (24 bits)
 *              erf(x)  = sign(x) * (c  + P1(s)/Q1(s))
 *              erfc(x) = (1-c)  - P1(s)/Q1(s) if x > 0
 *                        1+(c+P1(s)/Q1(s))    if x < 0
 *              |P1/Q1 - (erf(|x|)-c)| <= 2**-59.06
 *         Remark: here we use the taylor series expansion at x=1.
 *              erf(1+s) = erf(1) + s*Poly(s)
 *                       = 0.845.. + P1(s)/Q1(s)
 *         That is, we use rational approximation to approximate
 *                      erf(1+s) - (c = (single)0.84506291151)
 *         Note that |P1/Q1|< 0.078 for x in [0.84375,1.25]
 *         where
 *              P1(s) = degree 6 poly in s
 *              Q1(s) = degree 6 poly in s
 *
 *      3. For x in [1.25,1/0.35(~2.857143)],
 *              erfc(x) = (1/x)*exp(-x*x-0.5625+R1/S1)
 *              erf(x)  = 1 - erfc(x)
 *         where
 *              R1(z) = degree 7 poly in z, (z=1/x^2)
 *              S1(z) = degree 8 poly in z
 *
 *      4. For x in [1/0.35,28]
 *              erfc(x) = (1/x)*exp(-x*x-0.5625+R2/S2) if x > 0
 *                      = 2.0 - (1/x)*exp(-x*x-0.5625+R2/S2) if -6<x<0
 *                      = 2.0 - tiny            (if x <= -6)
 *              erf(x)  = sign(x)*(1.0 - erfc(x)) if x < 6, else
 *              erf(x)  = sign(x)*(1.0 - tiny)
 *         where
 *              R2(z) = degree 6 poly in z, (z=1/x^2)
 *              S2(z) = degree 7 poly in z
 *
 *      Note1:
 *         To compute exp(-x*x-0.5625+R/S), let s be a single
 *         precision number and s := x; then
 *              -x*x = -s*s + (s-x)*(s+x)
 *              exp(-x*x-0.5626+R/S) =
 *                      exp(-s*s-0.5625)*exp((s-x)*(s+x)+R/S);
 *      Note2:
 *         Here 4 and 5 make use of the asymptotic series
 *                        exp(-x*x)
 *              erfc(x) ~ ---------- * ( 1 + Poly(1/x^2) )
 *                        x*sqrt(pi)
 *         We use rational approximation to approximate
 *              g(s)=f(1/x^2) = log(erfc(x)*x) - x*x + 0.5625
 *         Here is the error bound for R1/S1 and R2/S2
 *              |R1/S1 - f(x)|  < 2**(-62.57)
 *              |R2/S2 - f(x)|  < 2**(-61.52)
 *
 *      5. For inf > x >= 28
 *              erf(x)  = sign(x) *(1 - tiny)  (raise inexact)
 *              erfc(x) = tiny*tiny (raise underflow) if x > 0
 *                      = 2 - tiny if x<0
 *
 *      7. Special case:
 *              erf(0)  = 0, erf(inf)  = 1, erf(-inf) = -1,
 *              erfc(0) = 1, erfc(inf) = 0, erfc(-inf) = 2,
 *              erfc/erf(NaN) is NaN
 */

//#include "libm.h"
/* Set the less significant 32 bits of a double from an int.  */
#define SET_LOW_WORD(d,lo)                        \
do {                                              \
  union {REAL f; uint64_t i;} __u;              \
  __u.f = (d);                                    \
  __u.i &= 0xffffffff00000000ull;                 \
  __u.i |= (uint32_t)(lo);                        \
  (d) = __u.f;                                    \
} while (0)

/* Get the more significant 32 bit int from a double.  */
#define GET_HIGH_WORD(hi,d)                       \
do {                                              \
  union {REAL f; uint64_t i;} __u;              \
  __u.f = (d);                                    \
  (hi) = __u.i >> 32;                             \
} while (0)

static const double
erx  = 8.45062911510467529297e-01;

//erx = 27691, /* 0x3FEB0AC1, 0x60000000 */
static const s1615
/*
 * Coefficients for approximation to  erf on [0,0.84375]
 */
efx8 = 33654, /* 0x3FF06EBA, 0x8214DB69 */
pp0 = 4207, /* 0x3FC06EBA, 0x8214DB68 */
pp1 = -10651, /* 0xBFD4CD7D, 0x691CB913 */
pp2 = -933, /* 0xBF9D2A51, 0xDBD7194F */
pp3 = -189, /* 0xBF77A291, 0x236668E4 */
pp4 = -1, /* 0xBEF8EAD6, 0x120016AC */
qq1 = 13039, /* 0x3FD97779, 0xCDDADC09 */
qq2 = 2131, /* 0x3FB0A54C, 0x5536CEBA */
qq3 = 167, /* 0x3F74D022, 0xC4D36B0F */
qq4 = 4, /* 0x3F215DC9, 0x221C1A10 */
qq5 = 0, /* 0xBED09C43, 0x42A26120 */
/*
 * Coefficients for approximation to  erf  in [0.84375,1.25]
 */
pa0 = -77, /* 0xBF6359B8, 0xBEF77538 */
pa1 = 13594, /* 0x3FDA8D00, 0xAD92B34D */
pa2 = -12197, /* 0xBFD7D240, 0xFBB8C3F1 */
pa3 = 10432, /* 0x3FD45FCA, 0x805120E4 */
pa4 = -3634, /* 0xBFBC6398, 0x3D3E28EC */
pa5 = 1163, /* 0x3FA22A36, 0x599795EB */
pa6 = -71, /* 0xBF61BF38, 0x0A96073F */
qa1 = 3487, /* 0x3FBB3E66, 0x18EEE323 */
qa2 = 17708, /* 0x3FE14AF0, 0x92EB6F33 */
qa3 = 2354, /* 0x3FB2635C, 0xD99FE9A7 */
qa4 = 4134, /* 0x3FC02660, 0xE763351F */
qa5 = 447, /* 0x3F8BEDC2, 0x6B51DD1C */
qa6 = 393, /* 0x3F888B54, 0x5735151D */

//static const double
/*
 * Coefficients for approximation to  erfc in [1.25,1/0.35]
 */
ra0 = -323, /* 0xBF843412, 0x600D6435 */
ra1 = -22736, /* 0xBFE63416, 0xE4BA7360 */
ra2 = -345985, /* 0xC0251E04, 0x41B0E726 */
ra3 = -2043915, /* 0xC04F300A, 0xE4CBA38D */
ra4 = -5321414, /* 0xC0644CB1, 0x84282266 */
ra5 = -6049140, /* 0xC067135C, 0xEBCCABB2 */
ra6 = -2663627, /* 0xC0545265, 0x57E4D2F2 */
ra7 = -321596, /* 0xC023A0EF, 0xC69AC25C */
sa1 = 643933, /* 0x4033A6B9, 0xBD707687 */
sa2 = 4510769, /* 0x4061350C, 0x526AE721 */
sa3 = 14239855, /* 0x407B290D, 0xD58A1A71 */
sa4 = 21148050, /* 0x40842B19, 0x21EC2868 */
sa5 = 14057739, /* 0x407AD021, 0x57700314 */
sa6 = 3559752, /* 0x405B28A3, 0xEE48AE2C */
sa7 = 215294, /* 0x401A47EF, 0x8E484A93 */
sa8 = -1980, /* 0xBFAEEFF2, 0xEE749A62 */
/*
 * Coefficients for approximation to  erfc in [1/.35,28]
 */
rb0 = -323, /* 0xBF843412, 0x39E86F4A */
rb1 = -26191, /* 0xBFE993BA, 0x70C285DE */
rb2 = -581893, /* 0xC031C209, 0x555F995A */
rb3 = -5263733, /* 0xC064145D, 0x43C5ED98 */
rb4 = -20891777, /* 0xC083EC88, 0x1375F228 */
rb5 = -33590317, /* 0xC0900461, 0x6A2E5992 */
rb6 = -15843957, /* 0xC07E384E, 0x9BDC383F */
sb1 = 994118, /* 0x403E568B, 0x261D5190 */
sb2 = 10675569, /* 0x40745CAE, 0x221B9F0A */
sb3 = 50355555, /* 0x409802EB, 0x189D5118 */
sb4 = 104852954, /* 0x40A8FFB7, 0x688C246A */
sb5 = 83658356, /* 0x40A3F219, 0xCEDF3BE6 */
sb6 = 15549351, /* 0x407DA874, 0xE79FE763 */
sb7 = -735345; /* 0xC03670E2, 0x42712D62 */

static union {
    s1615 as_s1615;
    REAL as_real;
} number;


/*
static inline REAL abs_changer(REAL x)
{
    
    number.as_real = x;//  *exc_syn_values;//
    s1615 x_s1615 = number.as_s1615; 
    
    number.as_s1615 = absk(x_s1615);
    input_t result = number.as_real;

    return result;
    
        
} 
*/

 //fait gagner vite fait 20 bytes
static accum fonc_abs(s1615 x){
    if(x<0){
        x=-1*x;
    }
    
    return x;
};

static double erfc1(REAL x)
{
	s1615 s,P,Q;
    
    number.as_real = x;
    s1615 x_s1615 = number.as_s1615;
    
	s = fonc_abs(x_s1615) - 32768;//1; //fabs(x) - 1; 

    
	P = pa0+s*(pa1+s*(pa2+s*(pa3+s*(pa4+s*(pa5+s*pa6)))));
	Q = 1+s*(qa1+s*(qa2+s*(qa3+s*(qa4+s*(qa5+s*qa6)))));
    /*
    number.as_s1615 = P;
    input_t P_real = number.as_real;
    number.as_s1615 = Q;
    input_t Q_real = number.as_real;
    */

	return 1 - erx - P/Q;
}

static double erfc2(uint32_t ix, REAL x)
{
	s1615 R,S;
	REAL z,s;

	if (ix < 0x3ff40000)  /* |x| < 1.25 */
		return erfc1(x);

	x = fonc_abs(x); //fabs(x);
	s = 1/(x*x);
    
    number.as_real = s;
    s1615 s_s1615 = number.as_s1615;
    
	if (ix < 0x4006db6d) {  /* |x| < 1/.35 ~ 2.85714 */
		R = ra0+s_s1615*(ra1+s_s1615*(ra2+s_s1615*(ra3+s_s1615*(ra4+s_s1615*(
		     ra5+s_s1615*(ra6+s_s1615*ra7))))));
		S = 1.0+s_s1615*(sa1+s_s1615*(sa2+s_s1615*(sa3+s_s1615*(sa4+s_s1615*(
		     sa5+s_s1615*(sa6+s_s1615*(sa7+s_s1615*sa8)))))));
	} else {                /* |x| > 1/.35 */
		R = rb0+s_s1615*(rb1+s_s1615*(rb2+s_s1615*(rb3+s_s1615*(rb4+s_s1615*(
		     rb5+s_s1615*rb6)))));
		S = 1.0+s_s1615*(sb1+s_s1615*(sb2+s_s1615*(sb3+s_s1615*(sb4+s_s1615*(
		     sb5+s_s1615*(sb6+s_s1615*sb7))))));
	}
	z = x;
	SET_LOW_WORD(z,0);
	return expk(-z*z-0.5625)*expk((z-x)*(z+x)+R/S)/x; //EXP
}
/*
double erf(double x)
{
	REAL r,s,z,y;
	uint32_t ix;
	int sign;

	GET_HIGH_WORD(ix, x);
	sign = ix>>31;
	ix &= 0x7fffffff;
	if (ix >= 0x7ff00000) {
		// erf(nan)=nan, erf(+-inf)=+-1 
		return 1-2*sign + 1/x;
	}
	if (ix < 0x3feb0000) {  /* |x| < 0.84375 */
//		if (ix < 0x3e300000) {  /* |x| < 2**-28 */
			/* avoid underflow */
//			return 0.125*(8*x + efx8*x);
//		}
//		z = x*x;
//		r = pp0+z*(pp1+z*(pp2+z*(pp3+z*pp4)));
//		s = 1.0+z*(qq1+z*(qq2+z*(qq3+z*(qq4+z*qq5))));
//		y = r/s;
//		return x + x*y;
//	}
//	if (ix < 0x40180000)  /* 0.84375 <= |x| < 6 */
//		y = 1 - erfc2(ix,x);
//	else
//		y = 1 - 0x1p-1022;
//	return sign ? -y : y;
//}

double erfc(REAL x)
{
	REAL r,s,z,y;
	uint32_t ix;
	int sign;

	GET_HIGH_WORD(ix, x);
	sign = ix>>31;
	ix &= 0x7fffffff;
	if (ix >= 0x7ff00000) {
		/* erfc(nan)=nan, erfc(+-inf)=0,2 */
		return 2*sign + 1/x;
	}
	if (ix < 0x3feb0000) {  /* |x| < 0.84375 */
		if (ix < 0x3c700000)  /* |x| < 2**-56 */
			return 1.0 - x;
		z = x*x;
		r = pp0+z*(pp1+z*(pp2+z*(pp3+z*pp4)));
		s = 1.0+z*(qq1+z*(qq2+z*(qq3+z*(qq4+z*qq5))));
		y = r/s;
		if (sign || ix < 0x3fd00000) {  /* x < 1/4 */
			return 1.0 - (x+x*y);
		}
		return 0.5 - (x - 0.5 + x*y);
	}
	if (ix < 0x403c0000) {  /* 0.84375 <= |x| < 28 */
		return sign ? 2 - erfc2(ix,x) : erfc2(ix,x);
	}


	return sign ? 2 - 0x1p-1022 : 0x1p-1022*0x1p-1022;
}



