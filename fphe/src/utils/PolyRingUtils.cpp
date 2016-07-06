#include "PolyRingUtils.h"

#include <cassert>
#include <cmath>

using namespace std;
using namespace NTL;

void PolyRingUtils::truncate(ZZ& c, const long& logMod) {
	c = sign(c) * trunc_ZZ(c, logMod);
}

void PolyRingUtils::addPolyRing(ZZX& res, ZZX& poly1, ZZX& poly2, const ZZ& mod, const long& deg) {
	long i;
	ZZX add;
	ZZ c;
	add.SetLength(deg);
	for (i = 0; i < deg; ++i) {
		c = (coeff(poly1, i) + coeff(poly2, i)) % mod;
		SetCoeff(add, i, c);
	}
	add.normalize();
	res = add;
}

void PolyRingUtils::addPolyRing2(ZZX& res, ZZX& poly1, ZZX& poly2, const long& logMod, const long& deg) {
	long i;
	ZZX add;
	ZZ c;
	add.SetLength(deg);
	for (i = 0; i < deg; ++i) {
		c = coeff(poly1, i) + coeff(poly2, i);
		truncate(c, logMod);
		SetCoeff(add, i, c);
	}
	add.normalize();
	res = add;
}

void PolyRingUtils::subPolyRing(ZZX& res, ZZX& poly1, ZZX& poly2, const ZZ& mod, const long& deg) {
	long i;
	ZZX sub;
	ZZ c;
	sub.SetLength(deg);
	for (i = 0; i < deg; ++i) {
		c = (coeff(poly1, i) - coeff(poly2, i)) % mod;
		SetCoeff(sub, i, c);
	}
	sub.normalize();
	res = sub;
}

void PolyRingUtils::subPolyRing2(ZZX& res, ZZX& poly1, ZZX& poly2, const long& logMod, const long& deg) {
	long i;
	ZZX add;
	ZZ c;
	add.SetLength(deg);
	for (i = 0; i < deg; ++i) {
		c = coeff(poly1, i) - coeff(poly2, i);
		truncate(c, logMod);
		SetCoeff(add, i, c);
	}
	add.normalize();
	res = add;
}

void PolyRingUtils::mulPolyRing(ZZX& res, ZZX& poly1, ZZX& poly2, const ZZ& mod, const long& d) {
	long i;
	ZZ tmp;
	ZZX poly;
	mul(poly, poly1, poly2);

	for (i = 0; i < deg(poly); ++i) {
		tmp = coeff(poly, i) % mod;
		SetCoeff(poly, i, tmp);
	}

	for (i = 0; i < d; ++i) {
		tmp = (coeff(poly, i) - coeff(poly, i + d)) % mod;
		SetCoeff(res, i, tmp);
	}
	res.normalize();
}

void PolyRingUtils::mulPolyRing2(ZZX& res, ZZX& poly1, ZZX& poly2, const long& logMod, const long& d) {
	long i;
	ZZ tmp;
	ZZX poly;
	mul(poly, poly1, poly2);

	for (i = 0; i < deg(poly); ++i) {
		tmp = coeff(poly, i);
		truncate(tmp, logMod);
		SetCoeff(poly, i, tmp);
	}

	poly.normalize();

	for (i = 0; i < d; ++i) {
		tmp = coeff(poly, i) - coeff(poly, i + d);
		truncate(tmp, logMod);
		SetCoeff(res, i, tmp);
	}
	res.normalize();
}

void PolyRingUtils::mulPolyByConstantRing(ZZX& res, ZZX& poly, const ZZ& cnst, const ZZ& mod, const long& deg) {
	long i;
	ZZX mul;
	ZZ c;
	mul.SetLength(deg);
	for (i = 0; i < deg; ++i) {
		c = (coeff(poly, i) * cnst) % mod;
		SetCoeff(mul, i, c);
	}
	mul.normalize();
	res = mul;
}

void PolyRingUtils::mulPolyByConstantRing2(ZZX& res, ZZX& poly, const ZZ& cnst, const long& logMod, const long& deg) {
	long i;
	ZZX mul;
	ZZ c;
	mul.SetLength(deg);
	for (i = 0; i < deg; ++i) {
		c = coeff(poly, i) * cnst;
		truncate(c, logMod);
		SetCoeff(mul, i, c);
	}
	mul.normalize();
	res = mul;
}

void PolyRingUtils::rightShiftPolyRing(ZZX& res, ZZX& poly, const long& bits, const ZZ& mod, const long& deg) {
	long i;
	ZZX mul;
	ZZ c;
	mul.SetLength(deg);
	for (i = 0; i < deg; ++i) {
		c = (coeff(poly, i) >> bits) % mod;
		SetCoeff(mul, i, c);
	}
	mul.normalize();
	res = mul;
}

void PolyRingUtils::rightShiftPolyRing2(ZZX& res, ZZX& poly, const long& bits, const long& logMod, const long& deg) {
	long i;
	ZZX mul;
	ZZ c;
	mul.SetLength(deg);
	for (i = 0; i < deg; ++i) {
		c = coeff(poly, i) >> bits;
		truncate(c, logMod);
		SetCoeff(mul, i, c);
	}
	mul.normalize();
	res = mul;
}

void PolyRingUtils::leftShiftPolyRing(ZZX& res, ZZX& poly, const long& bits, const ZZ& mod, const long& deg) {
	long i;
	ZZX mul;
	ZZ c;
	mul.SetLength(deg);
	for (i = 0; i < deg; ++i) {
		c = (coeff(poly, i) << bits) % mod;
		SetCoeff(mul, i, c);
	}
	mul.normalize();
	res = mul;
}

void PolyRingUtils::leftShiftPolyRing2(ZZX& res, ZZX& poly, const long& bits, const long& logMod, const long& deg) {
	long i;
	ZZX mul;
	ZZ c;
	mul.SetLength(deg);
	for (i = 0; i < deg; ++i) {
		c = coeff(poly, i) << bits;
		truncate(c, logMod);
		SetCoeff(mul, i, c);
	}
	mul.normalize();
	res = mul;
}

long PolyRingUtils::mobius(long n) {
  long p,e,arity=0;
  PrimeSeq s;
  while (n!=1)
    { p=s.next();
      e=0;
      while ((n%p==0)) { n=n/p; e++; }
      if (e>1) { return 0; }
      if (e!=0) { arity^=1; }
    }
  if (arity==0) { return 1; }
  return -1;
}

ZZX PolyRingUtils::Cyclotomic(long N) {
  ZZX Num,Den,G,F;
  set(Num); set(Den);
  long m,d;
  for (d=1; d<=N; d++)
    { if ((N%d)==0)
         { clear(G);
           SetCoeff(G,N/d,1); SetCoeff(G,0,-1);
           m=PolyRingUtils::mobius(d);
           if (m==1)       { Num*=G; }
           else if (m==-1) { Den*=G; }
         }
    }
  F=Num/Den;
  return F;
}

void PolyRingUtils::sampleGaussian(ZZX &res, long d, double stdev) {
  static double const Pi=4.0*atan(1.0); // Pi=3.1415..
  static long const bignum = 0xfffffff;
  // THREADS: C++11 guarantees these are initialized only once

  if (d<=0) d=deg(res)+1; if (d<=0) return;
  res.SetMaxLength(d); // allocate space for degree-(n-1) polynomial
  for (long i=0; i<d; i++) SetCoeff(res, i, ZZ::zero());

  // Uses the Box-Muller method to get two Normal(0,stdev^2) variables
  for (long i=0; i<d; i+=2) {
    double r1 = (1+RandomBnd(bignum))/((double)bignum+1);
    double r2 = (1+RandomBnd(bignum))/((double)bignum+1);
    double theta=2*Pi*r1;
    double rr= sqrt(-2.0*log(r2))*stdev;

    assert(rr < 8*stdev); // sanity-check, no more than 8 standard deviations

    // Generate two Gaussians RV's, rounded to integers
    long x = (long) floor(rr*cos(theta) +0.5);
    SetCoeff(res, i, x);
    if (i+1 < d) {
      x = (long) floor(rr*sin(theta) +0.5);
      SetCoeff(res, i+1, x);
    }
  }
  res.normalize(); // need to call this after we work on the coeffs
}

void PolyRingUtils::sampleUniform(ZZX& res, ZZ& B, long d) {

	if (d<=0) d=deg(res)+1; if (d<=0) return;
	if (B <= 0) {
		clear(res);
		return;
	}

	res.SetMaxLength(d); // allocate space for degree-(n-1) polynomial

	ZZ UB, tmp;

	UB =  2*B + 1;

	for (long i = 0; i < d; i++) {
		RandomBnd(tmp, UB);
		tmp -= B;
		SetCoeff(res, i, tmp);
	}
	res.normalize();
}

void PolyRingUtils::sampleUniform2(ZZX& res, long& logB, long d) {

	if (d<=0) d=deg(res)+1; if (d<=0) return;

	res.SetMaxLength(d); // allocate space for degree-(n-1) polynomial

	ZZ tmp;

	for (long i = 0; i < d; i++) {
		RandomBits(tmp, logB);
		SetCoeff(res, i, tmp);
	}
	res.normalize();
}
