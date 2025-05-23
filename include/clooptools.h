#ifndef FTYPES_H
#define FTYPES_H

#if 0
#define FORTRAN(s) s
#else
#define FORTRAN(s) s##_
#endif

#ifndef CQUADSIZE
#define CQUADSIZE QUADSIZE
#endif

#if QUAD && CQUADSIZE == 16

#include <quadmath.h>

#define _prec(f) f##q
typedef __float128 RealType;
typedef __float128 REAL;
#ifndef __cplusplus
typedef __complex128 ComplexType;
#endif

#define ToReal(r) (RealType)(r)
#define ToREAL(r) (REAL)(r)

#elif QUAD && CQUADSIZE == 10

#define _prec(f) f##l
#define RealType long double

#if QUADSIZE == CQUADSIZE

typedef RealType REAL;

#define ToReal(r) (r)
#define ToREAL(r) (r)

#else

#pragma pack(push, 1)
typedef struct {
  unsigned long long frac;
  unsigned short exp;
} REAL10;
typedef struct {
  char zero[6];
  unsigned long long frac;
  unsigned short exp;
} REAL16;
#pragma pack(pop)

typedef union {
  long double r10;
  REAL10 i10;
  REAL16 i16;
  unsigned long long i8[2];
} REAL;

static inline REAL ToREAL(const RealType r) {
  REAL n;
  n.i8[0] = 0;
  n.i16.frac = ((REAL *)&r)->i10.frac << 1;
  n.i16.exp = ((REAL *)&r)->i10.exp;
  return n;
}

static inline RealType ToReal(const REAL r) {
  REAL n;
  const long long z = r.i16.frac | (r.i16.exp & 0x7fff);
  n.i10.frac = (r.i16.frac >> 1) | ((z | -z) & 0x8000000000000000LL);
  n.i10.exp = r.i16.exp;
  return n.r10;
}

static inline void ToRealArray(RealType *out, const REAL *in, const int n) {
  int i;
  for( i = 0; i < n; ++i ) out[i] = ToReal(in[i]);
}

static inline void ToREALArray(REAL *out, const RealType *in, const int n) {
  int i;
  for( i = 0; i < n; ++i ) out[i] = ToREAL(in[i]);
}

#endif

#else

#define _prec(f) f
#define RealType double
typedef double REAL;

#define ToReal(r) (r)
#define ToREAL(r) (r)

#endif

typedef int INTEGER;
typedef const INTEGER CINTEGER;
typedef long long int INTEGER8;
typedef const REAL CREAL;
typedef struct { REAL re, im; } COMPLEX;
typedef const COMPLEX CCOMPLEX;
typedef char CHARACTER;
typedef const CHARACTER CCHARACTER;

#ifdef __cplusplus

#include <complex>
typedef std::complex<RealType> ComplexType;
#define ToComplex(c) ComplexType(ToReal((c).re), ToReal((c).im))
#define ToComplex2(r,i) ComplexType(r, i)
#define Re(x) std::real(x)
#define Im(x) std::imag(x)
#define Conjugate(x) std::conj(x)

#elif __STDC_VERSION__ >= 199901L

#include <complex.h>
#ifdef RealType
typedef RealType complex ComplexType;
#endif
#define ToComplex(c) (ToReal((c).re) + I*ToReal((c).im))
#define ToComplex2(r,i) (r + I*(i))
#define Re(c) _prec(creal)(c)
#define Im(c) _prec(cimag)(c)
#define Conjugate(c) _prec(conj)(c)

#else

#ifdef RealType
typedef struct { RealType re, im; } ComplexType;
#endif
#define ToComplex(c) (ComplexType){ToReal((c).re), ToReal((c).im)}
#define ToComplex2(r,i) (ComplexType){r, i}
#define Re(x) (x).re
#define Im(x) (x).im
#define Conjugate(x) (ComplexType){(x).re, -(x).im}

#endif

typedef const RealType cRealType;
typedef const ComplexType cComplexType;

#endif

/*
	clooptools.h
		the C/C++ header file with all definitions for LoopTools
		this file is part of LoopTools
		last modified 24 Mar 22 th
*/


#ifndef CLOOPTOOLS_H
#define CLOOPTOOLS_H

#define aa0 0
#define aa00 3
#define Naa 6

#define bb0 0
#define bb1 3
#define bb00 6
#define bb11 9
#define bb001 12
#define bb111 15
#define dbb0 18
#define dbb1 21
#define dbb00 24
#define dbb11 27
#define dbb001 30
#define Nbb 33

#define cc0 0
#define cc1 3
#define cc2 6
#define cc00 9
#define cc11 12
#define cc12 15
#define cc22 18
#define cc001 21
#define cc002 24
#define cc111 27
#define cc112 30
#define cc122 33
#define cc222 36
#define cc0000 39
#define cc0011 42
#define cc0012 45
#define cc0022 48
#define cc1111 51
#define cc1112 54
#define cc1122 57
#define cc1222 60
#define cc2222 63
#define Ncc 66

#define dd0 0
#define dd1 3
#define dd2 6
#define dd3 9
#define dd00 12
#define dd11 15
#define dd12 18
#define dd13 21
#define dd22 24
#define dd23 27
#define dd33 30
#define dd001 33
#define dd002 36
#define dd003 39
#define dd111 42
#define dd112 45
#define dd113 48
#define dd122 51
#define dd123 54
#define dd133 57
#define dd222 60
#define dd223 63
#define dd233 66
#define dd333 69
#define dd0000 72
#define dd0011 75
#define dd0012 78
#define dd0013 81
#define dd0022 84
#define dd0023 87
#define dd0033 90
#define dd1111 93
#define dd1112 96
#define dd1113 99
#define dd1122 102
#define dd1123 105
#define dd1133 108
#define dd1222 111
#define dd1223 114
#define dd1233 117
#define dd1333 120
#define dd2222 123
#define dd2223 126
#define dd2233 129
#define dd2333 132
#define dd3333 135
#define dd00001 138
#define dd00002 141
#define dd00003 144
#define dd00111 147
#define dd00112 150
#define dd00113 153
#define dd00122 156
#define dd00123 159
#define dd00133 162
#define dd00222 165
#define dd00223 168
#define dd00233 171
#define dd00333 174
#define dd11111 177
#define dd11112 180
#define dd11113 183
#define dd11122 186
#define dd11123 189
#define dd11133 192
#define dd11222 195
#define dd11223 198
#define dd11233 201
#define dd11333 204
#define dd12222 207
#define dd12223 210
#define dd12233 213
#define dd12333 216
#define dd13333 219
#define dd22222 222
#define dd22223 225
#define dd22233 228
#define dd22333 231
#define dd23333 234
#define dd33333 237
#define Ndd 240

#define ee0 0
#define ee1 3
#define ee2 6
#define ee3 9
#define ee4 12
#define ee00 15
#define ee11 18
#define ee12 21
#define ee13 24
#define ee14 27
#define ee22 30
#define ee23 33
#define ee24 36
#define ee33 39
#define ee34 42
#define ee44 45
#define ee001 48
#define ee002 51
#define ee003 54
#define ee004 57
#define ee111 60
#define ee112 63
#define ee113 66
#define ee114 69
#define ee122 72
#define ee123 75
#define ee124 78
#define ee133 81
#define ee134 84
#define ee144 87
#define ee222 90
#define ee223 93
#define ee224 96
#define ee233 99
#define ee234 102
#define ee244 105
#define ee333 108
#define ee334 111
#define ee344 114
#define ee444 117
#define ee0000 120
#define ee0011 123
#define ee0012 126
#define ee0013 129
#define ee0014 132
#define ee0022 135
#define ee0023 138
#define ee0024 141
#define ee0033 144
#define ee0034 147
#define ee0044 150
#define ee1111 153
#define ee1112 156
#define ee1113 159
#define ee1114 162
#define ee1122 165
#define ee1123 168
#define ee1124 171
#define ee1133 174
#define ee1134 177
#define ee1144 180
#define ee1222 183
#define ee1223 186
#define ee1224 189
#define ee1233 192
#define ee1234 195
#define ee1244 198
#define ee1333 201
#define ee1334 204
#define ee1344 207
#define ee1444 210
#define ee2222 213
#define ee2223 216
#define ee2224 219
#define ee2233 222
#define ee2234 225
#define ee2244 228
#define ee2333 231
#define ee2334 234
#define ee2344 237
#define ee2444 240
#define ee3333 243
#define ee3334 246
#define ee3344 249
#define ee3444 252
#define ee4444 255
#define Nee 258

enum {
  KeyA0 = 1,
  KeyBget = 1<<2,
  KeyC0 = 1<<4,
  KeyD0 = 1<<6,
  KeyE0 = 1<<8,
  KeyEget = 1<<10,
  KeyEgetC = 1<<12,
  KeyAll = KeyA0 + KeyBget + KeyC0 + KeyD0 + KeyE0 + KeyEget + KeyEgetC
};

enum {
  DebugA = 1,
  DebugB = 1<<1,
  DebugC = 1<<2,
  DebugD = 1<<3,
  DebugE = 1<<4,
  DebugAll = DebugA + DebugB + DebugC + DebugD + DebugE
};

typedef long long int memindex;

/****************************************************************/

#ifdef __cplusplus
extern "C" {
#endif

#define CACHEPTR(n,i) &FORTRAN(ltvars).cache[n][i]

#define EPSINDEX(i) i+FORTRAN(ltvars).epsi

enum { ltvars_ncache = 10 };

extern struct {		/* MUST match common block ltvars in lt.h! */
  COMPLEX cache[ltvars_ncache][2];
  INTEGER8 savedptr[ltvars_ncache][2];
  REAL maxdev;
  INTEGER epsi, warndigits, errdigits;
  INTEGER serial, versionkey;
  INTEGER debugkey, debugfrom, debugto;
} FORTRAN(ltvars);

#define _lt_Cr_(v) cRealType v
#define _lt_Cc_(v) cComplexType v
#define _lt_Fr_(v) CREAL *v
#define _lt_Fc_(v) CCOMPLEX *v
#define _lt_Id_(v) v

#if QUAD && CQUADSIZE != QUADSIZE
#define _lt_CFr_(v) v##_ = ToREAL(v)
#define _lt_CFc_(v) v##_ = {ToREAL(Re(v)), ToREAL(Im(v))}
#define _lt_Frp_(v) &v##_
#define _lt_Fcp_(v) &v##_
#define _lt_Fap_(v) v##_
#define _lt_Frd_(f) CREAL f(_lt_CFr_);
#define _lt_Fcd_(f) CCOMPLEX f(_lt_CFc_);
#define _lt_Fad_(v,n) COMPLEX v##_[n];
#define _lt_Fax_(v,n) ToRealArray((RealType *)v, (REAL *)v##_, 2*n);
#else
#define _lt_Frp_(v) &v
#define _lt_Fcp_(v) (CCOMPLEX *)&v
#define _lt_Fap_(v) (COMPLEX *)v
#define _lt_Frd_(f)
#define _lt_Fcd_(f)
#define _lt_Fad_(v,n)
#define _lt_Fax_(v,n)
#endif

/****************************************************************/

#define AARGS(t) t(m)

static inline memindex Aget(AARGS(_lt_Cr_)) {
  _lt_Frd_(AARGS)
  extern memindex FORTRAN(aget)(AARGS(_lt_Fr_));
  return FORTRAN(aget)(AARGS(_lt_Frp_));
}

static inline memindex AgetC(AARGS(_lt_Cc_)) {
  _lt_Fcd_(AARGS)
  extern memindex FORTRAN(agetc)(AARGS(_lt_Fc_));
  return FORTRAN(agetc)(AARGS(_lt_Fcp_));
}

static inline void Aput(ComplexType *res, AARGS(_lt_Cr_)) {
  _lt_Frd_(AARGS)
  _lt_Fad_(res, Naa)
  extern void FORTRAN(aput)(COMPLEX *res, AARGS(_lt_Fr_));
  FORTRAN(aput)(_lt_Fap_(res), AARGS(_lt_Frp_));
  _lt_Fax_(res, Naa)
}

static inline void AputC(ComplexType *res, AARGS(_lt_Cc_)) {
  _lt_Fcd_(AARGS)
  _lt_Fad_(res, Naa)
  extern void FORTRAN(aputc)(COMPLEX *res, AARGS(_lt_Fc_));
  FORTRAN(aputc)(_lt_Fap_(res), AARGS(_lt_Fcp_));
  _lt_Fax_(res, Naa)
}

static inline void Aputnocache(ComplexType *res, AARGS(_lt_Cr_)) {
  _lt_Frd_(AARGS)
  _lt_Fad_(res, Naa)
  extern void FORTRAN(aputnocache)(COMPLEX *res, AARGS(_lt_Fr_));
  FORTRAN(aputnocache)(_lt_Fap_(res), AARGS(_lt_Frp_));
  _lt_Fax_(res, Naa)
}

static inline void AputnocacheC(ComplexType *res, AARGS(_lt_Cc_)) {
  _lt_Fcd_(AARGS)
  _lt_Fad_(res, Naa)
  extern void FORTRAN(aputnocachec)(COMPLEX *res, AARGS(_lt_Fc_));
  FORTRAN(aputnocachec)(_lt_Fap_(res), AARGS(_lt_Fcp_));
  _lt_Fax_(res, Naa)
}

static inline COMPLEX *Acache(const memindex integral)
  { return CACHEPTR(0,integral); }

static inline COMPLEX *AcacheC(const memindex integral)
  { return CACHEPTR(1,integral); }

static inline ComplexType Aval(const int i, const memindex integral)
  { return ToComplex(Acache(integral)[i]); }

static inline ComplexType AvalC(const int i, const memindex integral)
  { return ToComplex(AcacheC(integral)[i]); }

static inline ComplexType A0i(const int i, AARGS(_lt_Cr_))
  { return Aval(EPSINDEX(i), Aget(AARGS(_lt_Id_))); }

static inline ComplexType A0iC(const int i, AARGS(_lt_Cc_))
  { return AvalC(EPSINDEX(i), AgetC(AARGS(_lt_Id_))); }

static inline ComplexType A0(AARGS(_lt_Cr_))
  { return A0i(aa0, AARGS(_lt_Id_)); }
static inline ComplexType A00(AARGS(_lt_Cr_))
  { return A0i(aa00, AARGS(_lt_Id_)); }

static inline ComplexType A0C(AARGS(_lt_Cc_))
  { return A0iC(aa0, AARGS(_lt_Id_)); }
static inline ComplexType A00C(AARGS(_lt_Cc_))
  { return A0iC(aa00, AARGS(_lt_Id_)); }

/****************************************************************/

#define BARGS(t) t(p), t(m1), t(m2)

static inline memindex Bget(BARGS(_lt_Cr_)) {
  _lt_Frd_(BARGS)
  extern memindex FORTRAN(bget)(BARGS(_lt_Fr_));
  return FORTRAN(bget)(BARGS(_lt_Frp_));
}

static inline memindex BgetC(BARGS(_lt_Cc_)) {
  _lt_Fcd_(BARGS)
  extern memindex FORTRAN(bgetc)(BARGS(_lt_Fc_));
  return FORTRAN(bgetc)(BARGS(_lt_Fcp_));
}

static inline void Bput(ComplexType *res, BARGS(_lt_Cr_)) {
  _lt_Frd_(BARGS)
  _lt_Fad_(res, Nbb)
  extern void FORTRAN(bput)(COMPLEX *res, BARGS(_lt_Fr_));
  FORTRAN(bput)(_lt_Fap_(res), BARGS(_lt_Frp_));
  _lt_Fax_(res, Nbb)
}

static inline void BputC(ComplexType *res, BARGS(_lt_Cc_)) {
  _lt_Fcd_(BARGS)
  _lt_Fad_(res, Nbb)
  extern void FORTRAN(bputc)(COMPLEX *res, BARGS(_lt_Fc_));
  FORTRAN(bputc)(_lt_Fap_(res), BARGS(_lt_Fcp_));
  _lt_Fax_(res, Nbb)
}

static inline void Bputnocache(ComplexType *res, BARGS(_lt_Cr_)) {
  _lt_Frd_(BARGS)
  _lt_Fad_(res, Nbb)
  extern void FORTRAN(bputnocache)(COMPLEX *res, BARGS(_lt_Fr_));
  FORTRAN(bputnocache)(_lt_Fap_(res), BARGS(_lt_Frp_));
  _lt_Fax_(res, Nbb)
}

static inline void BputnocacheC(ComplexType *res, BARGS(_lt_Cc_)) {
  _lt_Fcd_(BARGS)
  _lt_Fad_(res, Nbb)
  extern void FORTRAN(bputnocachec)(COMPLEX *res, BARGS(_lt_Fc_));
  FORTRAN(bputnocachec)(_lt_Fap_(res), BARGS(_lt_Fcp_));
  _lt_Fax_(res, Nbb)
}

static inline COMPLEX *Bcache(const memindex integral)
  { return CACHEPTR(2,integral); }

static inline COMPLEX *BcacheC(const memindex integral)
  { return CACHEPTR(3,integral); }

static inline ComplexType Bval(const int i, const memindex integral)
  { return ToComplex(Bcache(integral)[i]); }

static inline ComplexType BvalC(const int i, const memindex integral)
  { return ToComplex(BcacheC(integral)[i]); }

static inline ComplexType B0i(const int i, BARGS(_lt_Cr_))
  { return Bval(EPSINDEX(i), Bget(BARGS(_lt_Id_))); }

static inline ComplexType B0iC(const int i, BARGS(_lt_Cc_))
  { return BvalC(EPSINDEX(i), BgetC(BARGS(_lt_Id_))); }

static inline ComplexType B0(BARGS(_lt_Cr_))
  { return B0i(bb0, BARGS(_lt_Id_)); }
static inline ComplexType B1(BARGS(_lt_Cr_))
  { return B0i(bb1, BARGS(_lt_Id_)); }
static inline ComplexType B00(BARGS(_lt_Cr_))
  { return B0i(bb00, BARGS(_lt_Id_)); }
static inline ComplexType B11(BARGS(_lt_Cr_))
  { return B0i(bb11, BARGS(_lt_Id_)); }
static inline ComplexType B001(BARGS(_lt_Cr_))
  { return B0i(bb001, BARGS(_lt_Id_)); }
static inline ComplexType B111(BARGS(_lt_Cr_))
  { return B0i(bb111, BARGS(_lt_Id_)); }
static inline ComplexType DB0(BARGS(_lt_Cr_))
  { return B0i(dbb0, BARGS(_lt_Id_)); }
static inline ComplexType DB1(BARGS(_lt_Cr_))
  { return B0i(dbb1, BARGS(_lt_Id_)); }
static inline ComplexType DB00(BARGS(_lt_Cr_))
  { return B0i(dbb00, BARGS(_lt_Id_)); }
static inline ComplexType DB11(BARGS(_lt_Cr_))
  { return B0i(dbb11, BARGS(_lt_Id_)); }

static inline ComplexType B0C(BARGS(_lt_Cc_))
  { return B0iC(bb0, BARGS(_lt_Id_)); }
static inline ComplexType B1C(BARGS(_lt_Cc_))
  { return B0iC(bb1, BARGS(_lt_Id_)); }
static inline ComplexType B00C(BARGS(_lt_Cc_))
  { return B0iC(bb00, BARGS(_lt_Id_)); }
static inline ComplexType B11C(BARGS(_lt_Cc_))
  { return B0iC(bb11, BARGS(_lt_Id_)); }
static inline ComplexType B001C(BARGS(_lt_Cc_))
  { return B0iC(bb001, BARGS(_lt_Id_)); }
static inline ComplexType B111C(BARGS(_lt_Cc_))
  { return B0iC(bb111, BARGS(_lt_Id_)); }
static inline ComplexType DB0C(BARGS(_lt_Cc_))
  { return B0iC(dbb0, BARGS(_lt_Id_)); }
static inline ComplexType DB1C(BARGS(_lt_Cc_))
  { return B0iC(dbb1, BARGS(_lt_Id_)); }
static inline ComplexType DB00C(BARGS(_lt_Cc_))
  { return B0iC(dbb00, BARGS(_lt_Id_)); }
static inline ComplexType DB11C(BARGS(_lt_Cc_))
  { return B0iC(dbb11, BARGS(_lt_Id_)); }

/****************************************************************/

#define CARGS(t) t(p1), t(p2), t(p1p2), t(m1), t(m2), t(m3)

static inline memindex Cget(CARGS(_lt_Cr_)) {
  _lt_Frd_(CARGS)
  extern memindex FORTRAN(cget)(CARGS(_lt_Fr_));
  return FORTRAN(cget)(CARGS(_lt_Frp_));
}

static inline memindex CgetC(CARGS(_lt_Cc_)) {
  _lt_Fcd_(CARGS)
  extern memindex FORTRAN(cgetc)(CARGS(_lt_Fc_));
  return FORTRAN(cgetc)(CARGS(_lt_Fcp_));
}

static inline void Cput(ComplexType *res, CARGS(_lt_Cr_)) {
  _lt_Frd_(CARGS)
  _lt_Fad_(res, Ncc)
  extern void FORTRAN(cput)(COMPLEX *res, CARGS(_lt_Fr_));
  FORTRAN(cput)(_lt_Fap_(res), CARGS(_lt_Frp_));
  _lt_Fax_(res, Ncc)
}

static inline void CputC(ComplexType *res, CARGS(_lt_Cc_)) {
  _lt_Fcd_(CARGS)
  _lt_Fad_(res, Ncc)
  extern void FORTRAN(cputc)(COMPLEX *res, CARGS(_lt_Fc_));
  FORTRAN(cputc)(_lt_Fap_(res), CARGS(_lt_Fcp_));
  _lt_Fax_(res, Ncc)
}

static inline void C0nocache(ComplexType *res, CARGS(_lt_Cr_)) {
  _lt_Frd_(CARGS)
  _lt_Fad_(res, 3)
  extern void FORTRAN(c0nocache)(COMPLEX *res, CARGS(_lt_Fr_));
  FORTRAN(c0nocache)(_lt_Fap_(res), CARGS(_lt_Frp_));
  _lt_Fax_(res, 3)
}

static inline void C0nocacheC(ComplexType *res, CARGS(_lt_Cc_)) {
  _lt_Fcd_(CARGS)
  _lt_Fad_(res, 3)
  extern void FORTRAN(c0nocachec)(COMPLEX *res, CARGS(_lt_Fc_));
  FORTRAN(c0nocachec)(_lt_Fap_(res), CARGS(_lt_Fcp_));
  _lt_Fax_(res, 3)
}

static inline COMPLEX *Ccache(const memindex integral)
  { return CACHEPTR(4,integral); }

static inline COMPLEX *CcacheC(const memindex integral)
  { return CACHEPTR(5,integral); }

static inline ComplexType Cval(const int i, const memindex integral)
  { return ToComplex(Ccache(integral)[i]); }

static inline ComplexType CvalC(const int i, const memindex integral)
  { return ToComplex(CcacheC(integral)[i]); }

static inline ComplexType C0i(const int i, CARGS(_lt_Cr_))
  { return Cval(EPSINDEX(i), Cget(CARGS(_lt_Id_))); }

static inline ComplexType C0iC(const int i, CARGS(_lt_Cc_))
  { return CvalC(EPSINDEX(i), CgetC(CARGS(_lt_Id_))); }

static inline ComplexType C0(CARGS(_lt_Cr_))
  { return C0i(cc0, CARGS(_lt_Id_)); }

static inline ComplexType C0C(CARGS(_lt_Cc_))
  { return C0iC(cc0, CARGS(_lt_Id_)); }

/****************************************************************/

#define DARGS(t) t(p1), t(p2), t(p3), t(p4), t(p1p2), t(p2p3), \
  t(m1), t(m2), t(m3), t(m4)

static inline memindex Dget(DARGS(_lt_Cr_)) {
  _lt_Frd_(DARGS)
  extern memindex FORTRAN(dget)(DARGS(_lt_Fr_));
  return FORTRAN(dget)(DARGS(_lt_Frp_));
}

static inline memindex DgetC(DARGS(_lt_Cc_)) {
  _lt_Fcd_(DARGS)
  extern memindex FORTRAN(dgetc)(DARGS(_lt_Fc_));
  return FORTRAN(dgetc)(DARGS(_lt_Fcp_));
}

static inline void Dput(ComplexType *res, DARGS(_lt_Cr_)) {
  _lt_Frd_(DARGS)
  _lt_Fad_(res, Ndd)
  extern void FORTRAN(dput)(COMPLEX *res, DARGS(_lt_Fr_));
  FORTRAN(dput)(_lt_Fap_(res), DARGS(_lt_Frp_));
  _lt_Fax_(res, Ndd)
}

static inline void DputC(ComplexType *res, DARGS(_lt_Cc_)) {
  _lt_Fcd_(DARGS)
  _lt_Fad_(res, Ndd)
  extern void FORTRAN(dputc)(COMPLEX *res, DARGS(_lt_Fc_));
  FORTRAN(dputc)(_lt_Fap_(res), DARGS(_lt_Fcp_));
  _lt_Fax_(res, Ndd)
}

static inline void D0nocache(ComplexType *res, DARGS(_lt_Cr_)) {
  _lt_Frd_(DARGS)
  _lt_Fad_(res, 3)
  extern void FORTRAN(d0nocache)(COMPLEX *res, DARGS(_lt_Fr_));
  FORTRAN(d0nocache)(_lt_Fap_(res), DARGS(_lt_Frp_));
  _lt_Fax_(res, 3)
}

static inline void D0nocacheC(ComplexType *res, DARGS(_lt_Cc_)) {
  _lt_Fcd_(DARGS)
  _lt_Fad_(res, 3)
  extern void FORTRAN(d0nocachec)(COMPLEX *res, DARGS(_lt_Fc_));
  FORTRAN(d0nocachec)(_lt_Fap_(res), DARGS(_lt_Fcp_));
  _lt_Fax_(res, 3)
}

static inline COMPLEX *Dcache(const memindex integral)
  { return CACHEPTR(6,integral); }

static inline COMPLEX *DcacheC(const memindex integral)
  { return CACHEPTR(7,integral); }

static inline ComplexType Dval(const int i, const memindex integral)
  { return ToComplex(Dcache(integral)[i]); }

static inline ComplexType DvalC(const int i, const memindex integral)
  { return ToComplex(DcacheC(integral)[i]); }

static inline ComplexType D0i(const int i, DARGS(_lt_Cr_))
  { return Dval(EPSINDEX(i), Dget(DARGS(_lt_Id_))); }

static inline ComplexType D0iC(const int i, DARGS(_lt_Cc_))
  { return DvalC(EPSINDEX(i), DgetC(DARGS(_lt_Id_))); }

static inline ComplexType D0(DARGS(_lt_Cr_))
  { return D0i(dd0, DARGS(_lt_Id_)); }

static inline ComplexType D0C(DARGS(_lt_Cc_))
  { return D0iC(dd0, DARGS(_lt_Id_)); }

/****************************************************************/

#define EARGS(t) t(p1), t(p2), t(p3), t(p4), t(p5), \
  t(p1p2), t(p2p3), t(p3p4), t(p4p5), t(p5p1), \
  t(m1), t(m2), t(m3), t(m4), t(m5)

static inline memindex Eget(EARGS(_lt_Cr_)) {
  _lt_Frd_(EARGS)
  extern memindex FORTRAN(eget)(EARGS(_lt_Fr_));
  return FORTRAN(eget)(EARGS(_lt_Frp_));
}

static inline memindex EgetC(EARGS(_lt_Cc_)) {
  _lt_Fcd_(EARGS)
  extern memindex FORTRAN(egetc)(EARGS(_lt_Fc_));
  return FORTRAN(egetc)(EARGS(_lt_Fcp_));
}

static inline void Eput(ComplexType *res, EARGS(_lt_Cr_)) {
  _lt_Frd_(EARGS)
  _lt_Fad_(res, Nee)
  extern void FORTRAN(eput)(COMPLEX *res, EARGS(_lt_Fr_));
  FORTRAN(eput)(_lt_Fap_(res), EARGS(_lt_Frp_));
  _lt_Fax_(res, Nee)
}

static inline void EputC(ComplexType *res, EARGS(_lt_Cc_)) {
  _lt_Fcd_(EARGS)
  _lt_Fad_(res, Nee)
  extern void FORTRAN(eputc)(COMPLEX *res, EARGS(_lt_Fc_));
  FORTRAN(eputc)(_lt_Fap_(res), EARGS(_lt_Fcp_));
  _lt_Fax_(res, Nee)
}

static inline void E0nocache(ComplexType *res, EARGS(_lt_Cr_)) {
  _lt_Frd_(EARGS)
  _lt_Fad_(res, 3)
  extern void FORTRAN(e0nocache)(COMPLEX *res, EARGS(_lt_Fr_));
  FORTRAN(e0nocache)(_lt_Fap_(res), EARGS(_lt_Frp_));
  _lt_Fax_(res, 3)
}

static inline void E0nocacheC(ComplexType *res, EARGS(_lt_Cc_)) {
  _lt_Fcd_(EARGS)
  _lt_Fad_(res, 3)
  extern void FORTRAN(e0nocachec)(COMPLEX *res, EARGS(_lt_Fc_));
  FORTRAN(e0nocachec)(_lt_Fap_(res), EARGS(_lt_Fcp_));
  _lt_Fax_(res, 3)
}

static inline COMPLEX *Ecache(const memindex integral)
  { return CACHEPTR(8,integral); }

static inline COMPLEX *EcacheC(const memindex integral)
  { return CACHEPTR(9,integral); }

static inline ComplexType Eval(const int i, const memindex integral)
  { return ToComplex(Ecache(integral)[i]); }

static inline ComplexType EvalC(const int i, const memindex integral)
  { return ToComplex(EcacheC(integral)[i]); }

static inline ComplexType E0i(const int i, EARGS(_lt_Cr_))
  { return Eval(EPSINDEX(i), Eget(EARGS(_lt_Id_))); }

static inline ComplexType E0iC(const int i, EARGS(_lt_Cc_))
  { return EvalC(EPSINDEX(i), EgetC(EARGS(_lt_Id_))); }

static inline ComplexType E0(EARGS(_lt_Cr_))
  { return E0i(ee0, EARGS(_lt_Id_)); }

static inline ComplexType E0C(EARGS(_lt_Cc_))
  { return E0iC(ee0, EARGS(_lt_Id_)); }

/****************************************************************/

#define XARGS(t) t(x)

static inline ComplexType Li2(XARGS(_lt_Cr_)) {
  _lt_Frd_(XARGS)
  COMPLEX res;
  extern void FORTRAN(li2sub)(COMPLEX *res, XARGS(_lt_Fr_));
  FORTRAN(li2sub)(&res, XARGS(_lt_Frp_));
  return ToComplex(res);
}

static inline ComplexType Li2C(XARGS(_lt_Cc_)) {
  _lt_Fcd_(XARGS)
  COMPLEX res;
  extern void FORTRAN(li2csub)(COMPLEX *res, XARGS(_lt_Fc_));
  FORTRAN(li2csub)(&res, XARGS(_lt_Fcp_));
  return ToComplex(res);
}

static inline ComplexType Li2omx(XARGS(_lt_Cr_)) {
  _lt_Frd_(XARGS)
  COMPLEX res;
  extern void FORTRAN(li2omxsub)(COMPLEX *res, XARGS(_lt_Fr_));
  FORTRAN(li2omxsub)(&res, XARGS(_lt_Frp_));
  return ToComplex(res);
}

static inline ComplexType Li2omxC(XARGS(_lt_Cc_)) {
  _lt_Fcd_(XARGS)
  COMPLEX res;
  extern void FORTRAN(li2omxcsub)(COMPLEX *res, XARGS(_lt_Fc_));
  FORTRAN(li2omxcsub)(&res, XARGS(_lt_Fcp_));
  return ToComplex(res);
}

/****************************************************************/

#define _lt_setreal_(f) static inline void Set##f(XARGS(_lt_Cr_)) { \
  _lt_Frd_(XARGS) \
  extern void FORTRAN(set##f)(CREAL *); \
  FORTRAN(set##f)(XARGS(_lt_Frp_)); \
} 
#define _lt_getreal_(f) static inline RealType Get##f() { \
  extern REAL FORTRAN(get##f)(void); \
  return ToReal(FORTRAN(get##f)()); \
}

_lt_setreal_(mudim)
#define setmudim Setmudim
_lt_getreal_(mudim)
#define getmudim Getmudim

_lt_setreal_(delta)
#define setdelta Setdelta
_lt_getreal_(delta)
#define getdelta Getdelta

_lt_setreal_(uvdiv)
#define setuvdiv Setuvdiv
_lt_getreal_(uvdiv)
#define getuvdiv Getuvdiv

_lt_setreal_(lambda)
#define setlambda Setlambda
_lt_getreal_(lambda)
#define getlambda Getlambda

_lt_setreal_(minmass)
#define setminmass Setminmass
_lt_getreal_(minmass)
#define getminmass Getminmass

_lt_setreal_(maxdev)
#define setmaxdev Setmaxdev
_lt_getreal_(maxdev)
#define getmaxdev Getmaxdev

_lt_setreal_(diffeps)
#define setdiffeps Setdiffeps
_lt_getreal_(diffeps)
#define getdiffeps Getdiffeps

_lt_setreal_(zeroeps)
#define setzeroeps Setzeroeps
_lt_getreal_(zeroeps)
#define getzeroeps Getzeroeps

#define _lt_setint_(f) static inline void Set##f(const int x) { \
  extern void FORTRAN(set##f)(CINTEGER *); \
  FORTRAN(set##f)(&x); \
}
#define _lt_getint_(f) static inline int Get##f() { \
  extern INTEGER FORTRAN(get##f)(void); \
  return FORTRAN(get##f)(); \
}

_lt_setint_(warndigits)
#define setwarndigits Setwarndigits
_lt_getint_(warndigits)
#define getwarndigits Getwarndigits

_lt_setint_(errdigits)
#define seterrdigits Seterrdigits
_lt_getint_(errdigits)
#define geterrdigits Geterrdigits

_lt_setint_(versionkey)
#define setversionkey Setversionkey
_lt_getint_(versionkey)
#define getversionkey Getversionkey

_lt_setint_(debugkey)
#define setdebugkey Setdebugkey
_lt_getint_(debugkey)
#define getdebugkey Getdebugkey

_lt_setint_(cmpbits)
#define setcmpbits Setcmpbits
_lt_getint_(cmpbits)
#define getcmpbits Getcmpbits

_lt_getint_(epsi)
#define getepsi Getepsi

static inline void Setdebugrange(const int f, const int t) {
  extern void FORTRAN(setdebugrange)(CINTEGER *, CINTEGER *);
  FORTRAN(setdebugrange)(&f, &t);
}
#define setdebugrange Setdebugrange

/****************************************************************/

extern void FORTRAN(clearcache)(void);
extern void FORTRAN(markcache)(void);
extern void FORTRAN(restorecache)(void);

extern void FORTRAN(ltini)(void);
extern void FORTRAN(ltexi)(void);

#if !0
#define clearcache FORTRAN(clearcache)
#define markcache FORTRAN(markcache)
#define restorecache FORTRAN(restorecache)
#define ltini FORTRAN(ltini)
#define ltexi FORTRAN(ltexi)
#endif

/****************************************************************/

#ifdef __cplusplus
}
#endif

#endif

