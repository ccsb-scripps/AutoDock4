/* Prototypes for all user accessible RANLIB routines */
#ifndef _RANLIB_H
#define _RANLIB_H

#include "structs.h"

extern void advnst(FourByteLong k);
extern float genbet(float aa,float bb);
extern float genchi(float df);
extern float genexp(float av);
extern float genf(float dfn, float dfd);
extern float gengam(float a,float r);
extern void genmn(float *parm,float *x,float *work);
extern void genmul(FourByteLong n,float *p,FourByteLong ncat,FourByteLong *ix);
extern float gennch(float df,float xnonc);
extern float gennf(float dfn, float dfd, float xnonc);
extern float gennor(float av,float sd);
extern void genprm(FourByteLong *iarray,int larray);
extern float genunf(float low,float high);
extern void getsd(FourByteLong *iseed1,FourByteLong *iseed2);
extern void gscgn(FourByteLong getset,FourByteLong *g);
extern FourByteLong ignbin(FourByteLong n,float pp);
extern FourByteLong ignnbn(FourByteLong n,float p);
extern FourByteLong ignlgi(void);
extern FourByteLong ignpoi(float mu);
extern FourByteLong ignuin(FourByteLong low,FourByteLong high);
extern void initgn(FourByteLong isdtyp);
extern FourByteLong mltmod(FourByteLong a,FourByteLong s,FourByteLong m);
extern void phrtsd(char* phrase,FourByteLong* seed1,FourByteLong* seed2);
extern float ranf(void);
extern void setall(FourByteLong iseed1,FourByteLong iseed2);
extern void setant(FourByteLong qvalue);
extern void setgmn(float *meanv,float *covm,FourByteLong p,float *parm);
extern void setsd(FourByteLong iseed1,FourByteLong iseed2);
extern float sexpo(void);
extern float sgamma(float a);
extern float snorm(void);
extern float rcauchy(float, float);
extern float scauchy1(void);

#endif
