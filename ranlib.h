/*

 $Id: ranlib.h,v 1.8 2010/10/01 22:51:40 mp Exp $

 AutoDock 

Copyright (C) 2009 The Scripps Research Institute. All rights reserved.

 AutoDock is a Trade Mark of The Scripps Research Institute.

 This program is free software; you can redistribute it and/or
 modify it under the terms of the GNU General Public License
 as published by the Free Software Foundation; either version 2
 of the License, or (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

 */

/* Prototypes for all user accessible RANLIB routines */

#ifndef _RANLIB_H
#define _RANLIB_H

#include "typedefs.h"

extern void advnst(const FourByteLong k);
extern Real genbet(ConstReal  aa,ConstReal  bb);
extern Real genchi(ConstReal  df);
extern Real genexp(ConstReal  av);
extern Real genf(ConstReal  dfn, ConstReal  dfd);
extern Real gengam(ConstReal  a,ConstReal  r);
extern void genmn(const Real *const parm, /* not const */ Real *const x, /* not const */ Real *const work);
extern void genmul(const FourByteLong n,const Real *const p,const FourByteLong ncat,/* not const */ FourByteLong *const ix);
extern Real gennch(ConstReal  df,ConstReal  xnonc);
extern Real gennf(ConstReal  dfn, ConstReal  dfd, ConstReal  xnonc);
extern Real gennor(ConstReal  av,ConstReal  sd);
extern void genprm(/* not const */FourByteLong *const iarray,const int larray);
extern Real genunf(ConstReal  low,ConstReal  high);
extern void getsd(FourByteLong *const iseed1,FourByteLong *const iseed2);
extern void gscgn(const FourByteLong getset,const FourByteLong *g);
extern FourByteLong ignbin(const FourByteLong n,ConstReal  pp);
extern FourByteLong ignnbn(const FourByteLong n,ConstReal  p);
extern FourByteLong ignlgi(void);
extern FourByteLong ignpoi(ConstReal  mu);
extern FourByteLong ignuin(const FourByteLong low,const FourByteLong high);
extern void initgn(const FourByteLong isdtyp);
extern FourByteLong mltmod(const FourByteLong a,const FourByteLong s,const FourByteLong m);
extern void phrtsd(const char *const phrase,FourByteLong *const seed1,FourByteLong *const seed2);
extern Real ranf(void);
extern void setall(const FourByteLong iseed1,const FourByteLong iseed2);
extern void setant(const FourByteLong qvalue);
extern void setgmn(const Real *const meanv,Real *const covm,const FourByteLong p,Real *const parm);
extern void setsd(const FourByteLong iseed1,const FourByteLong iseed2);
extern Real sexpo(void);
extern Real sgamma(ConstReal a);
extern Real snorm(void);
extern Real rcauchy(ConstReal , ConstReal );
extern Real scauchy1(void);
extern Real scauchy2(void);

#endif
