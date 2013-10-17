/*

 $Id: ranlib.h,v 1.11 2013/10/17 23:39:06 mp Exp $

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

/* Prototypes for all user accessible RANLIB routines: 
   note only a few are used in AutoDock 
   Most are functions in ranlib.cc, a few are from the underlying com.cc
 */

#ifndef _RANLIB_H
#define _RANLIB_H

#include "typedefs.h"

//extern void advnst(const FourByteLong k);
extern Real genbet(ConstReal  aa,ConstReal  bb);
//extern Real genchi(ConstReal  df);
extern Real genexp(ConstReal  av);
extern Real genf(ConstReal  dfn, ConstReal  dfd);
//extern Real gengam(ConstReal  a,ConstReal  r);
extern void genmn(const Real *const parm, /* not const */ Real *const x, /* not const */ Real *const work);
extern void genmul(const FourByteLong n,const Real *const p,const FourByteLong ncat,/* not const */ FourByteLong *const ix);
//extern Real gennch(ConstReal  df,ConstReal  xnonc);
extern Real gennf(ConstReal  dfn, ConstReal  dfd, ConstReal  xnonc);
extern Real gennor(ConstReal  av,ConstReal  sd); // referenced by ls.h
extern void genprm(/* not const */FourByteLong *const iarray,const int larray);
extern Real genunf(ConstReal  low,ConstReal  high); // referenced by ls.h
//extern void getsd(FourByteLong *const iseed1,FourByteLong *const iseed2);
//extern FourByteLong ignbin(const FourByteLong n,ConstReal  pp);
extern FourByteLong ignnbn(const FourByteLong n,ConstReal  p);
extern FourByteLong ignlgi(void); // referenced by gs.cc
//extern FourByteLong ignpoi(ConstReal  mu);
extern FourByteLong ignuin(const FourByteLong low,const FourByteLong high); // referenced by gs.cc
extern FourByteLong mltmod(const FourByteLong a,const FourByteLong s,const FourByteLong m); // referenced by com.cc
//extern void phrtsd(const char *const phrase,FourByteLong *const seed1,FourByteLong *const seed2);

/* ranf() changed to macro to avoid extra function call, M Pique 2013, see ranlib.cc */
//extern Real ranf(void); // referenced by gs.cc
/*
     4.656613057E-10 is 1/M1  M1 is set in a data statement in IGNLGI
      and is currently 2147483563. If M1 changes, change this also.
*/
#define ranf() ((Real)(ignlgi()*4.656613057E-10))

extern void setall(const FourByteLong iseed1,const FourByteLong iseed2); // referenced by main.cc
//extern void setgmn(const Real *const meanv,Real *const covm,const FourByteLong p,Real *const parm);
extern void setsd(const FourByteLong iseed1,const FourByteLong iseed2); // reference by com.cc
//extern Real sexpo(void);
//extern Real sgamma(ConstReal a);
//extern Real snorm(void);
//extern Real rcauchy(ConstReal , ConstReal );
//extern Real scauchy1(void); // referenced by gencau.cc
extern Real scauchy2(void); // referenced by gencau.cc

#endif
