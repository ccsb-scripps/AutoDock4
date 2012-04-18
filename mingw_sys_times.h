#ifndef	_SYS_TIMES_H
#ifdef __cplusplus
extern "C" {
#endif
#define	_SYS_TIMES_H
//#include <_ansi.h>
//#include <machine/types.h>
#include <time.h>   // will set or define "clock_t CLK_TCK"
//#ifndef __clock_t_defined
//typedef _CLOCK_T_ clock_t;
//typedef clock_t _CLOCK_T_; // clock_t;
//#define __clock_t_defined
//#endif
/*  Get Process Times, P1003.1b-1993, p. 92 */
struct tms {
	clock_t	tms_utime;		/* user time */
	clock_t	tms_stime;		/* system time */
	clock_t	tms_cutime;		/* user time, children */
	clock_t	tms_cstime;		/* system time, children */
};
//clock_t _EXFUN(times,(struct tms *));
 __CRT_INLINE clock_t __cdecl times( struct tms *buffer )
{
// this version: supply totally bogus values for user and system CPU times:
buffer->tms_utime         = CLK_TCK;
buffer->tms_stime         = clock();
buffer->tms_cutime = clock();
buffer->tms_cstime = clock();
return clock();
}
#ifdef __cplusplus
}
#endif
#endif	/* !_SYS_TIMES_H */
