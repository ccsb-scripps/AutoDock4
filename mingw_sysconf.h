/* AutoDock subset from Mingw32 - for sysconf() and gethostname()
 * $Id: mingw_sysconf.h,v 1.1 2012/04/17 22:56:54 mp Exp $
 */

/*
 * This file is part of the Mingw32 package.
 *
 * unistd.h maps (roughly) to io.h
 */
#include <io.h>

#include <sys/locking.h>

 

#include <time.h>

#define _SC_CLK_TCK 1

 

__CRT_INLINE long sysconf(int name)

{

if(name == _SC_CLK_TCK)

{

return CLK_TCK;

}

else

{

//printf("Only _SC_CLK_TCK can be used in simulated sysconf");

//exit(-1);
	return (long) 0;

}

}

 

#ifdef __MINGW32__

#define HAVE_GETHOSTNAME

#include <windows.h>

//int gethostname_mingw (char *, size_t);

 

__CRT_INLINE int  __cdecl gethostname_mingw (char *name, size_t len)

{

  DWORD dlen = len;

  return (GetComputerName (name, &dlen) ? 0 : -1);

}

#define gethostname gethostname_mingw

#endif

