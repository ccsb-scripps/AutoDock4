//
// coliny.h
//
// Interface to Coliny optimizers
//

#ifndef __coliny_h
#define __coliny_h

#if defined(USING_COLINY)

#include <utilib/BasicArray.h>
#include <utilib/CommonIO.h>

using utilib::Flush;

//
// Initialize the 'algname' coliny optimizer over 'domain'
//
void coliny_init(char* algname, char* domain,
			utilib::BasicArray<double>& initpt);

//
// Perform minimization with a given seed and initial point. Return
// summary statistics
//
void coliny_minimize(int seed, utilib::BasicArray<double>& initpt,
				utilib::BasicArray<double>& finalpt,
				int& neval, int& niters);

#endif

#endif
