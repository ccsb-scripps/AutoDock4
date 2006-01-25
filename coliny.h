//
// coliny.h
//
// Interface to Coliny optimizers
//

#ifndef __coliny_h
#define __coliny_h

#if defined(USING_COLINY)

#include <vector>

//
// Initialize the 'algname' coliny optimizer over 'domain'
//
void coliny_init(char* algname, char* domain);

//
// Perform minimization with a given seed and initial point. Return
// summary statistics
//
void coliny_minimize(int seed, std::vector<double>& initpt,
				std::vector<double>& finalpt,
				int& neval, int& niters);

#endif

#endif
