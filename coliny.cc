//
// coliny.cc
//
// This file provides the low-level hooks to the Coliny optimizers.
// Coliny optimizers can be accessed within AutoDock using the command
// line:
//
//	coliny <algname> <numtrials>
//
// If the "help" algname is used, this code prints the list of 
// solvers that are supported within the current Coliny configuration.
//

#if USING_COLINY

#include "coliny.h"
#include <coliny/coliny.h>

//
// The AutoDock 'objective function' used within Coliny
//
double ADEvalFn(double* x, int n);

//
// Global COLIN problem
//
colin::OptProblem<BasicArray<double>,colin::AppResponse_Utilib> coliny_problem;
//
// Global Coliny solver
//
coliny::ColinySolver<colin::OptProblem<BasicArray<double>, colin::AppResponse_Utilib>,BasicArray<double> > coliny_solver;




////
//// Initialize the "algname" optimizer over the given domain.  An initial
//// point is generate as the midpoint over the domain.
////
void coliny_init(char* algname, char* domain)
{
//
// If 'algname' equals "help", then return after calling 
// ColinySolver::initialize
//
if (strcmp(algname,"help")==0) {
   coliny_solver.initialize(algname);
   return;
   }
//
// If 'domain' equals "help", then return after print the 
// solver options
//
if (strcmp(domain,"help")==0) {
   coliny_solver.initialize(algname);
   coliny_solver.help_parameters(cout);
   cout << flush;
   return;
   }
//
// Setup the OptProblem object
//
colin::OptSetup(coliny_problem, &ADEvalFn, domain);
//
// Initialize the OptSolver object
//
coliny_solver.initialize(algname);
//
// Read in the OptSolver parameters
//
ifstream ifstr;
char fname[256];
sprintf(fname,"%s.in",algname);
ifstr.open(fname);
if (ifstr)
   coliny_solver.read_parameter_values(ifstr);
//
// Create a default 'initial point'
//
/*
BasicArray<double> lower,upper;
coliny_problem.get_real_bounds(lower,upper);
initpt.resize(lower.size());
for (size_type i=0; i<lower.size(); i++)
  initpt[i] = (upper[i]+lower[i])/2.0;
*/
}


////
//// Perform minimize using a generic Coliny minimization utility
////
//// To turn on "full debugging", set the false flag to true.
////
void coliny_minimize(int seed, utilib::BasicArray<double>& initpt,
				utilib::BasicArray<double>& finalpt,
				int& neval, int& niters)
{
colin::real best_value;
coliny_solver.minimize(coliny_problem, initpt, seed, false, false, finalpt, best_value);
}


#endif

