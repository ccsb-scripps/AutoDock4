/*

 $Id: printEnergies.cc,v 1.2.2.1 2005/03/02 20:46:01 gillet Exp $

*/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

/* printEnergies.cc */

#include <math.h>

    #include <stdlib.h>
    #include <stdio.h>
    #include <string.h>
    #include "printEnergies.h"

extern FILE *logFile;
extern FILE *stateFile;
extern int write_stateFile;

#define print1000(file,x) pr(file,  ((fabs((x)) >= 0.001) && ((fabs(x)) <= 1000.)) ? "%+7.2f" : "%+11.2e" , (x));

void printEnergies( FloatOrDouble einter, FloatOrDouble eintra, FloatOrDouble torsFreeEnergy, char  *prefixString, int ligand_is_inhibitor )
{
    FloatOrDouble deltaG = 0.0;
    FloatOrDouble Ki = 1.0;
    FloatOrDouble edocked=0.0;
    // FloatOrDouble RJ = 8.31441;  // in J/K/mol, Gas Constant, Atkins Phys.Chem., 2/e
    FloatOrDouble Rcal = 1.9871917; // in cal/K/mol, Gas Constant, RJ/4.184
    FloatOrDouble TK = 298.15;      // Room temperature, in K

    edocked = einter + eintra;
    // equilibrium:   E  +  I  <=>    EI
    // binding:       E  +  I   ->    EI         K(binding),      Kb
    // dissociation:     EI     ->  E  +  I      K(dissociation), Kd
    //
    //                            1
    //         K(binding) = ---------------
    //                      K(dissociation)
    // so:
    //      ln K(binding) = -ln K(dissociation)
    //              ln Kb = -ln Kd
    // Ki = dissociation constant of the enzyme-inhibitor complex = Kd
    //      [E][I]
    // Ki = ------
    //       [EI]
    // so:
    //              ln Kb = -ln Ki
    // deltaG(binding)    = -R*T*ln Kb
    // deltaG(inhibition) =  R*T*ln Ki
    //
    // Binding and Inhibition occur in opposite directions, so we 
    // lose the minus-sign:  deltaG = R*T*lnKi,  _not_ -R*T*lnKi
    // => deltaG/(R*T) = lnKi
    // => Ki = exp(deltaG/(R*T))
    deltaG = einter + torsFreeEnergy;
    if (deltaG < 0.0) {
        Ki = exp((deltaG*1000.)/(Rcal*TK));
    }

    pr( logFile, "%sEstimated Free Energy of Binding    = ", prefixString);
    print1000(logFile,deltaG);
    pr( logFile, " kcal/mol  [=(1)+(3)]\n");

    if (deltaG < 0.0) {
        if (ligand_is_inhibitor == 1) {
            pr( logFile, "%sEstimated Inhibition Constant, Ki   = ", prefixString);
        } else {
            pr( logFile, "%sEstimated Dissociation Constant, Kd = ", prefixString);
        }
        print1000(logFile,Ki);
        pr( logFile, "       [Temperature = %.2f K]\n", TK);
    }

    pr( logFile, "%s\n", prefixString);

    pr( logFile, "%sFinal Docked Energy                 = ", prefixString);
    print1000(logFile,edocked);
    pr( logFile, " kcal/mol  [=(1)+(2)]\n");

    pr( logFile, "%s\n", prefixString);

    pr( logFile, "%s(1) Final Intermolecular Energy     = ", prefixString);
    print1000(logFile,einter);
    pr( logFile, " kcal/mol\n");

    pr( logFile, "%s(2) Final Internal Energy of Ligand = ", prefixString);
    print1000(logFile,eintra);
    pr( logFile, " kcal/mol\n");

    pr( logFile, "%s(3) Torsional Free Energy           = ", prefixString);
    print1000(logFile,torsFreeEnergy);
    pr( logFile, " kcal/mol\n");

    pr( logFile, "%s\n", prefixString);
}

void printStateEnergies( FloatOrDouble einter, FloatOrDouble eintra, FloatOrDouble torsFreeEnergy, char  *prefixString, int ligand_is_inhibitor )
{
    FloatOrDouble deltaG = 0.0;
    FloatOrDouble Ki = 1.0;
    FloatOrDouble edocked=0.0;
    // FloatOrDouble RJ = 8.31441;  // in J/K/mol, Gas Constant, Atkins Phys.Chem., 2/e
    FloatOrDouble Rcal = 1.9871917; // in cal/K/mol, Gas Constant, RJ/4.184
    FloatOrDouble TK = 298.15;      // Room temperature, in K

    edocked = einter + eintra;
    // equilibrium:   E  +  I  <=>    EI
    // binding:       E  +  I   ->    EI         K(binding),      Kb
    // dissociation:     EI     ->  E  +  I      K(dissociation), Kd
    //
    //                            1
    //         K(binding) = ---------------
    //                      K(dissociation)
    // so:
    //      ln K(binding) = -ln K(dissociation)
    //              ln Kb = -ln Kd
    // Ki = dissociation constant of the enzyme-inhibitor complex = Kd
    //      [E][I]
    // Ki = ------
    //       [EI]
    // so:
    //              ln Kb = -ln Ki
    // deltaG(binding)    = -R*T*ln Kb
    // deltaG(inhibition) =  R*T*ln Ki
    //
    // Binding and Inhibition occur in opposite directions, so we 
    // lose the minus-sign:  deltaG = R*T*lnKi,  _not_ -R*T*lnKi
    // => deltaG/(R*T) = lnKi
    // => Ki = exp(deltaG/(R*T))
    deltaG = einter + torsFreeEnergy;
    if (deltaG < 0.0) {
        Ki = exp((deltaG*1000.)/(Rcal*TK));
    }

    pr(stateFile,"\t\t<free_NRG_binding>");
    print1000(stateFile,deltaG);
    pr(stateFile,"</free_NRG_binding>\n");
    if (deltaG < 0.0) {
      if (ligand_is_inhibitor == 1) {
	pr(stateFile,"\t\t<Ki>");
	print1000(stateFile,Ki);
	pr(stateFile,"</Ki>\n");
      } else {
	pr(stateFile,"\t\t<Kd>");
	print1000(stateFile,Ki);
	pr(stateFile,"</Kd>\n");
      }
      pr(stateFile,"\t\t<Temp>%.2f</Temp>\n",TK); //temperature in K
    } 
      pr(stateFile,"\t\t<final_dock_NRG>");
      print1000(stateFile,edocked);
      pr(stateFile,"</final_dock_NRG>\n");
      
      pr(stateFile,"\t\t<final_intermol_NRG>");
      print1000(stateFile,einter);
      pr(stateFile,"</final_intermol_NRG>\n");

      pr(stateFile,"\t\t<internal_ligand_NRG>");
      print1000(stateFile,eintra);
      pr(stateFile,"</internal_ligand_NRG>\n");
      
      pr(stateFile,"\t\t<torsonial_free_NRG>");
      print1000(stateFile,torsFreeEnergy);
      pr(stateFile,"</torsonial_free_NRG>\n"); 
      
}
