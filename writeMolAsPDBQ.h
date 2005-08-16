#ifndef _WRITEMOLASPDBQ
#define _WRITEMOLASPDBQ

extern FILE *stateFile;
extern int write_stateFile;

void writeMolAsPDBQ(Molecule *mol, FILE *output);

#endif /* WRITEPDBQSTATE */
