/* print_avsfld.cc */


    #include <stdio.h>
    #include "print_avsfld.h"
    #include "autocomm.h"


void print_avsfld( FILE *logFile,
		   int veclen,
		   int natom,
		   int nframe,
		   int offset[VECLENMAX],
		   int stride,
		   char label[MAX_CHARS],
		   char filename[MAX_CHARS] )
{
    int i;

    fprintf( logFile, "AVSFLD: # AVS field file\n" );
    fprintf( logFile, "AVSFLD: #\n" );
    fprintf( logFile, "AVSFLD: # Created by AutoDock\n" );
    fprintf( logFile, "AVSFLD: #\n" );
    fprintf( logFile, "AVSFLD: ndim=2           # number of dimensions in the field\n" );
    fprintf( logFile, "AVSFLD: nspace=1         # number of physical coordinates\n" );
    fprintf( logFile, "AVSFLD: veclen=%-9d # vector size\n", veclen );
    fprintf( logFile, "AVSFLD: dim1=%-11d # atoms\n", natom );
    fprintf( logFile, "AVSFLD: dim2=%-11d # conformations\n", nframe );
    fprintf( logFile, "AVSFLD: data=float       # data type (byte,integer,float,double)\n" );
    fprintf( logFile, "AVSFLD: field=uniform    # field coordinate layout\n" );
    fprintf( logFile, "AVSFLD: label= %s\n", label );
    for (i=0; i<veclen; i++) {
        fprintf( logFile, "AVSFLD: variable %d file = %s filetype = ascii offset = %d stride = %d\n", 
					   i+1,    filename,                   offset[i],    stride );
    }
    fprintf( logFile, "AVSFLD: # end of file\n\n" );
    fflush( logFile );
}
/* EOF */
