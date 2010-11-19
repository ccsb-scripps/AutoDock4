/* dimLibrary.cc */

#include <math.h>

#include <stdio.h>
#include "dimLibrary.h"
#include "constants.h"
#include "ranlib.h"
#include "qmultiply.h"

extern FILE *logFile;
extern int nlig;  // value assigned in readPDBQT

/**
 *  The Position of N demension for a particle is the representaion of a ligand's state
 *  P.x[1..Nd] has following sectors:
 *  	x([0..2])n is the translation
 *  	x([3..6])n is rotation
 *  	x([7..ntor])n is the torsion section
 */
void make_state_from_position( State *S, Position P)
{
    register int i, j;
    int idx;
    
    //handle multi-ligand  -Huameng
    for(i = 0; i < nlig; i++) {
    	idx = 7*i;
    	//translation
	    S->T[i].x = P.x[idx + 0];
	    S->T[i].y = P.x[idx + 1];
	    S->T[i].z = P.x[idx + 2];
	    // quarterion
	    /*    	    
	    S->Q[i].x = P.x[idx + 3];
	    S->Q[i].y = P.x[idx + 4];
	    S->Q[i].z = P.x[idx + 5];
	    S->Q[i].w = P.x[idx + 6];
	    */  	
	   // rotation	 
	   S->Q[i].nx = P.x[idx + 3];
	   S->Q[i].ny = P.x[idx + 4];
	   S->Q[i].nz = P.x[idx + 5];
	   S->Q[i].ang = P.x[idx + 6];		   
	      	         		   
    }
    
    // all torsion part
    for(j = 7*nlig, i = 0; i < S->ntor; i++, j++) {
		S->tor[i] = P.x[j];
	}	
}


void copy_state_to_position(Position *R , State S)
{
    register int i, j, idx;
    //handle multi-ligand  -Huameng
    for(i = 0; i < nlig; i++) {
    	idx = 7*i;  	
    		 	
		R->x[idx + 0] = S.T[i].x;
		R->x[idx + 1] = S.T[i].y;
	    R->x[idx + 2] = S.T[i].z;
	        
	    // quarterion
	    /*	   
	    R->x[idx + 3] = S.Q[i].x;
	    R->x[idx + 4] = S.Q[i].y;
	    R->x[idx + 5] = S.Q[i].z;
	    R->x[idx + 6] = S.Q[i].w;
	   	*/
	   	// rotation
	    R->x[idx + 3] = S.Q[i].nx;
	    R->x[idx + 4] = S.Q[i].ny;
	    R->x[idx + 5] = S.Q[i].nz;
	    R->x[idx + 6] = S.Q[i].ang;
	       	 	        
    }
    
    for(j= 7*nlig, i=0; i < S.ntor; i++, j++) {
		R->x[j] = S.tor[i];
    }
}

/**
 * 
 * 
 * 
 */
void initializeStateDimensions(GridMapSetInfo *info, double *xmin, double *xmax, double pso_qstep, double pso_dstep, int d_freedom) 
{
	int i, d, offset;
	//handle multi-ligand -Huameng
	for(i = 0; i < nlig; i++)
	{								
		// 0, 1, 2 (x,y,z) for translation
		offset = 7*i;
		xmin[offset + 0] = info->lo[X];
    	xmax[offset + 0] = info->hi[X];
    						
		xmin[offset + 1] = info->lo[Y];
    	xmax[offset + 1] = info->hi[Y];	
    			
		xmin[offset + 2] = info->lo[Z];
    	xmax[offset + 2] = info->hi[Z];
			
		fprintf(logFile, "In initialise X: d=%d  xmin=%.3f	xmax=%.3f\n", offset + 0, xmin[offset + 0], xmax[offset + 0]);
		fprintf(logFile, "In initialise Y: d=%d  xmin=%.3f	xmax=%.3f\n", offset + 1, xmin[offset + 1], xmax[offset + 1]);		
		fprintf(logFile, "In initialise Z: d=%d  xmin=%.3f	xmax=%.3f\n", offset + 2, xmin[offset + 2], xmax[offset + 2]);
		fflush(logFile);
			
		// 7*i + 3,4,5,6 are for orientation
		xmin[offset + 3] = -1.0;
    	xmax[offset + 3] = 1.0;					
		xmin[offset + 4] = -1.0;
    	xmax[offset + 4] = 1.0;			
		xmin[offset + 5] = -1.0;
    	xmax[offset + 5] = 1.0;		    	
    	xmin[offset + 6] = -1.0;
    	xmax[offset + 6] = 1.0;
	    	
	    fprintf(logFile, "In initialise Rotation: d=%d  xmin=%.3f	xmax=%.3f\n", offset + 6, xmin[offset + 6], xmax[offset + 6]);
		fflush(logFile);
					    												 		    
	}
	
	// the rest are for torsion
	for(d = 7*nlig; d < d_freedom; d++)
	{
		xmin[d] = -PI;
	    xmax[d] = PI;	    
	    fprintf(logFile, "In initialise Torsion: xmin[%d]=%.3f	xmax[%d]=%.3f\n", d, xmin[d], d, xmax[d]);	    		
	}
	
	fprintf(logFile, "End initializeStateDimensions() of PSO \n");
	fflush(logFile);
}




/**
 * 
 */
void initialiseDimension(GridMapSetInfo *info, double *xmin, double *xmax, int D)
{
	int i, d, offset;
			
	//handle multi-ligand -Huameng
	for(i = 0; i < nlig; i++)
	{								
		// 0, 1, 2 (x,y,z) for translation
		offset = 7*i;
		xmin[offset + 0] = info->lo[X];
    	xmax[offset + 0] = info->hi[X];
    						
		xmin[offset + 1] = info->lo[Y];
    	xmax[offset + 1] = info->hi[Y];	
    			
		xmin[offset + 2] = info->lo[Z];
    	xmax[offset + 2] = info->hi[Z];
//#ifdef DEBUG			
		fprintf(logFile, "In initialise X: d=%d  xmin=%.3f	xmax=%.3f\n", offset + 0, xmin[offset + 0], xmax[offset + 0]);
		fprintf(logFile, "In initialise Y: d=%d  xmin=%.3f	xmax=%.3f\n", offset + 1, xmin[offset + 1], xmax[offset + 1]);		
		fprintf(logFile, "In initialise Z: d=%d  xmin=%.3f	xmax=%.3f\n", offset + 2, xmin[offset + 2], xmax[offset + 2]);
		fflush(logFile);
//#endif			
		// 7*i + 3,4,5,6 are for orientation
		xmin[offset + 3] = 0.0;
    	xmax[offset + 3] = 1.0;					
		xmin[offset + 4] = 0;
    	xmax[offset + 4] = 1.0;			
		xmin[offset + 5] = 0;
    	xmax[offset + 5] = 1.0;		    	
    	xmin[offset + 6] = -PI;
    	xmax[offset + 6] = PI;
//#ifdef DEBUG
		fprintf(logFile, "In initialise Q: d=%d  xmin=%.3f	xmax=%.3f\n", offset + 3, xmin[offset + 3], xmax[offset + 3]);
		fprintf(logFile, "In initialise Q: d=%d  xmin=%.3f	xmax=%.3f\n", offset + 4, xmin[offset + 4], xmax[offset + 4]);		
		fprintf(logFile, "In initialise Q: d=%d  xmin=%.3f	xmax=%.3f\n", offset + 5, xmin[offset + 5], xmax[offset + 5]);			    	
	    fprintf(logFile, "In initialise Q: d=%d  xmin=%.3f	xmax=%.3f\n", offset + 6, xmin[offset + 6], xmax[offset + 6]);
		fflush(logFile);
//#endif					    												 		    
	}
	
	// the rest are for torsion
	for(d = 7*nlig; d < D; d++)
	{
		xmin[d] = -PI;
	    xmax[d] = PI;
//#ifdef DEBUG	    
	    fprintf(logFile, "In initialise Torsion: xmin[%d]=%.3f	xmax[%d]=%.3f\n", d, xmin[d], d, xmax[d]);
//#endif	    		
	}
	
	fprintf(logFile, "End initialiseDimension() of PSO \n");
	fflush(logFile);
	/*
    for (d=0; d < D ; d++)
    {
    	switch(d)
        {
	        case 0:
				xmin[d] = xlo;
	            xmax[d] = xhi;
				break;
			case 1:
				xmin[d] = ylo;
				xmax[d] = yhi;
				break;
			case 2:
				xmin[d] = zlo;
				xmax[d] = zhi;
				break;
			case 3: // For Quaternion angles
				xmin[d] = 0;
				xmax[d] = 1;
				break;
			case 4:
				xmin[d] = 0;
				xmax[d] = 1;
				break;
			case 5:
				xmin[d] = 0;
				xmax[d] = 1;
				break;				
			case 6:
				xmin[d] = -PI;
				xmax[d] = PI;
				break;
			default: // For Torsional Angles
				xmin[d] = -PI;
				xmax[d] = PI;
				break;		
		}
	}
	*/
}


void initialiseParticle(int s, int D, Position *Xi, Velocity *Vi, double *xmin, double *xmax, double *Vmin, double *Vmax)
{
    int d = 0, i, j;
	double temp;
	
	// handle multi-ligand -Huameng	02/24/08
	for(i = 0; i < nlig; i++) {
		for(j = 0; j < 7; j++) {			
			if( j <= 2 ) {
				temp = xmin[d] - xmax[d];
				Vmax[d] = fabs(temp) / 2;
				Vmin[d] = -Vmax[d];
			} else if( j > 2 && j < 6) {
				temp = xmin[d] - xmax[d];
				Vmax[d] = fabs(temp);
				Vmin[d] = -Vmax[d];				
			} else {				
				Vmin[d] = -PI;
				Vmax[d] = PI;				
			}
			Xi[s].x[d] = random_range(xmin[d], xmax[d]);
			Xi[s].prev_x[d] = Xi[s].x[d];
			Vi[s].v[d] = random_range(Vmin[d], Vmax[d] ); 
			d++;
#ifdef DEBUG			
			pr(logFile, "Particle first 7 Dim: Xi[%d].x[%d]=%.3f\n", s, d, Xi[s].x[d]);	
			fflush(logFile);			
#endif																
		}						
	}
		
	// torsion part
	for( d = 7*nlig; d < D; d++ )
	{
		Vmin[d] = -PI;
		Vmax[d] = PI;
				
		Xi[s].x[d] = random_range(xmin[d], xmax[d]);
		Xi[s].prev_x[d] = Xi[s].x[d];
		Vi[s].v[d] = random_range(Vmin[d], Vmax[d] ); 							
	}
	//pr(logFile, "initialise torsions  Xi[%d].x[%d]=%.3f\n", s, d-1, Xi[s].x[d-1]);	
		
	
	/*	
	for (d=0; d<D; d++)
	{
		if(d > 2 && d < 6)
		{
			temp = xmin[d] - xmax[d];
			Vmax[d] = fabs(temp);
			Vmin[d] = -Vmax[d];
		} 
		else if (d > 5)
		{
			Vmax[d] = PI;
			Vmin[d] = -PI;
		} 
		else
		{
			temp = xmin[d] - xmax[d];
			Vmax[d] = fabs(temp) / 2;
			Vmin[d] = -Vmax[d];
		}
		Xi[s].x[d] = random_range(xmin[d], xmax[d]);
		Xi[s].prev_x[d] = Xi[s].x[d];
		Vi[s].v[d] = random_range(Vmin[d], Vmax[d] ); 
	}
	*/
}
			
void swarmActivity(int S, int D, Position *Xi, int nb_eval, int outlev)
{
	int s, d;
	double swarm_activity;
	double position_dist[S_max];    // Euclidian distance of position and prev position
	double diff;                    // just a helper variable for calculating the distance
	// RG Swarm activity begin
    swarm_activity = 0;
    for(s=0; s < S; s++)
	{
    	for(d = 0; d < D; d++)
        {
			diff = Xi[s].x[d] - Xi[s].prev_x[d];
			//pr(logFile, "swarm_move %d particle= %d dim= %d diff= %f\n", nb_eval, s, d, diff);
			// Differences of angles
			if (d > 5)
			{
				if (diff < -PI) 
					diff = diff + PI;
				if (diff > PI)  
					diff = diff - PI;
				
			}
			position_dist[s] += (diff * diff);
		}
		position_dist[s] = sqrt(position_dist[s]);
		swarm_activity += position_dist[s];
	}
	
	if (outlev >1)
	{
		swarm_activity = (swarm_activity) / (S * D);
		pr(logFile, "swarm_move %d SA= %f\n", nb_eval+1, swarm_activity);
	}
	// RG Swarm activity end
}

/* EOF */

			
