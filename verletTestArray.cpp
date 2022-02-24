//============================================================================
// Name        : verletTestArray.cpp
// Author      : Ivo Roghair
// Version     :
// Copyright   : Open
// Description : A small benchmark for moving Lagrangian points in parallel
//                using a Verlet integration scheme.
//============================================================================

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

#define PA 73
#define PB 51
#define DT 0.001

#define pi M_PI

typedef double vec3[2];

void initVelocities(vec3 pos, vec3 vel) {
	double x = pos[0];
	double y = pos[1];
	vel[0] = -2*pow(sin(pi*x),2.0)*sin(pi*y)*cos(pi*y);
	vel[1] =  2*pow(sin(pi*y),2.0)*sin(pi*x)*cos(pi*x);
}

void getAcceleration(vec3 pos, vec3 vel, vec3 acc) {
	double ux = vel[0];
	double uy = vel[1];
	double x = pos[0];
	double y = pos[1];
	acc[0] = ux * (-4*pi*cos(pi*x)*cos(pi*y)*sin(pi*x)*sin(pi*y)) +
			uy * (-2*pi*cos(2*pi*y)*pow(sin(pi*x),2));
	acc[1] =  ux*(2*pi*cos(2*pi*x)*pow(sin(pi*y),2)) +
			uy*(4*pi*cos(pi*x)*cos(pi*y)*sin(pi*x)*sin(pi*y));
}

int main(int argc, char *argv[]) {
	int p,t;
	vec3 *pos, *vel, *acc, *acc_old, *acc_old_old;
	char fname[80];
	FILE *fp = 0;
	int NP = 30000;
    int NSTEP = 10000;
	int ND = 100;


	if (argc < 5) {
		printf("Error - needs 3 or 4 arguments (./program [NParticles] [NTimeSteps] [NFileDumps] [NThreads]\n");
		exit(1);
	}
	else {
		NP = atoi(argv[1]);
		NSTEP = atoi(argv[2]);
		ND = atoi(argv[3]);
		if (argc == 5) {
			omp_set_num_threads(atoi(argv[4]));
		}
	}

	pos = (vec3*) malloc(NP*sizeof(vec3));
	vel = (vec3*) malloc(NP*sizeof(vec3));
	acc = (vec3*) malloc(NP*sizeof(vec3));
	acc_old = (vec3*) malloc(NP*sizeof(vec3));
	acc_old_old = (vec3*) malloc(NP*sizeof(vec3));

	// Create pseudo-random positions for particles
#pragma omp parallel for
	for (p=0; p<NP; p++) {
		pos[p][0] = p % PA / (PA+1.0);
		pos[p][1] = p % PB / (PB+1.0);
	}

	// Initialize particle velocities and old acceleration arrays
#pragma omp parallel for
	for (p=0; p<NP; p++) {
		initVelocities(pos[p],vel[p]);
		getAcceleration(pos[p], vel[p], acc_old[p]);
		getAcceleration(pos[p], vel[p], acc_old_old[p]);
	}

for (t=0; t<NSTEP; t++) {
		// If appropriate, open a file for output
		if (t % ND == 0) {
			sprintf(fname, "Pos%05d.csv",t);
			fp = fopen(fname,"w");
			fprintf(fp, "ID,X,Y,UX,UY,AX,AY\n");
		}

		#pragma omp parallel for
		for(p=0; p<NP; p++) {
			// Get local acceleration
			getAcceleration(pos[p], vel[p], acc[p]);
			// Euler
//			vel[p][0] += acc[p][0] * DT;
//			vel[p][1] += acc[p][1] * DT;
//			pos[p][0] += vel[p][0] * DT;
//			pos[p][1] += vel[p][1] * DT;

			// Verlet integration
			pos[p][0] += vel[p][0] * DT + ( 2.0*acc[p][0]/3.0 - acc_old[p][0]/6.0 ) * DT * DT;
			vel[p][0] += (acc[p][0]/3.0 + 5.0*acc_old[p][0]/6.0 - acc_old_old[p][0]/6.0 ) * DT;
			pos[p][1] += vel[p][1] * DT + ( 2.0*acc[p][1]/3.0 - acc_old[p][1]/6.0 ) * DT * DT;
			vel[p][1] += (acc[p][1]/3.0 + 5.0*acc_old[p][1]/6.0 - acc_old_old[p][1]/6.0 ) * DT;

			// Write output
			if (t % ND == 0)
#pragma omp critical
			{
				fprintf(fp, "%d,%1.5e,%1.5e,%1.5e,%1.5e,%1.5e,%1.5e\n",
					p, pos[p][0],pos[p][1],
					vel[p][0],vel[p][1],
					acc[p][0],acc[p][1]);
			}

			// Update old and older arrays
			acc_old_old[p][0] = acc_old[p][0];
			acc_old_old[p][1] = acc_old[p][1];
			acc_old[p][0] = acc[p][0];
			acc_old[p][1] = acc[p][1];
		}
		if (t % ND == 0) {
			fclose(fp);
			printf("t = %1.4e\n", t*DT);
		}
	}


	return EXIT_SUCCESS;


}
