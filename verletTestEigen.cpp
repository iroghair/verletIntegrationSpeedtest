//============================================================================
// Name        : verletTestArray.cpp
// Author      : Ivo Roghair
// Version     :
// Copyright   : Open
// Description : Hello World in C, Ansi-style
//============================================================================

#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <math.h>
#include <omp.h>
#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>

using namespace Eigen;
#define EIGEN_RUNTIME_NO_MALLOC
#define EIGEN_NO_DEBUG

#define PA 731
#define PB 237
#define DT 0.001

#define pi M_PI

typedef double vec3[2];

void initVelocities(MatrixXd pos, MatrixXd* vel) {
//	vel[0] = -2*pow(sin(pi*x),2.0)*sin(pi*y)*cos(pi*y);
//	vel[1] =  2*pow(sin(pi*y),2.0)*sin(pi*x)*cos(pi*x);
	vel->col(0) = -2*(pi*pos.col(0)).array().sin()*(pi*pos.col(0)).array().sin() * (pi*pos.col(1)).array().sin() * (pi*pos.col(1)).array().cos();
	vel->col(1) =  2*(pi*pos.col(1)).array().sin()*(pi*pos.col(1)).array().sin() * (pi*pos.col(0)).array().sin() * (pi*pos.col(0)).array().cos();
}

void getAcceleration(MatrixXd pos, MatrixXd vel, MatrixXd* acc) {
	acc->col(0) = vel.col(0).array() * (-4.0) * pi *(
					(pi*pos.col(0)).array().cos() *
					(pi*pos.col(1)).array().cos() *
					(pi*pos.col(0)).array().sin() *
					(pi*pos.col(1)).array().sin() ) +
					vel.col(1).array() * (-2.0*pi*(pi*2.0*pos.col(1)).array().cos() *
					(pi*pos.col(0)).array().sin() * (pi*pos.col(0)).array().sin());
	acc->col(1) = vel.col(0).array() * (2.0*pi*(2.0*pi*pos.col(0)).array().cos() *
					(pi*pos.col(1)).array().sin() * (pi*pos.col(1)).array().sin()) +
					vel.col(1).array() * (4.0 * pi *
					(pi*pos.col(0)).array().cos() *
					(pi*pos.col(1)).array().cos() *
					(pi*pos.col(0)).array().sin() *
					(pi*pos.col(1)).array().sin());
//
////	acc[0] = ux * (-4*pi*cos(pi*x)*cos(pi*y)*sin(pi*x)*sin(pi*y)) +
////			uy * (-2*pi*cos(2*pi*y)*pow(sin(pi*x),2));
//	acc[1] =  ux*(2*pi*cos(2*pi*x)*pow(sin(pi*y),2)) +
//			uy*(4*pi*cos(pi*x)*cos(pi*y)*sin(pi*x)*sin(pi*y));
}

int main(int argc, char *argv[]) {
	int p,t=0;

	char fname[80];
	FILE *fp;
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

	MatrixXd pos(NP,2);
	MatrixXd vel(NP,2);
	MatrixXd acc(NP,2);
	MatrixXd acc_old(NP,2);
	MatrixXd acc_old_old(NP,2);

	// Create pseudo-random positions for particles
	for (p=0; p<NP; p++) {
		pos(p,0) = p % PA / (PA+1.0);
		pos(p,1) = p % PB / (PB+1.0);
	}

	// Initialize particle velocities and old acceleration arrays
	initVelocities(pos,&vel);
	getAcceleration(pos,vel,&acc);
	getAcceleration(pos,vel,&acc_old);
	getAcceleration(pos,vel,&acc_old_old);

	for (t=0; t<NSTEP; t++) {

		getAcceleration(pos,vel,&acc);
		pos = pos + vel*DT + ( 2.0*acc/3.0 - acc_old/6.0 ) * DT * DT;
		vel = vel + (acc/3.0 + 5.0*acc_old/6.0 - acc_old_old/6.0 ) * DT;

		// Write output
		if (t % ND == 0)
		{
			sprintf(fname, "Pos%05d.csv",t);
			fp = fopen(fname,"w");
			fprintf(fp, "ID,X,Y,UX,UY,AX,AY\n");

			for (p = 0; p < NP; p++)
			  fprintf(fp, "%i,%1.5e,%1.5e,%1.5e,%1.5e,%1.5e,%1.5e\n",
				p,pos(p,0),pos(p,1),
				vel(p,0),vel(p,1),
				acc(p,0),acc(p,1));
		}

		// Update old and older arrays
		acc_old_old = acc_old;
		acc_old = acc;
	}

	return EXIT_SUCCESS;
}
