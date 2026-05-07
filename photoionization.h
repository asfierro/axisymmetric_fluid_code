#ifndef _PHOTOIONIZATION_H
#define _PHOTOIONIZATION_H

#include <vector>

#include "field_solver.h"
#include "mesh_data.h"
#include "species.h"
#include "utilities.h"
#include "run_parameters.h"
#include "reaction.h"
#include <fftw3-mpi.h>
#include <complex>

#define photo_skip 8

struct absorption_data
{
	std::vector<double> total_x;
	std::vector<double> total_y;
	std::vector<double> pi_x;
	std::vector<double> pi_y;
};

struct spectral_data
{
	std::vector<double> wavelength;
	std::vector<double> intensity;
	std::vector<double> pi_coefficient;
	std::vector<double> total_coefficient;
	std::vector<double> intensity_times_pi_coefficient;
};

struct fft_photoionization_model_data
{
	double l_1;
	double l_2;
	double l_3;
	double A_1;
	double A_2;
	double A_3;
	double p_o2;
	double q_factor;

	double *rhs_1;			// first term
	double *rhs_2;			// second term
	double *rhs_3;			// third term

  double *lhs_1;			// first term
	double *lhs_2;			// second term
	double *lhs_3;			// third term

	double *rhs_output_1;			// first term
	double *rhs_output_2;			// second term
	double *rhs_output_3;			// third term

	double *freq_solution_1;			// solution term 1
	double *freq_solution_2;			// solution term 2
	double *freq_solution_3;			// solution term 3

	// forward plans of RHS
	fftw_plan fwd_plan_1;
	fftw_plan fwd_plan_2;
	fftw_plan fwd_plan_3;

	// backward plans for solution
	fftw_plan bwd_plan_1;
	fftw_plan bwd_plan_2;
	fftw_plan bwd_plan_3;

	double **intensity_term;

	// needed to find the intensity based upon ionzation rate
	Reaction *ionization_reaction;
	int electron_index;		// electron index inside of species_list vector
};

void semiempirical_photoionization_model_initialize(mesh_data *, spectral_data *, std::vector<double> &, RunParameters *);
void semiempirical_photoionization_model(mesh_data *, std::vector<Species *> &, RunParameters *, double **, FieldSolver *, mesh_information *, spectral_data *, std::vector<double> &);
void simplified_photoionization_model(mesh_data *, std::vector<Species *> &, RunParameters * , double **, FieldSolver *);
void simplified_photoionization_model_allprocs(mesh_data *, std::vector<Species *> &, double **, FieldSolver *, RunParameters *, int);
void fft_photoionization_model_initialize(mesh_data *, std::vector<Species *> &, fft_photoionization_model_data *, RunParameters *);
void fft_photoionization_model(mesh_data *, std::vector<Species *> &, fft_photoionization_model_data *, double **, FieldSolver *, RunParameters *, int);

#endif
