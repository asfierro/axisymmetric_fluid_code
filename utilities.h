#ifndef UTILITIES_H
#define UTILITIES_H

#include "field_solver.h"
#include "mesh_data.h"
#include "species.h"
#include "reaction.h"
#include "run_parameters.h"

void transfer_data(RunParameters *,double **, int, int, int);
void transfer_data_two_cell(RunParameters *, double **, double *, double *, int, int, int);
void output_nodal_data(mesh_data *, FieldSolver *, RunParameters *, std::vector<Species *>&, int, std::string, std::string prefix="");
void output_species_data(mesh_data *, FieldSolver *, RunParameters *, std::vector<Species *>&, int, std::string, std::string prefix="");
void output_charge_density(mesh_data *, FieldSolver *, RunParameters *, int, std::string, std::string prefix="");
void output_intensity(mesh_data *, RunParameters *, std::vector<Species *> &, FieldSolver *, std::vector<double> &, int, std::string);
void output_volumetric_source_term(mesh_data *, RunParameters *, double **, int, std::string, std::string);
void populate_vectors_from_file(std::vector<double> &, std::vector<double> &, std::string);
void populate_vectors_from_file(std::vector<double> &, std::vector<double> &, std::string, double, double);
std::string remove_white_space(std::string);
double interpolate_data(std::vector<double> &, std::vector<double> &, double);
int get_current_processor_from_global_index(int, mesh_information *, int);
bool load_mesh_data(mesh_data *, FieldSolver *, RunParameters *, std::vector<Species *> &, int, std::string);
void check_large_density(mesh_data *, FieldSolver *, RunParameters *, std::vector<Species *> &, double);
void check_negative_density(mesh_data *, FieldSolver *, RunParameters *, std::vector<Species *> &, double);
void scale_volumetric_by_inverse_volume(mesh_data *, RunParameters *,double **);
void scale_volumetric_by_volume(mesh_data *, RunParameters *, double **);
void print_reactions(mesh_data *, std::vector<Species *> &);
double global_cfl(mesh_data *, FieldSolver *, RunParameters *, std::vector<Species *> &, double);
void output_directory_creation(mesh_data *, RunParameters *);
int get_species_index(std::vector<Species *> &, std::string);

#endif

