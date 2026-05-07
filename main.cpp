#include <mpi.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include "constants.h"
#include <cmath>
#include <Teuchos_DefaultMpiComm.hpp>
#include <Tpetra_Version.hpp>
#include <Tpetra_Core.hpp>
#include <Tpetra_Vector.hpp>
#include <Teuchos_Comm.hpp>
#include <Tpetra_Map.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_Array.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Amesos2.hpp>
#include <Amesos2_Solver.hpp>
#include <fftw3-mpi.h>

#include "input_deck.h"
#include "mesh_data.h"
#include "utilities.h"
#include "transport_data.h"
#include "species.h"
#include "field_solver.h"
#include "reaction.h"
#include "photon_solver.h"
#include <algorithm>
#include <chrono>
#include "photoionization.h"


// Teunissen dissertation
double psi(double a, double b)
{
	//return std::max(0, std::min(1.0,((2.0+x)/6.0), x));
	double aa = a * a;
	double ab = a * b;

	double psi_val;
	if(ab <= 0)
		psi_val = 0;
	else if(aa <= (0.25 * ab))
		psi_val = a;
	else if(aa <= (2.5 * ab))
		psi_val = (1.0/6.0) * (b + 2.0 * a);
	else
		psi_val = b;

	return psi_val;
}

void update_species_densities_finite_volume(mesh_data *my_data, FieldSolver *p, RunParameters *rp, std::vector<Species *> &species_list, double current_timestep)
{
	int end_y_val = my_data->my_y_size+1;
	if(my_data->my_rank == (my_data->world_size-1))
		end_y_val--;

	for(int j = 1; j < end_y_val; j++)
	{
		double r = (j+my_data->beg_y-1)*SPACE_STEP; // beginning of cell
		double r1 = (j+my_data->beg_y)*SPACE_STEP;		// end of cell
		double x_area = PI * r1 * r1 - PI * r * r;
		double y_area_plus_one = 2 * PI * r1 * SPACE_STEP;
		double y_area = 2 * PI * r * SPACE_STEP;
		double cell_volume = (PI * r1 * r1 * SPACE_STEP) - (PI * r * r * SPACE_STEP);

		for(int i = 0; i < (rp->x_size-1); i++)
		{
			double field_x_plus_surf = (p->GetElectricFieldXAt(i+1,j-1) + p->GetElectricFieldXAt(i+1,j)) / 2.0;
			double field_x_minus_surf = (p->GetElectricFieldXAt(i,j-1) + p->GetElectricFieldXAt(i,j)) / 2.0;
			double field_y_plus_surf = (p->GetElectricFieldYAt(i,j) + p->GetElectricFieldYAt(i+1,j)) / 2.0;
			double field_y_minus_surf = (p->GetElectricFieldYAt(i,j-1) + p->GetElectricFieldYAt(i+1,j-1)) / 2.0;

			double field_x_middle = (p->GetElectricFieldXAt(i,j-1) + p->GetElectricFieldXAt(i+1,j-1) +
															p->GetElectricFieldXAt(i,j) + p->GetElectricFieldXAt(i+1,j)) / 4.0;

			double field_y_middle = (p->GetElectricFieldYAt(i,j-1) + p->GetElectricFieldYAt(i+1,j-1) +
															p->GetElectricFieldYAt(i,j) + p->GetElectricFieldYAt(i+1,j)) / 4.0;
															
			double mag_e_field = sqrt(pow(field_x_middle,2) + pow(field_y_middle,2));
			double reduced_e_field = mag_e_field / N_B / 1e-21;
			//if(my_data->my_rank == 0 && j <= 1)
				//std::cout << i << "," << j << "," << reduced_e_field << std::endl;

			double e_x_t = 0.0;
			double e_y_t = 0.0;
			double mag_field_x_plus_one = 0.0;
			if(i == (rp->x_size-2))
			{
				e_x_t = (p->GetElectricFieldXAt(i+1,j-1) + p->GetElectricFieldXAt(i+1,j)) / 2.0;
				e_y_t = (p->GetElectricFieldYAt(i+1,j-1) + p->GetElectricFieldYAt(i+1,j)) / 2.0;
				mag_field_x_plus_one = sqrt(pow(e_x_t,2) + pow(e_y_t,2));
				mag_field_x_plus_one = 0.0;
			}
			else
			{
				e_x_t = (p->GetElectricFieldXAt(i+1,j-1) + p->GetElectricFieldXAt(i+2,j-1) +
															p->GetElectricFieldXAt(i+1,j) + p->GetElectricFieldXAt(i+2,j)) / 4.0;
				e_y_t = (p->GetElectricFieldYAt(i+1,j-1) + p->GetElectricFieldYAt(i+2,j-1) +
															p->GetElectricFieldYAt(i+1,j) + p->GetElectricFieldYAt(i+2,j)) / 4.0;
				mag_field_x_plus_one = (sqrt(pow(e_x_t,2) + pow(e_y_t,2)) + mag_e_field)/2.0;
			}
		
			double mag_field_x_minus_one = 0.0;
			if(i == 0)
			{
				e_x_t = (p->GetElectricFieldXAt(i,j-1) + p->GetElectricFieldXAt(i,j)) / 2.0;
				e_y_t = (p->GetElectricFieldYAt(i,j-1) + p->GetElectricFieldYAt(i,j)) / 2.0;
				mag_field_x_minus_one = sqrt(pow(e_x_t,2) + pow(e_y_t,2));
				mag_field_x_minus_one = 0.0;
			}
			else
			{
				e_x_t = (p->GetElectricFieldXAt(i-1,j-1) + p->GetElectricFieldXAt(i,j-1) +
															p->GetElectricFieldXAt(i-1,j) + p->GetElectricFieldXAt(i,j)) / 4.0;
				e_y_t = (p->GetElectricFieldYAt(i-1,j-1) + p->GetElectricFieldYAt(i,j-1) +
															p->GetElectricFieldYAt(i-1,j) + p->GetElectricFieldYAt(i,j)) / 4.0;
				mag_field_x_minus_one = (sqrt(pow(e_x_t,2) + pow(e_y_t,2)) + mag_e_field)/2.0;	
			}

			e_x_t = (p->GetElectricFieldXAt(i,j-1) + p->GetElectricFieldXAt(i+1,j-1) +
															p->GetElectricFieldXAt(i,j-2) + p->GetElectricFieldXAt(i+1,j-2)) / 4.0;
			e_y_t = (p->GetElectricFieldYAt(i,j-1) + p->GetElectricFieldYAt(i+1,j-1) +
															p->GetElectricFieldYAt(i,j-2) + p->GetElectricFieldYAt(i+1,j-2)) / 4.0;
			double mag_field_y_minus_one = (sqrt(pow(e_x_t,2) + pow(e_y_t,2)) + mag_e_field)/2.0;	
			
				
			e_x_t = (p->GetElectricFieldXAt(i,j) + p->GetElectricFieldXAt(i+1,j) +
															p->GetElectricFieldXAt(i,j+1) + p->GetElectricFieldXAt(i+1,j+1)) / 4.0;
			e_y_t = (p->GetElectricFieldYAt(i,j) + p->GetElectricFieldYAt(i+1,j) +
															p->GetElectricFieldYAt(i,j+1) + p->GetElectricFieldYAt(i+1,j+1)) / 4.0;
			double mag_field_y_plus_one = (sqrt(pow(e_x_t,2) + pow(e_y_t,2)) + mag_e_field)/2.0;

			for(int q = 0; q < species_list.size(); q++)
			{
				TransportData *tp = species_list[q]->GetTransportData();
				std::vector<Reaction *> g = species_list[q]->GetGainReactions();
				std::vector<Reaction *> l = species_list[q]->GetLossReactions();

				double v_x_plus = field_x_plus_surf * species_list[q]->GetChargeSign() * tp->GetMobility(mag_field_x_plus_one/N_B/1e-21)/N_B;
				double v_y_plus = field_y_plus_surf * species_list[q]->GetChargeSign() * tp->GetMobility(mag_field_y_plus_one/N_B/1e-21)/N_B;
				double v_x_minus = field_x_minus_surf * species_list[q]->GetChargeSign() * tp->GetMobility(mag_field_x_minus_one/N_B/1e-21)/N_B;
				double v_y_minus = field_y_minus_surf * species_list[q]->GetChargeSign() * tp->GetMobility(mag_field_y_minus_one/N_B/1e-21)/N_B;

				double d_m_x_plus = tp->GetDiffusionCoefficient(mag_field_x_plus_one/N_B/1e-21) / N_B;
				double d_m_y_plus = tp->GetDiffusionCoefficient(mag_field_y_plus_one/N_B/1e-21) / N_B;
				double d_m_x_minus = tp->GetDiffusionCoefficient(mag_field_x_minus_one/N_B/1e-21) / N_B;
				double d_m_y_minus = tp->GetDiffusionCoefficient(mag_field_y_minus_one/N_B/1e-21) / N_B;

				double mob_1 = 0.0;
				double mob_2 = 0.0;
				double mob_3 = 0.0;
				double mob_4 = 0.0;

				if(v_y_minus < 0)
				{
					// KOREN FLUX LIMITER
					double a = species_list[q]->GetOldDensity(i,j-1) - species_list[q]->GetOldDensity(i,j);
					double b = species_list[q]->GetOldDensity(i,j) - species_list[q]->GetOldDensity(i,j+1);
					double c = species_list[q]->GetOldDensity(i,j);
					mob_1 = v_y_minus * (c + psi(a,b)) * y_area;
				}
				else
				{
					// KOREN FLUX LIMITER
					double a = species_list[q]->GetOldDensity(i,j) - species_list[q]->GetOldDensity(i,j-1);
					double b = species_list[q]->GetOldDensity(i,j-1) - species_list[q]->GetOldDensity(i,j-2);
					double c = species_list[q]->GetOldDensity(i,j-1);
					mob_1 = v_y_minus * (c + psi(a,b)) * y_area;
				}
				if(v_x_plus < 0)
				{
					if(i == (rp->x_size-2))
					{
						mob_2 = 0.0;
					}
					else
					{
						// KOREN FLUX LIMITER
						double a = species_list[q]->GetOldDensity(i,j) - species_list[q]->GetOldDensity(i+1,j);
						double b = species_list[q]->GetOldDensity(i+1,j) - species_list[q]->GetOldDensity(i+2,j);
						double c = species_list[q]->GetOldDensity(i+1,j);
						mob_2 = -v_x_plus * (c + psi(a,b)) * x_area;
					}
				}
				else
				{
					if(i == (rp->x_size-2))
					{
						mob_2 = 0.0;
						mob_2 = -v_x_plus * species_list[q]->GetOldDensity(i,j) * x_area;
					}
					else if(i == 0)
					{
						mob_2 = -v_x_plus * species_list[q]->GetOldDensity(i,j) * x_area;
					}
					else
					{
						// KOREN FLUX LIMITER
						double a = species_list[q]->GetOldDensity(i+1,j) - species_list[q]->GetOldDensity(i,j);
						double b = species_list[q]->GetOldDensity(i,j) - species_list[q]->GetOldDensity(i-1,j);
						double c = species_list[q]->GetOldDensity(i,j);
						mob_2 = -v_x_plus * (c + psi(a,b)) * x_area;
					}
				}
				if(v_y_plus < 0)
				{
					// KOREN FLUX LIMITER
					double a = species_list[q]->GetOldDensity(i,j) - species_list[q]->GetOldDensity(i,j+1);
					double b = species_list[q]->GetOldDensity(i,j+1) - species_list[q]->GetOldDensity(i,j+2);
					double c = species_list[q]->GetOldDensity(i,j+1);
					mob_3 = -v_y_plus * (c + psi(a,b)) * y_area_plus_one;
				}
				else
				{
					// KOREN FLUX LIMITER
					double a = species_list[q]->GetOldDensity(i,j+1) - species_list[q]->GetOldDensity(i,j);
					double b = species_list[q]->GetOldDensity(i,j) - species_list[q]->GetOldDensity(i,j-1);
					double c = species_list[q]->GetOldDensity(i,j);
					mob_3 = -v_y_plus * (c + psi(a,b)) * y_area_plus_one;
				}
				if(v_x_minus < 0)
				{
					if(i == 0)
					{
						mob_4 = 0.0;
						mob_4 = v_x_minus * species_list[q]->GetOldDensity(i,j) * x_area;
					}
					else if(i == (rp->x_size-2))
					{
						mob_4 = v_x_minus * species_list[q]->GetOldDensity(i,j) * x_area;
					}
					else
					{
						// KOREN FLUX LIMITER
						double a = species_list[q]->GetOldDensity(i-1,j) - species_list[q]->GetOldDensity(i,j);
						double b = species_list[q]->GetOldDensity(i,j) - species_list[q]->GetOldDensity(i+1,j);
						double c = species_list[q]->GetOldDensity(i,j);
						mob_4 = v_x_minus * (c + psi(a,b)) * x_area;
					}
				}
				else
				{
					if(i == 0)
					{
						mob_4 = 0.0;
					}
					else if(i == 1)
					{
						mob_4 = v_x_minus * species_list[q]->GetOldDensity(i-1,j) * x_area;
					}
					else
					{
						// KOREN FLUX LIMITER
						double a = species_list[q]->GetOldDensity(i,j) - species_list[q]->GetOldDensity(i-1,j);
						double b = species_list[q]->GetOldDensity(i-1,j) - species_list[q]->GetOldDensity(i-2,j);
						double c = species_list[q]->GetOldDensity(i-1,j);
						mob_4 = v_x_minus * (c + psi(a,b)) * x_area;
					}
				}
				if(my_data->my_rank == 0 && j == 1)
				{
					mob_1 = 0.0;
				}
				if(my_data->my_rank == (my_data->world_size-1) && j == (my_data->my_y_size-1))
					mob_3 = 0.0;
				
				double diff_1 = 0.0;
				if(my_data->my_rank == 0 && j == 1)
					diff_1 = 0.0;
				else
					diff_1 = -d_m_y_minus * (species_list[q]->GetOldDensity(i,j) - species_list[q]->GetOldDensity(i,j-1)) / (1.0 * SPACE_STEP) * y_area;

				double diff_3 = 0.0;
				if(my_data->my_rank == (my_data->world_size-1) && j == my_data->my_y_size)
					diff_3 = 0.0;
				else
					diff_3 = d_m_y_plus * (species_list[q]->GetOldDensity(i,j+1) - species_list[q]->GetOldDensity(i,j)) / (1.0 * SPACE_STEP) * y_area_plus_one;

				double diff_2 = 0.0;
				if(i == (rp->x_size-2))
					diff_2 = 0.0;
					//diff_2 = d_m_x_plus * (0.0 - species_list[q]->GetOldDensity(i,j)) / (1.0 * SPACE_STEP) * x_area;
				else
					diff_2 = d_m_x_plus * (species_list[q]->GetOldDensity(i+1,j) - species_list[q]->GetOldDensity(i,j)) / (1.0 * SPACE_STEP) * x_area;
				
				double diff_4 = 0.0;
				if(i == 0)
					diff_4 = 0.0;
					//diff_4 = -d_m_x_minus * (species_list[q]->GetOldDensity(i,j) - 0.0) / (1.0 * SPACE_STEP) * x_area;
				else
					diff_4 = -d_m_x_minus * (species_list[q]->GetOldDensity(i,j) - species_list[q]->GetOldDensity(i-1,j)) / (1.0 * SPACE_STEP) * x_area;

				double gain_term = 0.0;
				std::vector<Species *> cur_species;
				for(int s = 0; s < g.size(); s++)														// gain processes
				{
					cur_species = g[s]->GetReactants();	
					double react_term = g[s]->GetReactionRate(reduced_e_field);
					for(int z = 0; z < cur_species.size(); z++)
					{
						react_term = react_term * cur_species[z]->GetOldDensity(i,j);
					}
					gain_term += react_term;
				}
				double loss_term = 0.0;
				for(int s = 0; s < l.size(); s++)														// loss processes
				{
					cur_species = l[s]->GetReactants();
					double react_term = l[s]->GetReactionRate(reduced_e_field);
					for(int z = 0; z < cur_species.size(); z++)
					{
						react_term = react_term * cur_species[z]->GetOldDensity(i,j);
					}
					loss_term += react_term;
				}
				///////////////////////////////////
				species_list[q]->SetDensity(i, j, species_list[q]->GetOldDensity(i,j) + ((mob_2 + mob_4 + mob_1 + mob_3)/cell_volume +
																																								(diff_1 + diff_2 + diff_3 + diff_4)/cell_volume +
																																							  (gain_term - loss_term)) * current_timestep); 
			}
		}
	}
}

void semi_implicit_finite_volume(mesh_data *my_data, FieldSolver *p, std::vector<Species *> &species_list, unsigned int time_int, double &current_timestep, RunParameters *rp)
{
	unsigned int max_iters = 0;
	for(int q = 0; q < species_list.size(); q++)
	{
		for(int i = 0; i < my_data->my_y_size; i++)
		{
			//species_list[q]->SetOldDensity(0,i,species_list[q]->GetOldDensity(1,i));
			//species_list[q]->SetDensity(0,i,species_list[q]->GetDensity(1,i));
			//species_list[q]->SetOldDensity(0,i,0.0);
			//species_list[q]->SetDensity(0,i,0.0);
			//species_list[q]->SetOldDensity(XSIZE-1,i,0.0);
			//species_list[q]->SetDensity(XSIZE-1,i,0.0);
		}
	}
	if(my_data->my_rank == 0)
	{
		for(int q = 0; q < species_list.size(); q++)
		{
			for(int i = 0; i < rp->x_size; i++)
			{
				species_list[q]->SetOldDensity(i,0,species_list[q]->GetOldDensity(i,1));
				species_list[q]->SetDensity(i,0,species_list[q]->GetDensity(i,1));
			}
		}
	}
	if(my_data->my_rank == (my_data->world_size-1))
	{
		for(int q = 0; q < species_list.size(); q++)
		{
			for(int i = 0; i < rp->x_size; i++)
			{
				species_list[q]->SetOldDensity(i,my_data->my_y_size,species_list[q]->GetOldDensity(i,my_data->my_y_size-1));
				species_list[q]->SetDensity(i,my_data->my_y_size,species_list[q]->GetDensity(i,my_data->my_y_size-1));
			}
		}
	}
	for(int q = 0; q < species_list.size(); q++)
	{
		species_list[q]->CopyDataToOld();
		species_list[q]->TransferData(my_data);
	}

	auto start = std::chrono::high_resolution_clock::now();
	for(int i = 0; i < 2; i++)
	{
		auto start_species = std::chrono::high_resolution_clock::now();
		update_species_densities_finite_volume(my_data, p, rp, species_list, current_timestep);
		auto end_species = std::chrono::high_resolution_clock::now();
		
		for(int q = 0; q < species_list.size(); q++)
			species_list[q]->TransferData(my_data);

		auto start_charge = std::chrono::high_resolution_clock::now();
		p->ComputeChargeDensity(species_list, my_data->beg_y);
		p->UpdateCharge();
		auto end_charge = std::chrono::high_resolution_clock::now();

		auto start_field = std::chrono::high_resolution_clock::now();
		if(rp->voltage_mode == "rf")
			p->UpdateBC(rp->voltage*sin(2.0*3.1415*rp->voltage_frequency*(time_int*TIME_STEP)));
		p->Solve();
		auto end_field = std::chrono::high_resolution_clock::now();

		std::chrono::duration<double> diff_field = end_field - start_field;
		std::chrono::duration<double> diff_species = end_species - start_species;
		std::chrono::duration<double> diff_charge = end_charge - start_charge;
		//max_iters += p->getNumberIterationsForSolve();
		//
		if((time_int % rp->status_stride) == 0 && my_data->my_rank == 0)
		{
			std::cout << "iteration = " << i << ", field solve iterations : " << p->getNumberIterationsForSolve() << ", field solve time = " << diff_field.count() << " seconds" <<  " , fvm time = " << diff_species.count() << " seconds" << ", charge density compute time = " << diff_charge.count() << " seconds" << std::endl;
		}
	}
	auto end = std::chrono::high_resolution_clock::now();

	if(my_data->my_rank == 0 && (time_int % rp->status_stride) == 0)
	{
		std::chrono::duration<double> diff = end - start;
		std::cout << "fvm + field solve time = " << diff.count() << " seconds" << std::endl;
	}

	if(max_iters > 25)
	{
		if(current_timestep/2.0 > 1e-16)
		{
			current_timestep = current_timestep / 2.0;
			if(my_data->my_rank == 0)
				std::cout << "DECREASING TIMESTEP to " << current_timestep << std::endl;
		}
		else
		{
			current_timestep = 1e-16;
			if(my_data->my_rank == 0)
				std::cout << "MINIMUM TIMESTEP REACHED " << current_timestep << std::endl;
		}
	}

	if(max_iters < 5)
	{
		if(current_timestep*2.0 < TIME_STEP)
		{
			current_timestep = current_timestep * 2.0;
			if(my_data->my_rank == 0)
				std::cout << "INCREASING TIMESTEP to " << current_timestep << std::endl;
		}
	}

	for(int q = 0; q < species_list.size(); q++)
		species_list[q]->CopyDataToOld();	
}


void populate_from_volumetric_source_term(mesh_data *my_data, std::vector<Species *> &species_list, RunParameters *rp, double **volumetric_source_term, int species_index, double current_timestep)
{
	for(int j = 1; j < my_data->my_y_size+1; j++)
	{
		for(int i = 0; i < rp->x_size; i++)
		{
			species_list[species_index]->SetDensity(i, j, species_list[species_index]->GetOldDensity(i,j) + volumetric_source_term[j-1][i] * current_timestep); 
		}
	}

	for(int q = 0; q < species_list.size(); q++)
		species_list[q]->CopyDataToOld();	
}

int main(int argc, char *argv[])
{

	MPI_Init(&argc, &argv);
	fftw_mpi_init();


	int my_rank;
	int world_size;

	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);

	mesh_data *data = new mesh_data;
	data->my_rank = my_rank;
	data->world_size = world_size;
	{

		Teuchos::RCP<const Teuchos::Comm<int> > t_comm (new Teuchos::MpiComm<int> (MPI_COMM_WORLD));
		/* INITIAL FIELD SOLVE 
		 * Field solver determines the alloaction
		 * of the problem to different processors
		 * so this must be called first 
		 * boundary conditions are set in here */
		if(argc < 2)
		{
			std::cout<<"Not enough arguments"<<std::endl;
			return -1;
		}
		else
		{
			std::ifstream testIn;
			
			testIn.open(argv[1], std::ios::in);
			if(!testIn.is_open())
			{
				std::cout<<"File cannot be opened"<<std::endl;
				return -1;

			}
				testIn.close();

		}
		InputDeck *id = new InputDeck(std::string(argv[1]));
		RunParameters *rp = new RunParameters;
		std::vector<Species *> species_list;
		id->ReadInputDeck(rp, data, species_list);

		FieldSolver *field_solver = new FieldSolver(t_comm, rp->voltage, BELOS, rp);
		field_solver->InitializeSolver();
		if(rp->voltage_mode == "rf")
			field_solver->UpdateBC(0.0);
		else
			field_solver->UpdateBC(rp->voltage);
	  field_solver->ConfigureInitialSolve();

		data->my_y_size = field_solver->GetYSize();		/* Get allocation of data to processors */
		data->beg_y = my_rank * (int)(rp->y_size / world_size);
		id->Finish(data, species_list, rp); // called last for any input that needs size of mesh set up

		/////////////////////////////////////////////
		MPI_Barrier(MPI_COMM_WORLD);
		if(data->my_rank == 0)
			std::cout << "creating volumetric source term" << std::endl;
		//////////////////////////////////////////////
		double **volumetric_source_term = new double*[data->my_y_size];
		double **volumetric_source_term_bourdon = new double*[data->my_y_size];
		double **volumetric_source_term_simplified = new double*[data->my_y_size];
		double **volumetric_source_term_fft = new double*[data->my_y_size];

		for(int i = 0; i < data->my_y_size; i++)
		{
			volumetric_source_term[i] = new double[rp->x_size];
			volumetric_source_term_bourdon[i] = new double[rp->x_size];
			volumetric_source_term_simplified[i] = new double[rp->x_size];
			volumetric_source_term_fft[i] = new double[rp->x_size];
		}

		for(int j = 0; j < data->my_y_size; j++)
		{
			for(int i = 0; i < rp->x_size; i++)
			{
				volumetric_source_term[j][i] = 0.0;
				volumetric_source_term_bourdon[j][i] = 0.0;
				volumetric_source_term_simplified[j][i] = 0.0;
				volumetric_source_term_fft[j][i] = 0.0;
			}
		}
		MPI_Barrier(MPI_COMM_WORLD);

		if(data->my_rank == 0)
		{
			std::cout << "assemble global mesh info" << std::endl;
		}
	
		mesh_information *global_mesh_info = new mesh_information[data->world_size];	
		// compute mesh information
		for(int i = 0; i < data->world_size; i++)
		{
			if(i == data->my_rank)
			{
				global_mesh_info[i].processor = data->my_rank;
				global_mesh_info[i].beg_y = data->beg_y;
				global_mesh_info[i].y_size = data->my_y_size;
			}
			MPI_Bcast(&global_mesh_info[i].processor, 1, MPI_INT, i, MPI_COMM_WORLD);
			MPI_Bcast(&global_mesh_info[i].beg_y, 1, MPI_INT, i, MPI_COMM_WORLD);
			MPI_Bcast(&global_mesh_info[i].y_size, 1, MPI_INT, i, MPI_COMM_WORLD);
		}

		if(data->my_rank == 0)
		{
			std::cout << "mesh info gathered" << std::endl;
			for(int i = 0; i < data->world_size; i++)
			{
				std::cout << global_mesh_info[i].processor << "\t" << global_mesh_info[i].beg_y << "\t" << global_mesh_info[i].y_size << std::endl;
			}
		}

		///////////////////////////////
		// PHOTOIONIZATION MODEL STUFF 
		spectral_data *spectral_info;
		std::vector<double> phi_array;
		fft_photoionization_model_data *fft_data;
		PhotonSolver *ps;

		if(rp->photoionization_method != "none")
		{
			// Uses experimental data and calculates integral model directly
			// initialize appropriate data
			if(rp->photoionization_method == "semi-empirical")
			{
				if(data->my_rank == 0)
					std::cout << "enabling semi-imperical integral photoionization method" << std::endl;
				spectral_info = new spectral_data;
				semiempirical_photoionization_model_initialize(data, spectral_info, phi_array, rp);
			}
			// FFT spectral method, uses fftw3 and intensity from Bourdon
			else if(rp->photoionization_method == "fft")
			{
				if(data->my_rank == 0)
					std::cout << "enabling FFT photoionization method" << std::endl;
				fft_data = new fft_photoionization_model_data;
				fft_photoionization_model_initialize(data, species_list, fft_data, rp);
				if(data->my_rank == 0)
					std::cout << "FFT initialize complete" << std::endl;
				if(data->my_rank == 0)
					std::cout << "enabling Helmholtz photoionization method" << std::endl;

				ps = new PhotonSolver(t_comm,rp);
				ps->InitializeSolver();
				ps->ConfigureSolver();
				if(fft_data->electron_index == -1)
				{
					rp->photoionization_method = "none";
					if(data->my_rank == 0)
						std::cout << "cannot finding species named e-, turning off" << std::endl;
				}
				if(fft_data->ionization_reaction == NULL)
				{
					rp->photoionization_method = "none";
					if(data->my_rank == 0)
						std::cout << "cannot find ionization reaction named n2_ionization, turning off" << std::endl;
				}
			}
			// 3-term helmholtz photoionization model
			else if(rp->photoionization_method == "helmholtz")
			{
				if(data->my_rank == 0)
					std::cout << "enabling Helmholtz photoionization method" << std::endl;

				ps = new PhotonSolver(t_comm,rp);
				ps->InitializeSolver();
				ps->ConfigureSolver();
				if(data->my_rank == 0)
					std::cout << "initial photon solve complete" << std::endl;
			}
			else
			{
				if(data->my_rank == 0)
					std::cout << "photoionization method not recognized, turning off" << std::endl;
				rp->photoionization_method = "none";
			}
		}

		// END PHOTOIONIZATION MODEL STUFF
		///////////////////////////////

		double current_time = 0.0;
		unsigned int time_int = 0;

		////////////////
		/// RESTART ///
		////////////////
		bool restart_enable = false;
		if(rp->restart != -1)
		{
			unsigned int res = rp->restart;
			std::cout << "attempting restart" << std::endl;
			if(!load_mesh_data(data, field_solver, rp, species_list, res, rp->output_directory))
			{
				std::cout << "error reading restart" << std::endl;
				MPI_Finalize();
				return 0;
			}
			else
			{
				current_time = (res+1) * TIME_STEP;
				time_int = res+1;
				restart_enable = true;
				std::cout << "restart complete, continuing sim" << std::endl;
			}
		}

	
		if(restart_enable == false)
		{
			if(data->my_rank == 0)
				std::cout << "perform initial solve" << std::endl;
			/***** don't touch these lines  *********/
			field_solver->PerformInitialSolve();
			if(data->my_rank == 0)
				std::cout << "initial solve complete" << std::endl;
		}

		for(int i = 0; i < species_list.size(); i++)
		{
			species_list[i]->TransferData(data);
			species_list[i]->CopyDataToOld();
		}
		field_solver->ComputeChargeDensity(species_list, data->beg_y);
		field_solver->UpdateCharge();
		if(rp->voltage_mode == "rf")
			field_solver->UpdateBC(0.0);
		field_solver->Solve();
		output_nodal_data(data, field_solver, rp, species_list, -1, rp->output_directory);
		/***** don't touch these lines  *********/
		MPI_Barrier(MPI_COMM_WORLD);

		///////////////////////
		// Set final simulation time here
		// timestep information is located in constants.h
		double t_final = rp->final_time;
		double current_timestep = TIME_STEP;
		////////////////////////
		print_reactions(data, species_list);
		if(rp->voltage_mode == "rf")
		{
			std::cout << "RF VOLTAGE MODE SELECTED" << std::endl;
			std::cout << "peak voltage = " << rp->voltage << std::endl;
			std::cout << "frequency = " << rp->voltage_frequency << std::endl;
		}
		//////////////////////////
		// MAIN TIME LOOP
		/////////////////////////
		while(current_time < t_final)
		{
			semi_implicit_finite_volume(data, field_solver, species_list, time_int, current_timestep, rp);

			if(rp->photoionization_method != "none")
			{
				if(rp->photoionization_method == "semi-empirical" && (time_int % 10) == 0)
				{
					semiempirical_photoionization_model(data, species_list, rp, volumetric_source_term, field_solver, global_mesh_info, spectral_info, phi_array);
				}
				else if(rp->photoionization_method == "fft" && (time_int % 1) == 0)
				{
					fft_photoionization_model(data, species_list, fft_data, volumetric_source_term_fft, field_solver, rp, time_int);		
				}
				else if(rp->photoionization_method == "helmholtz" && (time_int % 10) == 0)
				{
					ps->UpdateIonizationRate(species_list, field_solver);
					ps->Solve();
					for(int j = 0; j < data->my_y_size; j++)
						for(int i = 0; i < rp->x_size; i++)
							volumetric_source_term_bourdon[j][i] = ps->GetPhotonRateAt(i,j);
				}
				populate_from_volumetric_source_term(data, species_list, rp, volumetric_source_term_fft, 0, current_timestep); // populate electrons
				populate_from_volumetric_source_term(data, species_list, rp, volumetric_source_term_fft, 1, current_timestep); // populate O2+
			}
			/////////////////////////
			// simplified photoionization model 
			//simplified_photoionization_model_allprocs(data, species_list, volumetric_source_term_simplified, field_solver, rp, time_int);
			//populate_from_volumetric_source_term(data, species_list, volumetric_source_term_simplified, 0, current_timestep); // populate electrons
			//populate_from_volumetric_source_term(data, species_list, volumetric_source_term_simplified, 1, current_timestep); // populate O2+
			//

			// output stuff
			if(time_int % rp->output_stride == 0)
			{
				output_nodal_data(data, field_solver, rp, species_list, time_int, rp->output_directory);
				output_species_data(data, field_solver, rp, species_list, time_int, rp->output_directory);
				output_volumetric_source_term(data, rp, volumetric_source_term_fft, time_int, "fft", rp->output_directory);
				//output_volumetric_source_term(data, rp, volumetric_source_term_bourdon, time_int, "bourdon", rp->output_directory);
				//output_intensity(data, species_list, field_solver, phi_array, time_int, rp->output_directory);
				//output_volumetric_source_term(data, volumetric_source_term_bourdon, time_int, "bourdon", rp->output_directory);
			}

			double max_cfl = global_cfl(data, field_solver, rp, species_list, current_timestep);
			if(data->my_rank == 0 & (time_int % rp->status_stride) == 0)
			{
				std::cout << "timestep # = " <<  time_int << " current time = " << current_time << " max cfl = " << max_cfl << std::endl;
			}
			////////////////////////
			MPI_Barrier(MPI_COMM_WORLD);
			current_time += current_timestep;
			time_int++;
		}
		
	}
	MPI_Barrier(MPI_COMM_WORLD);

	MPI_Finalize();

	return 0;
}


