#include "utilities.h"

#include <iostream>
#include <ostream>
#include <fstream>
#include <vector>
#include <cstring>
#include <cstdlib>
#include <iomanip>
#include <string>
#include <sys/stat.h>
#include <sys/types.h>


void transfer_data(RunParameters *rp, double **data, int my_rank, int world_size, int my_y_size)
{
	/* send data to next processor */
	if((my_rank%2) == 0)
	{
		MPI_Send(data[my_y_size], rp->x_size, MPI_DOUBLE, my_rank+1, 0, MPI_COMM_WORLD);
	}
	else
	{	
		MPI_Recv(data[0], rp->x_size, MPI_DOUBLE, my_rank-1, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}

	if((my_rank%2) == 1 && my_rank != (world_size-1))
	{
		MPI_Send(data[my_y_size], rp->x_size, MPI_DOUBLE, my_rank+1, 0, MPI_COMM_WORLD);
	}
	else if(my_rank != 0 && my_rank != (world_size-1))
	{
		MPI_Recv(data[0], rp->x_size, MPI_DOUBLE, my_rank-1, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}

	/* send data to previous processor  */
	if((my_rank%2) == 0 && my_rank != 0)
	{
		MPI_Send(data[1], rp->x_size, MPI_DOUBLE, my_rank-1, 0, MPI_COMM_WORLD);
	}
	else if(my_rank != 0 && my_rank != (world_size-1))
	{	
		MPI_Recv(data[my_y_size+1], rp->x_size, MPI_DOUBLE, my_rank+1, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}

	if((my_rank%2) == 1)
	{
		MPI_Send(data[1], rp->x_size, MPI_DOUBLE, my_rank-1, 0, MPI_COMM_WORLD);
	}
	else 
	{
		MPI_Recv(data[my_y_size+1], rp->x_size, MPI_DOUBLE, my_rank+1, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}
}

void transfer_data_two_cell(RunParameters *rp, double **data_send, double *data_after, double *data_before, int my_rank, int world_size, int my_y_size)
{
	/* send data to next processor */
	if((my_rank%2) == 0)
	{
		MPI_Send(data_send[my_y_size-1], rp->x_size, MPI_DOUBLE, my_rank+1, 0, MPI_COMM_WORLD);
	}
	else
	{	
		MPI_Recv(data_before, rp->x_size, MPI_DOUBLE, my_rank-1, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}

	if((my_rank%2) == 1 && my_rank != (world_size-1))
	{
		MPI_Send(data_send[my_y_size-1], rp->x_size, MPI_DOUBLE, my_rank+1, 0, MPI_COMM_WORLD);
	}
	else if(my_rank != 0 && my_rank != (world_size-1))
	{
		MPI_Recv(data_before, rp->x_size, MPI_DOUBLE, my_rank-1, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}

	/* send data to previous processor  */
	if((my_rank%2) == 0 && my_rank != 0)
	{
		MPI_Send(data_send[2], rp->x_size, MPI_DOUBLE, my_rank-1, 0, MPI_COMM_WORLD);
	}
	else if(my_rank != 0 and my_rank != (world_size-1))
	{	
		MPI_Recv(data_after, rp->x_size, MPI_DOUBLE, my_rank+1, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}

	if((my_rank%2) == 1)
	{
		MPI_Send(data_send[2], rp->x_size, MPI_DOUBLE, my_rank-1, 0, MPI_COMM_WORLD);
	}
	else 
	{
		MPI_Recv(data_after, rp->x_size, MPI_DOUBLE, my_rank+1, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}
}

bool load_mesh_data(mesh_data *my_data, FieldSolver *p, RunParameters *rp,std::vector<Species *> &species_list, int time_int, std::string od)
{
	std::stringstream ss;
	ss << my_data->my_rank;
	std::string my_rank_str = ss.str();

	std::stringstream time_ss;
	time_ss << time_int;
	std::string time_int_str = time_ss.str();

	bool success = true;
	
	std::string file = od + "mesh_data/data_out_" + time_int_str + "_" + my_rank_str + ".csv";
	std::ifstream infile;
	infile.open(file,std::ios::in);

	std::string line;
	long long int line_ct = 0;

	if(infile.is_open())
	{
		std::getline(infile, line); // read single header line
		while(infile)
		{
			long long int i = line_ct % rp->x_size;
			long long int j = line_ct / rp->x_size;
			std::getline(infile, line);
		
			if(!line.empty())
			{
				char *cstr = new char[line.length()+1];
				std::strcpy(cstr, line.c_str());

				char *line_list = std::strtok(cstr,","); // read x val
				line_list = std::strtok(nullptr,",");			 // read y val

				line_list = std::strtok(nullptr,",");
				p->SetPotentialAt(i,j,atof(line_list));

				line_list = std::strtok(nullptr,",");
				p->SetElectricFieldXAt(i,j, atof(line_list));
			
				line_list = std::strtok(nullptr,",");
				p->SetElectricFieldYAt(i,j, atof(line_list));

				//for(int q = 0; q < species_list.size(); q++)
				//{
					//line_list = std::strtok(nullptr,",");
					//species_list[q]->SetDensity(i, j+1, atof(line_list));
				//}
				//p->SetPotentialAt(i,j) = std::string::atod(std::string(line_list[2]));
				//p->SetElectriclFieldXAt(i,j) = std::string::atod(std::string(line_list[3]));
				//p->SetElectriclFieldYAt(i,j) = std::string::atod(std::string(line_list[4]));
				//for(int q = 0; q < species_list.size(); q++)
					//species_list[q]->SetDensity(i,j+1) = std::string::atod(std::string(line_list[5+q]));
				line_ct++;
				delete[] cstr;
			}
		}		 
		infile.close();
		//std::cout << "LINE COUNT = " << line_ct << std::endl;
		//for(int q = 0; q < species_list.size(); q++)
			//species_list[q]->CopyDataToOld();
	}
	else
	{
		std::cout << "could not open restart file: " << file << std::endl;
		success = false;
	}

	file = od + "species_data/species_out_" + time_int_str + "_" + my_rank_str + ".csv";
	infile.open(file,std::ios::in);
	line_ct = 0;
	if(success && infile.is_open())
	{
		std::getline(infile, line); // read single header line
		while(infile)
		{
			long long int i = line_ct % (rp->x_size-1);
			long long int j = line_ct / (rp->x_size-1);
			std::getline(infile, line);
			if(!line.empty())
			{
				char *cstr = new char[line.length()+1];
				std::strcpy(cstr, line.c_str());
				char *line_list = std::strtok(cstr,","); // read x val
				line_list = std::strtok(nullptr,","); // read y val
				line_list = std::strtok(nullptr,","); // read z val

				for(int q = 0; q < species_list.size(); q++)
				{
					line_list = std::strtok(nullptr,",");

					species_list[q]->SetDensity(i, j+1, atof(line_list));
				}

				line_ct++;
				delete [] cstr;
			}
		}
		for(int q = 0; q < species_list.size(); q++)
			species_list[q]->CopyDataToOld();
		infile.close();
	}
	else
	{
		if(success == true)
		{
			std::cout << "could not open restart file: " << file << std::endl;
		}
		success = false;
	}

	return success;
}

void output_intensity(mesh_data *my_data, RunParameters *rp, std::vector<Species *> &species_list, FieldSolver *p, std::vector<double> &phi_array, int time_int, std::string od)
{
	std::ofstream outfile;
	std::stringstream ss;
	ss << my_data->my_rank;
	std::string my_rank_str = ss.str();

	std::stringstream time_ss;
	time_ss << time_int;
	std::string time_int_str = time_ss.str();

	//std::string file = "data/data_out_" + my_rank_str + "_" + time_int_str + ".dat";
	std::string file = od + "intensity_data/data_out_" + my_rank_str + "_" + time_int_str + ".dat";
	outfile.open(file.c_str(), std::ios::out);

	TransportData *tp;
	int electron_index;
	// find electron's
	for(int q = 0; q < species_list.size(); q++)
	{
		if(species_list[q]->GetSpeciesType() == ELECTRON)
		{
			electron_index = q;
			break;
		}
	}
	tp = species_list[electron_index]->GetTransportData();

	for(int j = 0; j < my_data->my_y_size; j++)		
	{
		for(int i = 0; i < rp->x_size; i++)
		{
			double mag_e_field = sqrt(pow(p->GetElectricFieldXAt(i,j),2) + pow(p->GetElectricFieldYAt(i,j),2));
			double t_electron_density = species_list[electron_index]->GetDensity(i, j+1);		// add one to get into species coordinates for y direction
			double reduced_e_field = mag_e_field / N_B / 1e-21;
			double u_m = tp->GetMobility(reduced_e_field) / N_B;
			
			double cell_power = pow(t_electron_density * u_m * Q_E * pow(mag_e_field,2),1.0) / (1.2 / 1.57e-8);

			outfile << i*SPACE_STEP << "\t" << (j+my_data->beg_y)*SPACE_STEP << "\t" << cell_power << "\t" << mag_e_field << "\t" << t_electron_density << "\t" << u_m;
			outfile << std::endl;
		}
	}
	outfile.close();
}

void output_charge_density(mesh_data *my_data, FieldSolver *p, RunParameters *rp, int time_int, std::string od, std::string prefix)
{
	std::ofstream outfile;
	std::stringstream ss;
	ss << my_data->my_rank;
	std::string my_rank_str = ss.str();

	std::stringstream time_ss;
	time_ss << time_int;
	std::string time_int_str = time_ss.str();

	//std::string file = "data/data_out_" + my_rank_str + "_" + time_int_str + ".dat";
	std::string file = od + "mesh_data/charge_density_out_" + time_int_str + "_" + my_rank_str + prefix + ".csv";
	outfile.open(file.c_str(), std::ios::out);

	outfile<<std::setprecision(10);
	outfile<<"x,y,z,charge_density,";
	outfile << std::endl;
	//std::cout << "my rank = " << my_data->my_rank << " file = " << file << " y_size = " << my_data->my_y_size << std::endl;
	for(int j = 0; j < my_data->my_y_size; j++)
	{
		for(int i = 0; i < rp->x_size; i++)
		{
			outfile << i*SPACE_STEP << "," << (j+my_data->beg_y)*SPACE_STEP << "," << "0," << p->GetChargeDensity(i,j);
			outfile << std::endl;

		}
	}
	//std::cout << "my rank = " << my_data->my_rank << " output complete " << std::endl;
	outfile.close();
}

void output_nodal_data(mesh_data *my_data, FieldSolver *p, RunParameters *rp, std::vector<Species *> &species_list, int time_int, std::string od, std::string prefix)
{
	std::ofstream outfile;
	std::stringstream ss;
	ss << my_data->my_rank;
	std::string my_rank_str = ss.str();

	std::stringstream time_ss;
	time_ss << time_int;
	std::string time_int_str = time_ss.str();

	//std::string file = "data/data_out_" + my_rank_str + "_" + time_int_str + ".dat";
	std::string file = od + "mesh_data/data_out_" + time_int_str + "_" + my_rank_str + prefix + ".csv";
	outfile.open(file.c_str(), std::ios::out);

	outfile<<std::setprecision(10);
	outfile<<"x,y,z,potential,e_x,e_y,charge_density";
	outfile << std::endl;
	for(int j = 0; j < my_data->my_y_size; j++)
	{
		for(int i = 0; i < rp->x_size; i++)
		{
			outfile << i*SPACE_STEP << "," << (j+my_data->beg_y)*SPACE_STEP << "," << "0," << p->GetPotentialAt(i,j) << 
							"," << p->GetElectricFieldXAt(i,j) << "," << p->GetElectricFieldYAt(i,j) << "," << p->GetChargeDensity(i,j) << std::endl;
		}
	}
	//std::cout << "my rank = " << my_data->my_rank << " output complete " << std::endl;
	outfile.close();
}

void output_species_data(mesh_data *my_data, FieldSolver *p, RunParameters *rp,std::vector<Species *> &species_list, int time_int, std::string od, std::string prefix)
{
	std::ofstream outfile;
	std::stringstream ss;
	ss << my_data->my_rank;
	std::string my_rank_str = ss.str();

	std::stringstream time_ss;
	time_ss << time_int;
	std::string time_int_str = time_ss.str();

	//std::string file = "data/data_out_" + my_rank_str + "_" + time_int_str + ".dat";
	std::string file = od + "species_data/species_out_" + time_int_str + "_" + my_rank_str + prefix + ".csv";
	outfile.open(file.c_str(), std::ios::out);

	outfile<<std::setprecision(10);
	outfile<<"x,y,z,";
	for(int q = 0; q < species_list.size(); q++)
	{
		if(q == (species_list.size() - 1))
			outfile << species_list[q]->GetName();
		else
			outfile << species_list[q]->GetName() << ",";
	}
	outfile << std::endl;
	//std::cout << "my rank = " << my_data->my_rank << " file = " << file << " y_size = " << my_data->my_y_size << std::endl;
	//
	int end_y_val = my_data->my_y_size;
	if(my_data->my_rank == (my_data->world_size-1))
		end_y_val--;
	for(int j = 0; j < end_y_val; j++)
	{
		for(int i = 0; i < rp->x_size-1; i++)
		{
			outfile << i*SPACE_STEP << "," << (j+my_data->beg_y)*SPACE_STEP << "," << "0,";
			for(int q = 0; q < species_list.size(); q++)
			{
				if(q == (species_list.size() - 1))
					outfile << species_list[q]->GetDensity(i,j+1);
				else
					outfile << species_list[q]->GetDensity(i,j+1) << ",";
			}
			outfile << std::endl;

		}
	}
	outfile.close();
}

void output_volumetric_source_term(mesh_data *my_data, RunParameters *rp, double **volumetric_source_term, int time_int, std::string pref, std::string od)
{
	std::ofstream outfile;
	std::stringstream ss;
	ss << my_data->my_rank;
	std::string my_rank_str = ss.str();

	std::stringstream time_ss;
	time_ss << time_int;
	std::string time_int_str = time_ss.str();

	std::string file = od + "source_term_data/" + pref + "_source_term_out_" + my_rank_str + "_" + time_int_str + ".dat";
	//std::string file = "data/source_term_out_" + my_rank_str + "_" + time_int_str + ".dat";
	outfile.open(file.c_str(), std::ios::out);

	for(int j = 0; j < my_data->my_y_size; j++)
	{
		for(int i = 0; i < rp->x_size; i++)
		{
			outfile << i*SPACE_STEP << "\t" << (j+my_data->beg_y)*SPACE_STEP << "\t" << volumetric_source_term[j][i];
			outfile << std::endl;
		}
	}
	outfile.close();
}

double interpolate_data(std::vector<double> &x, std::vector<double> &y, double value)
{
	if(value >= x[x.size()-2])
		return (y[y.size()-1]);
	else if(value <= x[1])
		return y[0];
	else
	{
		int before_index = 0;
		int after_index = 0;
		for(int i = 0; i < x.size(); i++)
		{
			if(x[i] > value)
			{
				before_index = i-1;
				after_index = i;
				break;
			}
		}

		return (y[before_index] + (value - x[before_index]) * (y[after_index] - y[before_index] ) 
															/ (x[after_index] - x[before_index] ) );
	}	
}

void populate_vectors_from_file(std::vector<double> &x, std::vector<double> &y, std::string filename)
{
	std::ifstream infile;
	infile.open(filename.c_str(), std::ios::in);
	std::string line;

	if(infile.is_open())
	{
		while(std::getline(infile, line))
		{
			if(!line.empty())
			{
				int pos = line.find("\t");
				std::string temp = line.substr(0,pos);
				x.push_back(std::stod(temp.c_str()));

				temp = line.substr(pos+1,line.size());
				y.push_back(std::stod(temp.c_str()));
			}
		}
		infile.close();
	}
	else
	{
		std::cout << "ERROR: could not open " << filename << std::endl;
	}
}

void populate_vectors_from_file(std::vector<double> &x, std::vector<double> &y, std::string filename, double x_scale, double y_scale)
{
	std::ifstream infile;
	infile.open(filename.c_str(), std::ios::in);
	std::string line;

	if(infile.is_open())
	{
		while(std::getline(infile, line))
		{
			if(!line.empty())
			{
				int pos = line.find("\t");
				std::string temp = line.substr(0,pos);
				x.push_back(atof(temp.c_str())*x_scale);

				temp = line.substr(pos+1,line.size());
				y.push_back(atof(temp.c_str())*y_scale);
			}
		}
		infile.close();
	}
	else
	{
		std::cout << "ERROR: could not open " << filename << std::endl;
	}
}

int get_current_processor_from_global_index(int j, mesh_information *global_mesh_info, int world_size)
{
	for(int i = 0; i < world_size; i++)
	{
		if(j >= global_mesh_info[i].beg_y && j < (global_mesh_info[i].beg_y + global_mesh_info[i].y_size) )
			return global_mesh_info[i].processor;
	}
	return -1;
}

std::string remove_white_space(std::string *in)
{
	std::string temp(in->data());
	std::string::iterator end_pos = std::remove(temp.begin(), temp.end(), ' ');
	temp.erase(end_pos, temp.end());

	end_pos = std::remove(temp.begin(), temp.end(), '\t');
	temp.erase(end_pos, temp.end());
	return temp;
}

void check_large_density(mesh_data *my_data, FieldSolver *p, RunParameters *rp, std::vector<Species *> &species_list, double current_timestep)
{
	for(int q = 0; q < species_list.size(); q++)
	{
		TransportData *tp = species_list[q]->GetTransportData();
		for(int j = 1; j < my_data->my_y_size+1; j++)
		{
			for(int i = 1; i < (rp->x_size-1); i++)
			{
				double mag_e_field = sqrt(pow(p->GetElectricFieldXAt(i,j-1),2) + pow(p->GetElectricFieldYAt(i,j-1),2));
				double reduced_e_field = mag_e_field / N_B / 1e-21;

				double u_m = species_list[q]->GetChargeSign() * tp->GetMobility(reduced_e_field) / N_B;
				double d_m = tp->GetDiffusionCoefficient(reduced_e_field) / N_B;

				
				//if(species_list[q]->GetDensity(i,j) > 1e21 && (species_list[q]->GetSpeciesType() == ELECTRON || species_list[q]->GetSpeciesType() == ION))
				if(mag_e_field > 1e9 && (species_list[q]->GetSpeciesType() == ELECTRON || species_list[q]->GetSpeciesType() == ION))
				{
					int flip = 0;
					if(species_list[q]->GetSpeciesType() == ELECTRON)
						flip = 1;
					std::cout << "LARGE DENSITY FOUND on processor " << my_data->my_rank << " type = " << species_list[q]->GetSpeciesType() << std::endl;
					std::cout << " > density = " << species_list[q]->GetDensity(i,j) << std::endl;
					std::cout << " > other density = " << species_list[flip]->GetDensity(i,j) << std::endl;
					std::cout << " > i,j = " << i << "," << j << std::endl;
					std::cout << " > e_field = " << mag_e_field << std::endl;
					std::cout << " > potential = " << p->GetPotentialAt(i,j-1) << std::endl;
					std::cout << " > reduced_e_field = " << reduced_e_field << std::endl;
					std::cout << " > mobility = " << u_m << std::endl;
					std::cout << " > diffusion = " << d_m << std::endl;
					p->PrintMatrixB();
				}
			}
		}
	}
}

void check_negative_density(mesh_data *my_data, FieldSolver *p, RunParameters *rp, std::vector<Species *> &species_list, double current_timestep)
{
	for(int q = 0; q < species_list.size(); q++)
	{
		TransportData *tp = species_list[q]->GetTransportData();
		for(int j = 1; j < my_data->my_y_size+1; j++)
		{
			for(int i = 1; i < (rp->x_size-1); i++)
			{
				double mag_e_field = sqrt(pow(p->GetElectricFieldXAt(i,j-1),2) + pow(p->GetElectricFieldYAt(i,j-1),2));
				double reduced_e_field = mag_e_field / N_B / 1e-21;

				double u_m = species_list[q]->GetChargeSign() * tp->GetMobility(reduced_e_field) / N_B;
				double d_m = tp->GetDiffusionCoefficient(reduced_e_field) / N_B;

				if(species_list[q]->GetDensity(i,j) < 0.0)
				{
					std::cout << "NEGATIVE DENSITY FOUND on processor " << my_data->my_rank << std::endl;
					std::cout << " > e_field = " << mag_e_field << std::endl;
					std::cout << " > reduced_e_field = " << reduced_e_field << std::endl;
					std::cout << " > mobility = " << u_m << std::endl;
					std::cout << " > diffusion = " << d_m << std::endl;
				}
			}
		}
	}
}

void scale_volumetric_by_inverse_volume(mesh_data *my_data, RunParameters *rp, double **volumetric_source_term)
{
	for(int j = 1; j < my_data->my_y_size+1; j++)
	{
		double r = (j+my_data->beg_y-1)*SPACE_STEP; // beginning of cell
		double r1 = (j+my_data->beg_y)*SPACE_STEP;		// end of cell
		double cell_volume = (PI * r1 * r1 * SPACE_STEP) - (PI * r * r * SPACE_STEP);

		for(int i = 0; i < rp->x_size; i++)
		{
			volumetric_source_term[j-1][i] = volumetric_source_term[j-1][i] / cell_volume;
		}
	}
}

void scale_volumetric_by_volume(mesh_data *my_data, RunParameters *rp, double **volumetric_source_term)
{
	for(int j = 1; j < my_data->my_y_size+1; j++)
	{
		double r = (j+my_data->beg_y-1)*SPACE_STEP; // beginning of cell
		double r1 = (j+my_data->beg_y)*SPACE_STEP;		// end of cell
		double cell_volume = (PI * r1 * r1 * SPACE_STEP) - (PI * r * r * SPACE_STEP);

		for(int i = 0; i < rp->x_size; i++)
		{
			volumetric_source_term[j-1][i] = volumetric_source_term[j-1][i] * cell_volume;
		}
	}
}

void print_reactions(mesh_data *my_data, std::vector<Species *> &species_list)
{
	if(my_data->my_rank == 0)
	{
		std::cout << "#####################" << std::endl;
		std::cout << "## SPECIES LIST    ##" << std::endl;
		std::cout << "#####################" << std::endl;
		for(int q = 0; q < species_list.size(); q++)
		{
			std::cout << "species name = " << species_list[q]->GetName() << " charge = " << species_list[q]->GetSpeciesCharge() << std::endl;
		}
		std::cout << "#####################" << std::endl;
		std::cout << "#####################" << std::endl;
		std::cout << std::endl;
		std::cout << "#####################" << std::endl;
		std::cout << "## REACTIONS       ##" << std::endl;
		std::cout << "#####################" << std::endl;
		for(int q = 0; q < species_list.size(); q++)
		{
			std::cout << "<<< species name = " << species_list[q]->GetName() << " >>>" << std::endl;
			std::cout << "GAIN REACTIONS" << std::endl;
			std::vector<Reaction *> g = species_list[q]->GetGainReactions();
			std::vector<Reaction *> l = species_list[q]->GetLossReactions();

			for(int s = 0; s < g.size(); s++)														// gain processes
			{
				std::vector<Species *> rs = g[s]->GetReactants();	
				std::cout << "reaction name = " << g[s]->GetReactionName() << std::endl;
				if(g[s]->GetReactionType() == FIXED_REACTION_RATE)
					std::cout << " > fixed reaction rate selected with rate = " << g[s]->GetReactionRate(0.0) << std::endl;
				else
					std::cout << " > file reaction rate selected" << std::endl;
				std::cout << " >> ";
				for(int z = 0; z < rs.size(); z++)
				{
					if(z != (rs.size() - 1))
						std::cout << rs[z]->GetName() << " + " ;
					else
						std::cout << rs[z]->GetName();
				}
				rs = g[s]->GetProducts();
				std::cout << " -> ";
				for(int z = 0; z < rs.size(); z++)
				{
					if(z != (rs.size() - 1))
						std::cout << rs[z]->GetName() << " + " ;
					else
						std::cout << rs[z]->GetName();
				}
				std::cout << std::endl;
			}
			std::cout << "LOSS REACTIONS" << std::endl;
			for(int s = 0; s < l.size(); s++)														// loss processes
			{
				std::vector<Species *> ls = l[s]->GetReactants();
				std::cout << "reaction name = " << l[s]->GetReactionName() << std::endl;
				if(l[s]->GetReactionType() == FIXED_REACTION_RATE)
					std::cout << " > fixed reaction rate selected with rate = " << l[s]->GetReactionRate(0.0) << std::endl;
				else
					std::cout << " > file reaction rate selected" << std::endl;
				std::cout << " >> ";
				for(int z = 0; z < ls.size(); z++)
				{
					if(z != (ls.size() - 1))
						std::cout << ls[z]->GetName() << " + ";
					else
						std::cout << ls[z]->GetName();
				}
				ls = l[s]->GetProducts();
				std::cout << " -> ";
				for(int z = 0; z < ls.size(); z++)
				{
					if(z != (ls.size() - 1))
						std::cout << ls[z]->GetName() << " + ";
					else
						std::cout << ls[z]->GetName();
				}
				std::cout << std::endl;
			}
			std::cout << std::endl;
		}
	}
}

int get_species_index(std::vector<Species *> &species_list, std::string to_find)
{
	int species_index = -1;
	for(int q = 0; q < species_list.size(); q++)
	{
		if(species_list[q]->GetName() == to_find)
		{
			species_index = q;
			break;
		}
	}
	return species_index;
}

double global_cfl(mesh_data *my_data, FieldSolver *p, RunParameters *rp,std::vector<Species *> &species_list, double current_timestep)
{
	double max_cfl = 0.0;
	int electron_index = 0;
	for(int q = 0; q < species_list.size(); q++)
	{
		if(species_list[q]->GetSpeciesType() == ELECTRON)
		{
			electron_index = q;
			break;
		}
	}
	int q = electron_index;
	TransportData *tp = species_list[q]->GetTransportData();

	for(int j = 1; j < my_data->my_y_size+1; j++)
	{
		for(int i = 1; i < (rp->x_size-1); i++)
		{
			double field_x_plus_surf = (p->GetElectricFieldXAt(i+1,j-1) + p->GetElectricFieldXAt(i+1,j)) / 2.0;
			double field_x_minus_surf = (p->GetElectricFieldXAt(i,j-1) + p->GetElectricFieldXAt(i,j)) / 2.0;
			double field_y_plus_surf = (p->GetElectricFieldYAt(i,j) + p->GetElectricFieldYAt(i+1,j)) / 2.0;
			double field_y_minus_surf = (p->GetElectricFieldYAt(i,j-1) + p->GetElectricFieldYAt(i+1,j-1)) / 2.0;

			double field_x_middle = (p->GetElectricFieldXAt(i,j-1) + p->GetElectricFieldXAt(i+1,j-1) +
															p->GetElectricFieldXAt(i,j) + p->GetElectricFieldXAt(i+1,j)) / 4.0;

			double field_y_middle = (p->GetElectricFieldYAt(i,j-1) + p->GetElectricFieldYAt(i+1,j-1) +
															p->GetElectricFieldYAt(i,j) + p->GetElectricFieldYAt(i+1,j)) / 4.0;

			//double mag_e_field = sqrt(pow(p->GetElectricFieldXAt(i,j-1),2) + pow(p->GetElectricFieldYAt(i,j-1),2));
			double mag_e_field = sqrt(pow(field_x_middle,2) + pow(field_y_middle,2));
			double reduced_e_field = mag_e_field / N_B / 1e-21;

			double u_m = species_list[q]->GetChargeSign() * tp->GetMobility(reduced_e_field) / N_B;
			double d_m = tp->GetDiffusionCoefficient(reduced_e_field) / N_B;
			double v_x = field_x_middle * u_m;
			double v_y = field_y_middle * u_m;

			double cfl = current_timestep * (sqrt(pow(v_x,2)+pow(v_y,2))) / SPACE_STEP;
			if(cfl > max_cfl)
				max_cfl = cfl;
		}
	}

	double max_recv = 0.0;
	MPI_Reduce(&max_cfl, &max_recv, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	return max_recv;
}

void output_directory_creation(mesh_data *my_data, RunParameters *rp)
{
	if(my_data->my_rank == 0)
	{
		std::string out_dir = rp->output_directory;
		std::vector<std::string> directories;
		directories.push_back("/intensity_data");
		directories.push_back("/mesh_data");
		directories.push_back("/source_term_data");
		directories.push_back("/species_data");

		if(mkdir(out_dir.c_str(), 0755) == -1)
		{
			if(errno == EEXIST)
			{
				std::cout << "Directory already exists: " << out_dir << std::endl;
			}
			else
			{
				std::cerr << "Error: " << strerror(errno) << std::endl;
			}

		}
		else
		{
			std::cout << "Directory create successfully: " << out_dir << std::endl;
		}
		for(int i = 0; i < directories.size(); i++)
		{
			std::string sub_dir = out_dir + directories[i];

			if(mkdir(sub_dir.c_str(), 0755) == -1)
			{
				if(errno == EEXIST)
				{
					std::cout<< "Directory already exists: " << sub_dir << std::endl;
				}
				else
				{
					std::cerr << "Error: "<< strerror(errno) << std::endl;
				}
			}
			else
			{
				std::cout << "Directory created successfully: " << sub_dir << std::endl;
			}
		}
	}

}
