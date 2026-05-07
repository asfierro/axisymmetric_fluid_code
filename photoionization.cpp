#include "photoionization.h"
#include "transport_data.h"

// initialize semi-empircal integral model
void semiempirical_photoionization_model_initialize(mesh_data *my_data, spectral_data *spectral_info, std::vector<double> &phi_array, RunParameters *rp)
{
	// read absorption data - raw data is in Angstrom vs cm^2, convert to nm and m^2
	absorption_data *o2_absorption = new absorption_data;
	populate_vectors_from_file(o2_absorption->total_x, o2_absorption->total_y, "absorption_data/fennelly_total_abs.dat",1,1);
	populate_vectors_from_file(o2_absorption->pi_x, o2_absorption->pi_y, "absorption_data/fennelly_o2_ion.dat",1,1);
	//populate_vectors_from_file(o2_absorption->total_x, o2_absorption->total_y, "absorption_data/o2_total_abs.txt",0.1,1e-4);
	//populate_vectors_from_file(o2_absorption->pi_x, o2_absorption->pi_y, "absorption_data/o2_photoionization.txt",0.1,1e-4);

	// read spectral data
	populate_vectors_from_file(spectral_info->wavelength, spectral_info->intensity, "spectra/N2_100mTorr.txt");
	for(int i = 0; i < spectral_info->wavelength.size(); i++)
	{
		spectral_info->pi_coefficient.push_back( N_B * 0.22 * interpolate_data(o2_absorption->pi_x, o2_absorption->pi_y, spectral_info->wavelength[i]));
		spectral_info->total_coefficient.push_back( N_B * 0.22 * interpolate_data(o2_absorption->total_x, o2_absorption->total_y, spectral_info->wavelength[i]));
	}

	double t_integral = 0;
	for(int i = 0; i < spectral_info->intensity.size(); i++)
	{
		spectral_info->intensity[i] *= 16693.9067 * 5e-4;                     // constant to calibrate by to get W/m
		double photon_energy = H_P * c_light / (spectral_info->wavelength[i] * 1e-9);
		spectral_info->intensity[i] = spectral_info->intensity[i] / photon_energy;          // Convert to photons/m/s
		// assuming volume of experimental plasma is something like 5 mm long by 1mm in radius - cylinder
		spectral_info->intensity[i] = spectral_info->intensity[i] / 1.57e-8;      // Convert to photons/m/s/m^3
		spectral_info->intensity_times_pi_coefficient.push_back (spectral_info->intensity[i] * spectral_info->pi_coefficient[i]);// now in photons/m/s/m^4
		//t_integral += spectral_info->intensity[i];
		//if(data->my_rank == 0)
			//std::cout << i << "\t" << spectral_info->wavelength[i] << "\t" <<  spectral_info->pi_coefficient[i] << std::endl;
	}

	if(my_data->my_rank == 0)
		std::cout << "photon tabulation" << std::endl;
	for(int i = 0; i < (5*std::max(rp->x_size,rp->y_size)); i++)
	{
		phi_array.push_back(0.0);
		for(int q = 0; q < spectral_info->wavelength.size(); q++)
		{
			phi_array[i] += spectral_info->intensity_times_pi_coefficient[q] * exp(-spectral_info->total_coefficient[q] * i*SPACE_STEP) *
				(spectral_info->wavelength[1] - spectral_info->wavelength[0])*1e-9; // multiply by dLambda to get photons/s/m^4
		}
		//if(data->my_rank == 0)
		//{
			//std::cout << "d = " << i * SPACE_STEP << "\t" << phi_array[i] << std::endl;
		//}
	}
	if(my_data->my_rank == 0)
	{
		for(int i = 0; i < phi_array.size(); i++)
		{
			//std::cout << "phi_array[" << i << "] = " << phi_array[i] << " d = " << i*SPACE_STEP << std::endl;
		}
	//std::cout << "phi_array at phi[0] = " << phi_array[0] << std::endl;
	//std::cout << "phi_array at phi[1] = " << phi_array[1] << std::endl;
	}
	MPI_Barrier(MPI_COMM_WORLD);
	if(my_data->my_rank == 0)
	{
		std::cout << "spectra integral = " << (t_integral * (spectral_info->wavelength[1] - spectral_info->wavelength[0]) * 1e-9) << std::endl;
		std::cout << "o2 absorption[0] = " << o2_absorption->pi_x[0] << "," << o2_absorption->pi_y[0] << std::endl;
		std::cout << "o2 total[0] = " << o2_absorption->total_x[0] << "," << o2_absorption->total_y[0] << std::endl;
		std::cout << "total coefficient[0] = " << spectral_info->wavelength[0] << "," << spectral_info->total_coefficient[0] << "," <<
                                                  interpolate_data(o2_absorption->total_x, o2_absorption->total_y, spectral_info->wavelength[0]) << std::endl;

  }
  MPI_Barrier(MPI_COMM_WORLD);
}


/* solves the Zheleznyak model integral directly */
void semiempirical_photoionization_model(mesh_data *my_data, std::vector<Species *> &species_list, RunParameters *rp, double **volumetric_source_term, FieldSolver *p, mesh_information *global_mesh_info,
													spectral_data *spectral_info, std::vector<double> &phi_array)
{
	double current_proc = 0;
	int my_data_counter = 0;
	double my_phi;
	int electron_index = 0;
	double mobility;
	TransportData *tp;
	int theta_division = 10;
	int z_division = 1;
	int r_division = 1;
	double d_phi = 2.0*PI/(double)theta_division;
	double local_dz = (SPACE_STEP*photo_skip) / (double)z_division;
	double local_dr = (SPACE_STEP*photo_skip) / (double)r_division;

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

	for(int j = 0; j < rp->y_size-1; j+=photo_skip)
	{
		current_proc = get_current_processor_from_global_index(j, global_mesh_info, my_data->world_size);
		//std::cout << "current proc = " << current_proc << ", j = " << j << ", my_proc = " << my_data->my_rank << ", my_data_counter = " << my_data_counter << std::endl;
		for(int i = 0; i < rp->x_size-1; i+=photo_skip)
		{
			double global_x = (i) * SPACE_STEP;
			double global_y = (j) * SPACE_STEP;
			double phi = 0.0;
			double r_global = (j)*SPACE_STEP; // beginning of cell
			double r1_global = (j+1)*SPACE_STEP;		// end of cell

			for(int local_j = 0; local_j < my_data->my_y_size; local_j+=photo_skip)		
			{
				double r_local = (local_j+my_data->beg_y)*SPACE_STEP; // beginning of cell
				double r1_local = (local_j+photo_skip+my_data->beg_y)*SPACE_STEP;		// end of cell
				
				// integration over theta, small piece of rotational volume
				//double cell_volume = ((PI * r1 * r1 * SPACE_STEP) - (PI * r * r * SPACE_STEP)) / (double)theta_division;
				for(int local_i = 0; local_i < rp->x_size-1; local_i+=photo_skip)
				{
					//if(local_i != i && (local_j+my_data->beg_y) != j)
					{
						double local_to_global_y = (local_j + my_data->beg_y) * SPACE_STEP;
						double local_to_global_x = (local_i) * SPACE_STEP;
						double field_x_middle = (p->GetElectricFieldXAt(local_i,local_j) + p->GetElectricFieldXAt(local_i+1,local_j) +
																p->GetElectricFieldXAt(local_i,local_j+1) + p->GetElectricFieldXAt(local_i+1,local_j+1)) / 4.0;

						double field_y_middle = (p->GetElectricFieldYAt(local_i,local_j) + p->GetElectricFieldYAt(local_i+1,local_j) +
																p->GetElectricFieldYAt(local_i,local_j+1) + p->GetElectricFieldYAt(local_i+1,local_j+1)) / 4.0;
					
						double mag_e_field = sqrt(pow(field_x_middle,2) + pow(field_y_middle,2));
						double t_electron_density = species_list[electron_index]->GetDensity(local_i, local_j+1);		// add one to get into species coordinates for y direction
						double reduced_e_field = mag_e_field / N_B / 1e-21;
						double u_m = tp->GetMobility(reduced_e_field) / N_B;

						// cell power is normalized by 1.2 W as was done in experiment
						// J dot E = sigma * E^2, sigma = n_e * q * u_m
						// assuming volume of experimental plasma is something like 5 mm long by 1mm in radius - cylinder
						double cell_power = pow(t_electron_density * u_m * Q_E * pow(mag_e_field,2),1.0) / (1.2 / 1.57e-8);
						//double reduction = 4.0 * PI * pow(theta_r, 2);		// account for distance loss

						if(cell_power < 0)
							cell_power = 0.0;

						for(int z_b = 0; z_b < z_division; z_b++)
						{
						for(int r_b = 0; r_b < r_division; r_b++)
						{
						// need to break up phi in both source and target volumes into an integration
						// Integration over theta direction in cylindrical coordinates 
						// lets call v1 the source
						for(int v1_theta = 0; v1_theta < theta_division; v1_theta++)
						{
							// lets call v2 the target
							//for(int v2_theta = 0; v2_theta < theta_division; v2_theta++)
							//{
								double a_theta = (v1_theta * d_phi - 0.0); //v2_theta * d_phi);
								double theta_r = sqrt(pow(local_to_global_y + r_b*local_dr,2) + pow(global_y,2) + 
										pow((local_to_global_x+z_b*local_dz) - global_x,2) - 2.0*global_y*(local_to_global_y+r_b*local_dr)*cos(a_theta));
								//double theta_r = sqrt(pow(local_to_global_y,2) + pow(global_y,2) + 
										//pow(local_to_global_x - (global_x),2) - 2.0*(global_y)*local_to_global_y*cos(a_theta));

								int phi_index = (int) (theta_r / SPACE_STEP);

								if(theta_r < 1e-9)
								{
									theta_r = SPACE_STEP/4.0;
									phi_index = 0;
								}

								//reduction = cell_power / reduction;

								// linearly interpolate phi
								if(theta_r > SPACE_STEP/2.0)
								{
									double phi_interp = phi_array[phi_index] + (theta_r - phi_index*SPACE_STEP) * (phi_array[phi_index+1] - phi_array[phi_index]) / 
												((phi_index+1)*SPACE_STEP - phi_index*SPACE_STEP);
 									phi += cell_power * phi_interp / (4.0 * PI) / (pow(theta_r,2)) * (local_to_global_y+r_b*local_dr) * d_phi * local_dz * local_dr;
									/*if(phi > 1e26)
									{
										std::cout << "my_rank = " << my_data->my_rank << "\tphi = " << phi << "\tbeg_y = " << my_data->beg_y << 
													"\tlocal_i = " << local_i << "\tlocal_j = " << local_j << "\tphi_interp = " << phi_interp << 
													"\tglobal_i = " << i  << "\tglobal_j = " << j << "\tmag_field = " << mag_e_field << "\tphi_array = " << phi_array[phi_index] << 
													"\tphi_inndex = " << phi_index << std::endl;
									}*/
 								//phi += (local_to_global_y+r_b*local_dr) * d_phi * local_dz * local_dr;
 								//phi += cell_power * phi_interp / (8.0 * PI) * (SPACE_STEP*photo_skip) * (pow(r1_local,2) - pow(r_local,2)) / (pow(theta_r,2)) *  d_phi;
								}
								else
								{
								//double phi_interp = phi_array[phi_index];
								//phi += cell_power * phi_interp / (8.0 * PI) * (SPACE_STEP*photo_skip) * (pow(r1_local,2) - pow(r_local,2)) / (pow(theta_r,2)) * (global_y+r_b*local_dr) * d_phi * local_dz * local_dr;
								}
						//}
							}
							}
							}
						}
				}
			}
			// Account for V1 (volume 1) in equation 1 of Zheleznyak
			// once again assumes that the global_cell_volume does not change appreciably over an increase in r, photo skip must be small enough
			// 12/7/22: 
			// stop multiplying by global_cell_volume because it wasn't dividing out in populate_from_volumetric_source term properly?  not sure what was happening
			//phi *= ((spectral_info->wavelength[1] - spectral_info->wavelength[0]) * 1e-9);// * global_cell_volume; //* pow(SPACE_STEP,2);
			MPI_Reduce(&phi, &my_phi, 1, MPI_DOUBLE, MPI_SUM, current_proc, MPI_COMM_WORLD);
				
			if(current_proc == my_data->my_rank)
			{
				//std::cout << "current proc = " << current_proc << std::endl;
				for(int skip_i = 0; skip_i < photo_skip; skip_i++)
				{
					for(int skip_j = 0; skip_j < photo_skip; skip_j++)
					{
						// once again assumes that the global_cell_volume does not change appreciably over an increase in r, photo skip must be small enough
						if( (i+skip_i) < rp->x_size && (my_data_counter+skip_j) < my_data->my_y_size)
							volumetric_source_term[my_data_counter+skip_j][i+skip_i] = my_phi;
					}
				}
			}
		}

	  if(current_proc == my_data->my_rank)
		{
			my_data_counter = my_data_counter + photo_skip;
			if(my_data_counter == my_data->my_y_size)
				my_data_counter = 0;
		}
	}
}

void fft_photoionization_model_initialize(mesh_data *my_data, std::vector<Species *> &species_list, fft_photoionization_model_data *fft_data, RunParameters *rp)
{
	//fft_data->l_1 = 0.0553*1e2;
	//fft_data->l_2 = 0.1460*1e2;
	fft_data->l_1 = 0.0974*1e2;
	fft_data->l_2 = 0.5877*1e2;
	fft_data->l_3 = 0.8900*1e2;

	//fft_data->A_1 = 1.986e-4*1e4;
	//fft_data->A_2 = 0.0051*1e4;
	fft_data->A_1 = 0.0021*1e4;
	fft_data->A_2 = 0.1775*1e4;
	fft_data->A_3 = 0.4886*1e4;

	fft_data->p_o2 = 150.0;

	// quenching factor + photoinization efficiency
	// squiggle * pq / (pq + p)
	fft_data->q_factor = 0.06 * 30.0 / (30.0 + 760.0);

	// Initialize RHS vector sizes, y = YSIZE-1 does not need to be solved for
	// or x = 0 or XSIZE-1, these are the BC's
	const ptrdiff_t size_array[2] = {rp->y_size-1, rp->x_size-2};
	ptrdiff_t local_y_size = my_data->my_y_size; // rp->y_size / my_data->world_size;

	if(my_data->my_rank == (my_data->world_size-1))
		local_y_size--;
	ptrdiff_t alloc_local = 0;
	ptrdiff_t my_fft_size = 0;
	ptrdiff_t my_fft_start = 0;

	//alloc_local = fftw_mpi_local_size_many(2, size_array, 1, local_y_size, MPI_COMM_WORLD, &my_fft_size, &my_fft_start);
	alloc_local = fftw_mpi_local_size_2d(size_array[0], size_array[1], MPI_COMM_WORLD, &my_fft_size, &my_fft_start);
	if(my_data->my_rank == 0)
		std::cout << "FFT PHOTOIONZATION RANK 0 size = " << alloc_local << "," << my_fft_size << " requested size = " << local_y_size << "," << my_fft_start << std::endl;
	if(my_data->my_rank == 1)
		std::cout << "FFT PHOTOIONZATION RANK 1 size = " << alloc_local << "," << my_fft_size << " requested size = " << local_y_size << "," << my_fft_start << std::endl;
	if(my_data->my_rank == (my_data->world_size-1))
		std::cout << "FFT PHOTOIONZATION RANK LAST = " << alloc_local << "," << my_fft_size << " requested size = " << local_y_size << "," << my_fft_start << std::endl;

	fft_data->rhs_1 = fftw_alloc_real(alloc_local);
	fft_data->rhs_2 = fftw_alloc_real(alloc_local);
	fft_data->rhs_3 = fftw_alloc_real(alloc_local);

	fft_data->lhs_1 = fftw_alloc_real(alloc_local);
	fft_data->lhs_2 = fftw_alloc_real(alloc_local);
	fft_data->lhs_3 = fftw_alloc_real(alloc_local);

	fft_data->rhs_output_1 = fftw_alloc_real(alloc_local);
	fft_data->rhs_output_2 = fftw_alloc_real(alloc_local);
	fft_data->rhs_output_3 = fftw_alloc_real(alloc_local);

	fft_data->freq_solution_1 = fftw_alloc_real(alloc_local);
	fft_data->freq_solution_2 = fftw_alloc_real(alloc_local);
	fft_data->freq_solution_3 = fftw_alloc_real(alloc_local);

	fft_data->fwd_plan_1 = fftw_mpi_plan_r2r_2d(rp->y_size-1, rp->x_size-2, fft_data->rhs_1, 
			fft_data->rhs_output_1, MPI_COMM_WORLD, FFTW_REDFT01, FFTW_RODFT00, FFTW_MEASURE);
	fft_data->fwd_plan_2 = fftw_mpi_plan_r2r_2d(rp->y_size-1, rp->x_size-2, fft_data->rhs_2, 
			fft_data->rhs_output_2, MPI_COMM_WORLD, FFTW_REDFT01, FFTW_RODFT00, FFTW_MEASURE);
	fft_data->fwd_plan_3 = fftw_mpi_plan_r2r_2d(rp->y_size-1, rp->x_size-2, fft_data->rhs_3, 
			fft_data->rhs_output_3, MPI_COMM_WORLD, FFTW_REDFT01, FFTW_RODFT00, FFTW_MEASURE);

	fft_data->bwd_plan_1 = fftw_mpi_plan_r2r_2d(rp->y_size-1, rp->x_size-2, fft_data->lhs_1,
			fft_data->freq_solution_1, MPI_COMM_WORLD, FFTW_REDFT10, FFTW_RODFT00, FFTW_ESTIMATE);
	fft_data->bwd_plan_2 = fftw_mpi_plan_r2r_2d(rp->y_size-1, rp->x_size-2, fft_data->lhs_2,
			fft_data->freq_solution_2, MPI_COMM_WORLD, FFTW_REDFT10, FFTW_RODFT00, FFTW_MEASURE);
	fft_data->bwd_plan_3 = fftw_mpi_plan_r2r_2d(rp->y_size-1, rp->x_size-2, fft_data->lhs_3,
			fft_data->freq_solution_3, MPI_COMM_WORLD, FFTW_REDFT10, FFTW_RODFT00, FFTW_MEASURE);

	fft_data->electron_index = get_species_index(species_list, "e-");
	std::vector<Reaction *> g = species_list[fft_data->electron_index]->GetGainReactions();
	fft_data->ionization_reaction = NULL;
	for(int q = 0; q < g.size(); q++)
	{
		if(g[q]->GetReactionName() == "n2_ionization")
		{
			fft_data->ionization_reaction = g[q];
			break;
		}
	}
}

void fft_photoionization_model(mesh_data *my_data, std::vector<Species *> &species_list, fft_photoionization_model_data *fft_data, double **volumetric_source_term, FieldSolver *p, RunParameters *rp, int ts_int)
{
  //update intensity everywhere
	int use_y_size = my_data->my_y_size;
	if(my_data->my_rank == (my_data->world_size-1))
		use_y_size--;

	for(int j = 0; j < use_y_size; j++)
	{
		for(int i = 1; i < rp->x_size-1; i++)
		{
			int index = j * (rp->x_size-2) + (i-1);
			double field_x_middle = (p->GetElectricFieldXAt(i,j) + p->GetElectricFieldXAt(i+1,j) +
                                p->GetElectricFieldXAt(i,j+1) + p->GetElectricFieldXAt(i+1,j+1)) / 4.0;

      double field_y_middle = (p->GetElectricFieldYAt(i,j) + p->GetElectricFieldYAt(i+1,j) +
                                p->GetElectricFieldYAt(i,j+1) + p->GetElectricFieldYAt(i+1,j+1)) / 4.0;

      double mag_e = sqrt(pow(field_x_middle,2) + pow(field_y_middle,2));
      double E_n = mag_e / N_B / 1e-21;
			
			/*
			double z_r = (double)i * SPACE_STEP;
			double r_r = ((double)j +(double) my_data->beg_y) * SPACE_STEP;
			double z_0 = (double)rp->x_size/2.0 * SPACE_STEP;
			double SD = 1e-8;
			double intensity = 3.5e28 * exp(-pow(z_r - z_0,2) / SD - pow(r_r,2) / SD);*/

			fft_data->rhs_1[index] = -1.0 * fft_data->A_1 * fft_data->p_o2 * fft_data->p_o2 * fft_data->q_factor *
										fft_data->ionization_reaction->GetReactionRate(E_n) * 
										species_list[fft_data->electron_index]->GetDensity(i, j+1) * N_B;

			fft_data->rhs_2[index] = -1.0 * fft_data->A_2 * fft_data->p_o2 * fft_data->p_o2 * fft_data->q_factor *
										fft_data->ionization_reaction->GetReactionRate(E_n) * 
										species_list[fft_data->electron_index]->GetDensity(i, j+1) * N_B;

			fft_data->rhs_3[index] = -1.0 * fft_data->A_3 * fft_data->p_o2 * fft_data->p_o2 * fft_data->q_factor *
										fft_data->ionization_reaction->GetReactionRate(E_n) * 
										species_list[fft_data->electron_index]->GetDensity(i, j+1) * N_B;

			//if(my_data->my_rank == 0 && j == 0)
				//std::cout << intensity << std::endl;
			//fft_data->rhs_1[index] = -1.0 * fft_data->A_1 * pow(fft_data->p_o2,2) * intensity;
			//fft_data->rhs_2[index] = -1.0 * fft_data->A_2 * pow(fft_data->p_o2,2) * intensity;
			//fft_data->rhs_3[index] = -1.0 * fft_data->A_3 * pow(fft_data->p_o2,2) * intensity;
		}
	}
	fftw_execute(fft_data->fwd_plan_1); 
	fftw_execute(fft_data->fwd_plan_2); 
	fftw_execute(fft_data->fwd_plan_3); 

	double L_z = (rp->x_size-1) * SPACE_STEP;
	double L_r = (rp->y_size-1) * SPACE_STEP;

	for(int j = 0; j < use_y_size; j++)
	{
		for(int i = 1; i < (rp->x_size-1); i++)
		{
			int index = j * (rp->x_size-2) + (i-1);

			double k_r = PI * (j+my_data->beg_y+0.5) / L_r;
			double k_z = PI * (i+1) / L_z;
			double denom1 = (-1.0) * (2.0*pow(k_r,2) + pow(k_z,2) + pow(fft_data->l_1 * fft_data->p_o2, 2));
			double denom2 = (-1.0) * (2.0*pow(k_r,2) + pow(k_z,2) + pow(fft_data->l_2 * fft_data->p_o2, 2));
			double denom3 = (-1.0) * (2.0*pow(k_r,2) + pow(k_z,2) + pow(fft_data->l_3 * fft_data->p_o2, 2));

			fft_data->lhs_1[index] = fft_data->rhs_output_1[index] / denom1;
			fft_data->lhs_2[index] = fft_data->rhs_output_2[index] / denom2;
			fft_data->lhs_3[index] = fft_data->rhs_output_3[index] / denom3;
		}
	}

	fftw_execute(fft_data->bwd_plan_1);
	fftw_execute(fft_data->bwd_plan_2);
	fftw_execute(fft_data->bwd_plan_3);

	double norm_factor = 4.0 * (rp->x_size-2) * (rp->y_size-1);
	for(int j = 0; j < use_y_size; j++)
	{
		for(int i = 1; i < (rp->x_size-1); i++)
		{
			int flat_index = j * (rp->x_size-2) + (i-1);
			//int full_index = j * rp->x_size + i;

			volumetric_source_term[j][i] = (fft_data->freq_solution_1[flat_index] + 
																			fft_data->freq_solution_2[flat_index] +
																			fft_data->freq_solution_3[flat_index]) / norm_factor;
			//volumetric_source_term[j][i] = fft_data->rhs_1[flat_index];
		}
	}
	// apply BC's
	for(int j = 0; j < my_data->my_y_size; j++)
	{
		volumetric_source_term[j][0] = 0.0;
		volumetric_source_term[j][rp->x_size-1] = 0.0;
	}

	if(my_data->my_rank == 0)
	{
		for(int i = 0; i < rp->x_size; i++)
		{
			volumetric_source_term[0][i] = volumetric_source_term[1][i];
		}
	}
	if(my_data->my_rank == (my_data->world_size-1))
	{
		for(int i = 0; i < rp->x_size; i++)
		{
			volumetric_source_term[my_data->my_y_size-1][i] = 0.0;
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
}


void simplified_photoionization_model(mesh_data *my_data, std::vector<Species *> &species_list, RunParameters *rp, double **volumetric_source_term, FieldSolver *p)
{
	// -------------- Find electron species for transport data
	TransportData *tp;
	// find electron's
	int electron_index = 0;
	for(int q = 0; q < species_list.size(); q++)
	{
		if(species_list[q]->GetSpeciesType() == ELECTRON)
		{
			electron_index = q;
			break;
		}
	}
	tp = species_list[electron_index]->GetTransportData();

	// ------------ Find peak electric fields and intensity on axis and broadcast
	int peak_z_index = 0;
	double I_0 = 0.0;
	std::vector<int> threshold_indexes;
	if(my_data->my_rank == 0)
	{
		for(int i = 0; i < rp->x_size-1; i++)	
		{
			double field_x_middle = (p->GetElectricFieldXAt(i,0) + p->GetElectricFieldXAt(i+1,0) +
	                   p->GetElectricFieldXAt(i,1) + p->GetElectricFieldXAt(i+1,1)) / 4.0;
			double field_y_middle = (p->GetElectricFieldYAt(peak_z_index,0) + p->GetElectricFieldYAt(peak_z_index+1,0) +
	                   p->GetElectricFieldYAt(peak_z_index,1) + p->GetElectricFieldYAt(peak_z_index+1,1)) / 4.0;

			double mag_e = sqrt(pow(field_x_middle,2) + pow(field_y_middle,2));

			if(mag_e  > 60e5)
			{
				threshold_indexes.push_back(i);
			}
		}
	}
	// ----------- Find n2 ionization reaction
	std::vector<Reaction *> g = species_list[electron_index]->GetGainReactions();
	Reaction *ionization_reaction = NULL;
	for(int q = 0; q < g.size(); q++)
	{
		if(g[q]->GetReactionName() == "n2_ionization")
			ionization_reaction = g[q];
	}

	if(ionization_reaction == NULL)
		return;

	/*
	if(my_data->my_rank == 0)
	{
		std::cout << "vector size = " << threshold_indexes.size() << " photo params: "; //peak_z_index << "\t" << I_0 << std::endl;
		for(int i = 0; i < threshold_indexes.size(); i++)
		{
			std::cout << threshold_indexes[i] << std::endl;
		}
	}*/

	// ------ get Ionization rate at peak field location
	std::vector<double> threshold_intensities;
	if(my_data->my_rank == 0)
	{
		for(int i = 0; i < threshold_indexes.size(); i++)
		{
			peak_z_index = threshold_indexes[i];
			double field_x_middle = (p->GetElectricFieldXAt(peak_z_index,0) + p->GetElectricFieldXAt(peak_z_index+1,0) +
	                   p->GetElectricFieldXAt(peak_z_index,1) + p->GetElectricFieldXAt(peak_z_index+1,1)) / 4.0;
	
			double field_y_middle = (p->GetElectricFieldYAt(peak_z_index,0) + p->GetElectricFieldYAt(peak_z_index+1,0) +
	                   p->GetElectricFieldYAt(peak_z_index,1) + p->GetElectricFieldYAt(peak_z_index+1,1)) / 4.0;
			double mag_e = sqrt(pow(field_x_middle,2) + pow(field_y_middle,2));
			double E_n = mag_e / N_B / 1e-21;

			// ------ intensity peak given by Bourdon, equation (2) assuming squiggly ~0.06 and pq = 30 torr and p = 760 torr
			I_0 = 0.06 * (0.037974) * ionization_reaction->GetReactionRate(E_n) * species_list[electron_index]->GetDensity(peak_z_index, 1) * N_B;
			threshold_intensities.push_back(I_0);
		}
	}
	
	/*
	if(my_data->my_rank == 0)
	{
		std::cout << "photo params: "; //peak_z_index << "\t" << I_0 << std::endl;
		for(int i = 0; i < threshold_indexes.size(); i++)
		{
			std::cout << threshold_intensities[i] << std::endl;
		}
		std::cout << std::endl;
	}*/

	double p_o2 = 168.0;
	int vector_size = threshold_indexes.size();
	// ----------- Broadcast z-indexes to all processors
	MPI_Bcast(&vector_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
	// ----------- Broadcast peak intensity to all processors
	int *peak_indexes = new int[vector_size];
	double *peak_intensities = new double[vector_size];
	if(my_data->my_rank == 0)
	{
		for(int q = 0; q < threshold_indexes.size(); q++)
		{
			peak_indexes[q] = threshold_indexes[q];
			peak_intensities[q] = threshold_intensities[q];
		}
	}
	MPI_Bcast(peak_indexes, vector_size, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(peak_intensities, vector_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);


	// Three Term Model from Bourdon, lambda and A given in paper
	// ----------- K coefficients - sqrt ((lambda * p_o2)^2 + (pi / l)^2 )
	double k1 = sqrt(pow(5.53*p_o2,2) + pow(PI/(rp->x_size*SPACE_STEP),2));
	double k2 = sqrt(pow(14.16*p_o2,2) + pow(PI/(rp->x_size*SPACE_STEP),2));
	double k3 = sqrt(pow(89.0*p_o2,2) + pow(PI/(rp->x_size*SPACE_STEP),2));


	// ----------- Calculate photoionization source term everywhere on my proc
	for(int i = 0; i < rp->x_size-1; i++)
	{
		for(int j = 0; j < my_data->my_y_size; j++)		
		{
			volumetric_source_term[j][i] = 0.0;
			double z_r = i * SPACE_STEP;
			double r_r = (j+my_data->beg_y) * SPACE_STEP;
			for(int q = 0; q < vector_size; q++)
			{
				peak_z_index = peak_indexes[q];
				I_0 = peak_intensities[q];

				// ----------- C Coefficients - A * p_o2^2 * I_0
				double C1 = 1.986*pow(p_o2,2)*I_0;
				double C2 = 51.0*pow(p_o2,2)*I_0;
				double C3 = 4886.0*pow(p_o2,2)*I_0;

				double i_term = exp( -(pow(z_r - (peak_z_index*SPACE_STEP),2) + pow(r_r,2)) / 5e-8);
				double s_ph1 = (pow(k1,2)*exp(k1*r_r) + pow(k1,2)*exp(-k1*r_r) + C1*i_term) / (pow(5.53*p_o2,2));
				double s_ph2 = (pow(k2,2)*exp(k2*r_r) + pow(k2,2)*exp(-k2*r_r) + C2*i_term) / (pow(14.16*p_o2,2));
				double s_ph3 = (pow(k3,2)*exp(k3*r_r) + pow(k3,2)*exp(-k3*r_r) + C3*i_term) / (pow(89.0*p_o2,2));
				volumetric_source_term[j][i] += s_ph1 + s_ph2 + s_ph3;
			}
		}
	}
	if(my_data->my_rank != 0)
	{
		if(vector_size == 1)
		{
			delete peak_indexes;
			delete peak_intensities;
		}
		else if(vector_size > 1)
		{
			delete [] peak_indexes;
			delete [] peak_intensities;
		}
	}
}

void simplified_photoionization_model_allprocs(mesh_data *my_data, std::vector<Species *> &species_list, double **volumetric_source_term, FieldSolver *p, RunParameters *rp, int ts_int)
{
	// -------------- Find electron species for transport data
	TransportData *tp;
	// find electron's
	int electron_index = 0;
	for(int q = 0; q < species_list.size(); q++)
	{
		if(species_list[q]->GetSpeciesType() == ELECTRON)
		{
			electron_index = q;
			break;
		}
	}
	tp = species_list[electron_index]->GetTransportData();

	// ----------- Find n2 ionization reaction
	std::vector<Reaction *> g = species_list[electron_index]->GetGainReactions();
	Reaction *ionization_reaction = NULL;
	for(int q = 0; q < g.size(); q++)
	{
		if(g[q]->GetReactionName() == "n2_ionization")
			ionization_reaction = g[q];
	}

	if(ionization_reaction == NULL)
		return;

	// ------------ Find electric fields that meet threshold
	double I_0 = 0.0;
	struct xy
	{
		int z;
		int r;
		double intensity;
	};

	std::vector<xy> threshold_indexes;
	if(my_data->my_rank < 10000)
	{
	for(int i = 0; i < rp->x_size-1; i++)	
	{
		for(int j = 0; j < my_data->my_y_size; j++)
		{
			xy temp_coord;
			temp_coord.r = j + my_data->beg_y; // store global coordinate
			temp_coord.z = i;

			double field_x_middle = (p->GetElectricFieldXAt(i,j) + p->GetElectricFieldXAt(i+1,j) +
	                  p->GetElectricFieldXAt(i,j+1) + p->GetElectricFieldXAt(i+1,j+1)) / 4.0;
			double field_y_middle = (p->GetElectricFieldYAt(i,j) + p->GetElectricFieldYAt(i+1,j) +
	                  p->GetElectricFieldYAt(i,j+1) + p->GetElectricFieldYAt(i+1,j+1)) / 4.0;

			double En = sqrt(pow(field_x_middle,2) + pow(field_y_middle,2)) / N_B / 1e-21;

			temp_coord.intensity = 0.06 * (0.037974) * ionization_reaction->GetReactionRate(En) * species_list[electron_index]->GetDensity(temp_coord.z, temp_coord.r+1-my_data->beg_y) * N_B;

			if(temp_coord.intensity  > 1e23)	 // approx 50 when collisional ionization rate > 1e-21 m^3/s
			{
				threshold_indexes.push_back(temp_coord);
			}
		}
	}
	}
	/*
	// ------ get Ionization rate at peak field location
	for(int i = 0; i < threshold_indexes.size(); i++)
	{
		int peak_z_index = threshold_indexes[i].z;
		int peak_r_index = threshold_indexes[i].r - my_data->beg_y;	 // convert back to local coordinates
		double field_x_middle = (p->GetElectricFieldXAt(peak_z_index,peak_r_index) + p->GetElectricFieldXAt(peak_z_index+1,peak_r_index) +
	                  p->GetElectricFieldXAt(peak_z_index,peak_r_index+1) + p->GetElectricFieldXAt(peak_z_index+1,peak_r_index+1)) / 4.0;

		double field_y_middle = (p->GetElectricFieldYAt(peak_z_index,peak_r_index) + p->GetElectricFieldYAt(peak_z_index+1,peak_r_index) +
                   p->GetElectricFieldYAt(peak_z_index,peak_r_index+1) + p->GetElectricFieldYAt(peak_z_index+1,peak_r_index+1)) / 4.0;
		double mag_e = sqrt(pow(field_x_middle,2) + pow(field_y_middle,2));
		double E_n = mag_e / N_B / 1e-21;

		// ------ intensity peak given by Bourdon, equation (2) assuming squiggly ~0.06 and pq = 30 torr and p = 760 torr
		I_0 = 0.06 * (0.037974) * ionization_reaction->GetReactionRate(E_n) * species_list[electron_index]->GetDensity(peak_z_index, peak_r_index+1) * N_B;
		threshold_indexes[i].intensity = I_0;
	}*/

	double p_o2 = 168.0;
	int vector_size = threshold_indexes.size();
	// ----------- Broadcast number of indexes that meet threshold from all processors to all processors
	int *processor_threshold = new int[my_data->world_size];
	//std::cout << my_data->my_rank << "\t" << vector_size << std::endl;
	for(int i = 0; i < my_data->world_size; i++)
		processor_threshold[i] = 0;

	MPI_Allgather(&vector_size, 1, MPI_INT, processor_threshold, 1, MPI_INT, MPI_COMM_WORLD);
	int *displacements = new int[my_data->world_size];
	int number_of_elements = 0;
	for(int i = 0; i < my_data->world_size; i++)
	{
		number_of_elements += processor_threshold[i];
		displacements[i] = 0;
	}
	for(int i = 1; i < my_data->world_size; i++)
	{
		for(int j = 0; j < i; j++)
			displacements[i] += processor_threshold[j];
	}
	if((ts_int % rp->status_stride) == 0 && my_data->my_rank == 0)
		std::cout << "elements contributing to photoionization = " << number_of_elements << std::endl;

	int *all_indexes_r;
	int *all_indexes_z;
	double *all_intensities;
	if(number_of_elements > 0)
	{
		all_indexes_r = new int[number_of_elements];
		all_indexes_z = new int[number_of_elements];
		all_intensities = new double[number_of_elements];
	}

	// ----------- Broadcast peak intensity to from all processors to all processors
	int *peak_indexes_r;
	int *peak_indexes_z;
	double *peak_intensities;
	if(vector_size > 0)	
	{
		peak_indexes_r = new int[vector_size];
		peak_indexes_z = new int[vector_size];
		peak_intensities = new double[vector_size];
	}
	for(int q = 0; q < vector_size; q++)
	{
		peak_indexes_r[q] = threshold_indexes[q].r;
		peak_indexes_z[q] = threshold_indexes[q].z;
		peak_intensities[q] = threshold_indexes[q].intensity;
	}
	//std::cout << "commiting all gather" << std::endl;
	if(number_of_elements > 0)
	{
		MPI_Allgatherv(peak_indexes_r, vector_size, MPI_INT, all_indexes_r, processor_threshold, displacements, MPI_INT, MPI_COMM_WORLD);
		MPI_Allgatherv(peak_indexes_z, vector_size, MPI_INT, all_indexes_z, processor_threshold, displacements, MPI_INT, MPI_COMM_WORLD);
		MPI_Allgatherv(peak_intensities, vector_size, MPI_DOUBLE, all_intensities, processor_threshold, displacements, MPI_DOUBLE, MPI_COMM_WORLD);


		//std::cout << "photo calculation" << std::endl;
		// Three Term Model from Bourdon, lambda and A given in paper
		// ----------- K coefficients - sqrt ((lambda * p_o2)^2 + (pi / l)^2 )
		double k1 = sqrt(pow(5.53*p_o2,2) + pow(PI/(rp->x_size*SPACE_STEP),2));
		double k2 = sqrt(pow(14.16*p_o2,2) + pow(PI/(rp->x_size*SPACE_STEP),2));
		double k3 = sqrt(pow(89.0*p_o2,2) + pow(PI/(rp->x_size*SPACE_STEP),2));

		// ----------- Calculate photoionization source term everywhere on my proc
		double exp_const1 = 1.0 / (exp(k1*(rp->y_size*SPACE_STEP)) + exp(-k1*(rp->y_size*SPACE_STEP)));
		double exp_const2 = 1.0 / (exp(k2*(rp->y_size*SPACE_STEP)) + exp(-k2*(rp->y_size*SPACE_STEP)));
		double exp_const3 = 1.0 / (exp(k3*(rp->y_size*SPACE_STEP)) + exp(-k3*(rp->y_size*SPACE_STEP)));
		for(int i = 0; i < rp->y_size-1; i++)
		{
			for(int j = 0; j < my_data->my_y_size; j++)		
			{
				volumetric_source_term[j][i] = 0.0;
				double z_r = i * SPACE_STEP;
				double r_r = (j+my_data->beg_y) * SPACE_STEP;
				double i_term = 0.0;

				for(int q = 0; q < number_of_elements; q++)
				{
					int peak_z_index = all_indexes_z[q];
					int peak_r_index = all_indexes_r[q];
					I_0 = all_intensities[q];

					// ----------- C Coefficients - A * p_o2^2 * I_0
					//double C1 = 1.986*p_o2*p_o2*I_0;
					//double C2 = 51.0*p_o2*p_o2*I_0;
					//double C3 = 4886.0*p_o2*p_o2*I_0;
					

					// Gaussian source
					//i_term += I_0 * exp( -((z_r-(peak_z_index*SPACE_STEP))*(z_r-(peak_z_index*SPACE_STEP)) 
															 //+ (r_r-(peak_r_index*SPACE_STEP))*(r_r-(peak_r_index*SPACE_STEP))) / 5e-8);
					// Laplace Source
					//i_term += I_0 * exp( -fabs(z_r-(peak_z_index*SPACE_STEP))/8e-6 - fabs(r_r-(peak_r_index*SPACE_STEP))/8e-6 ) ;

					// Cauchy (Lorentz) Source
					i_term += I_0 / (1.0 + pow((z_r-(peak_z_index*SPACE_STEP))/2e-6,2) + pow((r_r-(peak_r_index*SPACE_STEP))/2e-6,2)) ;
					
					
															 /*
					double s_ph1 = (k1*k1*exp(k1*r_r)*exp_const1 + k1*k1*exp(-k1*r_r)*exp_const1 + C1*i_term) / ((5.53*5.53*p_o2*p_o2));
					double s_ph2 = (k2*k2*exp(k2*r_r)*exp_const2 + k2*k2*exp(-k2*r_r)*exp_const2 + C2*i_term) / ((14.16*14.16*p_o2*p_o2));
					double s_ph3 = (k3*k3*exp(k3*r_r)*exp_const3 + k3*k3*exp(-k3*r_r)*exp_const3 + C3*i_term) / ((89.0*89.0*p_o2*p_o2));
					volumetric_source_term[j][i] += s_ph1 + s_ph2 + s_ph3;*/
				}
				double s_ph1 = (k1*k1*exp(k1*r_r)*exp_const1 + k1*k1*exp(-k1*r_r)*exp_const1 + 1.986*p_o2*p_o2*i_term) / ((5.53*5.53*p_o2*p_o2));
				double s_ph2 = (k2*k2*exp(k2*r_r)*exp_const2 + k2*k2*exp(-k2*r_r)*exp_const2 + 51.0*p_o2*p_o2*i_term) / ((14.16*14.16*p_o2*p_o2));
				double s_ph3 = (k3*k3*exp(k3*r_r)*exp_const3 + k3*k3*exp(-k3*r_r)*exp_const3 + 4886.0*p_o2*p_o2*i_term) / ((89.0*89.0*p_o2*p_o2));
				volumetric_source_term[j][i] += s_ph1 + s_ph2 + s_ph3;
			}
		}
	}
	//std::cout << "deleting memory" << std::endl;
	if(vector_size == 1)
	{
		delete peak_indexes_r;
		delete peak_indexes_z;
		delete peak_intensities;
	}
	else if(vector_size > 1)
	{
		delete [] peak_indexes_r;
		delete [] peak_indexes_z;
		delete [] peak_intensities;
	}
	if(number_of_elements == 1)
	{
		delete all_indexes_r;
		delete all_indexes_z;
		delete all_intensities;
	}
	else if(number_of_elements > 1)
	{
		delete [] all_indexes_r;
		delete [] all_indexes_z;
		delete [] all_intensities;
	}

	delete [] processor_threshold;
}




