#include "species.h"

#include "utilities.h"
#include "reaction.h"

Species::Species(std::string in_name, int x, int y, int i_type, RunParameters *new_rp)
{
	name = in_name;
	x_size = x;
	y_size = y;
	rp = new_rp;
	transport = new TransportData();

	density = new double*[y_size+2];
	old_density = new double*[y_size+2];
	old_old_density = new double*[y_size+2];
	
	charge = 0.0;
	charge_sign = 1.0;
	for(int i = 0; i < y_size+2; i++)
	{
		density[i] = new double[x_size];
		old_density[i] = new double[x_size];
		old_old_density[i] = new double[x_size];
	}

	for(int j = 0; j < y_size+2; j++)
	{
		for(int i = 0; i < x_size; i++)
		{
			density[j][i] = 0.0;
			old_density[j][i] = 0.0;
			old_old_density[j][i] = 0.0;
		}
	}
	
	species_type = i_type;

	if(species_type == ELECTRON || species_type == NEGATIVE_ION)
	{
		charge = -Q_E;
		charge_sign = -1.0;
	}
	else if(species_type == ION)
		charge = Q_E;

	density_after = new double[x_size];
	density_before = new double[x_size];
	old_density_after = new double[x_size];
	old_density_before = new double[x_size];
}

bool Species::isCharged()
{
	if(charge != 0.0)
		return true;
	else 
		return false;
}

void Species::CopyDataToOld()
{
	for(int j = 0; j < y_size+2; j++)
	{
		for(int i = 0; i < x_size; i++)
		{
			SetOldDensity(i,j,density[j][i]);
		}
	}

	for(int i = 0; i < x_size; i++)
	{
		old_density_before[i] = density_before[i];
		old_density_after[i] = density_after[i];
	}
}

void Species::TransferData(mesh_data *data)
{
	transfer_data(rp, old_density, data->my_rank, data->world_size, y_size);
	transfer_data(rp, density, data->my_rank, data->world_size, y_size);
	transfer_data_two_cell(rp, density, density_after, density_before, data->my_rank, data->world_size, y_size);
}

void Species::SetTransportDataFixed(double d_constant, double u_constant)
{
	transport->SetFixedQuantities(d_constant, u_constant);
}

void Species::SetTransportDataFile(std::string d_f, std::string m_f)
{
	transport->ReadDiffusionFile(d_f);
	transport->ReadMobilityFile(m_f);
}

void Species::StoreDensity()
{
	for(int j = 0; j < y_size+2; j++)
	{
		for(int i = 0; i < x_size; i++)
		{
			old_old_density[j][i] = density[j][i];
		}
	}
}

void Species::AverageDensity()
{
	for(int j = 0; j < y_size+2; j++)
	{
		for(int i = 0; i < x_size; i++)
		{
			density[j][i] = (density[j][i] + old_old_density[j][i]) / 2.0;
		}
	}
}

void Species::SetDensity(int x, int y, double value)
{
	density[y][x] = value;
}

void Species::SetOldDensity(int x, int y, double value)
{
	old_density[y][x] = value;
}

double Species::GetDensity(int x, int y)
{
	if(y == (y_size+2))
		return density_after[x];
	else if(y == -1)
		return density_before[x];
	else
		return density[y][x];
}

double Species::GetOldDensity(int x, int y)
{
	if(y == (y_size+2))
		return old_density_after[x];
	else if(y == -1)
		return old_density_before[x];
	else
		return old_density[y][x];
}



