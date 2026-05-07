#include "transport_data.h"

#include "utilities.h"

TransportData::TransportData()
{
	mobility_quantity_type = 0;
	diffusion_quantity_type = 0;
}

void TransportData::ReadDiffusionFile(std::string d_f)
{
	populate_vectors_from_file(reduced_field_diffusion, species_diffusion, d_f);
	diffusion_quantity_type = 1;
}

void TransportData::ReadMobilityFile(std::string m_f)
{
	populate_vectors_from_file(reduced_field_mobility, species_mobility, m_f);
	mobility_quantity_type = 1;
}

void TransportData::SetFixedQuantities(double d, double m)
{
	diffusion_fixed = d;
	mobility_fixed = m;

	diffusion_quantity_type = 2;
	mobility_quantity_type = 2;
}

double TransportData::GetDiffusionCoefficient(double in_field)
{
	if(diffusion_quantity_type == 2)
		return diffusion_fixed;

	////////////////////////////
	// check out of bounds for our data
	if(in_field < reduced_field_diffusion[0])
		return species_diffusion[0];
	if(in_field > reduced_field_diffusion[reduced_field_diffusion.size()-1])
		return species_diffusion[reduced_field_diffusion.size() - 1];
	////////////////////////////
	
	////////////////////////////
	// look up data
	return interpolate_data(reduced_field_diffusion, species_diffusion, in_field);
	/////////////////////////////
}

double TransportData::GetMobility(double in_field)
{
	if(mobility_quantity_type == 2)
		return mobility_fixed;

	////////////////////////////
	// check out of bounds for our data
	if(in_field < reduced_field_mobility[0])
		return species_mobility[0];
	if(in_field > reduced_field_mobility[reduced_field_mobility.size()-1])
		return species_mobility[reduced_field_mobility.size() - 1];
	////////////////////////////
	
	////////////////////////////
	// look up data
	return interpolate_data(reduced_field_mobility, species_mobility, in_field);
	/////////////////////////////
}

