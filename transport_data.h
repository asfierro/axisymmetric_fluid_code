#ifndef TRANSPORT_DATA_H
#define TRANSPORT_DATA_H

#include <ostream>
#include <fstream>
#include <vector>
#include <string>

class TransportData
{
	public:
		TransportData();
		void ReadDiffusionFile(std::string);
		void ReadMobilityFile(std::string);
		void SetFixedQuantities(double, double);
		double GetDiffusionCoefficient(double);
		double GetMobility(double);

	private:
		std::vector<double> reduced_field_diffusion;
		std::vector<double> species_diffusion;

		std::vector<double> reduced_field_mobility;
		std::vector<double> species_mobility;

		int diffusion_quantity_type;
		int mobility_quantity_type;
	
		double diffusion_fixed;
		double mobility_fixed;

};

#endif

