#ifndef SPECIES_H
#define SPECIES_H

#include <string>
#include <vector>

#include "transport_data.h"
#include "constants.h"
#include "mesh_data.h"
#include "run_parameters.h"

class Reaction;

class Species
{
	public:
		Species(std::string, int, int, int, RunParameters *);
		double GetDensity(int, int);
		double GetDensity(int, int, mesh_data *);
		double GetOldDensity(int, int);
		void SetTransportDataFile(std::string, std::string);
		void SetTransportDataFixed(double, double);
		TransportData *GetTransportData() { return transport; }
		void SetDensity(int, int, double);
		void SetOldDensity(int, int, double);
		int GetSpeciesType() { return species_type; }
		double GetSpeciesCharge() { return charge; }
		void CopyDataToOld();
		void TransferData(mesh_data *);
		std::string GetName() { return name; }
		void AverageDensity();
		void StoreDensity();
		double GetBeginningTimestepDensity(int i, int j) { return old_old_density[j][i]; }

		double GetCharge() { return charge; }
		double GetChargeSign() { return charge_sign; }

		bool isCharged();

		void AddGainReaction(Reaction *g) { gain_reactions.push_back(g); }
		void AddLossReaction(Reaction *l) { loss_reactions.push_back(l); }

		std::vector<Reaction *> GetGainReactions() { return gain_reactions; }
		std::vector<Reaction *> GetLossReactions() { return loss_reactions; }

	private:
		double **old_old_density;
		double **old_density;
		double **density;
		std::string name;
		TransportData *transport;
		RunParameters *rp;
		
		double *density_after;
		double *density_before;
		double *old_density_after;
		double *old_density_before;
		std::vector<Reaction *> gain_reactions;
		std::vector<Reaction *> loss_reactions;

		int x_size;
		int y_size;
		
		double charge;
		double charge_sign;
		int species_type;
};

#endif

