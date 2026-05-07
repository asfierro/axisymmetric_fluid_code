#ifndef REACTION_H
#define REACTION_H

#include <vector>
#include <string>

#include "constants.h"

class Species;

class Reaction
{
	public:
		Reaction(std::string, std::vector<Species *> &, std::vector<Species *> &, int);
		double GetReactionRate(double);
		void SetFixedReactionRate(double);
		void SetReactionDataFile(std::string);
		std::vector<Species *> GetReactants() { return reactants; }
		std::vector<Species *> GetProducts() { return products; }
		std::string GetReactionName() { return reaction_name; }
		int GetReactionType() { return reaction_type; }

	private:
		std::vector<Species *> reactants;
		std::vector<Species *> products;
	
		std::vector<double> reaction_rate_x;
		std::vector<double> reaction_rate_y;

		double fixed_reaction_rate;

		int reaction_type;
		std::string reaction_name;
};

#endif

