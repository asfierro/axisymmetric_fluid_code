#include "reaction.h"

#include "utilities.h"
#include "species.h"

Reaction::Reaction(std::string in_name, std::vector<Species *> &r, std::vector<Species *> &p, int r_type)
{
	reaction_name = in_name;

	reactants = r;
	products = p;

	reaction_type = r_type;
}

void Reaction::SetFixedReactionRate(double in_rate)
{
	reaction_type = FIXED_REACTION_RATE;
	fixed_reaction_rate = in_rate;
}

void Reaction::SetReactionDataFile(std::string filename)
{
	reaction_type = FILE_REACTION_RATE;
	populate_vectors_from_file(reaction_rate_x, reaction_rate_y, filename);
}

double Reaction::GetReactionRate(double in_field)
{
	if(reaction_type == FIXED_REACTION_RATE)
		return fixed_reaction_rate;
	
	//////////////////////////////////
	// check out of bounds for our data
	if(in_field < reaction_rate_x[0])
		return reaction_rate_y[0];
	if(in_field > reaction_rate_x[reaction_rate_x.size() - 1])
		return reaction_rate_y[reaction_rate_x.size() - 1];

	//////////////////////////////////
	// look up data
	return interpolate_data(reaction_rate_x, reaction_rate_y, in_field);	
	///////////////////////////////////
}

