#include "input_deck.h"

#include "utilities.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <cmath>

bool InputDeck::OpenFile(std::ifstream &files)
{
	files.open(filename,std::ios::in);
	if(files.is_open())
		return true;
	else
		return false;
}

void InputDeck::CloseFile(std::ifstream &files)
{
	files.close();
}

std::string InputDeck::GetToken(std::string line, std::string token)
{
	int find = line.find(token);
	int find_next = line.find("=",find+token.length());
	int find_last = line.find(" ",find_next+2);
	
	return line.substr(find_next+2, find_last-find_next-2);
}

void InputDeck::ReadInputDeck(RunParameters *rp, mesh_data *my_data, std::vector<Species *> &species_list)
{
	if(my_data->my_rank == 0)
		std::cout << "reading input deck" << std::endl;
	std::ifstream infile;
	if(OpenFile(infile))
	{
		std::string line;
		while(std::getline(infile,line))
		{
			input_lines.push_back(line);	
			input_lines_read.push_back(0);
		}

		rp->voltage = ReadVoltage(my_data);
		rp->x_size = ReadXSize(my_data);
		rp->y_size = ReadYSize(my_data);
		ReadVoltageMode(my_data, rp);
		rp->output_directory = ReadOutputDirectory(my_data);
		rp->photoionization_method = ReadPhotoionizationMethod(my_data);
		rp->restart = ReadRestart(my_data);
		rp->output_stride = ReadOutputStride(my_data);
		rp->status_stride = ReadStatusStride(my_data);
		rp->final_time = ReadFinalTime(my_data);
		CloseFile(infile);
		output_directory_creation(my_data, rp);
	}
	else
	{
		std::cout << "could not open input deck file" << std::endl;
	}
}

void InputDeck::Finish(mesh_data *my_data, std::vector<Species *> &species_list, RunParameters *rp)
{
	ReadSpeciesInput(my_data, species_list, rp);	
	ReadSpeciesInteractions(my_data, species_list);
	ReadSpeciesInitialDensity(my_data, rp, species_list);
	CheckNonReadLines(my_data);
	std::cout << "outputting to directory " << rp->output_directory << std::endl;
	std::cout << "restart status = " << rp->restart << std::endl;
}

void InputDeck::CheckNonReadLines(mesh_data *my_data)
{
	if(my_data->my_rank == 0)
	{
		std::cout << "### CHECKING INPUT DECK FOR NON PARSED LINES ###" << std::endl;
		for(int i = 0; i < input_lines.size(); i++)
		{
			if(input_lines_read[i] == 0)
			{
				std::cout << input_lines[i] << std::endl;
			}
		}
		std::cout << "### END ###" << std::endl;
	}
}

int InputDeck::ReadStatusStride(mesh_data *my_data)
{
	if(my_data->my_rank == 0)
		std::cout << "reading status stride" << std::endl;
	int status_stride = 100;
	for(int i = 0; i < input_lines.size(); i++)
	{
		std::string line = input_lines[i];
		if(line.find("status stride") != std::string::npos)
		{
			status_stride = std::stoi(GetToken(line,"stride"));
			input_lines_read[i] = 1;
		}
	}
	return status_stride;
}

int InputDeck::ReadOutputStride(mesh_data *my_data)
{
	if(my_data->my_rank == 0)
		std::cout << "reading output stride" << std::endl;
	int output_stride = 1000;
	for(int i = 0; i < input_lines.size(); i++)
	{
		std::string line = input_lines[i];
		if(line.find("output stride") != std::string::npos)
		{
			output_stride = std::stoi(GetToken(line,"stride"));
			input_lines_read[i] = 1;
		}
	}
	return output_stride;
}

double InputDeck::ReadFinalTime(mesh_data *my_data)
{
	if(my_data->my_rank == 0)
		std::cout << "reading sim time" << std::endl;
	double final_time = 1e-9;
	for(int i = 0; i < input_lines.size(); i++)
	{
		std::string line = input_lines[i];
		if(line.find("final simulation time") != std::string::npos)
		{
			final_time = std::stod(GetToken(line,"time"));
			input_lines_read[i] = 1;
		}
	}
	return final_time;
}


int InputDeck::ReadRestart(mesh_data *my_data)
{
	if(my_data->my_rank == 0)
		std::cout << "reading restart" << std::endl;
	int restart_file = -1;
	for(int i = 0; i < input_lines.size(); i++)
	{
		std::string line = input_lines[i];
		if(line.find("restart step") != std::string::npos)
		{
			std::string restart_file_t = GetToken(line, "step");
			restart_file = std::stoi(restart_file_t);
			input_lines_read[i] = 1;
		}
	}
	return restart_file;
}

std::string InputDeck::ReadOutputDirectory(mesh_data *my_data)
{
	if(my_data->my_rank == 0)
		std::cout << "reading output directory" << std::endl;
	std::string out_dir;
	for(int i = 0; i < input_lines.size(); i++)
	{
		std::string line = input_lines[i];
		if(line.find("define output_directory") != std::string::npos)
		{
			out_dir = GetToken(line, "output_directory");
			input_lines_read[i] = 1;
		}
	}
	return out_dir;	
}

std::string InputDeck::ReadPhotoionizationMethod(mesh_data *my_data)
{
	if(my_data->my_rank == 0)
		std::cout << "reading photoionization method selection" << std::endl;
	std::string method;
	for(int i = 0; i < input_lines.size(); i++)
	{
		std::string line = input_lines[i];
		if(line.find("define photoionization_method") != std::string::npos)
		{
			method = GetToken(line, "photoionization_method");
			input_lines_read[i] = 1;
		}
	}
	return method;	
}

void InputDeck::ReadSpeciesInitialDensity(mesh_data *my_data, RunParameters *rp ,std::vector<Species *> &species_list)
{
	if(my_data->my_rank == 0)
		std::cout << "reading species initial densities" << std::endl;
	for(int i = 0; i < input_lines.size(); i++)
	{
		std::string line = input_lines[i];
		if(line.find("initial density") != std::string::npos)
		{
			std::string species_name = GetToken(line, "species");	
			std::string shape = GetToken(line, "shape");

			int si = -1;
			for(int q = 0; q < species_list.size(); q++)
			{
				if(species_list[q]->GetName() == species_name)
					si = q;
			}
			if(si != -1)
			{
				if(shape == "uniform")
				{
					double den = std::stod(GetToken(line,"value"));
					for(int i_c = 0; i_c < rp->x_size; i_c++)
					{
						for(int j_c = 0; j_c < my_data->my_y_size+1; j_c++)
						{
							species_list[si]->SetDensity(i_c, j_c, den);
						}
					}
				}
				else if(shape == "flux_uniform")
				{
					double den = std::stod(GetToken(line,"value"));

					double area_one = PI*SPACE_STEP*SPACE_STEP;    // end of cell

					for(int i_c = 0; i_c < rp->x_size; i_c++)
					{
						for(int j_c = 1; j_c < my_data->my_y_size+1; j_c++)
						{
							double r1 = (j_c+my_data->beg_y)*SPACE_STEP;    // end of cell
							double r = (j_c+my_data->beg_y-1)*SPACE_STEP; // beginning of cell
							double x_area = PI * r1 * r1 - PI * r * r;
							species_list[si]->SetDensity(i_c, j_c, den*(area_one/x_area));
						}
					}
				}
				else if(shape == "gaussian")
				{
					double peak_den = std::stod(GetToken(line,"peak"));
					double background_den = std::stod(GetToken(line,"background"));
					double x_cen = std::stod(GetToken(line,"x_center"));	
					double y_cen = std::stod(GetToken(line,"y_center"));	
					double sigma = std::stod(GetToken(line,"sigma"));	
					
					for(int i_c = 0; i_c < rp->x_size; i_c++)
					{
						for(int j_c = 1; j_c < my_data->my_y_size+1; j_c++)
						{
							double x_pos = i_c * SPACE_STEP;
							double y_pos = (j_c+my_data->beg_y) * SPACE_STEP;
							double val = background_den + peak_den*exp(-1.0*pow(x_pos-x_cen,2)/sigma - pow(y_pos-y_cen,2)/sigma);
							species_list[si]->SetDensity(i_c, j_c, val);
						}
					}
				}
				input_lines_read[i] = 1;
			}
		}
	}
}
int InputDeck::ReadXSize(mesh_data *my_data)
{
	if(my_data->my_rank == 0)
		std::cout<<"reading X size" << std::endl;
	int in_x_size = 0;
	for(int i = 0; i < input_lines.size(); i++)
	{
		std::string line = input_lines[i];
		if(line.find("define x_size") != std::string::npos)
		{
			std::string str_x_size = GetToken(line, "x_size");
			in_x_size = std::stoi(str_x_size);
			input_lines_read[i] = 1;
		}
	}	
	return in_x_size;	
}
int InputDeck::ReadYSize(mesh_data *my_data)
{
	if(my_data->my_rank == 0)
		std::cout<<"reading Y size"<< std::endl;
	int in_y_size = 0;
	for(int i = 0; i < input_lines.size(); i++)
	{
		std::string line = input_lines[i];
		if(line.find("define y_size") != std::string::npos)
		{
			std::string str_y_size = GetToken(line, "y_size");
			in_y_size = std::stoi(str_y_size);
			input_lines_read[i] = 1;
		}
	}
	return in_y_size;
}

double InputDeck::ReadVoltage(mesh_data *my_data)
{
	if(my_data->my_rank == 0)
		std::cout << "reading input voltage" << std::endl;
	double in_voltage = 0.0;
	for(int i = 0; i < input_lines.size(); i++)
	{
		std::string line = input_lines[i];
		if(line.find("define voltage") != std::string::npos)
		{
			std::string str_voltage = GetToken(line, "voltage");		
			in_voltage = std::stod(str_voltage);
			input_lines_read[i] = 1;	
		}
	}
	return in_voltage;
}

void InputDeck::ReadVoltageMode(mesh_data *my_data, RunParameters *rp)
{
	if(my_data->my_rank == 0)
		std::cout << "reading voltage mode" << std::endl;
	for(int i = 0; i < input_lines.size(); i++)
	{
		std::string line = input_lines[i];
		if(line.find("voltage mode") != std::string::npos)
		{
			rp->voltage_mode = GetToken(line, "mode");		
			if(rp->voltage_mode == "rf")
				rp->voltage_frequency = std::stod(GetToken(line, "frequency"));

			input_lines_read[i] = 1;	
		}
	}
}

void InputDeck::ReadSpeciesInteractions(mesh_data *my_data, std::vector<Species *> &species_list)
{
	if(my_data->my_rank == 0)
		std::cout << "reading species interactions" << std::endl;
	for(int i = 0; i < input_lines.size(); i++)
	{
		std::string line = input_lines[i];
		if(line.find("define gain_reaction") != std::string::npos || line.find("define loss_reaction") != std::string::npos)
		{
			std::string reactant_string = GetToken(line, "reactants");
			std::string product_string = GetToken(line, "products");
			std::string reaction_name = GetToken(line, "name");
			std::string rd_type = GetToken(line, "reaction_data");

			std::vector<Species *> reactants;
			std::vector<Species *> products;

			std::stringstream r_parse(reactant_string);
			std::stringstream p_parse(product_string);
			while(r_parse.good())
			{
				std::string substr;
				getline(r_parse, substr, ',');
				for(int q = 0; q < species_list.size(); q++)
				{
					if(substr == species_list[q]->GetName())
						reactants.push_back(species_list[q]);
				}		
			}
			while(p_parse.good())
			{
				std::string substr;
				getline(p_parse, substr, ',');
				for(int q = 0; q < species_list.size(); q++)
				{
					if(substr == species_list[q]->GetName())
						products.push_back(species_list[q]);
				}		
			}

			Reaction *r;
			if(rd_type == "file")
			{
				r = new Reaction(reaction_name, reactants, products, FILE_REACTION_RATE);	
				r->SetReactionDataFile(GetToken(line, "reaction_file"));
			}
			else if(rd_type == "fixed")
			{
				r = new Reaction(reaction_name, reactants, products, FIXED_REACTION_RATE);
				r->SetFixedReactionRate(std::stod(GetToken(line, "reaction_rate")));
			}

			int gain_or_loss_reaction = 1; // defaults to gain reaction
			if(line.find("define loss_reaction") != std::string::npos)
				gain_or_loss_reaction = 0;		// 0 means loss reaction
			for(int q = 0; q < species_list.size(); q++)
			{
				for(int z = 0; z < products.size(); z++)
				{
					if(species_list[q]->GetName() == products[z]->GetName())
					{
						if(gain_or_loss_reaction == 1)
							species_list[q]->AddGainReaction(r);
						else
							species_list[q]->AddLossReaction(r);
					}
				}
			}
			/*
			if(my_data->my_rank == 0)
			{
				std::cout << reactant_string << std::endl;
				for(int q = 0; q < reactants.size(); q++)
					std::cout << reactants[q]->GetName() << " - ";
				std::cout << std::endl;
				std::cout << product_string << std::endl;
				for(int q = 0; q < products.size(); q++)
					std::cout << products[q]->GetName() << " - ";
				std::cout << std::endl;
				std::cout << reaction_name << std::endl;
				std::cout << rd_type << std::endl;
			}*/
			input_lines_read[i] = 1;	
		}
	}
}

void InputDeck::ReadSpeciesInput(mesh_data *my_data, std::vector<Species *> &species_list, RunParameters *rp)
{
	if(my_data->my_rank == 0)
		std::cout << "reading species input" << std::endl;
	for(int i = 0; i < input_lines.size(); i++)
	{
		std::string line = input_lines[i];
		if(line.find("define species") != std::string::npos)
		{
			int species_length = species_list.size();
			std::string name = GetToken(line, "name");
			std::string type = GetToken(line, "type");
			std::string charge = GetToken(line, "charge");
			std::string tp_type = GetToken(line, "transport_data");
			int int_species_type = 0;
			if(type == "electron")
				int_species_type = ELECTRON;
			else if(type == "ion")
				int_species_type = ION;
			else if(type == "neutral")
				int_species_type = NEUTRAL;
			else if(type == "negative_ion")
				int_species_type = NEGATIVE_ION;
			std::string mobility = GetToken(line, "mobility");
			std::string diffusion = GetToken(line, "diffusion");
			species_list.push_back(new Species(name, rp->x_size, my_data->my_y_size, int_species_type, rp));
			if(tp_type == "fixed")
			{
				species_list[species_length]->SetTransportDataFixed(std::stod(diffusion), std::stod(mobility));
			}
			else if(tp_type == "file")
			{
				species_list[species_length]->SetTransportDataFile(diffusion, mobility);
			}
			input_lines_read[i] = 1;	
		}
	}
}

