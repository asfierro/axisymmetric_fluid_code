#ifndef RUN_PARAMETERS_H
#define RUN_PARAMETERS_H

#include <string>

struct RunParameters
{
	double voltage;
	int x_size;
	int y_size;
	std::string voltage_mode;
	double voltage_frequency;
	std::string output_directory;
	int restart;
	int output_stride;
	int status_stride;
	double final_time;
	std::string photoionization_method;
};

#endif
