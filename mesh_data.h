#ifndef _MESH_DATA_H
#define _MESH_DATA_H

struct mesh_data
{
	int my_rank;
	int my_y_size;
	int beg_y;
	int world_size;
};	

struct mesh_information
{
	int processor;
	int beg_y;
	int y_size;
};

#endif

