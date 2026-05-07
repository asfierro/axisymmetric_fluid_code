#ifndef _CONSTANTS_H
#define _CONSTANTS_H

//#define XSIZE	 	1200
//#define YSIZE	 	1024

#define N_B					2.5e25					// background gas density
#define SPACE_STEP	5e-6
#define TIME_STEP		1e-13
#define Q_E					1.602e-19
#define EPSILON_0		8.854e-12
#define PI		3.14159
#define H_P					6.626e-34
#define c_light			3.0e8

#define 	AMESOS		0
#define		BELOS			1

#define NEUTRAL			0
#define ELECTRON		1
#define ION					2
#define NEGATIVE_ION		3

// REACTION TYPES
#define FIXED_REACTION_RATE		1
#define FILE_REACTION_RATE		2

#endif

