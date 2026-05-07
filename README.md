# 2d fluid axisymmetric

contact: andrew.fierro@nmt.edu

Build Environment Requirements:
	- MPI (parallelization)
	- FFTW (fast fourier transforms)
	- Trilinos (linear solver for potential)

see example build script

This code has been published in 

A. Fierro, A. Alibalazadeh, J. Stephens, C. Moore, "Massively parallel axisymmetric fluid model for streamer discharges,"
Comp. Phys. Comm., 2025.

and

A. Alibalazadeh, A. Fierro, M. Gilmore, 
"Enhancing Photoionization Rate Calculations in Low-Temperature Plasmas Using Spectral Methods", Comp. Phys. Comm., 2026.

run with mpirun -np <procs> ./main air.in


mesh structure is different between quantities in the field solver and
the mesh quantities in the mesh_data.h, row vs. column ordered for data processing

in mesh_data.h, each quantity is made with +2 size to allow for data transfers between
previous and next processor.  For example, electron density may look like this
if the my_y_size is equal to 4, the actual number of rows will be 6. 

```		
[ 0 		0 		0			0			0			0				< -------- buffer row (y = 5)				^
  1e6		1e6		1e6		1e6		1e6		1e6			(y = 4)															|
  1e7		1e7		1e7		1e7		1e7		1e7			(y = 3)															|
  1e7		1e7		1e7		1e7		1e7		1e7			(y = 2)															|	y-direction
  1e7		1e7		1e7		1e7		1e7		1e7			(y = 1)															|
  0 		0 		0			0			0			0 ]			< -------- buffer row (y = 0)				|

	x direction ------------->
```
rows 0 and 5 are data from the previous and next processors.  Rows 1 and 4
are passed to previous and next processors in their corresponding buffer rows.
to access an element, you would use [y_index][x_index]

Field solver class stores data differently, to access an element you use [x_index][y_index] but 
does not contain buffer elements directly.  Field solver transfers buffer rows and stores them
separately, so if you request +1 my_y_size outside of a processors domain, that information is 
available to the current processor and the request is handled appropriately.


