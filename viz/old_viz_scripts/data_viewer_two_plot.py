import math
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import ticker

num_procs = 128
x_size = 0
y_size = 0
space_step = 0
skip = 1
time_slices = 10000
file_begin = 50000
file_stride = 10000

constant_file = open("constants.h","r")
for line in constant_file:
	line = line.strip().split()
	
	if(len(line) > 0 and len(line) == 3):
		if(line[1] == "XSIZE"):
			x_size = int(line[2]);
		if(line[1] == "YSIZE"):
			y_size = int(line[2]);
		if(line[1] == "SPACE_STEP"):
			space_step = float(line[2]);
	
print("x_size = " + str(x_size))
print("y_size = " + str(y_size))
print("space_step = " + str(space_step))
constant_file.close()

x = np.linspace(0,x_size*space_step,x_size)
y = np.linspace(0,y_size*space_step,y_size)
Z_potential = np.ones((y_size,x_size))
Z_field_x = np.ones((y_size,x_size))
Z_field_y = np.ones((y_size,x_size))
Z_field_mag = np.ones((y_size,x_size))
Z_electron_density = np.ones((y_size,x_size))
Z_ion_density = np.ones((y_size,x_size))
volumetric_term = np.ones((y_size,x_size))

X_g, Y_g = np.meshgrid(x, y)

for i in range(time_slices):
	c_time_slice = i*file_stride + file_begin;
	proc_y_size = int(y_size / num_procs)
	y_size_add = y_size % num_procs

	cur_y_pos = 0
	for n in range(num_procs):
		if n == (num_procs-1):
			y_size = y_size + y_size_add
		infile = open("/users/asfierro/wheeler-scratch/2d_fluid_out/data/data_out_" + str(n) + "_" + str(c_time_slice) + ".dat","r");
		#infile = open("data/data_out_" + str(n) + "_" + str(c_time_slice) + ".dat","r");
		print("opening file: " + str(n) + " with time plane = " + str(c_time_slice))
		cur_x_pos = 0
		for line in infile:
			line = line.strip().split()

			Z_potential[cur_y_pos][cur_x_pos] = float(line[2])
			Z_field_x[cur_y_pos][cur_x_pos] = float(line[3])
			Z_field_y[cur_y_pos][cur_x_pos] = float(line[4])
			Z_electron_density[cur_y_pos][cur_x_pos] = float(line[5])
			Z_ion_density[cur_y_pos][cur_x_pos] = float(line[6])

			Z_field_mag[cur_y_pos][cur_x_pos] = np.log10(math.sqrt(Z_field_x[cur_y_pos][cur_x_pos]**2 + Z_field_y[cur_y_pos][cur_x_pos]**2))
			if(Z_electron_density[cur_y_pos][cur_x_pos] > 0.0):
				Z_electron_density[cur_y_pos][cur_x_pos] = np.log10(Z_electron_density[cur_y_pos][cur_x_pos])
			if(Z_ion_density[cur_y_pos][cur_x_pos] > 0.0):
				Z_ion_density[cur_y_pos][cur_x_pos] = np.log10(Z_ion_density[cur_y_pos][cur_x_pos])

			if Z_potential[cur_y_pos][cur_x_pos] < 0:
				print(str(cur_x_pos) + "\t" + str(cur_y_pos) + "\t" + str(Z_potential[cur_y_pos][cur_x_pos]))
		
			cur_x_pos = cur_x_pos + 1
			if cur_x_pos == x_size:
				cur_x_pos = 0
				cur_y_pos = cur_y_pos + 1

		infile.close();
	cur_y_pos = 0
	for n in range(num_procs):
		if n == (num_procs-1):
			y_size = y_size + y_size_add
		infile = open("/users/asfierro/wheeler-scratch/2d_fluid_out/source_term_data/source_term_out_" + str(n) + "_" + str(c_time_slice) + ".dat","r");
		#infile = open("data/source_term_out_" + str(n) + "_" + str(c_time_slice) + ".dat","r");
		print("opening file: " + str(n) + " with time plane = " + str(c_time_slice))
		cur_x_pos = 0
		for line in infile:
			line = line.strip().split()
			volumetric_term[cur_y_pos][cur_x_pos] = float(line[2])
			if(volumetric_term[cur_y_pos][cur_x_pos] > 0.0):
				volumetric_term[cur_y_pos][cur_x_pos] = np.log10(volumetric_term[cur_y_pos][cur_x_pos])
			cur_x_pos = cur_x_pos + 1
			if cur_x_pos == x_size:
				cur_x_pos = 0
				cur_y_pos = cur_y_pos + 1

		infile.close()


	num_levels = np.linspace(0,2000,15)
	plt.figure(figsize=(12,3))
	#num_levels = np.linspace(1.0e8, 1.0e19, 50)
	#num_levels = [1e8, 1e9, 1e10, 1e11, 1e12, 1e13, 1e14, 1e15, 1e16, 1e17, 1e18]
	num_levels = np.linspace(8,20,20)
	plt.subplot(1,2,1)
	plt.contourf(X_g,Y_g,Z_electron_density, levels=num_levels, extend='both', cmap="plasma")
	plt.colorbar()

	plt.subplot(1,2,2)
	plt.contourf(X_g, Y_g, volumetric_term, extend='both')
	plt.colorbar()

	#infile = open("data/data_out_" + str(n) + ".dat","r");
	plt.savefig("images/two_plot_img_out_" + str(c_time_slice) + ".png", dpi=200)
	plt.clf()
	plt.close()

