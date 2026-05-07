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
time_slices = 160
file_begin = 0
file_stride = 10000
x_size = 2500
y_size = 1024
space_step = 5e-6
'''
constant_file2 = open("../constants.h","r")
for line in constant_file:
	line = line.strip().split()
	
	if(len(line) > 0 and len(line) == 3):
		if(line[1] == "XSIZE"):
			x_size = int(line[2]);
		if(line[1] == "YSIZE"):
			y_size = int(line[2]);
		if(line[1] == "SPACE_STEP"):
			space_step = float(line[2]);
	
constant_file.close()
'''
print("x_size = " + str(x_size))
print("y_size = " + str(y_size))
print("space_step = " + str(space_step))

x = np.linspace(space_step/2.0,(x_size-2)*space_step+space_step/2.0,x_size-1)
y = np.linspace(space_step/2.0,(y_size-2)*space_step+space_step/2.0,y_size-1)
print(x)
print(y)
Z_ne = np.ones((y_size-1,x_size-1))
Z_ion = np.ones((y_size-1,x_size-1))

X_g, Y_g = np.meshgrid(x, y)
X_g = X_g / 1e-3
Y_g = Y_g / 1e-3

dir_i = "/users/asfierro/carc-scratch/air_streamer_5um/"
for i in range(time_slices):
	c_time_slice = i*file_stride + file_begin;
	proc_y_size = int(y_size / num_procs)
	y_size_add = y_size % num_procs

	cur_y_pos = 0
	for n in range(num_procs):
		if n == (num_procs-1):
			y_size = y_size + y_size_add
		infile = open(dir_i + "species_data/species_out_" + str(c_time_slice) + "_" + str(n) + ".csv","r");
		#infile = open("data/data_out_" + str(n) + "_" + str(c_time_slice) + ".dat","r");
		print("opening file: " + str(n) + " with time plane = " + str(c_time_slice))
		cur_x_pos = 0
		infile.readline()
		for line in infile:
			line = line.strip().split(',')

			Z_ne[cur_y_pos][cur_x_pos] = float(line[3])
			Z_ion[cur_y_pos][cur_x_pos] = float(line[4])


			if(Z_ne[cur_y_pos][cur_x_pos] > 0.0):
				Z_ne[cur_y_pos][cur_x_pos] = np.log10(Z_ne[cur_y_pos][cur_x_pos])
			if(Z_ion[cur_y_pos][cur_x_pos] > 0.0):
				Z_ion[cur_y_pos][cur_x_pos] = np.log10(Z_ion[cur_y_pos][cur_x_pos])

		
			cur_x_pos = cur_x_pos + 1
			if cur_x_pos == (x_size-1):
				cur_x_pos = 0
				cur_y_pos = cur_y_pos + 1

		infile.close();
	cur_y_pos = 0
	'''
	for n in range(num_procs):
		if n == (num_procs-1):
			y_size = y_size + y_size_add
		infile = open(dir_i + "source_term_data/direct_source_term_out_" + str(n) + "_" + str(c_time_slice) + ".dat","r");
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
	'''

	num_levels = np.linspace(8,18,20)
	plt.figure(figsize=(12,12.0*(y_size/x_size)), dpi=500)
	plt.subplot(2,1,1)
	plt.contourf(X_g,Y_g,Z_ne,levels=num_levels,extend='both')
	plt.contourf(X_g,Y_g*-1,Z_ne,levels=num_levels,extend='both')
	plt.colorbar()

	plt.subplot(2,1,2)
	plt.contourf(X_g,Y_g,Z_ion,levels=num_levels,extend='both')
	plt.contourf(X_g,Y_g*-1,Z_ion,levels=num_levels,extend='both')
	plt.colorbar()


	#infile = open("data/data_out_" + str(n) + ".dat","r");
	plt.savefig("../images/species_out_" + str(c_time_slice) + ".png", dpi=200)
	plt.clf()
	plt.close()

