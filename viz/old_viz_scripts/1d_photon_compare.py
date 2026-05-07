import math
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import ticker

num_procs = 1
x_size = 0
y_size = 0
space_step = 0
skip = 1
time_slices = 10000
file_begin = 1010 
file_stride = 1000
constant_file = open("../constants.h","r")
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
Z_bourdon = np.ones((y_size,x_size))
Z_direct = np.ones((y_size,x_size))
Z_electron_density = np.ones((y_size,x_size))
Z_field_mag = np.ones((y_size,x_size))
Z_intensity = np.ones((y_size,x_size))

X_g, Y_g = np.meshgrid(x, y)
X_g = X_g / 1e-3
Y_g = Y_g / 1e-3

dir_i = "/users/asfierro/wheeler-scratch/2d_fluid_axisymmetric/data/"
for i in range(time_slices):
	c_time_slice = i*file_stride + file_begin;
	proc_y_size = int(y_size / num_procs)
	y_size_add = y_size % num_procs

	cur_y_pos = 0
	for n in range(num_procs):
		if n == (num_procs-1):
			y_size = y_size + y_size_add
		infile = open(dir_i + "mesh_data/data_out_" + str(n) + "_" + str(c_time_slice) + ".dat","r");
		#infile = open("data/data_out_" + str(n) + "_" + str(c_time_slice) + ".dat","r");
		print("opening file: " + str(n) + " with time plane = " + str(c_time_slice))
		cur_x_pos = 0
		infile.readline()
		for line in infile:
			line = line.strip().split()
			Z_electron_density[cur_y_pos][cur_x_pos] = float(line[5])

			if(Z_electron_density[cur_y_pos][cur_x_pos] > 0.0):
				Z_electron_density[cur_y_pos][cur_x_pos] = np.log10(Z_electron_density[cur_y_pos][cur_x_pos])

			x_field = float(line[3])
			y_field = float(line[4])
			Z_field_mag[cur_y_pos][cur_x_pos] = math.sqrt(x_field**2 + y_field**2)
			#if Z_potential[cur_y_pos][cur_x_pos] < 0:
				#print(str(cur_x_pos) + "\t" + str(cur_y_pos) + "\t" + str(Z_potential[cur_y_pos][cur_x_pos]))
		
			cur_x_pos = cur_x_pos + 1
			if cur_x_pos == x_size:
				cur_x_pos = 0
				cur_y_pos = cur_y_pos + 1
		infile.close()

	cur_y_pos = 0
	for n in range(num_procs):
		if n == (num_procs-1):
			y_size = y_size + y_size_add
		infile = open(dir_i + "source_term_data/direct_source_term_out_" + str(n) + "_" + str(c_time_slice) + ".dat","r");
		#infile = open("data/source_term_out_" + str(n) + "_" + str(c_time_slice) + ".dat","r");
		print("opening file: " + str(n) + " with time plane = " + str(c_time_slice))
		cur_x_pos = 0
		for line in infile:
			line = line.strip().split()
			Z_direct[cur_y_pos][cur_x_pos] = float(line[2])
			if(Z_direct[cur_y_pos][cur_x_pos] > 0.0):
				Z_direct[cur_y_pos][cur_x_pos] = np.log10(Z_direct[cur_y_pos][cur_x_pos])
			cur_x_pos = cur_x_pos + 1
			if cur_x_pos == x_size:
				cur_x_pos = 0
				cur_y_pos = cur_y_pos + 1

		infile.close()

	cur_y_pos = 0
	for n in range(num_procs):
		if n == (num_procs-1):
			y_size = y_size + y_size_add
		infile = open(dir_i + "source_term_data/bourdon_source_term_out_" + str(n) + "_" + str(c_time_slice) + ".dat","r");
		#infile = open("data/source_term_out_" + str(n) + "_" + str(c_time_slice) + ".dat","r");
		print("opening file: " + str(n) + " with time plane = " + str(c_time_slice))
		cur_x_pos = 0
		for line in infile:
			line = line.strip().split()
			Z_bourdon[cur_y_pos][cur_x_pos] = float(line[2])
			if(Z_bourdon[cur_y_pos][cur_x_pos] > 0.0):
				Z_bourdon[cur_y_pos][cur_x_pos] = np.log10(Z_bourdon[cur_y_pos][cur_x_pos])
			cur_x_pos = cur_x_pos + 1
			if cur_x_pos == x_size:
				cur_x_pos = 0
				cur_y_pos = cur_y_pos + 1
		infile.close()

	cur_y_pos = 0
	for n in range(num_procs):
		if n == (num_procs-1):
			y_size = y_size + y_size_add
		infile = open(dir_i + "intensity_data/data_out_" + str(n) + "_" + str(c_time_slice) + ".dat","r");
		#infile = open("data/source_term_out_" + str(n) + "_" + str(c_time_slice) + ".dat","r");
		print("opening file: " + str(n) + " with time plane = " + str(c_time_slice))
		cur_x_pos = 0
		for line in infile:
			line = line.strip().split()
			Z_intensity[cur_y_pos][cur_x_pos] = float(line[2])
			cur_x_pos = cur_x_pos + 1
			if cur_x_pos == x_size:
				cur_x_pos = 0
				cur_y_pos = cur_y_pos + 1
		infile.close()

	plt.figure(figsize=(5,10))
	plt.subplot(4,1,1)
	plt.plot(x,Z_electron_density[0,:])
	#plt.xlim(1e-3,4e-3)
	plt.ylim(12,22)

	plt.subplot(4,1,2)
	plt.plot(x,Z_direct[0,:])
	plt.plot(x,Z_bourdon[0,:])
	#plt.xlim(1e-3,4e-3)
	plt.ylim(15,30)

	#plt.subplot(2,3,3)
	#plt.plot(x,Z_bourdon[0,:])

	plt.subplot(4,1,3)
	plt.plot(x,Z_field_mag[0,:])
	#plt.xlim(1e-3,4e-3)

	plt.subplot(4,1,4)
	plt.plot(x,Z_intensity[0,:])
	#plt.xlim(1e-3,4e-3)

	plt.savefig("../images_1d/1d_axial_photon_out_" + str(c_time_slice) + ".png", dpi=200)
	plt.clf()
	plt.close()

