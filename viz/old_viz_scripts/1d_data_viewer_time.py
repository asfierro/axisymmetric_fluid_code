import math
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import ticker

num_procs = 32
x_size = 0
y_size = 0
space_step = 0
skip = 1
time_slices = 1
file_begin = 0
file_stride = 500

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
Z_potential = np.ones((y_size,x_size))
Z_field_x = np.ones((y_size,x_size))
Z_field_y = np.ones((y_size,x_size))
Z_field_mag = np.ones((y_size,x_size))
Z_electron_density = np.ones((y_size,x_size))
Z_ion_density = np.ones((y_size,x_size))
Z_total_density = np.ones((y_size,x_size))
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
		infile = open("/users/asfierro/wheeler-scratch/2d_fluid_axisymmetric/data_process/data_initial/data_out_" + str(n) + "_" + str(c_time_slice) + ".dat","r");
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

			#Z_field_mag[cur_y_pos][cur_x_pos] = np.log10(math.sqrt(Z_field_x[cur_y_pos][cur_x_pos]**2 + Z_field_y[cur_y_pos][cur_x_pos]**2))
			Z_field_mag[cur_y_pos][cur_x_pos] = math.sqrt(Z_field_x[cur_y_pos][cur_x_pos]**2 + Z_field_y[cur_y_pos][cur_x_pos]**2)
			
			Z_total_density
			if(Z_electron_density[cur_y_pos][cur_x_pos] > 0.0):
				Z_electron_density[cur_y_pos][cur_x_pos] = np.log10(Z_electron_density[cur_y_pos][cur_x_pos])
			if(Z_ion_density[cur_y_pos][cur_x_pos] > 0.0):
				Z_ion_density[cur_y_pos][cur_x_pos] = np.log10(Z_ion_density[cur_y_pos][cur_x_pos])

			#if Z_potential[cur_y_pos][cur_x_pos] < 0:
				#print(str(cur_x_pos) + "\t" + str(cur_y_pos) + "\t" + str(Z_potential[cur_y_pos][cur_x_pos]))
		
			cur_x_pos = cur_x_pos + 1
			if cur_x_pos == x_size:
				cur_x_pos = 0
				cur_y_pos = cur_y_pos + 1

		infile.close();
	'''
	cur_y_pos = 0
	for n in range(num_procs):
		if n == (num_procs-1):
			y_size = y_size + y_size_add
		infile = open("/users/asfierro/wheeler-scratch/2d_fluid_axisymmetric/source_term_data/direct_source_term_out_" + str(n) + "_" + str(c_time_slice) + ".dat","r");
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
	#t_mag = Z_field_mag[0,:]
	#for i in range(len(t_mag)):
		#print(str(x[i]) + "\t" + str(t_mag[i])) 

	file_out = open("/users/asfierro/wheeler-scratch/2d_fluid_axisymmetric/data_process/data_initial/1d_data_out_" + str(c_time_slice) + ".dat", "w")
	for z in range(x_size):
		file_out.write(str(x[i]) + "\t" + str(Z_potential[0,z]) + "\n")
	file_out.close()
	'''

	num_levels = np.linspace(0,1500,15)
	plt.figure(figsize=(12,7.5))
	plt.subplot(2,4,1)
	plt.plot(x,Z_potential[0,:])
	plt.ylim(-6000,20000)

	plt.subplot(2,4,2)
	plt.plot(x,Z_field_x[0,:])
	plt.ylim(0,5e7)

	plt.subplot(2,4,3)
	plt.plot(x,Z_field_y[0,:])
	plt.ylim(-1e6,5e6)

	plt.subplot(2,4,4)
	plt.plot(x,Z_field_mag[0,:])
	plt.ylim(0,1.8e7)

	plt.subplot(2,4,5)
	plt.plot(x,Z_electron_density[0,:])
	plt.ylim(11,22)
	#plt.contourf(X_g,Y_g,Z_electron_density, levels=num_levels, extend='both', cmap="plasma")

	plt.subplot(2,4,6)
	plt.plot(x,Z_ion_density[0,:])
	plt.ylim(11,22)
	#plt.contourf(X_g,Y_g,Z_ion_density, levels=num_levels, extend='both', cmap="plasma")

	plt.subplot(2,4,7)
	plt.plot(x,volumetric_term[0,:])
	plt.ylim(20,32)
	#plt.contourf(X_g, Y_g, volumetric_term, extend='both')
	#infile = open("data/data_out_" + str(n) + ".dat","r");
	plt.savefig("../images_1d/img_out_" + str(c_time_slice) + ".png", dpi=200)
	plt.clf()
	plt.close()
	'''
