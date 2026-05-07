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
file_begin = 5000
file_stride = 500

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
		infile = open("/users/asfierro/wheeler-scratch/2d_fluid_out/data2/data_out_" + str(n) + "_" + str(c_time_slice) + ".dat","r");
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
			
			Z_total_density[cur_y_pos][cur_x_pos] = Z_ion_density[cur_y_pos][cur_x_pos] - Z_electron_density[cur_y_pos][cur_x_pos]
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
	cur_y_pos = 0


	num_levels = np.linspace(0,1500,15)
	plt.figure(figsize=(12,7.5))
	plt.subplot(2,3,1)
	plt.plot(x,Z_potential[int(y_size/2),:])
	plt.ylim(-6000,10000)

	plt.subplot(2,3,2)
	plt.plot(x,Z_field_x[int(y_size/2),:])
	plt.ylim(-3e7,3e7)

	plt.subplot(2,3,3)
	plt.plot(x,Z_field_y[int(y_size/2),:])
	plt.ylim(-10000,10000)

	plt.subplot(2,3,4)
	plt.plot(x,Z_electron_density[int(y_size/2),:])
	plt.ylim(10,20)
	#plt.contourf(X_g,Y_g,Z_electron_density, levels=num_levels, extend='both', cmap="plasma")

	plt.subplot(2,3,5)
	plt.plot(x,Z_ion_density[int(y_size/2),:])
	plt.ylim(10,20)
	#plt.contourf(X_g,Y_g,Z_ion_density, levels=num_levels, extend='both', cmap="plasma")

	plt.subplot(2,3,6)
	plt.plot(x,Z_total_density[int(y_size/2),:])
	#plt.ylim(10,20)
	#plt.contourf(X_g, Y_g, volumetric_term, extend='both')
	#infile = open("data/data_out_" + str(n) + ".dat","r");
	plt.savefig("images_1d/img_out_" + str(c_time_slice) + ".png", dpi=200)
	plt.clf()
	plt.close()
