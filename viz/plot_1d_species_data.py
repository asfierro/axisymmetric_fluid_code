import math
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import ticker

num_procs = 1
x_size = 2500
y_size = 1024
space_step = 5e-6
skip = 1
time_slices = 16
file_begin = 0
file_stride = 10000
'''
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
	plt.plot(x,Z_ne[0,:], label=str(c_time_slice))
	#plt.plot(x,Z_ion[0,:], label=str(c_time_slice))
print(Z_ne[0][int(x_size/2)])
plt.legend()
plt.ylim(10,23)
plt.savefig("../images/1d_species_out.png")

#file_out = open(directory + "1d_data_out_" + str(c_time_slice) + ".dat", "w")
#for z in range(x_size):
#file_out.write(str(x[z]) + "\t" + str(Z_potential[0,z]) + "\t" + str(Z_field_x[0,z]) + "\t" + \
	#str(Z_field_y[0,z]) + "\t" + str(Z_field_mag[0,z]) + "\t" + str(Z_electron_density[0,z]) + "\t" + str(Z_ion_density[0,z]) + "\n")
#file_out.close()

