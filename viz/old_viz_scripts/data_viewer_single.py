import math
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import numpy as np

num_procs = 32
x_size = 0
y_size = 0
space_step = 0
skip = 1
time_plane = 0

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

x = np.linspace(0,x_size*space_step/1e-3,x_size)
y = np.linspace(0,y_size*space_step/1e-3,y_size)

X_g, Y_g = np.meshgrid(x, y)
Z_potential = np.ones((y_size,x_size))
Z_field_x = np.ones((y_size,x_size))
Z_field_y = np.ones((y_size,x_size))
Z_electron_density = np.ones((y_size,x_size))
Z_ion_density = np.ones((y_size,x_size))

proc_y_size = int(y_size / num_procs)
y_size_add = y_size % num_procs

cur_y_pos = 0
for n in range(num_procs):
	if n == (num_procs-1):
		y_size = y_size + y_size_add
	#Einfile = open("data/data_out_" + str(n) + "_" + str(time_plane) + ".dat","r");
	infile = open("/users/asfierro/wheeler-scratch/2d_fluid_axisymmetric/data/data_out_" + str(n) + "_" + str(time_plane) + ".dat","r");
	print("opening file: " + str(n))
	cur_x_pos = 0
	for line in infile:
		line = line.strip().split()

		Z_potential[cur_y_pos][cur_x_pos] = float(line[2])
		Z_field_x[cur_y_pos][cur_x_pos] = float(line[3])
		Z_field_y[cur_y_pos][cur_x_pos] = float(line[4])
		Z_electron_density[cur_y_pos][cur_x_pos] = float(line[5])
		Z_ion_density[cur_y_pos][cur_x_pos] = float(line[6])

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

plt.figure(figsize=(12,4))

num_levels = np.linspace(0,2500,50)
plt.subplot(2,3,1)
plt.contourf(X_g,Y_g,Z_potential,extend='both')
plt.colorbar()

num_levels = np.linspace(-1e8,1e8,50)
plt.subplot(2,3,2)
plt.contourf(X_g,Y_g,Z_field_x, extend='both')
plt.colorbar()

num_levels = np.linspace(-1e4,1e4,50)
plt.subplot(2,3,3)
plt.contourf(X_g,Y_g,Z_field_y, extend='both')
plt.colorbar()

num_levels = np.linspace(8,20,25)
plt.subplot(2,3,4)
plt.contourf(X_g,Y_g,Z_electron_density, extend='both')
plt.colorbar()

num_levels = np.linspace(8,20,25)
plt.subplot(2,3,5)
plt.contourf(X_g,Y_g,Z_ion_density, extend='both')
plt.colorbar()

plt.savefig("img_out_" + str(time_plane) + ".png", dpi=200)
plt.show()

