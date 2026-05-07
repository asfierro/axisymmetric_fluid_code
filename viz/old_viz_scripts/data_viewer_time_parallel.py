import math
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import ticker
from threading import Thread
import os

num_procs = 256
x_size = 0
y_size = 0
space_step = 0
skip = 1
time_slices = 1000
file_begin = 50000
file_stride = 5000

def process_data(thread_id, num_threads):
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
	print("thread_id = " + str(thread_id))
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

	fb = file_begin + thread_id*file_stride
	for i in range(time_slices):
		c_time_slice = fb + i*num_threads*file_stride;
		proc_y_size = int(y_size / num_procs)
		y_size_add = y_size % num_procs
		print("proc " +  str(thread_id) + " processing stride = " + str(c_time_slice))
		cur_y_pos = 0
		for n in range(num_procs):
			if n == (num_procs-1):
				y_size = y_size + y_size_add
			infile = open("/users/asfierro/wheeler-scratch/2d_fluid_axisymmetric/data/data_out_" + str(n) + "_" + str(c_time_slice) + ".dat","r");
			#infile = open("data/data_out_" + str(n) + "_" + str(c_time_slice) + ".dat","r");
			print("proc " + str(thread_id) + " opening file: " + str(n) + " with time plane = " + str(c_time_slice))
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
			infile = open("/users/asfierro/wheeler-scratch/2d_fluid_axisymmetric/source_term_data/source_term_out_" + str(n) + "_" + str(c_time_slice) + ".dat","r");
			#infile = open("data/source_term_out_" + str(n) + "_" + str(c_time_slice) + ".dat","r");
			print("proc " + str(thread_id) + " opening file: " + str(n) + " with time plane = " + str(c_time_slice))
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
		X_g = X_g / 1e-3
		num_levels = np.linspace(0,20000,15)
		plt.figure(figsize=(12,3))
		plt.subplot(2,3,2)
		plt.contourf(X_g,Y_g,Z_potential,levels=num_levels, extend='both')
		plt.contourf(X_g,Y_g*-1,Z_potential,levels=num_levels, extend='both')
		plt.colorbar()

		num_levels = np.linspace(-2e6,2e6,50)
		plt.subplot(2,3,2)
		plt.contourf(X_g,Y_g,Z_field_x, extend='both')
		plt.contourf(X_g,Y_g*-1,Z_field_x, extend='both')
		plt.colorbar()

		num_levels = np.linspace(-1e6,1e6,50)
		plt.subplot(2,3,3)
		plt.contourf(X_g,Y_g,Z_field_y, levels=num_levels, extend='both')
		plt.contourf(X_g,Y_g*-1,Z_field_y, levels=num_levels, extend='both')
		plt.colorbar()
		#num_levels = np.linspace(3,8,25)
		#plt.subplot(2,3,4)
		#plt.contourf(X_g,Y_g,Z_field_mag, levels=num_levels, extend='both')
		#plt.colorbar()

		#num_levels = np.linspace(1.0e8, 1.0e19, 50)
		#num_levels = [1e8, 1e9, 1e10, 1e11, 1e12, 1e13, 1e14, 1e15, 1e16, 1e17, 1e18]
		num_levels = np.linspace(13,19,20)
		plt.subplot(2,3,4)
		plt.contourf(X_g,Y_g/1e-3,Z_electron_density, levels=num_levels, extend='both', cmap="plasma")
		plt.contourf(X_g,Y_g*-1/1e-3,Z_electron_density, levels=num_levels, extend='both', cmap="plasma")
		plt.colorbar()

		num_levels = np.linspace(13,19,20)
		plt.subplot(2,3,5)
		plt.contourf(X_g,Y_g/1e-3,Z_ion_density, levels=num_levels, extend='both', cmap="plasma")
		plt.contourf(X_g,Y_g*-1/1e-3,Z_ion_density, levels=num_levels, extend='both', cmap="plasma")
		plt.colorbar()

		plt.subplot(2,3,6)
		num_levels = np.logspace(8,13,8)
		plt.contourf(X_g, Y_g/1e-3, volumetric_term, extend='both')
		plt.contourf(X_g, Y_g*-1/1e-3, volumetric_term, extend='both')
		plt.colorbar()

		#infile = open("data/data_out_" + str(n) + ".dat","r");
		plt.savefig("../images/img_out_" + str(c_time_slice) + ".png", dpi=200)
		plt.clf()
		plt.close()

if __name__ == "__main__":
	#num_threads = int(os.environ.get("SLURM_NTASKS"))	
	num_threads = 2
	print("num_threads = " + str(num_threads))
	if num_threads > 1:
		for t in range(num_threads):
			process_thread = Thread(target=process_data,args=(t,num_threads,))
			process_thread.start()
			
	else:
		print("num threads not detected or equal to 0")
'''

'''
