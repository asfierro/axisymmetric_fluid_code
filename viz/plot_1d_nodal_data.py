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

x = np.linspace(0,(x_size-1)*space_step,x_size)
y = np.linspace(0,(y_size-1)*space_step,y_size)
Z_potential = np.ones((y_size,x_size))
Z_field_x = np.ones((y_size,x_size))
Z_field_y = np.ones((y_size,x_size))
Z_field_mag = np.ones((y_size,x_size))
Z_cd = np.ones((y_size,x_size))
volumetric_term = np.ones((y_size,x_size))

X_g, Y_g = np.meshgrid(x, y)
X_g = X_g / 1e-3

dir_i = "/users/asfierro/carc-scratch/air_streamer_5um/"

plt.subplot(211)
for i in range(time_slices):
        c_time_slice = i*file_stride + file_begin;
        proc_y_size = int(y_size / num_procs)
        y_size_add = y_size % num_procs

        cur_y_pos = 0
        for n in range(num_procs):
                if n == (num_procs-1):
                        y_size = y_size + y_size_add
                infile = open(dir_i + "mesh_data/data_out_" + str(c_time_slice) + "_" + str(n) + ".csv","r");
                print("opening file: " + str(n) + " with time plane = " + str(c_time_slice))
                cur_x_pos = 0
                infile.readline()
                for line in infile:
                        line = line.strip().split(',')

                        Z_potential[cur_y_pos][cur_x_pos] = float(line[3])
                        Z_field_x[cur_y_pos][cur_x_pos] = float(line[4])
                        Z_field_y[cur_y_pos][cur_x_pos] = float(line[5])
                        Z_cd[cur_y_pos][cur_x_pos] = float(line[6])

                        Z_field_mag[cur_y_pos][cur_x_pos] = (math.sqrt(Z_field_x[cur_y_pos][cur_x_pos]**2 + Z_field_y[cur_y_pos][cur_x_pos]**2))
                        cur_x_pos = cur_x_pos + 1
                        if cur_x_pos == (x_size-1):
                                cur_x_pos = 0
                                cur_y_pos = cur_y_pos + 1

                infile.close();
        plt.plot(x,Z_field_mag[0,:], label=str(c_time_slice))
        print(Z_field_mag[0,:])
#plt.legend()
plt.subplot(212)

for i in range(time_slices):
        c_time_slice = i*file_stride + file_begin;
        proc_y_size = int(y_size / num_procs)
        y_size_add = y_size % num_procs

        cur_y_pos = 0
        for n in range(num_procs):
                if n == (num_procs-1):
                        y_size = y_size + y_size_add
                infile = open(dir_i + "source_term_data/fft_source_term_out_" + str(n) + "_" + str(c_time_slice) + ".dat","r");
                print("opening file: " + str(n) + " with time plane = " + str(c_time_slice))
                cur_x_pos = 0
                infile.readline()
                for line in infile:
                        line = line.strip().split()
                        volumetric_term[cur_y_pos][cur_x_pos] = float(line[2])
                        volumetric_term[cur_y_pos][cur_x_pos] = np.log10(volumetric_term[cur_y_pos][cur_x_pos])
                        cur_x_pos = cur_x_pos + 1
                        if cur_x_pos == x_size:
                                cur_x_pos = 0
                                cur_y_pos = cur_y_pos + 1
                infile.close();
        plt.plot(x,volumetric_term[0,:], label=str(c_time_slice))
'''

for i in range(time_slices):
        c_time_slice = i*file_stride + file_begin;
        proc_y_size = int(y_size / num_procs)
        y_size_add = y_size % num_procs

        cur_y_pos = 0
        for n in range(num_procs):
                if n == (num_procs-1):
                        y_size = y_size + y_size_add
                infile = open(dir_i + "source_term_data/bourdon_source_term_out_" + str(n) + "_" + str(c_time_slice) + ".dat","r");
                print("opening file: " + str(n) + " with time plane = " + str(c_time_slice))
                cur_x_pos = 0
                infile.readline()
                for line in infile:
                        line = line.strip().split()
                        volumetric_term[cur_y_pos][cur_x_pos] = float(line[2])
                        volumetric_term[cur_y_pos][cur_x_pos] = np.log10(volumetric_term[cur_y_pos][cur_x_pos])
                        cur_x_pos = cur_x_pos + 1
                        if cur_x_pos == x_size:
                                cur_x_pos = 0
                                cur_y_pos = cur_y_pos + 1
                infile.close();
        plt.plot(x,volumetric_term[0,:], label=str(c_time_slice))
'''
print(Z_field_mag[0][int(x_size/2)])
print(volumetric_term[0][int(x_size/2)])
plt.ylim(18,29)
plt.legend()
plt.savefig("../images/1d_field_out.png")
























