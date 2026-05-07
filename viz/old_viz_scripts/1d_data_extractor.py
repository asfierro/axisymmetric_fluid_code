import math
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import ticker

num_procs = 1
x_size = 4000
y_size = 3200
space_step = 3.125e-6
skip = 1
time_slices = 10000
file_begin = 8000
file_stride = 12000

directory = "/users/asfierro/wheeler-scratch/2d_fluid_axisymmetric/data_process/data_set_k/"

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
                infile = open(directory + "data_out_" + str(n) + "_" + str(c_time_slice) + ".dat","r");
                #infile = open("data/data_out_" + str(n) + "_" + str(c_time_slice) + ".dat","r");
                print("opening file: " + str(n) + " with time plane = " + str(c_time_slice))
                cur_x_pos = 0
                infile.readline()
                for line in infile:
                        line = line.strip().split()

                        Z_potential[cur_y_pos][cur_x_pos] = float(line[2])
                        Z_field_x[cur_y_pos][cur_x_pos] = float(line[3])
                        Z_field_y[cur_y_pos][cur_x_pos] = float(line[4])
                        Z_electron_density[cur_y_pos][cur_x_pos] = float(line[5])
                        Z_ion_density[cur_y_pos][cur_x_pos] = float(line[6])
                        Z_field_mag[cur_y_pos][cur_x_pos] = math.sqrt(Z_field_x[cur_y_pos][cur_x_pos]**2 + Z_field_y[cur_y_pos][cur_x_pos]**2)
                        
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
        file_out = open(directory + "1d_data_out_" + str(c_time_slice) + ".dat", "w")
        for z in range(x_size):
                file_out.write(str(x[z]) + "\t" + str(Z_potential[0,z]) + "\t" + str(Z_field_x[0,z]) + "\t" + \
                        str(Z_field_y[0,z]) + "\t" + str(Z_field_mag[0,z]) + "\t" + str(Z_electron_density[0,z]) + "\t" + str(Z_ion_density[0,z]) + "\n")
        file_out.close()

