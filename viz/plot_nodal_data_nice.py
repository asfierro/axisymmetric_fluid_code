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
time_slices = 1
file_begin = 150000
file_stride = 10000
x_size = 2500
y_size = 1024
space_step = 5e-6
dir_i = "/users/asfierro/carc-scratch/air_streamer_5um/"

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

x_species = np.linspace(space_step/2.0,(x_size-2)*space_step+space_step/2.0,x_size-1)
y_species = np.linspace(space_step/2.0,(y_size-2)*space_step+space_step/2.0,y_size-1)
Z_ne = np.ones((y_size-1,x_size-1))
Z_ion = np.ones((y_size-1,x_size-1))

X_g, Y_g = np.meshgrid(x, y)
X_g = X_g / 1e-3
Y_g = Y_g / 1e-3

X_g_s, Y_g_s = np.meshgrid(x_species, y_species)
X_g_s = X_g_s / 1e-3
Y_g_s = Y_g_s / 1e-3

for i in range(time_slices):
        c_time_slice = i*file_stride + file_begin;
        proc_y_size = int(y_size / num_procs)
        y_size_add = y_size % num_procs

        cur_y_pos = 0
        for n in range(num_procs):
                if n == (num_procs-1):
                        y_size = y_size + y_size_add
                infile = open(dir_i + "mesh_data/data_out_" + str(c_time_slice) + "_" + str(n) + ".csv","r");
                #infile = open("data/data_out_" + str(n) + "_" + str(c_time_slice) + ".dat","r");
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

                        #if Z_potential[cur_y_pos][cur_x_pos] < 0:
                                #print(str(cur_x_pos) + "\t" + str(cur_y_pos) + "\t" + str(Z_potential[cur_y_pos][cur_x_pos]))
                
                        cur_x_pos = cur_x_pos + 1
                        if cur_x_pos == x_size:
                                cur_x_pos = 0
                                cur_y_pos = cur_y_pos + 1

                infile.close();
        cur_y_pos = 0
        for n in range(num_procs):
                         
                if n == (num_procs-1):
                        y_size = y_size + y_size_add
                infile = open(dir_i + "source_term_data/fft_source_term_out_" + str(n) + "_" + str(c_time_slice) + ".dat","r");
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

        cur_y_pos = 0
        for n in range(num_procs):
                if n == (num_procs-1):
                        y_size = y_size + y_size_add
                infile = open(dir_i + "species_data/species_out_" + str(c_time_slice) + "_" + str(n) + ".csv","r")
                print("opening file: " + str(n) + " with time plane = " + str(c_time_slice))
                cur_x_pos = 0
                infile.readline()
                for line in infile:
                        line = line.strip().split(',')
                        Z_ne[cur_y_pos][cur_x_pos] = float(line[3])
                        Z_ion[cur_y_pos][cur_x_pos] = float(line[4])


                        if(Z_ne[cur_y_pos][cur_x_pos] > 0.0):
                                Z_ne[cur_y_pos][cur_x_pos] = np.log10(Z_ne[cur_y_pos][cur_x_pos])
                                Z_ion[cur_y_pos][cur_x_pos] = np.log10(Z_ion[cur_y_pos][cur_x_pos])
                        
                        cur_x_pos = cur_x_pos + 1
                        if(cur_x_pos == (x_size-1)):
                            cur_x_pos = 0
                            cur_y_pos = cur_y_pos + 1
                infile.close()
        cur_y_pos = 0


        plt.figure(figsize=(9.375,1.50))    
        
        num_levels = np.linspace(12,21,25)
        plt.subplot(1,3,1)
        plt.contourf(X_g_s,Y_g_s,Z_ne,levels=num_levels,extend='both', cmap="plasma")
        plt.contourf(X_g_s,Y_g_s*-1,Z_ne,levels=num_levels,extend='both', cmap="plasma")
        plt.ylim(-2,2)
        plt.colorbar()
        plt.xticks([])
        plt.yticks([])

        num_levels = np.linspace(0,1.3e7,25)
        plt.subplot(1,3,2)
        plt.contourf(X_g,Y_g,Z_field_mag,levels=num_levels,extend='both', cmap="plasma")
        plt.contourf(X_g,Y_g*-1,Z_field_mag,levels=num_levels,extend='both', cmap="plasma")
        plt.ylim(-2,2)
        plt.colorbar()
        plt.xticks([])
        plt.yticks([])

        plt.subplot(1,3,3)
        num_levels = np.linspace(22,27,25)
        plt.contourf(X_g, Y_g, volumetric_term, levels=num_levels,extend='both', cmap="plasma")
        plt.contourf(X_g, Y_g*-1, volumetric_term, levels=num_levels,extend='both', cmap="plasma")
        plt.ylim(-2,2)
        plt.colorbar()
        plt.xticks([])
        plt.yticks([])


        plt.savefig("../images/mesh_data_out_" + str(c_time_slice) + ".png", dpi=200)
        plt.clf()
        plt.close()

