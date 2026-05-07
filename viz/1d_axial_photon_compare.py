import math
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import ticker

num_procs = 1
x_size = 4096
y_size = 4096
space_step = 5e-6
skip = 1
time_slices = 1
file_begin = 0
file_stride = 1000

x = np.linspace(0,x_size*space_step,x_size)
y = np.linspace(0,y_size*space_step,y_size)
Z_bourdon = np.ones((y_size,x_size))
Z_fft = np.ones((y_size,x_size))

X_g, Y_g = np.meshgrid(x, y)
X_g = X_g / 1e-3
Y_g = Y_g / 1e-3

dir_i = "/users/asfierro/carc-scratch/air_streamer_compare/"
for i in range(time_slices):
        c_time_slice = i*file_stride + file_begin;
        proc_y_size = int(y_size / num_procs)
        y_size_add = y_size % num_procs

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
                        Z_fft[cur_y_pos][cur_x_pos] = float(line[2])
                        if(Z_fft[cur_y_pos][cur_x_pos] > 0.0):
                                Z_fft[cur_y_pos][cur_x_pos] = np.log10(Z_fft[cur_y_pos][cur_x_pos])
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

        plt.subplot(1,1,1)
        plt.plot(x,Z_bourdon[0,:], label="bourdon")
        plt.plot(x,Z_fft[0,:], label="fft")
        plt.ylim(15,30)
        plt.legend()

        plt.savefig("../images/1d_zdirection_photon_out_" + str(c_time_slice) + ".png", dpi=200)
        plt.clf()
        plt.close()

