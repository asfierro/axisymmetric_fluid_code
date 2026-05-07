import struct
import math

num_procs = 128
file_select = 120000
skip = 4
x_size = 2500
y_size = 1024

outfile = open("/users/asfierro/carc-scratch/air_streamer_5um/mesh_data_combined_" + str(file_select) + "_.bin","wb");
for n in range(num_procs):
        infile = open("/users/asfierro/carc-scratch/air_streamer_5um/mesh_data/data_out_" + str(file_select) + "_" + str(n) + ".csv","r");
        print("opening file: " + str(n) + " with time plane = " + str(file_select))

        ct = 0
        cur_x = 0
        cur_y = 0

        infile.readline()
        for line in infile:
                if (cur_x % skip) == 0 and (cur_y % skip) == 0:
                        line = line.strip().split(',')
                        x = float(line[0])
                        y = float(line[1])
                        v = float(line[2])
                        ex = float(line[3])
                        ey = float(line[4])
                        emag = math.sqrt(ex**2 + ey**2)

                        outfile.write(struct.pack('d',x))
                        outfile.write(struct.pack('d',y))
                        outfile.write(struct.pack('d',emag))
                        #outfile.write(struct.pack('d',v))
                        #outfile.write(struct.pack('d',ex))
                        #outfile.write(struct.pack('d',ey))
                        #outfile.write(struct.pack('d',id))
                
                cur_x = cur_x + 1
                if cur_x == x_size:
                        cur_x = 0
                        cur_y = cur_y + 1

        infile.close()
outfile.close() 

