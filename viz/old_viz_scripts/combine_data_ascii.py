import struct
import math

num_procs = 32
file_select = 2000
skip = 1
x_size = 4000
y_size = 3200

outfile = open("/users/asfierro/wheeler-scratch/2d_fluid_axisymmetric/data/mesh_data/data_combined_" + str(file_select) + ".dat","w");
for n in range(num_procs):
	infile = open("/users/asfierro/wheeler-scratch/2d_fluid_axisymmetric/data/mesh_data/data_out_" + str(n) + "_" + str(file_select) + ".dat","r");
	print("opening file: " + str(n) + " with time plane = " + str(file_select))

	ct = 0
	cur_x = 0
	cur_y = 0

	infile.readline()
	for line in infile:
		if (cur_x % skip) == 0 and (cur_y % skip) == 0:
			line = line.strip().split()
			x = float(line[0])
			y = float(line[1])
			v = float(line[2])
			ex = float(line[3])
			ey = float(line[4])
			ed = float(line[5])
			id = float(line[6])
			emag = math.sqrt(ex**2 + ey**2)

			outfile.write(str(x) + "\t" + str(y) + "\t" + str(ed) + "\t" + str(emag) + "\n")
		
		cur_x = cur_x + 1
		if cur_x == x_size:
			cur_x = 0
			cur_y = cur_y + 1

	infile.close()
outfile.close()	

