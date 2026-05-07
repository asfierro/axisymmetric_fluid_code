import math

infile = open("ionization_n2.dat","r")

mean_energy = []
x = []
y = [] 
bolsig_alpha = []

for line in infile:
	line = line.strip().split()
	x.append(float(line[0]))
	y.append(float(line[1]))

infile.close()
infile = open("mean_energy.dat","r")
for line in infile:
	line = line.strip().split()
	mean_energy.append(float(line[1]))
infile.close()

infile = open("bolsig_ionization_coefficient.dat","r")
for line in infile:
	line = line.strip().split()
	bolsig_alpha.append(float(line[1])*2.5e25)
infile.close()

for i in range(len(x)):
	vel = math.sqrt(2.0 * mean_energy[i] * 1.602e-19 / 9.11e-31)
	alpha = y[i] / vel * 2.5e25

	
	field = x[i] * 2.5e25 * 1e-21

	analytic = (1.1944e6 + 4.3666e26 / math.pow(field,3)) * math.exp(-2.73e7 / field)
	field = field / 1e5
	print("at E/n = " + str(x[i]) + ", field = " + str(field) + " kV/cm, velocity = " + str(vel) + ", alpha = " + str(alpha) 
			+ ", analytic = " + str(analytic) + ", bolsig = " + str(bolsig_alpha[i]))
