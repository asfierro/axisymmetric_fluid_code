import math
import numpy as np

e_n = np.linspace(0.1,5000,1000)
alpha = []
for i in range(len(e_n)):
	e_field = e_n[i] * 1e-21 * 2.414e25

	f = 1.1944e6 + 4.3666e26 / (math.pow(e_field,3))
	s = math.exp(-2.73e7 / e_field)

	alpha.append(f*s - 340.75)

ofile = open("alpha_air_2.dat","w")
for i in range(len(alpha)):
	ofile.write(str(e_n[i]) + "\t" + str(alpha[i]) + "\n")
ofile.close()	

