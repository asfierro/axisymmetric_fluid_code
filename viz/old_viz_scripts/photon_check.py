import math
from scipy import interpolate

n_b = 2.5e25
n_o2 = 2.5e25 * 0.22

spectral_intensity = open("spectra/N2_100mTorr.txt", "r")
o2_photoionization = open("absorption_data/o2_photoionization.txt", "r")
total_absorption = open("absorption_data/o2_total_abs.txt", "r")

wavelength = []
intensity = []

wavelength_iz = []
absorption_iz = []

wavelength_abs = []
absorption_abs =[]

for line in spectral_intensity:
	line = line.strip().split()
	wavelength.append(float(line[0])*1e-9)
	intensity.append(float(line[1]))

spectral_intensity.close()

for line in o2_photoionization:
	line = line.strip().split()
	wavelength_iz.append(float(line[0])*1e-10)
	absorption_iz.append(float(line[1])/1e4)

iz_interp = interpolate.interp1d(wavelength_iz, absorption_iz)
o2_photoionization.close()

for line in total_absorption:
	line = line.strip().split()
	wavelength_abs.append(float(line[0])*1e-10)
	absorption_abs.append(float(line[1])/1e4)
abs_interp = interpolate.interp1d(wavelength_abs, absorption_abs)

distance = 1e-3
phi = 0.0

for i in range(len(wavelength)):
	c_w = wavelength[i]
	c_e = 6.626e-34 * 3e8 / c_w
	c_i = intensity[i] * 16693.9067 * 1e-6 / c_e / 1.57e-8
	print(str(c_w) + "\t" + str(c_i) + "\t" + str(iz_interp(c_w)) + "\t" + str(abs_interp(c_w)))
	phi = phi + c_i * (iz_interp(c_w) * n_o2) * math.exp(-abs_interp(c_w) * n_o2 * distance)

phi = phi * (wavelength[1] - wavelength[0])

volume = 3.1415 * 2e-6 * 2e-6 * 2e-6
q_phi = phi * volume  / (4 * 3.1415 * distance * distance) 
print(phi)
print(q_phi)
