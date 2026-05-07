# for air
# from Kossyi, Plasma Sources Sci. Technol., 1, 207, 1992
import math

T_e = 11600
T_gas = 300

# e- + e- + N2+ -> e- + N2+
# e- + e- + O2+ -> e- + O2+
k_43 = 1e-19 * (300/T_e)**4.5 / 1e12

# e- + O2 + O2+ -> O2 + O2
# e- + N2 + O2+ -> N2 + O2
# e- + N2 + N2+ -> N2 + N2
# e- + O2 + N2+ -> O2 + N2
k_44 = 6e-27 * (300/T_e)**1.5 / 1e12

# e + o2 + o2 -> o2- + o2
k_45 = 1.4e-29 * (300/T_e) * math.exp(-600/T_gas) * math.exp((700 * (T_e - T_gas)) / (T_e * T_gas)) / 1e12

# e + O2 + n2 -> o2- + N2
k_46 = 1.07e-31 * (300/T_e)**2 * math.exp(-70/T_gas) * math.exp((1500 * (T_e - T_gas)) / (T_e * T_gas)) / 1e12


print("k_43 = " + str(k_43))
print("k_44 = " + str(k_44))
print("k_45 = " + str(k_45))
print("k_46 = " + str(k_46))
