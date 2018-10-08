import matplotlib.pyplot as plt
import numpy as np
import scipy
import math
import os
from scipy.interpolate import interp1d
from astropy.table import Table
from astropy.io import fits


#####################################################
#                     PROBLEM 2                     # 
#####################################################

#define some constants
c = 3.0e10
mag_0_fluxdensity = 3631 * 1e-23
'''
#pop data into arrays
A0V_spectrum = '/home/pbc463/Desktop/Fall18_Courses/Galaxies/HW1/tools/Spectrum_A0V.txt'
subaru_g = '/home/pbc463/Desktop/Fall18_Courses/Galaxies/HW1/tools/subaru_g'
subaru_i = '/home/pbc463/Desktop/Fall18_Courses/Galaxies/HW1/tools/subaru_i'
subaru_r = '/home/pbc463/Desktop/Fall18_Courses/Galaxies/HW1/tools/subaru_r'
subaru_y = '/home/pbc463/Desktop/Fall18_Courses/Galaxies/HW1/tools/subaru_y'
subaru_z = '/home/pbc463/Desktop/Fall18_Courses/Galaxies/HW1/tools/subaru_z'

A0V_spectrum_data = np.zeros((1221,2))
flat_spectrum_data = np.zeros((1221,2))
subaru_g_data = np.zeros((200,2))
subaru_i_data = np.zeros((200,2))
subaru_r_data = np.zeros((200,2))
subaru_y_data = np.zeros((501,2))
subaru_z_data = np.zeros((200,2))

#grab data from text files
with open(A0V_spectrum) as fp:
	for i, line in enumerate(fp):
		if i != 0:
			split = line.split()
			A0V_spectrum_data[i-1,0] = float(split[0])
			A0V_spectrum_data[i-1,1] = float(split[1])
			
for i in range(len(flat_spectrum_data)):
	flat_spectrum_data[i,0] = float(A0V_spectrum_data[i,0])
	flat_spectrum_data[i,1] = float(mag_0_fluxdensity)
			
with open(subaru_g) as fp:
	for i, line in enumerate(fp):
		if i != 0:
			split = line.split()
			subaru_g_data[i-1,0] = float(split[0])
			subaru_g_data[i-1,1] = float(split[1])
		
with open(subaru_i) as fp:
	for i, line in enumerate(fp):
		if i != 0:
			split = line.split()
			subaru_i_data[i-1,0] = float(split[0])
			subaru_i_data[i-1,1] = float(split[1])

with open(subaru_r) as fp:
	for i, line in enumerate(fp):
		if i != 0:
			split = line.split()
			subaru_r_data[i-1,0] = float(split[0])
			subaru_r_data[i-1,1] = float(split[1])
			
with open(subaru_y) as fp:
	for i, line in enumerate(fp):
		if i != 0:
			split = line.split()
			subaru_y_data[i-1,0] = float(split[0])
			subaru_y_data[i-1,1] = float(split[1])

with open(subaru_z) as fp:
	for i, line in enumerate(fp):
		if i != 0:
			split = line.split()
			subaru_z_data[i-1,0] = float(split[0])
			subaru_z_data[i-1,1] = float(split[1])


#plot spectrum of A0V and filters
fig = plt.figure()
plt.style.use('ggplot')
plt.xlabel(r'$\lambda$' + ' ' + '['r'$\mu$$m$'']')
plt.ylabel('$S_\\nu$' + ' ' + r'[$Jy$]')
plt.plot(A0V_spectrum_data[:800,0] * 1e-4, A0V_spectrum_data[:800,1] *(1e23), color = 'k', label = 'A0V', linewidth = 0.9)
plt.plot(flat_spectrum_data[:750,0] * 1e-4, flat_spectrum_data[:750,1] *(1e23),label = '3631 Jy', linewidth = 0.9)
plt.plot(subaru_g_data[:,0] * 1e-4, subaru_g_data[:,1] * 3631, label = 'g', linewidth = 0.75, color = 'g')
plt.plot(subaru_r_data[:,0] * 1e-4, subaru_r_data[:,1] * 3631, label = 'r', linewidth = 0.75, color = 'r')
plt.plot(subaru_i_data[:,0] * 1e-4, subaru_i_data[:,1] * 3631, label = 'i', linewidth = 0.75, color = 'b')
plt.plot(subaru_z_data[:,0] * 1e-4, subaru_z_data[:,1] * 3631, label = 'z', linewidth = 0.75, color = 'm')
plt.plot(subaru_y_data[:,0] * 1e-4, subaru_y_data[:,1] * 3631, label = 'y', linewidth = 0.75, color = 'y')
plt.legend(frameon=False, loc='upper right', ncol=1, fontsize = 'small')
plt.tight_layout()
#plt.show()

fig.savefig("A0V_spectrum.pdf", bbox_inches='tight')


#interpolate to get transmission values at wavelengths where A0V star spectrum was measured
ginterp = np.interp(A0V_spectrum_data[369:465, 0], subaru_g_data[:,0], subaru_g_data[:,1])
iinterp = np.interp(A0V_spectrum_data[513:640, 0], subaru_i_data[:,0], subaru_i_data[:,1])
rinterp = np.interp(A0V_spectrum_data[448:544, 0], subaru_r_data[:,0], subaru_r_data[:,1])
yinterp = np.interp(A0V_spectrum_data[632:711, 0], subaru_y_data[:,0], subaru_y_data[:,1])
zinterp = np.interp(A0V_spectrum_data[577:696, 0], subaru_z_data[:,0], subaru_z_data[:,1])

#calculate magnitudes in different bands

def FluxDensity_numerator (nu, F_nu, t_nu):
	value = ((1/nu) * F_nu * t_nu)
	return value
	
def FluxDensity_denominator (nu, t_nu):
	value = ((1/nu) * t_nu)
	return value

#Crude integration for flux densities
gFluxDensity_num = 0
rFluxDensity_num = 0
iFluxDensity_num = 0
zFluxDensity_num = 0
yFluxDensity_num = 0

gFluxDensity_den = 0
iFluxDensity_den = 0
rFluxDensity_den = 0
yFluxDensity_den = 0
zFluxDensity_den = 0

for i in range(len(ginterp)):
	gFluxDensity_num += FluxDensity_numerator(A0V_spectrum_data[i+369, 0], A0V_spectrum_data[i+369, 1], ginterp[i])
	gFluxDensity_den += FluxDensity_denominator(A0V_spectrum_data[i+369, 0], ginterp[i])
	
for i in range(len(iinterp)):
	iFluxDensity_num += FluxDensity_numerator(A0V_spectrum_data[i+513, 0], A0V_spectrum_data[i+513, 1], iinterp[i])
	iFluxDensity_den += FluxDensity_denominator(A0V_spectrum_data[i+513, 0], iinterp[i])
	
for i in range(len(rinterp)):
	rFluxDensity_num += FluxDensity_numerator(A0V_spectrum_data[i+448, 0], A0V_spectrum_data[i+448, 1], rinterp[i])
	rFluxDensity_den += FluxDensity_denominator(A0V_spectrum_data[i+448, 0], rinterp[i])	
	
for i in range(len(yinterp)):
	yFluxDensity_num += FluxDensity_numerator(A0V_spectrum_data[i+632, 0], A0V_spectrum_data[i+632, 1], yinterp[i])
	yFluxDensity_den += FluxDensity_denominator(A0V_spectrum_data[i+632, 0], yinterp[i])
	
for i in range(len(zinterp)):
	zFluxDensity_num += FluxDensity_numerator(A0V_spectrum_data[i+577, 0], A0V_spectrum_data[i+577, 1], zinterp[i])
	zFluxDensity_den += FluxDensity_denominator(A0V_spectrum_data[i+577, 0], zinterp[i])

print ('g: ' + str(1e23 * gFluxDensity_num / gFluxDensity_den))
print ('r: ' + str(1e23 * rFluxDensity_num / rFluxDensity_den))
print ('i: ' + str(1e23 * iFluxDensity_num / iFluxDensity_den))
print ('z: ' + str(1e23 * zFluxDensity_num / zFluxDensity_den))
print ('y: ' + str(1e23 * yFluxDensity_num / yFluxDensity_den))

print (mag_0_fluxdensity)

A0_gmag = -2.5 * math.log10((gFluxDensity_num / gFluxDensity_den) / mag_0_fluxdensity)
A0_imag = -2.5 * math.log10((iFluxDensity_num / iFluxDensity_den) / mag_0_fluxdensity)
A0_rmag = -2.5 * math.log10((rFluxDensity_num / rFluxDensity_den) / mag_0_fluxdensity)
A0_ymag = -2.5 * math.log10((yFluxDensity_num / yFluxDensity_den) / mag_0_fluxdensity)
A0_zmag = -2.5 * math.log10((zFluxDensity_num / zFluxDensity_den) / mag_0_fluxdensity)

print('gmag: ' + str(A0_gmag) + ' rmag: ' + str(A0_rmag) + ' imag: ' + str(A0_imag) + ' zmag: ' + str(A0_zmag) + ' ymag: ' + str(A0_ymag))

#calculate and record Luminosities
Luminosities = np.zeros((1221,3))
D_Vega = 7.68 * (3.086e18)
L_solar = 3.839e33


for i in range(len(A0V_spectrum_data)):
	lambd = A0V_spectrum_data[i,0] * 1e-8 #in cm
	F_nu = A0V_spectrum_data[i,1] #erg/s/cm2/Hz
	nu = 3e10 / lambd # in Hz
	F_lambda = (1/(3e10)) * math.pow(nu,2) * F_nu
	Luminosities[i,0] = lambd * 1e4 #convert values to microns for plotting purposes
	Luminosities[i,1] = nu * F_nu * (4*math.pi*math.pow(D_Vega,2))
	Luminosities[i,2] = lambd * F_lambda * (4*math.pi*math.pow(D_Vega,2)) + 1e35


fig = plt.figure()
plt.style.use('ggplot')
plt.xlabel(r'$\lambda\,$' '['r'$\mu$$m$'']')
plt.ylabel(r'$\nu$$L_{\nu}\,[erg\, s^{-1}\, Hz]$' + ' or ' + r'$\lambda$$L_{\lambda}\,[erg\, s^{-1}\, \AA^{-1}]$')
plt.plot(Luminosities[360:750,0], Luminosities[360:750,1], color = 'k', label = r'$\nu$$L_{\nu}$', linewidth = 0.9)
plt.plot(Luminosities[360:750,0], Luminosities[360:750,2], label = r'$\lambda$$L_{\lambda}$' + ' offset by +$1x10^{35}$', linewidth = 0.9)
plt.legend(frameon=False, loc='upper right', ncol=1, fontsize = 'small')
plt.tight_layout()
#plt.show()

fig.savefig("Luminosities.pdf", bbox_inches='tight')
'''

#####################################################
#                     PROBLEM 3                     # 
#####################################################

#Define constants
M_jup = 1.898e30#grams
M_solar = 1.989e33#grams
L_solar = 3.839e33#ergs/s
t_solar = 1e10
G = 6.67e-8#cm^3g^-1s^-2
dm = 1.98748e32
dl = 1e33
d = 7.68 * 3.086e18 #parsecs to cm

def cumulative_mass_fraction(low_mass_limit, high_mass_limit):

	#Define function for Salpeter IMF
	def IMF(m):
		return math.pow(m, -2.35)

	starting_mass = low_mass_limit
	mass_distribution = []
	n_distribution = []
	x = []

	#calculate stellar mass distribution
	#calculate stellar number distribution

	while starting_mass < high_mass_limit:
		x.append(starting_mass / M_solar)
		dn = IMF(starting_mass) * dm
		mass_distribution.append(dn * starting_mass)
		starting_mass += dm
		n_distribution.append(dn)

	total_n = sum(n_distribution)
	total_mass = sum(mass_distribution)
	cumulative_mass_above_m = np.zeros(len(mass_distribution))

	rolling_cumulative_mass = 0
	for i in range(len(mass_distribution)-1, 0, -1):
		rolling_cumulative_mass += mass_distribution[i]
		cumulative_mass_above_m[i] = rolling_cumulative_mass

	cumulative_mass_fraction = cumulative_mass_above_m / total_mass

	#calculate expectation value	
	expectation_m = 0
	for i in range(len(x)):
		expectation_m += x[i] * n_distribution[i] / total_n
	return(x, cumulative_mass_fraction, expectation_m, np.asarray(n_distribution) / total_n)


low_mass_limit = 80 * M_jup
high_mass_limit = 100 * M_solar
int_mass_limit = 2.82 * M_solar
second_int_mass_limit = 2.14 * M_solar
	
hund_solar = cumulative_mass_fraction(low_mass_limit, high_mass_limit)
int_solar = cumulative_mass_fraction(low_mass_limit, int_mass_limit)
second_int_solar = cumulative_mass_fraction(low_mass_limit, second_int_mass_limit)


'''
fig = plt.figure()
plt.style.use('ggplot')
plt.xlabel(r'log M$_{s}\,[M_{\odot}]$')
plt.ylabel(r'f(>M$_s$)')
plt.loglog(hund_solar[0][1:], hund_solar[1][1:])
plt.tight_layout()
#plt.show()

fig.savefig("CSMF.pdf", bbox_inches='tight')
'''

#calculating mass and age relationship

#M > 20 M_solar
t1 = t_solar / 3200

#2 M_solar < M < 20 M_solar
def t2(M):
	t2 = t_solar * ((M / M_solar)**-2.5) * (1/1.5)
	return t2

#0.43 M_solar < M < 2 M_solar
def t3(M):
	t3 = t_solar * ((M / M_solar)**-3)
	return t3

#M < 0.43 M_solar
def t4(M):
	t4 = t_solar * ((M / M_solar)**-1.3) * (1/.23)
	return t4

low_mass = 0.4 * M_solar
high_mass = 100 * M_solar

masses = np.linspace(low_mass, high_mass, 100000)
ages = np.zeros(len(masses))

for i in range(len(masses)):
	if masses[i] < 0.43 * M_solar:
		ages[i] = t4(masses[i])
	if masses[i] > 0.43 * M_solar and masses[i] < 2 * M_solar:
		ages[i] = t3(masses[i])
	if masses[i] > 2 * M_solar and masses[i] < 20 * M_solar:
		ages[i] = t2(masses[i])
	if masses[i] > 20 * M_solar:
		ages[i] = t1
		
ages = ages / 1e6
masses = masses / M_solar

'''
fig = plt.figure()
plt.style.use('ggplot')
plt.xlabel(r'mass $[M_{\odot}]$')
plt.ylabel(r'age [Myr]')
plt.loglog(masses, ages, color = 'k', linewidth = 0.9)
plt.tight_layout()
#plt.show()

fig.savefig("massvage.pdf", bbox_inches='tight')
'''

ages_in_question = [.01, .05, 1]
masses_in_question = []

def mass(t):
	mass = ((1.5 * (t/10)) ** (1/-2.5))
	return mass

for age in ages_in_question:
	masses_in_question.append(mass(age))
	
print(masses_in_question)


def cumulative_lum_fraction(low_lum_limit, high_lum_limit):

	#Define function for Salpeter IMF
	def IMF(l, exp):
		return math.pow(l, exp)

	starting_luminosity = low_lum_limit
	luminosity_distribution = []
	n_distribution = []
	luminosities = []
	x = []

	#calculate stellar lum distribution
	#calculate stellar number distribution
	
	salpeter = -2.35
	exponent = 0
	
	while starting_luminosity < high_lum_limit:
		if starting_luminosity < 3.3e-2 * L_solar:
			exp = (salpeter - 1.3) / 2.3
			x.append(((starting_luminosity / L_solar) / 0.23)**(1/2.3))
		if starting_luminosity > 3.3e-2 * L_solar and starting_luminosity < 16 * L_solar:
			exp = (salpeter - 3) / 4
			x.append((starting_luminosity / L_solar)**(1/4))
		if starting_luminosity > 16 * L_solar and starting_luminosity < 5.37e4 * L_solar:
			x.append(((starting_luminosity / L_solar) / 1.5)**(1/3.5))
		if starting_luminosity > 5.37e4 * L_solar:
			exp = salpeter
			x.append((starting_luminosity / L_solar) / 3200)
		
		dn = IMF(starting_luminosity, exponent) * dl
		luminosity_distribution.append(dn * starting_luminosity)
		starting_luminosity += dl
		n_distribution.append(dn)

	total_n = sum(n_distribution)
	total_luminosity = sum(luminosity_distribution)
	cumulative_luminosity_above_l = np.zeros(len(luminosity_distribution))

	rolling_cumulative_luminosity = 0
	for i in range(len(luminosity_distribution)-1, 0, -1):
		rolling_cumulative_luminosity += luminosity_distribution[i]
		cumulative_luminosity_above_l[i] = rolling_cumulative_luminosity

	cumulative_luminosity_fraction = cumulative_luminosity_above_l / total_luminosity

	return(x, cumulative_luminosity_fraction)

CSLF1 = cumulative_lum_fraction(6.13e-4 * L_solar, 3.2e5 * L_solar)
CSLF2 = cumulative_lum_fraction(6.13e-4 * L_solar, 1.42e3 * L_solar)
CSLF3 = cumulative_lum_fraction(6.13e-4 * L_solar, 21.5 * L_solar)

fig = plt.figure()
plt.style.use('ggplot')
plt.xlabel(r'log M$_{s}\,[M_{\odot}]$')
plt.ylabel(r'f_L(>M$_s$)')
plt.loglog(CSLF1[0][1:], CSLF1[1][1:], label = '0-age population', linewidth = 0.9)
plt.loglog(CSLF2[0][1:], CSLF2[1][1:], label = '500 Myr old population', linewidth = 0.9)
plt.loglog(CSLF3[0][1:], CSLF3[1][1:], label = '1 Gyr old population', linewidth = 0.9)
plt.legend(frameon=False, loc='lower center', ncol=1, fontsize = 'small')
plt.tight_layout()
#plt.show()

fig.savefig("CSLF.pdf", bbox_inches='tight')


#####################################################
#                     PROBLEM 4                     # 
#####################################################

#pop data into arrays
A0V_spectrum = 'tools/EEM_dwarf_UBVIJHK_colors_Teff'

SpT = []
Teff = []
LogL = []
Msun = []

with open(A0V_spectrum) as fp:
	for i, line in enumerate(fp):
		split = line.split()
		if i > 30 and i < 105:
			SpT.append(split[0])
			Teff.append(float(split[1]))
			LogL.append(float(split[5]))
			Msun.append(float(split[18]))

Teff.append(45000)
LogL.append(3.2e5)
Msun.append(100)

x = np.linspace(80*M_jup / M_solar, 100, 1000)
T = interp1d(Msun, Teff)
L = interp1d(Msun, LogL)

interpedTemps = T(x)
interpedLums = L(x)
		

'''
#plot Temp against Mass
fig = plt.figure()
plt.style.use('ggplot')
plt.ylabel(r'T$_eff$')
plt.xlabel(r'M [$M_\odot$]')
plt.plot(x, interpedTemps, color = 'k', linewidth = 0.9)
#plt.legend(frameon=False, loc='upper right', ncol=1, fontsize = 'small')
plt.tight_layout()
#plt.show()

fig.savefig("TeffvsM.pdf", bbox_inches='tight')			

#plot Mass against Luminosity
fig = plt.figure()
plt.style.use('ggplot')
plt.ylabel(r'logL [$L_\odot$]')
plt.xlabel(r'M [$M_\odot$]')
plt.plot(x, interpedLums, color = 'k', linewidth = 0.9)
#plt.legend(frameon=False, loc='upper right', ncol=1, fontsize = 'small')
plt.tight_layout()
#plt.show()

fig.savefig("LvsM.pdf", bbox_inches='tight')
'''

#part b

#read in list of files and grab temps corresponding to the file

file_temps = np.zeros(61)

with open('tools/kurucz93/kp00/files.lis') as fp:
	for i, line in enumerate(fp):
		right = line[:-6]
		left = right[5:]
		file_temps[i] = (left)

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]

#find nearest temps in Kurucz files corresponding to interpolated temperatures and respective masses
nearestTemps = []
for i in range(len(interpedTemps)):
	nearestTemps.append(int(find_nearest(file_temps, interpedTemps[i])))

n_distribution = hund_solar[3]
n_distribution_500 = int_solar[3]
n_distribution_1 = second_int_solar[3]

print(len(n_distribution_500), len(n_distribution_1))

LumbyTemp = np.zeros((1221, len(x)))
LumbyTemp_500 = np.zeros((1221, 28))
LumbyTemp_1 = np.zeros((1221,21))
wavelengths = np.zeros(1221)

FluxbyTemp = np.zeros((1221, len(x)))
FluxbyTemp_500 = np.zeros((1221, 28))
FluxbyTemp_1 = np.zeros((1221,21))

for i in range(len(x)):

	fits_table_filename = 'tools/kurucz93/kp00/kp00_' + str(nearestTemps[i]) + '.fits'

	hdul = fits.open(fits_table_filename)  # open a FITS file
	data = hdul[1].data  # assume the first extension is a table
	columns = data.columns

	g = 'g45'
	
	if nearestTemps[i] < 9000:
		g = 'g45'
	if nearestTemps[i] < 36000 and nearestTemps[i] >= 9000:
		g = 'g40'
	if nearestTemps[i] < 41000 and nearestTemps[i] >= 36000:
		g = 'g45'
	if nearestTemps[i] >= 41000:
		g = 'g50'
		
		
	R_star = (G * x[i] * M_solar / (10^(int(g[1:]))))**(1/2)
	if i ==1:
		print(R_star)
	# show the values in field "g45"
	wavelengths = data['WAVELENGTH']
	fluxes = data[g] * ((R_star/d)**2)
	LumbyTemp[:,i] = fluxes * 4 * math.pi * (R_star**2) * n_distribution[i]
	FluxbyTemp[:,i] = fluxes * n_distribution[i]
	if i < 28:
		LumbyTemp_500[:,i] = fluxes * 4 * math.pi * (R_star**2) * n_distribution_500[i]
		FluxbyTemp_500[:,i] = fluxes * n_distribution_500[i]
	if i < 21:
		LumbyTemp_1[:,i] = fluxes * 4 * math.pi * (R_star**2) * n_distribution_1[i]
		FluxbyTemp_1[:,i] = fluxes * n_distribution_1[i]
	#print(fluxes * n_distribution[i])
	hdul.close()


wavelengthsInMicrons = wavelengths * 1e-4
wavelengthsInCm = wavelengths * 1e-8
nu = c / wavelengthsInCm

cumulativeL_nu = np.zeros(len(LumbyTemp))
cumulativeF_nu = np.zeros(len(LumbyTemp))

for i in range(len(LumbyTemp)):
	cumLnuofnu = (sum(LumbyTemp[i,:]) * c) / (nu[i]**2)
	cumFluxofnu = (sum(FluxbyTemp[i,:]) * c) / (nu[i]**2)
	cumulativeL_nu[i] = cumLnuofnu
	cumulativeF_nu[i] = cumFluxofnu

'''
#plot spectra
fig = plt.figure()
plt.style.use('ggplot')
plt.xlabel(r'$\lambda$ [$\mu m$]')
plt.ylabel(r'$\nu$$L_{\nu}\,[erg\, s^{-1}\, Hz]$')
plt.semilogx(wavelengthsInMicrons, nu * cumulativeL_nu, color = 'k', linewidth = 0.9)
#plt.legend(frameon=False, loc='upper right', ncol=1, fontsize = 'small')
plt.tight_layout()
#plt.show()

fig.savefig("compositespectra.pdf", bbox_inches='tight')
'''

#part c
LumbyTemp_range1 = np.zeros((len(LumbyTemp), 5))
LumbyTemp_range2 = np.zeros((len(LumbyTemp), 10))
LumbyTemp_range3 = np.zeros((len(LumbyTemp), 99))
LumbyTemp_range4 = np.zeros((len(LumbyTemp), 886))

#80M_J - .5M_sun
#.5M_sun - 1M_sun
#1M_sun - 10M_sun
#10M_sun - 100M_sun

#create separate arrays for each mass range for clarity
counter = 0
for i in range(len(LumbyTemp_range1[0])):
	LumbyTemp_range1[:,i] = LumbyTemp[:, counter]
	counter += 1
for i in range(len(LumbyTemp_range2[0])):
	LumbyTemp_range2[:,i] = LumbyTemp[:, counter]
	counter += 1
for i in range(len(LumbyTemp_range3[0])):
	LumbyTemp_range3[:,i] = LumbyTemp[:, counter]
	counter += 1
for i in range(len(LumbyTemp_range4[0])):
	LumbyTemp_range4[:,i] = LumbyTemp[:, counter]
	counter += 1

cumulativeL_nu_range1 = np.zeros(len(LumbyTemp_range1))
cumulativeL_nu_range2 = np.zeros(len(LumbyTemp_range2))
cumulativeL_nu_range3 = np.zeros(len(LumbyTemp_range3))
cumulativeL_nu_range4 = np.zeros(len(LumbyTemp_range4))

#compute the sums of luminosities by wavelength
for i in range(len(LumbyTemp_range1)):
	cumLnuofnu = (sum(LumbyTemp_range1[i,:]) * c) / (nu[i]**2)
	cumulativeL_nu_range1[i] = cumLnuofnu
for i in range(len(LumbyTemp_range2)):
	cumLnuofnu = (sum(LumbyTemp_range2[i,:]) * c) / (nu[i]**2)
	cumulativeL_nu_range2[i] = cumLnuofnu
for i in range(len(LumbyTemp_range1)):
	cumLnuofnu = (sum(LumbyTemp_range3[i,:]) * c) / (nu[i]**2)
	cumulativeL_nu_range3[i] = cumLnuofnu
for i in range(len(LumbyTemp_range4)):
	cumLnuofnu = (sum(LumbyTemp_range4[i,:]) * c) / (nu[i]**2)
	cumulativeL_nu_range4[i] = cumLnuofnu

'''
#plot spectra
fig = plt.figure()
plt.style.use('ggplot')
plt.xlabel(r'$\lambda$ [$\mu m$]')
plt.ylabel(r'$\nu$$L_{\nu}\,[erg\, s^{-1}\, Hz]$')
plt.semilogx(wavelengthsInMicrons, nu * cumulativeL_nu_range1, color = 'k', linewidth = 0.9)
#plt.legend(frameon=False, loc='upper right', ncol=1, fontsize = 'small')
plt.tight_layout()
#plt.show()

fig.savefig("compositespectrarange1.pdf", bbox_inches='tight')

fig = plt.figure()
plt.style.use('ggplot')
plt.xlabel(r'$\lambda$ [$\mu m$]')
plt.ylabel(r'$\nu$$L_{\nu}\,[erg\, s^{-1}\, Hz]$')
plt.semilogx(wavelengthsInMicrons, nu * cumulativeL_nu_range2, color = 'k', linewidth = 0.9)
#plt.legend(frameon=False, loc='upper right', ncol=1, fontsize = 'small')
plt.tight_layout()
#plt.show()

fig.savefig("compositespectrarange2.pdf", bbox_inches='tight')

fig = plt.figure()
plt.style.use('ggplot')
plt.xlabel(r'$\lambda$ [$\mu m$]')
plt.ylabel(r'$\nu$$L_{\nu}\,[erg\, s^{-1}\, Hz]$')
plt.semilogx(wavelengthsInMicrons, nu * cumulativeL_nu_range3, color = 'k', linewidth = 0.9)
#plt.legend(frameon=False, loc='upper right', ncol=1, fontsize = 'small')
plt.tight_layout()
#plt.show()

fig.savefig("compositespectrarange3.pdf", bbox_inches='tight')

fig = plt.figure()
plt.style.use('ggplot')
plt.xlabel(r'$\lambda$ [$\mu m$]')
plt.ylabel(r'$\nu$$L_{\nu}\,[erg\, s^{-1}\, Hz]$')
plt.semilogx(wavelengthsInMicrons, nu * cumulativeL_nu_range4, color = 'k', linewidth = 0.9)
#plt.legend(frameon=False, loc='upper right', ncol=1, fontsize = 'small')
plt.tight_layout()
#plt.show()

fig.savefig("compositespectrarange4.pdf", bbox_inches='tight')
'''

#part d
#500 Myr

cumulativeL_nu_500 = np.zeros(len(LumbyTemp))
cumulativeF_nu_500 = np.zeros(len(LumbyTemp))

for i in range(len(LumbyTemp_500)):
	cumLnuofnu =(sum(LumbyTemp_500[i,:]) * c) / (nu[i]**2)
	cumFnuofnu =(sum(FluxbyTemp_500[i,:]) * c) / (nu[i]**2)
	cumulativeL_nu_500[i] = cumLnuofnu
	cumulativeF_nu_500[i] = cumFnuofnu

#part c
LumbyTemp_range1_500 = np.zeros((len(LumbyTemp), 5))
LumbyTemp_range2_500 = np.zeros((len(LumbyTemp), 9))
LumbyTemp_range3_500 = np.zeros((len(LumbyTemp), 12))

#80M_J - .5M_sun
#.5M_sun - 1M_sun
#1M_sun - 10M_sun
#10M_sun - 100M_sun

#create separate arrays for each mass range for clarity
counter = 0
for i in range(len(LumbyTemp_range1_500[0])):
	LumbyTemp_range1_500[:,i] = LumbyTemp_500[:, counter]
	counter += 1
for i in range(len(LumbyTemp_range2_500[0])):
	LumbyTemp_range2_500[:,i] = LumbyTemp_500[:, counter]
	counter += 1
for i in range(len(LumbyTemp_range3_500[0])):
	LumbyTemp_range3_500[:,i] = LumbyTemp_500[:, counter]
	counter += 1

cumulativeL_nu_range1_500 = np.zeros(len(LumbyTemp_range1_500))
cumulativeL_nu_range2_500 = np.zeros(len(LumbyTemp_range2_500))
cumulativeL_nu_range3_500 = np.zeros(len(LumbyTemp_range3_500))

#compute the sums of luminosities by wavelength
for i in range(len(LumbyTemp_range1_500)):
	cumLnuofnu = (sum(LumbyTemp_range1_500[i,:]) * c) / (nu[i]**2)
	cumulativeL_nu_range1_500[i] = cumLnuofnu
for i in range(len(LumbyTemp_range2_500)):
	cumLnuofnu = (sum(LumbyTemp_range2_500[i,:]) * c) / (nu[i]**2)
	cumulativeL_nu_range2_500[i] = cumLnuofnu
for i in range(len(LumbyTemp_range1_500)):
	cumLnuofnu = (sum(LumbyTemp_range3_500[i,:]) * c) / (nu[i]**2)
	cumulativeL_nu_range3_500[i] = cumLnuofnu

'''
#plot spectra
fig = plt.figure()
plt.style.use('ggplot')
plt.xlabel(r'$\lambda$ [$\mu m$]')
plt.ylabel(r'$\nu$$L_{\nu}\,[erg\, s^{-1}\, Hz]$')
plt.semilogx(wavelengthsInMicrons, nu * cumulativeL_nu_range1_500, color = 'k', linewidth = 0.9)
#plt.legend(frameon=False, loc='upper right', ncol=1, fontsize = 'small')
plt.tight_layout()
#plt.show()

fig.savefig("compositespectrarange1_500.pdf", bbox_inches='tight')

fig = plt.figure()
plt.style.use('ggplot')
plt.xlabel(r'$\lambda$ [$\mu m$]')
plt.ylabel(r'$\nu$$L_{\nu}\,[erg\, s^{-1}\, Hz]$')
plt.semilogx(wavelengthsInMicrons, nu * cumulativeL_nu_range2_500, color = 'k', linewidth = 0.9)
#plt.legend(frameon=False, loc='upper right', ncol=1, fontsize = 'small')
plt.tight_layout()
#plt.show()

fig.savefig("compositespectrarange2_500.pdf", bbox_inches='tight')

fig = plt.figure()
plt.style.use('ggplot')
plt.xlabel(r'$\lambda$ [$\mu m$]')
plt.ylabel(r'$\nu$$L_{\nu}\,[erg\, s^{-1}\, Hz]$')
plt.semilogx(wavelengthsInMicrons, nu * cumulativeL_nu_range3_500, color = 'k', linewidth = 0.9)
#plt.legend(frameon=False, loc='upper right', ncol=1, fontsize = 'small')
plt.tight_layout()
#plt.show()

fig.savefig("compositespectrarange3_500.pdf", bbox_inches='tight')
'''

#part d
#1 Gyr

cumulativeL_nu_1 = np.zeros(len(LumbyTemp))
cumulativeF_nu_1 = np.zeros(len(LumbyTemp))

for i in range(len(LumbyTemp_1)):
	cumLnuofnu =(sum(LumbyTemp_1[i,:]) * c) / (nu[i]**2)
	cumFnuofnu =(sum(FluxbyTemp_1[i,:]) * c) / (nu[i]**2)
	cumulativeL_nu_1[i] = cumLnuofnu
	cumulativeF_nu_1[i] = cumFnuofnu

#part c
LumbyTemp_range1_1 = np.zeros((len(LumbyTemp), 5))
LumbyTemp_range2_1 = np.zeros((len(LumbyTemp), 9))
LumbyTemp_range3_1 = np.zeros((len(LumbyTemp), 6))

#80M_J - .5M_sun
#.5M_sun - 1M_sun
#1M_sun - 10M_sun
#10M_sun - 100M_sun

#create separate arrays for each mass range for clarity
counter = 0
for i in range(len(LumbyTemp_range1_1[0])):
	LumbyTemp_range1_1[:,i] = LumbyTemp_1[:, counter]
	counter += 1
for i in range(len(LumbyTemp_range2_1[0])):
	LumbyTemp_range2_1[:,i] = LumbyTemp_1[:, counter]
	counter += 1
for i in range(len(LumbyTemp_range3_1[0])):
	LumbyTemp_range3_1[:,i] = LumbyTemp_1[:, counter]
	counter += 1


cumulativeL_nu_range1_1 = np.zeros(len(LumbyTemp_range1_1))
cumulativeL_nu_range2_1 = np.zeros(len(LumbyTemp_range2_1))
cumulativeL_nu_range3_1 = np.zeros(len(LumbyTemp_range3_1))

#compute the sums of luminosities by wavelength
for i in range(len(LumbyTemp_range1_1)):
	cumLnuofnu = (sum(LumbyTemp_range1_1[i,:]) * c) / (nu[i]**2)
	cumulativeL_nu_range1_1[i] = cumLnuofnu
for i in range(len(LumbyTemp_range2_1)):
	cumLnuofnu = (sum(LumbyTemp_range2_1[i,:]) * c) / (nu[i]**2)
	cumulativeL_nu_range2_1[i] = cumLnuofnu
for i in range(len(LumbyTemp_range1_1)):
	cumLnuofnu = (sum(LumbyTemp_range3_1[i,:]) * c) / (nu[i]**2)
	cumulativeL_nu_range3_1[i] = cumLnuofnu


'''
#plot spectra
fig = plt.figure()
plt.style.use('ggplot')
plt.xlabel(r'$\lambda$ [$\mu m$]')
plt.ylabel(r'$\nu$$L_{\nu}\,[erg\, s^{-1}\, Hz]$')
plt.semilogx(wavelengthsInMicrons, nu * cumulativeL_nu_range1_1, color = 'k', linewidth = 0.9)
#plt.legend(frameon=False, loc='upper right', ncol=1, fontsize = 'small')
plt.tight_layout()
#plt.show()

fig.savefig("compositespectrarange1_1.pdf", bbox_inches='tight')

fig = plt.figure()
plt.style.use('ggplot')
plt.xlabel(r'$\lambda$ [$\mu m$]')
plt.ylabel(r'$\nu$$L_{\nu}\,[erg\, s^{-1}\, Hz]$')
plt.semilogx(wavelengthsInMicrons, nu * cumulativeL_nu_range2_1, color = 'k', linewidth = 0.9)
#plt.legend(frameon=False, loc='upper right', ncol=1, fontsize = 'small')
plt.tight_layout()
#plt.show()

fig.savefig("compositespectrarange2_1.pdf", bbox_inches='tight')

fig = plt.figure()
plt.style.use('ggplot')
plt.xlabel(r'$\lambda$ [$\mu m$]')
plt.ylabel(r'$\nu$$L_{\nu}\,[erg\, s^{-1}\, Hz]$')
plt.semilogx(wavelengthsInMicrons, nu * cumulativeL_nu_range3_1, color = 'k', linewidth = 0.9)
#plt.legend(frameon=False, loc='upper right', ncol=1, fontsize = 'small')
plt.tight_layout()
#plt.show()

fig.savefig("compositespectrarange3_1.pdf", bbox_inches='tight')




#plot mass ranges together in same plot
fig = plt.figure()
plt.style.use('ggplot')
plt.xlabel(r'$\lambda$ [$\mu m$]')
plt.ylabel(r'$\nu$$L_{\nu}\,[erg\, s^{-1}\, Hz]$')
plt.loglog(wavelengthsInMicrons, nu * cumulativeL_nu_range1, color = 'k', linewidth = 0.9, label = r'80 M$_J$ $<$ M $\leq$0.5 M$_\odot$')
plt.loglog(wavelengthsInMicrons, nu * cumulativeL_nu_range2, linewidth = 0.9, label = r'0.5 M$_\odot$ $<$ M $\leq$ 1 M$_\odot$')
plt.loglog(wavelengthsInMicrons, nu * cumulativeL_nu_range3, linewidth = 0.9, label = r'1 M$_\odot$ $<$ M $\leq$10 M$_\odot$')
plt.loglog(wavelengthsInMicrons, nu * cumulativeL_nu_range4, linewidth = 0.9, label = r'10 M$_\odot$ $<$ M $\leq$100 M$_\odot$')
plt.legend(frameon=False, loc='lower center', ncol=1, fontsize = 'small')
plt.tight_layout()
#plt.show()

fig.savefig("composite0age.pdf", bbox_inches='tight')

fig = plt.figure()
plt.style.use('ggplot')
plt.xlabel(r'$\lambda$ [$\mu m$]')
plt.ylabel(r'$\nu$$L_{\nu}\,[erg\, s^{-1}\, Hz]$')
plt.loglog(wavelengthsInMicrons, nu * cumulativeL_nu_range1_500, color = 'k', linewidth = 0.9, label = r'80 M$_J$ $<$ M $\leq$0.5 M$_\odot$')
plt.loglog(wavelengthsInMicrons, nu * cumulativeL_nu_range2_500, linewidth = 0.9, label = r'0.5 M$_\odot$ $<$ M $\leq$ 1 M$_\odot$')
plt.loglog(wavelengthsInMicrons, nu * cumulativeL_nu_range3_500, linewidth = 0.9, label = r'1 M$_\odot$ $<$ M $\leq$10 M$_\odot$')
plt.legend(frameon=False, loc='lower center', ncol=1, fontsize = 'small')
plt.tight_layout()
#plt.show()

fig.savefig("composite500Myr.pdf", bbox_inches='tight')

fig = plt.figure()
plt.style.use('ggplot')
plt.xlabel(r'$\lambda$ [$\mu m$]')
plt.ylabel(r'$\nu$$L_{\nu}\,[erg\, s^{-1}\, Hz]$')
plt.loglog(wavelengthsInMicrons, nu * cumulativeL_nu_range1_1, color = 'k', linewidth = 0.9, label = r'80 M$_J$ $<$ M $\leq$0.5 M$_\odot$')
plt.loglog(wavelengthsInMicrons, nu * cumulativeL_nu_range2_1, linewidth = 0.9, label = r'0.5 M$_\odot$ $<$ M $\leq$ 1 M$_\odot$')
plt.loglog(wavelengthsInMicrons, nu * cumulativeL_nu_range3_1,linewidth = 0.9, label = r'1 M$_\odot$ $<$ M $\leq$10 M$_\odot$')
plt.legend(frameon=False, loc='lower center', ncol=1, fontsize = 'small')
plt.tight_layout()
#plt.show()

fig.savefig("composite1Gyr.pdf", bbox_inches='tight')
'''
np.savetxt('LumbyTemp.csv', cumulativeF_nu, delimiter = ',')
np.savetxt('LumbyTemp_500.csv', cumulativeF_nu_500, delimiter = ',')
np.savetxt('LumbyTemp_1.csv', cumulativeF_nu_1, delimiter = ',')
np.savetxt('wavelengthsinMicrons.csv', wavelengthsInMicrons, delimiter = ',')

