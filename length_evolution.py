from math import *
from matplotlib import pyplot as plt
import numpy as np
import random
import os.path
from scipy.optimize import curve_fit

force = 120

def read_data(force,k):
	data = open('Data/forcesdevel_pf'+str(force)+'_'+str(k)+'.dat','r') 
	t, l = np.loadtxt(data, usecols=(0,1), unpack=True)
	data.close()

	return t, l

t100_0, l100_0 = read_data(force,0)
t100_1, l100_1 = read_data(force,1)
t100_2, l100_2 = read_data(force,2)
t100_3, l100_3 = read_data(force,3)
t100_4, l100_4 = read_data(force,4)

def funct(t, tcut, a, tau, c):
		return a * np.exp(- (t-tcut)/tau ) + c

def fit_decay(tcut, times, lengths):
	
	# selct relevant part of data
	t_part = times[tcut:tcut+25]
	l_part = lengths[tcut:tcut+25]

	# fit function
	def func(t, a, tau, c):
		return a * np.exp(- (t-tcut)/tau ) + c

	# performe fit
	popt, pcov = curve_fit(func, t_part, l_part)

	return popt, pcov, t_part, l_part

f2, b, t_part2, l_part2 = fit_decay(200, t100_2, l100_2)
print f2, b
f3, b, t_part3, l_part3 = fit_decay(300, t100_3, l100_3)
print f3, b
f4, b, t_part4, l_part4 = fit_decay(400, t100_4, l100_4)
print f4, b

#PLOT
# plot evolution
fig = plt.figure(1)
ax = fig.add_subplot(111)
plt.gcf().subplots_adjust(top=0.85) # to be able to read the lables
plt.gcf().subplots_adjust(bottom=0.15)
plt.gcf().subplots_adjust(left=0.15)
axes = plt.gca()
axes.xaxis.set_tick_params(labelsize=15, width=1)
axes.yaxis.set_tick_params(labelsize=15, width=1)

plt.plot(t100_4, l100_4, label='t = 400')
plt.plot(t100_3, l100_3, label='t = 300')
plt.plot(t100_2, l100_2, label='t = 200')

plt.plot(t_part2, l_part2)
plt.plot(t100_1, l100_1, label='t = 100')
plt.plot(t100_0, l100_0, label='t = 50')

plt.xlabel('t [au]', fontsize = 16)
plt.ylabel('length [au]', fontsize = 16)
plt.legend(loc='upper right', numpoints = 1, fontsize = 12)


plt.figure(2)
plt.plot(t_part2, l_part2, 'bo', label="data")
plt.plot(t_part2, funct(t_part2, 200, *f2), 'b-', label="fit")
plt.plot(t_part2, l_part3, 'go', label="data")
plt.plot(t_part2, funct(t_part3, 300, *f3), 'g-', label="fit")
plt.plot(t_part2, l_part4, 'ro', label="data")
plt.plot(t_part2, funct(t_part4, 400, *f4), 'r-', label="fit")

plt.legend()

plt.show()