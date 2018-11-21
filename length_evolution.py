from math import *
from matplotlib import pyplot as plt
import numpy as np
import random
import os.path
from scipy.optimize import curve_fit

force = 100
xi = '0'

def read_data(force,k):
	#data = open('TestData/time_evolution_pf'+str(force)+'_'+str(k)+'.dat','r') 
	#t, l, cmx, cmy, cmz = np.loadtxt(data, usecols=(0, 2, 3, 4, 5), unpack=True)
	data = open('Data/time_evolution_pf'+str(force)+'_'+str(k)+'.dat','r') 
	t, l, cmx, cmy, cmz = np.loadtxt(data, usecols=(0, 1, 2, 3, 4), unpack=True)
	data.close()

	v = cacl_vel(cmx, cmy, cmz, t)

	return t, l, v

def read_data_individual(force,k, xi):
	#data = open('TestData/time_evolution_pf'+str(force)+'_'+str(k)+'.dat','r') 
	#t, l, cmx, cmy, cmz = np.loadtxt(data, usecols=(0, 2, 3, 4, 5), unpack=True)
	data = open('T1/xi'+str(xi)+'/job_'+str(k)+'/Data/time_evolution_pf'+str(force)+'_'+str(k)+'.dat','r') 
	t, l, cmx, cmy, cmz = np.loadtxt(data, usecols=(0, 1, 2, 3, 4), unpack=True)
	data.close()

	v = cacl_vel(cmx, cmy, cmz, t)

	return t, l, v

def average_data(force,xi):
	lav = [0.0 for i in range(101)]
	vav = [0.0 for i in range(100)]
	for i in range(1,11):
		t, l, v = read_data_individual(force,i,xi)
		lav = lav + l
		vav = vav + np.array(v)
	lav = lav/10.0
	vav = vav/10.0
	return lav, vav


def cacl_vel(cmx, cmy, cmz, t):
	# calulate center of mass
	cm = []
	for i in range(len(cmx)):
		average = np.sqrt(cmx[i]*cmx[i] + cmy[i]*cmy[i] + cmz[i]*cmz[i])
		cm.append(average)

	# calculate velocity
	v = [0]
	for i in range(1,len(cm)-1):
		vel = (cm[i+1] - cm[i-1])/(t[i+1] - t[i-1])
		v.append(vel)
	return v

lav0, vav0 = average_data(force,'0')
lav001, vav001 = average_data(force,'001')
lav01, vav01 = average_data(force,'01')
lav1, vav1 = average_data(force,'1')
lav10, vav10 = average_data(force,'10')


t1, l1, v1 = read_data_individual(force,1,xi)
#t2, l2, v2 = read_data_individual(force,2,xi)
#t3, l3, v3 = read_data_individual(force,3,xi)
#t4, l4, v4 = read_data_individual(force,4,xi)


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


#plt.plot(t1, l1, label='run 1')
#plt.plot(t2, l2, label='run 2')
#plt.plot(t3, l3, label='run 3')
#plt.plot(t4, l4, label='run 4')
plt.plot(t1, lav0, label=r'$\xi$ = 0.0')
plt.plot(t1, lav001, label=r'$\xi$ = 0.01')
plt.plot(t1, lav01, label=r'$\xi$ = 0.1')
plt.plot(t1, lav1, label=r'$\xi$ = 1.0')
plt.plot(t1, lav10, label=r'$\xi$ = 10.0')
#plt.plot(t1, l1, label='xi = 0.0001')
#plt.plot(t2, l2, label='xi = 0.01')
#plt.plot(t3, l3, label='xi = 1.0')
#plt.plot(t4, l4, label='xi = 10.0')
#plt.plot(t5, l5, label='xi = 50.0')

plt.xlabel('t [au]', fontsize = 16)
plt.ylabel('length [au]', fontsize = 16)
plt.legend(loc='lower right', numpoints = 1, fontsize = 12)

fig = plt.figure(2)
ax = fig.add_subplot(111)
plt.gcf().subplots_adjust(top=0.85) # to be able to read the lables
plt.gcf().subplots_adjust(bottom=0.15)
plt.gcf().subplots_adjust(left=0.15)
axes = plt.gca()
axes.xaxis.set_tick_params(labelsize=15, width=1)
axes.yaxis.set_tick_params(labelsize=15, width=1)

t1 = np.delete(t1, -1)
#t2 = np.delete(t2, -1)
#t3 = np.delete(t3, -1)
#t4 = np.delete(t4, -1)
#t5 = np.delete(t5, -1)


#plt.plot(t1, v1, label='run 1')
#plt.plot(t2, v2, label='run 2')
#plt.plot(t3, v3, label='run 3')
#plt.plot(t4, v4, label='run 4')
plt.plot(t1, vav0, label=r'$\xi$ = 0.0')
plt.plot(t1, vav001, label=r'$\xi$ = 0.01')
plt.plot(t1, vav01, label=r'$\xi$ = 0.1')
plt.plot(t1, vav1, label=r'$\xi$ = 1.0')
plt.plot(t1, vav10, label=r'$\xi$ = 10.0')

#plt.semilogy(t1, v1, label='xi = 0.0001')
#plt.semilogy(t2, v2, label='xi = 0.01')
#plt.semilogy(t3, v3, label='xi = 1.0')
#plt.semilogy(t4, v4, label='xi = 10.0')
#plt.semilogy(t5, v5, label='xi = 50.0')

plt.xlabel('t [au]', fontsize = 16)
plt.ylabel('v [au]', fontsize = 16)
plt.legend(loc='lower right', numpoints = 1, fontsize = 12)

plt.show()


# FIT
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

'''
f2, b, t_part2, l_part2 = fit_decay(200, t100_2, l100_2)
print f2, b
f3, b, t_part3, l_part3 = fit_decay(300, t100_3, l100_3)
print f3, b
f4, b, t_part4, l_part4 = fit_decay(400, t100_4, l100_4)
print f4, b
'''

'''
plt.figure(2)
plt.plot(t_part2, l_part2, 'bo', label="data")
plt.plot(t_part2, funct(t_part2, 200, *f2), 'b-', label="fit")
plt.plot(t_part2, l_part3, 'go', label="data")
plt.plot(t_part2, funct(t_part3, 300, *f3), 'g-', label="fit")
plt.plot(t_part2, l_part4, 'ro', label="data")
plt.plot(t_part2, funct(t_part4, 400, *f4), 'r-', label="fit")

plt.legend()
'''
