# -*- coding: utf-8 -*-
#############################   Code description   ############################


## This code defines a library with a 
## Created by: Joel Alpízar Castillo.
## TU Delft
## Version: 1.0

###############################################################################
###########################   Imported modules   ##############################

import time
from math import ceil
import csvreader


## The function ---------------------------------------------------------------

def mean_effective_T(depth):
    return 0.0318*depth + 8.01768


## The function ---------------------------------------------------------------

def T_sky(T_amb, RH = 70, N = 0):
    from math import exp
#    return 0.0552*(T_amb)**1.5 + 2.652*N
    return 0.62+0.056*(6.11*(RH/100)*exp(17.63*T_amb/(T_amb+243)))**0.5
#    return 0.7

## The function ---------------------------------------------------------------

def Wind_forced_convection(v):
    
    # Duffie and Beckman
#    return 2.8 + 3*v

    # McAdams correlation
    if v < 5:
        return 5.7 + 3.8*v
    else:   # 5 <= v < 10
        return 6.47 + v**0.78


## The function ---------------------------------------------------------------

def radiative_heat_transfer_coefficient(T_amb, T_s, epsilon):
    from scipy.constants import Stefan_Boltzmann as sigma    
#    h_r = 4*epsilon*sigma*(T_amb + 273.15)**3 
#    h_r = epsilon*sigma*((T_amb + 273.15)**2 + T_sky(T_amb)**2)*((T_amb + 273.15) + T_sky(T_amb))
#    h_r = epsilon*sigma*((T_s + 273.15)**2 + T_sky(T_amb)**2)*((T_s + 273.15) + T_sky(T_amb))
    
    return epsilon*sigma*((T_s + 273.15)**2 + (T_amb + 273.15)**2)*((T_s + 273.15) + (T_amb + 273.15))    
#    return epsilon*sigma*((T_s + 273.15)**2 + T_sky(T_amb)**2)*((T_s + 273.15) + T_sky(T_amb))

## The function ---------------------------------------------------------------

def heat_transfer_coefficient(T_amb, T_s, v, epsilon = 0.85):
    return Wind_forced_convection(v) + radiative_heat_transfer_coefficient(T_amb, T_s, epsilon)


## The function ---------------------------------------------------------------    

def thermal_radiation(T_amb):
    from scipy.constants import Stefan_Boltzmann as sigma
    return sigma*((T_amb + 273.15)**4 - (T_sky(T_amb))**4)


## The function ---------------------------------------------------------------
    
def effective_surface_T(T_amb, T_s, S, depth, alpha_0, epsilon, k, v = 5):

    delta_R = thermal_radiation(T_amb)
    h = heat_transfer_coefficient(T_amb, T_s, v, epsilon)
    
    T_e = T_amb + alpha_0*S/h - epsilon*delta_R/h
    
    b = k/(depth*h)
    
    return (b*mean_effective_T(depth) + T_e)/(1 + b)

def initial_temperature_gradient(T_end, T_surface, depth, y):
    return (T_end - T_surface)*y/depth + T_surface

## The function ---------------------------------------------------------------
    
def Thermal_diffusivity(k, density, c):
    
    return k/(density*c)

## The function ---------------------------------------------------------------

def estimate_soil_Temperature(T_amb, S, depth, start_day, end_day, k, density_soil, c_p, alpha_0, epsilon, dt, dy, T_surf_0 = False):

    from numpy import linspace, arange, empty, transpose
    
    # Boundary conditions
    T_surface = effective_surface_T(T_amb[0], T_amb[0], S[0], depth, alpha_0, epsilon, k)
    T_end = mean_effective_T(depth)
    
    
#    alpha = Thermal_diffusivity(k, density_soil, c_p)    # Soil's thermal diffusivity, dimensionless.
    r = Thermal_diffusivity(k, density_soil, c_p)*dt/dy**2 
    n = int(depth/dy)
    
    y = linspace(dy/2, depth-dy/2, n)
    t_final = int((end_day - start_day)*24*3600)         # In s
    t = arange(0, t_final, dt)
    
    
    if not T_surf_0:
        T_0 = [initial_temperature_gradient(T_end, T_surface, depth, i) for i in y] # Create a linear relationship between Tsurf and Tend
    else:
        T_0 = [initial_temperature_gradient(T_end, T_surf_0, depth, i) for i in y] # Create a linear relationship between Tsurf and Tend
    
    dTdt = empty(n)
    
    
    T = T_0
    T_stored = [T]
    
    counter = 1
    
    
    for j in range(1,len(t)):
    
        for i in range(1,n-1):
            dTdt[i] = r*(T[i-1] - 2*T[i] + T[i+1])
        dTdt[0] = r*(T_surface - 2*T[0] + T[1])         # Upper boundary condition
        dTdt[n-1] = r*(T[n-2] - 2*T[n-1] + T_end)       # Lower boundary condition
        
        T = T + dTdt    
        T_surface = effective_surface_T(T_amb[j], T[0], S[j], depth, alpha_0, epsilon, k)
        
        if counter == 15:
            counter = 1
            T_stored.append(T)
        else:
            counter += 1
        
    #    print('For the time ', round(j/60, 2), ', timestep ', j, ', the ambient temperature is ', round(T_amb[j], 2), ', the radiation is ', round(S[j], 2), ', and the surface temperature is ', round(T[0], 2))
    #    plt.clf()    
    #    plt.figure(1)
    #    plt.plot(y,T)
    #    plt.axis([0, depth, -10, 50])
    #    plt.xlabel('Depth [m]')
    #    plt.ylabel('Temperature [°C]')
    #    plt.show()
    #    plt.pause(0.01)
    
    
    base = empty(((end_day - start_day)*24*4,int(depth/dy)))
    
    
    for timestep in range(len(T_stored)):
        for y in range(len(T_stored[timestep])):
            base[timestep][y] = T_stored[timestep][y]
    
    base = transpose(base)
    
    return base


###############################################################################
################################   Simulations ################################
    

# Paremeters
    

dt = 60
dy = 0.01
depth = 6
T_surf_0 = 0.1
y_f = 2.5

start_day = 0
end_day = 365
horizon = 4*0

k = 1.19                # Soil's conductivity, in W/(m-K).
density_soil = 2029.80  # Soil's density, in kg/m3.
c_p = 2*756.108         # Soil's specific heat capacity, in J/(kg-K).


alpha_0 = 0.25          # Soil's thermal absorptivity, dimensionless.
epsilon = 0.25          # Soil's thermal emissivity, dimensionless.


CSVDataTamb = csvreader.read_data(csv='Tamb_15min.csv', address='')
CSVDataRad = csvreader.read_data(csv='Radiation_1min.csv', address='')
CSVDataTamb.data2array()
CSVDataRad.data2array()
T_amb2 = [i[0] for i in CSVDataTamb.ar[start_day*24*4:(end_day + ceil(horizon/(24*4)))*24*4]]

from numpy import ones
T_amb = []
for i in T_amb2:
    T_amb += list(ones(15)*i)

S = [i[0]/3600 for i in CSVDataRad.ar[start_day*24*60:(end_day + ceil(horizon/(24*60)))*24*60]]


# Simulation
start = time.time()    # The timer is initializad.
    
soil_grad = estimate_soil_Temperature(T_amb, S, depth, start_day, end_day, k, density_soil, c_p, alpha_0, epsilon, dt, dy, T_surf_0)

end = time.time()    # The timer is initializad.
totalelapsed = end - start  # The total time is calculated.


# Plotting


base = soil_grad[0:int(y_f/dy)][:]

import matplotlib.pyplot as plt
from numpy import arange, savetxt

plt.rcParams.update({
#    "text.usetex": True,
    "font.family": "Times New Roman",
    'font.size': 18
})

#savetxt('Soil_dy_10cm.csv', base, delimiter =";", fmt ='% s')

        
fig1, (ax1) = plt.subplots(1, 1, constrained_layout=True, sharey=True)
a1 = plt.imshow(base, interpolation='nearest', aspect='auto', origin='upper', vmin = -15, vmax = 35)
cbar1 = fig1.colorbar(a1, ax=ax1)
cbar1.set_ticks(arange(-15, 36, step=5), [round(i,0) for i in arange(-15, 36, step=5)])
cbar1.set_label('Temperature of the soil, T$_{s}$(y), [°C]')
#ax1.set_title('Year 2021')
#plt.xticks(arange(1, int((end_day - start_day)*24*4), step=24*31*4), ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])  # Set label locations. 
plt.xticks(arange(1, int((end_day - start_day)*24*4), step=24*31*8), ['Jan', 'Mar', 'May', 'Jul', 'Sep', 'Nov'])  # Set label locations. 
ax1.set_xlabel('Year 2021 [Months]')
plt.yticks(arange(0, int(y_f/dy), step=20), [round(i,1) for i in arange(0, y_f, step=0.2)])
ax1.set_ylabel('Depth, y, [m]')
            
#fig2 = plt.figure(2)
#plt.plot(base[0,:], 'g', label='Surface') 
#plt.plot(base[25,:], 'b', label='$y = 25$ cm') 
#plt.plot(base[50,:], 'r', label='$y = 50$ cm') 
#plt.plot(base[75,:], 'c', label='$y = 75$ cm') 
#plt.plot(base[100,:], 'y', label='$y = 100$ cm') 
#plt.plot(base[125,:], 'k', label='$y = 125$ cm') 
#plt.plot(base[149,:], 'm', label='$y = 150$ cm') 
#plt.legend(loc='upper right')
#plt.grid()
#plt.xlim([start_day*24*4, end_day*24*4])
#plt.xticks(arange(1, (end_day-start_day)*24*4, step=24*31*4), ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])  # Set label locations. 
##plt.xlabel('Time [h]')
#plt.ylabel('Soil temperature, $T_{s}$, [°C]')
#plt.show()    
#
#fig3, (axs3) = plt.subplots(4, sharex=True)
#axs3[0].plot(T_amb2[100*24*4:101*24*4+1], 'm', label='Ambient', linestyle='dotted') 
#axs3[0].plot(base[0,100*24*4:101*24*4+1], 'g', label='Surface') 
#axs3[0].plot(base[10,100*24*4:101*24*4+1], 'b', label='$y = 10$ cm') 
#axs3[0].plot(base[20,100*24*4:101*24*4+1], 'r', label='$y = 20$ cm') 
#axs3[0].plot(base[30,100*24*4:101*24*4+1], 'c', label='$y = 30$ cm') 
#axs3[0].plot(base[40,100*24*4:101*24*4+1], 'y', label='$y = 40$ cm') 
#axs3[0].plot(base[50,100*24*4:101*24*4+1], 'k', label='$y = 50$ cm') 
#axs3[0].set_ylim([-5, 21])
#axs3[0].set_yticks(arange(-5, 21, step=5), [round(i,0) for i in arange(-5, 21, step=5)])
#axs3[0].set_title('Spring (10 of April)')
#axs3[0].set_ylabel('T$_{s}$(y), [°C]')
##axs3[0].legend(loc='lower center', ncol=6)
#axs3[0].grid()
#
#axs3[1].plot(T_amb2[200*24*4:201*24*4+1], 'm', label='Ambient', linestyle='dotted') 
#axs3[1].plot(base[0,200*24*4:201*24*4+1], 'g', label='Surface') 
#axs3[1].plot(base[10,200*24*4:201*24*4+1], 'b', label='$y = 10$ cm') 
#axs3[1].plot(base[20,200*24*4:201*24*4+1], 'r', label='$y = 20$ cm') 
#axs3[1].plot(base[30,200*24*4:201*24*4+1], 'c', label='$y = 30$ cm') 
#axs3[1].plot(base[40,200*24*4:201*24*4+1], 'y', label='$y = 40$ cm') 
#axs3[1].plot(base[50,200*24*4:201*24*4+1], 'k', label='$y = 50$ cm') 
#axs3[1].set_ylim([0, 26])
#axs3[1].set_yticks(arange(0, 26, step=5), [round(i,0) for i in arange(0, 26, step=5)])
#axs3[1].set_title('Summer (19 of July)')
#axs3[1].set_ylabel('T$_{s}$(y), [°C]')
##axs3[1].legend(loc='lower center', ncol=6)
#axs3[1].grid()
#
#axs3[2].plot(T_amb2[300*24*4:301*24*4+1], 'm', label='Ambient', linestyle='dotted') 
#axs3[2].plot(base[0,300*24*4:301*24*4+1], 'g', label='Surface') 
#axs3[2].plot(base[10,300*24*4:301*24*4+1], 'b', label='$y = 10$ cm') 
#axs3[2].plot(base[20,300*24*4:301*24*4+1], 'r', label='$y = 20$ cm') 
#axs3[2].plot(base[30,300*24*4:301*24*4+1], 'c', label='$y = 30$ cm') 
#axs3[2].plot(base[40,300*24*4:301*24*4+1], 'y', label='$y = 40$ cm') 
#axs3[2].plot(base[50,300*24*4:301*24*4+1], 'k', label='$y = 50$ cm') 
#axs3[2].set_ylim([0, 16])
#axs3[2].set_yticks(arange(0, 16, step=5), [round(i,0) for i in arange(0, 16, step=5)])
#axs3[2].set_title('Autumn (27 of October)')
#axs3[2].set_ylabel('T$_{s}$(y), [°C]')
##axs3[2].legend(loc='lower center', ncol=6)
#axs3[2].grid()
#
#axs3[3].plot(T_amb2[1*24*4:2*24*4+1], 'm', label='Ambient', linestyle='dotted') 
#axs3[3].plot(base[0,1*24*4:2*24*4+1], 'g', label='Surface') 
#axs3[3].plot(base[10,1*24*4:2*24*4+1], 'b', label='$y = 10$ cm') 
#axs3[3].plot(base[20,1*24*4:2*24*4+1], 'r', label='$y = 20$ cm') 
#axs3[3].plot(base[30,1*24*4:2*24*4+1], 'c', label='$y = 30$ cm') 
#axs3[3].plot(base[40,1*24*4:2*24*4+1], 'y', label='$y = 40$ cm') 
#axs3[3].plot(base[50,1*24*4:2*24*4+1], 'k', label='$y = 50$ cm') 
#axs3[3].set_ylim([-10, 6])
#axs3[3].set_yticks(arange(-10, 6, step=5), [round(i,0) for i in arange(-10, 6, step=5)])
#axs3[3].set_title('Winter (2 of January)')
#axs3[3].set_ylabel('T$_{s}$(y), [°C]')
##axs3[3].legend(loc='lower center', ncol=6)
#axs3[3].grid()

#plt.xticks(arange(0, 24*4+1, step=4*2), [i for i in range(0,25,2)])  # Set label locations. 
#plt.xlim([0, 24*4])
#labels = ['Ambient', 'Surface', '$y = 10$ cm', '$y = 20$ cm', '$y = 30$ cm', '$y = 40$ cm', '$y = 50$ cm',]
#fig3.legend(labels, loc='lower center', ncol=7)
#plt.show()

fig4 = plt.figure(4)
for i in range(int(len(base[0])/(48*4))):
    plt.plot(base[:,i*12*4], 'g')
plt.xlim([0, 250])
plt.ylim([-15, 20])
plt.ylabel('Soil temperature, $T_{s}$, [°C]')
plt.xlabel('Depth, y, [m]')
plt.grid()
plt.show