# -*- coding: utf-8 -*-
"""
Created on Sat Nov 11 07:12:14 2023

@author: jalpizarcastil
"""

from numpy import arange
import matplotlib.pyplot as plt 

import csvreader
from MCES_library import update_TESS, HP_Power

############################### TESS Validation ###############################
###############################################################################

CSVDataTsoil = csvreader.read_data(csv='Soil_dy_10cm.csv', address='')
CSVDataTsoil.data2cols()
Tsoil = CSVDataTsoil.col


T_charge = 95 + 273
T_network = 40 + 273
T_amb = 25 + 273
T_max = 90 + 273
T_min = 50 + 273
T_0 = T_max

time_step = []


T_TESS = []
Qdot_charge = []
Qdot_TESS = []
Qdot_SD = []
TESS_state = [T_0, 0, 0]


point = 0
step = 0
start_discharge = False
start_charge = False

while point < 3:
    if point == 0:
        Qdot_HP_TESS = 0
        [TESS_state[0], TESS_state[1], TESS_state[2]] = update_TESS(False, T_network, TESS_state[0], T_soil = Tsoil[step], Qdot_HP = Qdot_HP_TESS)
        
        if TESS_state[0] <= T_min:
            point +=1
            start_charge = True

    elif point == 1:
        [l, Qdot_HP_TESS] = HP_Power(True, T_amb, T_0, T_sup = T_charge)
        [TESS_state[0], TESS_state[1], TESS_state[2]] = update_TESS(False, T_network, TESS_state[0], T_soil = Tsoil[step], Qdot_HP = Qdot_HP_TESS)
        
        if TESS_state[0] >= T_max:
            point +=1
            start_discharge = True            

    elif point == 2:
        Qdot_HP_TESS = 0
        [TESS_state[0], TESS_state[1], TESS_state[2]] = update_TESS(True, T_network, TESS_state[0], T_soil = Tsoil[step], Qdot_HP = Qdot_HP_TESS)
        
        if TESS_state[0] <= T_min:
            point +=1
    
    Qdot_charge.append(Qdot_HP_TESS)
    T_TESS.append(TESS_state[0])
    Qdot_TESS.append(TESS_state[1])
    Qdot_SD.append(TESS_state[2])
    
    time_step.append(step*0.25)
    if start_charge:
        start_charge_time = time_step[-1]
        start_charge = False
    if start_discharge:
        start_discharge_time = time_step[-1]   
        start_discharge = False            
        
    step +=1
    
###############################################################################    
    
plt.rcParams.update({
#    "text.usetex": True,
    "font.family": "Times New Roman",
    'font.size': 18
})    
    
f1 = plt.figure(1)
plt.plot(time_step, [i-273 for i in T_TESS], 'r', label='TESS Temperature, $T_{TESS}$')
#plt.legend(loc='lower center')
plt.grid()
plt.xlim([0, max(time_step)])
plt.xticks(arange(0, 1101, step=100))  # Set label locations. 
plt.xlabel('Time [h]')
plt.ylabel('TESS Temperature, $T_{TESS}$, [°C]')
#plt.title('Temperature inside the house')
plt.show() 

f2 = plt.figure(2)
plt.plot(time_step, [i/1000 for i in Qdot_charge], 'r', label='Charge')
plt.plot(time_step, [i/1000 for i in Qdot_TESS], 'b', label='Discharge')
plt.legend(loc='lower right', ncol = 2)
plt.grid()
plt.xlim([start_charge_time, max(time_step)])
plt.xticks(arange(start_charge_time, max(time_step), step=10))  # Set label locations. 
plt.xlabel('Time [h]')
plt.ylabel('Thermal power, $\dot{Q}_{TESS}$, [kW]')
#plt.title('Temperature inside the house')
plt.show() 

f3 = plt.figure(3)
plt.plot(time_step, [i/1000 for i in Qdot_SD], 'r', label='Self-Discharge')
#plt.legend(loc='lower center', ncol = 2)
plt.grid()
plt.xlim([0, max(time_step)])
plt.xticks(arange(0, 1101, step=100))  # Set label locations. 
plt.xlabel('Time [h]')
plt.ylabel('Thermal power transfered to the soil, $\dot{Q}_{TESS}^{SD}$, [kW]')
#plt.title('Temperature inside the house')
plt.show() 