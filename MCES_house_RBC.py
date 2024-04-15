# -*- coding: utf-8 -*-
#############################   Code description   ############################


## This code defines a library with a 
## Created by: Joel Alpízar Castillo.
## TU Delft
## Version: 2.0

###############################################################################
###########################   Imported modules   ##############################
import MCES_library
import csvreader
import time
from math import ceil

from numpy import arange, array, empty, zeros#, cos, pi, append


###############################################################################
##############################   Parameters   #################################

start_day = 0
end_day = 364
horizon = 4*0

#CSVDataPV = csvreader.read_data(csv='PV_15min.csv', address='')
CSVDataTamb = csvreader.read_data(csv='Tamb_15min.csv', address='')
CSVDataP_Load = csvreader.read_data(csv='Load_Profile_15min.csv', address='', delim=',')
CSVDataRad = csvreader.read_data(csv='Radiation_1min.csv', address='')
CSVDataTsoil = csvreader.read_data(csv='Soil_dy_10cm.csv', address='')
#CSVDataPV.data2array()
CSVDataTamb.data2array()
CSVDataP_Load.data2array()
CSVDataRad.data2array()
CSVDataTsoil.data2cols()
#P_PV_av = [i[0]*n_modules*module_power/module_power_ref/1000 for i in CSVDataPV.ar]
T_amb = [i[0]+273 for i in CSVDataTamb.ar[start_day*24*4:(end_day + ceil(horizon/(24*4)))*24*4]]
P_Load = [i for i in CSVDataP_Load.ar[0]]
a = arange(0,len(CSVDataRad.ar),15)
#G = array([CSVDataRad.ar[i][0] for i in a])
G = [CSVDataRad.ar[i][0] for i in a]
Tsoil = CSVDataTsoil.col

# Initial conditions
dt = 60*15                          # In s
t_final = int((end_day - start_day)*24*3600/dt)         # In s

# Thermal demand
T_0 = 20 + 273              # In K
T_set_day = [T_0 - 3]*int((6-0)*4) + [T_0]*int((22-6)*4)+ [T_0 - 3]*int((24-22)*4)
T_set = array(T_set_day*((end_day + ceil(horizon/(24*4)))-start_day))


# PV
#n_modules = 10
#module_power_ref = 0.315
#module_power = 0.400  
#P_PV = []
#CSVDataPV = csvreader.read_data(csv='PV_15min.csv', address='')
#CSVDataPV.data2array()
#P_PV_av = [i[0]*n_modules*module_power/module_power_ref/1000 for i in CSVDataPV.ar]

# PVT
PVT_active = False
n_modules_list = [0] # arange(0,11) # 
T_sup = 50 + 273
#T_glass_0 = 40
#T_PV_0 = 56
#T_a_0 = 55
#T_f_0 = 53.5
T_glass_0 = 20
T_PV_0 = 18
T_a_0 = 20
T_f_0 = 20
T_PVT_layers = [T_glass_0, T_PV_0, T_a_0, T_f_0]
T_tank_PVT = 50



# TESS
TESS_vol_list =  [0] # arange(0,11) #
depth = 0.2
dy = 0.1
Height_TESS = 1.8
#thickness_TESS = 0.4
TESS_active = False
TESS_mode = 'cycling' #  ['cycling', 'full']
T_TESS_0 = 75 + 273         # In K
Tsoil = [i[int(depth/dy):int((depth + Height_TESS)/dy)] for i in Tsoil]

PVT_2_TESS = False
HP_2_TESS = False

if TESS_mode == 'cycling':
    TESS_min_op_T = T_sup
    TESS_max_op_T = 90 + 273
elif TESS_mode == 'full':
    TESS_min_op_T = 90 + 273
    TESS_max_op_T = 95 + 273


# HP
HP_active = True

# Thermal Network
T4 = 40 + 273

# BESS
SoCmax = 0.9
SoCmin = 0.2
SoC_0_BESS = 0.5
BESS_capacity = 3.36
P_BESS_max = 0
P_Grid_max = 0
Capacity_BESS = 3.36
SoC_BESS = [SoC_0_BESS]                     # In %
charge_efficiency = 0.943
discharge_efficiency = 0.943

# Grid
CO2_factor = 0.523

#time_registry = []
elapsed_time_cycle_registry = empty(shape=(len(n_modules_list), len(TESS_vol_list)), dtype='float')
cold_days_registry = empty(shape=(len(n_modules_list), len(TESS_vol_list)), dtype='float')
Qdot_Losses_registry = empty(shape=(len(n_modules_list), len(TESS_vol_list)), dtype='float')
Qdot_PVT_Network_registry = empty(shape=(len(n_modules_list), len(TESS_vol_list)), dtype='float')
Qdot_SC_TESS_registry = empty(shape=(len(n_modules_list), len(TESS_vol_list)), dtype='float')
Qdot_TESS_registry = empty(shape=(len(n_modules_list), len(TESS_vol_list)), dtype='float')
Qdot_TESS_SD_registry = empty(shape=(len(n_modules_list), len(TESS_vol_list)), dtype='float')
Qdot_HP_registry = empty(shape=(len(n_modules_list), len(TESS_vol_list)), dtype='float')
Qdot_HP_TESS_registry = empty(shape=(len(n_modules_list), len(TESS_vol_list)), dtype='float')

T_TESS_registry = [[[] for i in range(len(TESS_vol_list))] for j in range(len(TESS_vol_list))]

P_PVT_registry = empty(shape=(len(n_modules_list), len(TESS_vol_list)), dtype='float')
P_BESS_stored_registry = empty(shape=(len(n_modules_list), len(TESS_vol_list)), dtype='float')
P_BESS_provided_registry = empty(shape=(len(n_modules_list), len(TESS_vol_list)), dtype='float')
P_HP_registry = empty(shape=(len(n_modules_list), len(TESS_vol_list)), dtype='float')
P_Grid_consumed_registry = empty(shape=(len(n_modules_list), len(TESS_vol_list)), dtype='float')
P_Grid_returned_registry = empty(shape=(len(n_modules_list), len(TESS_vol_list)), dtype='float')

SoC_BESS_registry = [[[] for i in range(len(TESS_vol_list))] for j in range(len(TESS_vol_list))]


start = time.time()

for TESS_vol in TESS_vol_list:

    if TESS_vol != 0 and TESS_active == True:
        TESS_active_iteration = True
    else:
        TESS_active_iteration = False   
    
    for n_modules in n_modules_list:
        
        if n_modules != 0 and PVT_active == True:
            PVT_active_iteration = True
        else:
            PVT_active_iteration = False
         

        ###############################################################################
        ##############################   Simulation   #################################
        
        Thermal_Components = [PVT_active_iteration, TESS_active_iteration, HP_active]     # PVT, TESS, HP
          
           
        start_cycle = time.time()    # The timer is initializad.
         
        [T_in, Qdot_PVT_Network, Qdot_SC_TESS, Qdot_TESS, Qdot_HP, Qdot_HP_TESS, T_registry, Qdot_Losses, T_TESS, Qdot_TESS_SD, fluid_registry, T_tank_PVT_registry, P_Grid, P_PVT, P_HP, P_BESS, SoC_BESS] = MCES_library.Thermal_Electrical_model(Thermal_Components, T_0, T4, T_TESS_0, T_tank_PVT, T_PVT_layers, TESS_min_op_T, TESS_max_op_T, SoC_0_BESS, T_amb, T_set, G, Tsoil, P_Load, t_final, n_modules, P_BESS_max, P_Grid_max, Capacity_BESS, PVT_2_TESS, HP_2_TESS, m = TESS_vol*1000)
                                                                                                                                                                                                                           
        end_cycle = time.time()    # The timer is ended.
            
        elapsed_time_cycle = (end_cycle - start_cycle)/60

        cold_days = 0
        for i in range(int(len(T_in[2:])/96)):
            if min(T_in[2+i*96:2+(i+1)*96]) < 16 + 273:
                cold_days += 1        
        
        ###############################################################################
        ###############################   CSV logging   ################################
        
#        file_name = 'PVT_' + str(n_modules) + '_TESS_' + str(int(TESS_vol)) + '.csv'
#        
#        results = zip(T_in[2:], Qdot_PVT_Network, Qdot_SC_TESS, Qdot_TESS, Qdot_HP, T_registry[0], T_registry[1], T_registry[2], T_registry[3], T_registry[4], Qdot_Losses, T_TESS[1:], Qdot_TESS_SD, fluid_registry[0], fluid_registry[1], fluid_registry[2], fluid_registry[3], fluid_registry[4][1:], T_tank_PVT_registry, P_Grid, P_PVT, P_HP, P_BESS, SoC_BESS[1:])
#        
#        
#        if 'file_name' in globals():
#            import csv
#            with open(file_name, "w") as f:
#                writer = csv.writer(f)
#                for row in results:
#                    writer.writerow(row)
#                f.close()

        ###############################################################################
        ################################   Logging   ##################################
        
        elapsed_time_cycle_registry[n_modules, TESS_vol] = elapsed_time_cycle
        cold_days_registry[n_modules, TESS_vol] = cold_days
        Qdot_Losses_registry[n_modules, TESS_vol] = -0.25*sum(Qdot_Losses)/1000
        Qdot_PVT_Network_registry[n_modules, TESS_vol] = 0.25*sum(Qdot_PVT_Network)/1000
        Qdot_SC_TESS_registry[n_modules, TESS_vol] = 0.25*sum(Qdot_SC_TESS)/1000
        Qdot_TESS_registry[n_modules, TESS_vol] = 0.25*sum(Qdot_TESS)/1000
        Qdot_TESS_SD_registry[n_modules, TESS_vol] = 0.25*sum(Qdot_TESS_SD)/1000
        Qdot_HP_registry[n_modules, TESS_vol] = 0.25*sum(Qdot_HP)/1000
        Qdot_HP_TESS_registry[n_modules, TESS_vol] = 0.25*sum(Qdot_HP_TESS)/1000
        P_PVT_registry[n_modules, TESS_vol] = 0.25*sum(P_PVT)
        P_BESS_stored_registry[n_modules, TESS_vol] = -0.25*sum([i for i in P_BESS if i < 0])
        P_BESS_provided_registry[n_modules, TESS_vol] = 0.25*sum([i for i in P_BESS if i > 0])
        P_HP_registry[n_modules, TESS_vol] = 0.25*sum(P_HP)
        P_Grid_consumed_registry[n_modules, TESS_vol] = 0.25*sum([i for i in P_Grid if i > 0])
        P_Grid_returned_registry[n_modules, TESS_vol] = -0.25*sum([i for i in P_Grid if i < 0])
        
        T_TESS_registry[n_modules][TESS_vol] = T_TESS
        SoC_BESS_registry[n_modules][TESS_vol] = SoC_BESS
        
        ###############################################################################
        ###############################   Printing   ##################################
        
        
        print('----------------------------------------------------------------------')
        print('')
        print('PVT_' + str(n_modules) + '_TESS_' + str(int(TESS_vol)))
        print('Elapsed time per run: ', elapsed_time_cycle, ' min')
        print('')
        print('Cold days: ', cold_days)
        print('Thermal energy losses: ',  -0.25*sum(Qdot_Losses)/1000, ' kWh')
        print('Thermal energy from the PVT to the thermal network: ',  0.25*sum(Qdot_PVT_Network)/1000, ' kWh')
        print('Thermal energy from the PVT to the TESS: ',  0.25*sum(Qdot_SC_TESS)/1000, ' kWh')
        print('Thermal energy from the TESS to the thermal network: ',  0.25*sum(Qdot_TESS)/1000, ' kWh')
        print('Thermal energy lost from the TESS to the soil: ',  0.25*sum(Qdot_TESS_SD)/1000, ' kWh')
        print('Thermal energy from the HP to the TESS: ',  0.25*sum(Qdot_HP_TESS)/1000, ' kWh')
        print('Thermal energy from the HP to the thermal network: ',  0.25*sum(Qdot_HP)/1000, ' kWh')
        print('Energy produced from the PVT: ', 0.25*sum(P_PVT), ' kWh')
        print('Energy stored in the BESS: ', -0.25*sum([i for i in P_BESS if i < 0]), ' kWh')
        print('Energy provided by the BESS: ', 0.25*sum([i for i in P_BESS if i > 0]), ' kWh')
        print('Energy consumed by the HP: ', 0.25*sum(P_HP), 'kWh')
        print('Energy consumed from the grid: ', 0.25*sum([i for i in P_Grid if i > 0]), ' kWh')
        print('Energy returned to the grid: ', -0.25*sum([i for i in P_Grid if i < 0]), ' kWh')
        print('')
    

COP = [[(Qdot_HP_registry[m][v] + Qdot_HP_TESS_registry[m][v])/(P_HP_registry[m][v]) for m in n_modules_list.tolist()] for v in TESS_vol_list.tolist()]

end = time.time()
elapsed_time = (end - start)/60

print('----------------------------------------------------------------------')
print('')
print('Total elapsed time: ', elapsed_time, ' min')
print('----------------------------------------------------------------------')

###############################################################################
################################   Ploting   ##################################

#

import matplotlib.pyplot as plt

plt.rcParams.update({
#    "text.usetex": True,
    "font.family": "Times New Roman",
    'font.size': 16
})

#fig10, (ax10) = plt.subplots(1, 1, constrained_layout=True, sharey=True)
#a10 = plt.imshow(elapsed_time_cycle_registry, cmap='plasma', interpolation='bicubic', origin='lower')   # , vmin = -15, vmax = 35, , interpolation='bicubic'
#cbar10 = fig10.colorbar(a10, ax=ax10)
##cbar10.set_ticks(arange(-15, 36, step=5), [round(i,0) for i in arange(-15, 36, step=5)])
#cbar10.set_label('Time per case, t, [s]')
##ax10.set_title('Year 2021')
##plt.xticks(arange(1, int((end_day - start_day)*24*4), step=24*31*4), ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])  # Set label locations. 
#ax10.set_xlabel('TESS volume, [m$^{3}$]')
##plt.yticks(arange(0, int(y_f/dy), step=10), [round(i,1) for i in arange(0, y_f, step=0.1)])
#ax10.set_ylabel('Number of PVT modules, [-]')
#
#
fig11, (ax11) = plt.subplots(1, 1, constrained_layout=True, sharey=True)
a11 = plt.imshow(cold_days_registry, cmap='plasma', interpolation='bicubic', origin='lower')   # , vmin = -15, vmax = 35
cbar11 = fig11.colorbar(a11, ax=ax11)
#cbar11.set_ticks(arange(-15, 36, step=5), [round(i,0) for i in arange(-15, 36, step=5)])
cbar11.set_label('Number of cold days, [-]')
#ax11.set_title('Year 2021')
#plt.xticks(arange(1, int((end_day - start_day)*24*4), step=24*31*4), ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])  # Set label locations. 
ax11.set_xlabel('TESS volume, [m$^{3}$]')
#plt.yticks(arange(0, int(y_f/dy), step=10), [round(i,1) for i in arange(0, y_f, step=0.1)])
ax11.set_ylabel('Number of PVT modules, [-]')


fig12, (ax12) = plt.subplots(1, 1, constrained_layout=True, sharey=True)
a12 = plt.imshow(Qdot_Losses_registry, cmap='plasma', interpolation='bicubic', origin='lower')   # , vmin = -15, vmax = 35
cbar12 = fig12.colorbar(a12, ax=ax12)
#cbar12.set_ticks(arange(-15, 36, step=5), [round(i,0) for i in arange(-15, 36, step=5)])
cbar12.set_label('Thermal losses of the house, Q$_{L}$, [kWh]')
#ax12.set_title('Year 2021')
#plt.xticks(arange(1, int((end_day - start_day)*24*4), step=24*31*4), ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])  # Set label locations. 
ax12.set_xlabel('TESS volume, [m$^{3}$]')
#plt.yticks(arange(0, int(y_f/dy), step=10), [round(i,1) for i in arange(0, y_f, step=0.1)])
ax12.set_ylabel('Number of PVT modules, [-]')


fig13, (ax13) = plt.subplots(1, 1, constrained_layout=True, sharey=True)
a13 = plt.imshow(Qdot_PVT_Network_registry, cmap='plasma', interpolation='bicubic', origin='lower')   # , vmin = -15, vmax = 35
cbar13 = fig13.colorbar(a13, ax=ax13)
#cbar13.set_ticks(arange(-15, 36, step=5), [round(i,0) for i in arange(-15, 36, step=5)])
cbar13.set_label('Thermal energy from the PVT to the network, Q$_{PVT}^{N}$, [kWh]')
#ax13.set_title('Year 2021')
#plt.xticks(arange(1, int((end_day - start_day)*24*4), step=24*31*4), ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])  # Set label locations. 
ax13.set_xlabel('TESS volume, [m$^{3}$]')
#plt.yticks(arange(0, int(y_f/dy), step=10), [round(i,1) for i in arange(0, y_f, step=0.1)])
ax13.set_ylabel('Number of PVT modules, [-]')

fig14, (ax14) = plt.subplots(1, 1, constrained_layout=True, sharey=True)
a14 = plt.imshow(Qdot_SC_TESS_registry, cmap='plasma', interpolation='bicubic', origin='lower')   # , vmin = -15, vmax = 35
cbar14 = fig14.colorbar(a14, ax=ax14)
#cbar14.set_ticks(arange(-15, 36, step=5), [round(i,0) for i in arange(-15, 36, step=5)])
cbar14.set_label('Thermal energy from the PVT to the TESS, Q$_{PVT}^{TESS}$, [kWh]')
#ax14.set_title('Year 2021')
#plt.xticks(arange(1, int((end_day - start_day)*24*4), step=24*31*4), ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])  # Set label locations. 
ax14.set_xlabel('TESS volume, [m$^{3}$]')
#plt.yticks(arange(0, int(y_f/dy), step=10), [round(i,1) for i in arange(0, y_f, step=0.1)])
ax14.set_ylabel('Number of PVT modules, [-]')


fig15, (ax15) = plt.subplots(1, 1, constrained_layout=True, sharey=True)
a15 = plt.imshow(Qdot_TESS_registry, cmap='plasma', interpolation='bicubic', origin='lower')   # , vmin = -15, vmax = 35
cbar15 = fig15.colorbar(a15, ax=ax15)
#cbar15.set_ticks(arange(-15, 36, step=5), [round(i,0) for i in arange(-15, 36, step=5)])
cbar15.set_label('Thermal energy from the TESS, Q$_{TESS}$, [kWh]')
#ax15.set_title('Year 2021')
#plt.xticks(arange(1, int((end_day - start_day)*24*4), step=24*31*4), ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])  # Set label locations. 
ax15.set_xlabel('TESS volume, [m$^{3}$]')
#plt.yticks(arange(0, int(y_f/dy), step=10), [round(i,1) for i in arange(0, y_f, step=0.1)])
ax15.set_ylabel('Number of PVT modules, [-]')

fig28, (ax28) = plt.subplots(1, 1, constrained_layout=True, sharey=True)
a28 = plt.imshow(Qdot_TESS_SD_registry, cmap='plasma', interpolation='bicubic', origin='lower')   # , vmin = -15, vmax = 35
cbar28 = fig28.colorbar(a28, ax=ax28)
#cbar28.set_ticks(arange(-15, 36, step=5), [round(i,0) for i in arange(-15, 36, step=5)])
cbar28.set_label('Thermal energy lost from the TESS, Q$_{TESS}^{SD}$, [kWh]')
#ax28.set_title('Year 2021')
#plt.xticks(arange(1, int((end_day - start_day)*24*4), step=24*31*4), ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])  # Set label locations. 
ax28.set_xlabel('TESS volume, [m$^{3}$]')
#plt.yticks(arange(0, int(y_f/dy), step=10), [round(i,1) for i in arange(0, y_f, step=0.1)])
ax28.set_ylabel('Number of PVT modules, [-]')


fig16, (ax16) = plt.subplots(1, 1, constrained_layout=True, sharey=True)
a16 = plt.imshow(Qdot_HP_registry, cmap='plasma', interpolation='bicubic', origin='lower')   # , vmin = -15, vmax = 35
cbar16 = fig16.colorbar(a16, ax=ax16)
#cbar16.set_ticks(arange(-15, 36, step=5), [round(i,0) for i in arange(-15, 36, step=5)])
cbar16.set_label('Thermal energy from the HP, Q$_{HP}$, [kWh]')
#ax16.set_title('Year 2021')
#plt.xticks(arange(1, int((end_day - start_day)*24*4), step=24*31*4), ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])  # Set label locations. 
ax16.set_xlabel('TESS volume, [m$^{3}$]')
#plt.yticks(arange(0, int(y_f/dy), step=10), [round(i,1) for i in arange(0, y_f, step=0.1)])
ax16.set_ylabel('Number of PVT modules, [-]')

fig27, (ax27) = plt.subplots(1, 1, constrained_layout=True, sharey=True)
a27 = plt.imshow(Qdot_HP_TESS_registry, cmap='plasma', interpolation='bicubic', origin='lower')   # , vmin = -15, vmax = 35
cbar27 = fig27.colorbar(a27, ax=ax27)
#cbar27.set_ticks(arange(-15, 36, step=5), [round(i,0) for i in arange(-15, 36, step=5)])
cbar27.set_label('Thermal energy from the HP to the TESS, Q$_{HP}^{TESS}$, [kWh]')
#ax27.set_title('Year 2021')
#plt.xticks(arange(1, int((end_day - start_day)*24*4), step=24*31*4), ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])  # Set label locations. 
ax27.set_xlabel('TESS volume, [m$^{3}$]')
#plt.yticks(arange(0, int(y_f/dy), step=10), [round(i,1) for i in arange(0, y_f, step=0.1)])
ax27.set_ylabel('Number of PVT modules, [-]')


fig17, (ax17) = plt.subplots(1, 1, constrained_layout=True, sharey=True)
a17 = plt.imshow(P_PVT_registry, cmap='plasma', interpolation='bicubic', origin='lower')   # , vmin = -15, vmax = 35
cbar17 = fig17.colorbar(a17, ax=ax17)
#cbar17.set_ticks(arange(-15, 36, step=5), [round(i,0) for i in arange(-15, 36, step=5)])
cbar17.set_label('Electric energy produced by the PVT, E$_{PVT}$, [kWh]')
#ax17.set_title('Year 2021')
#plt.xticks(arange(1, int((end_day - start_day)*24*4), step=24*31*4), ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])  # Set label locations. 
ax17.set_xlabel('TESS volume, [m$^{3}$]')
#plt.yticks(arange(0, int(y_f/dy), step=10), [round(i,1) for i in arange(0, y_f, step=0.1)])
ax17.set_ylabel('Number of PVT modules, [-]')


fig18, (ax18) = plt.subplots(1, 1, constrained_layout=True, sharey=True)
a18 = plt.imshow(P_BESS_stored_registry, cmap='plasma', interpolation='bicubic', origin='lower')   # , vmin = -15, vmax = 35
cbar18 = fig18.colorbar(a18, ax=ax18)
#cbar18.set_ticks(arange(-15, 36, step=5), [round(i,0) for i in arange(-15, 36, step=5)])
cbar18.set_label('Electric energy stored in the BESS, E$_{BESS}^{in}$, [kWh]')
#ax18.set_title('Year 2021')
#plt.xticks(arange(1, int((end_day - start_day)*24*4), step=24*31*4), ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])  # Set label locations. 
ax18.set_xlabel('TESS volume, [m$^{3}$]')
#plt.yticks(arange(0, int(y_f/dy), step=10), [round(i,1) for i in arange(0, y_f, step=0.1)])
ax18.set_ylabel('Number of PVT modules, [-]')


fig19, (ax19) = plt.subplots(1, 1, constrained_layout=True, sharey=True)
a19 = plt.imshow(P_BESS_provided_registry, cmap='plasma', interpolation='bicubic', origin='lower')   # , vmin = -15, vmax = 35
cbar19 = fig19.colorbar(a19, ax=ax19)
#cbar19.set_ticks(arange(-15, 36, step=5), [round(i,0) for i in arange(-15, 36, step=5)])
cbar19.set_label('Electric energy supplied by the BESS, E$_{BESS}^{out}$, [kWh]')
#ax19.set_title('Year 2021')
#plt.xticks(arange(1, int((end_day - start_day)*24*4), step=24*31*4), ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])  # Set label locations. 
ax19.set_xlabel('TESS volume, [m$^{3}$]')
#plt.yticks(arange(0, int(y_f/dy), step=10), [round(i,1) for i in arange(0, y_f, step=0.1)])
ax19.set_ylabel('Number of PVT modules, [-]')


fig20, (ax20) = plt.subplots(1, 1, constrained_layout=True, sharey=True)
a20 = plt.imshow(P_HP_registry, cmap='plasma', interpolation='bicubic', origin='lower')   # , vmin = -15, vmax = 35
cbar20 = fig20.colorbar(a20, ax=ax20)
#cbar20.set_ticks(arange(-15, 36, step=5), [round(i,0) for i in arange(-15, 36, step=5)])
cbar20.set_label('Electric energy consumed by the HP, E$_{HP}$, [kWh]')
#ax20.set_title('Year 2021')
#plt.xticks(arange(1, int((end_day - start_day)*24*4), step=24*31*4), ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])  # Set label locations. 
ax20.set_xlabel('TESS volume, [m$^{3}$]')
#plt.yticks(arange(0, int(y_f/dy), step=10), [round(i,1) for i in arange(0, y_f, step=0.1)])
ax20.set_ylabel('Number of PVT modules, [-]')


fig21, (ax21) = plt.subplots(1, 1, constrained_layout=True, sharey=True)
a21 = plt.imshow(P_Grid_consumed_registry, cmap='plasma', interpolation='bicubic', origin='lower')   # , vmin = -15, vmax = 35
cbar21 = fig21.colorbar(a21, ax=ax21)
#cbar21.set_ticks(arange(-15, 36, step=5), [round(i,0) for i in arange(-15, 36, step=5)])
cbar21.set_label('Electric energy consumed from the grid, E$_{G}^{in}$, [kWh]')
#ax21.set_title('Year 2021')
#plt.xticks(arange(1, int((end_day - start_day)*24*4), step=24*31*4), ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])  # Set label locations. 
ax21.set_xlabel('TESS volume, [m$^{3}$]')
#plt.yticks(arange(0, int(y_f/dy), step=10), [round(i,1) for i in arange(0, y_f, step=0.1)])
ax21.set_ylabel('Number of PVT modules, [-]')


fig22, (ax22) = plt.subplots(1, 1, constrained_layout=True, sharey=True)
a22 = plt.imshow(P_Grid_returned_registry, cmap='plasma', interpolation='bicubic', origin='lower')   # , vmin = -15, vmax = 35
cbar22 = fig22.colorbar(a22, ax=ax22)
#cbar22.set_ticks(arange(-15, 36, step=5), [round(i,0) for i in arange(-15, 36, step=5)])
cbar22.set_label('Electric energy sent back to the grid, E$_{G}^{in}$, [kWh]')
#ax22.set_title('Year 2021')
#plt.xticks(arange(1, int((end_day - start_day)*24*4), step=24*31*4), ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])  # Set label locations. 
ax22.set_xlabel('TESS volume, [m$^{3}$]')
#plt.yticks(arange(0, int(y_f/dy), step=10), [round(i,1) for i in arange(0, y_f, step=0.1)])
ax22.set_ylabel('Number of PVT modules, [-]')        
        

f23 = plt.figure()
plt.plot([i-273 for i in T_TESS_registry[1][1]], 'g', label='1 PVT - 1 TESS')
plt.plot([i-273 for i in T_TESS_registry[1][10]], 'b', label='1 PVT - 10 TESS')
plt.plot([i-273 for i in T_TESS_registry[5][5]], 'r', label='5 PVT - 5 TESS')
plt.plot([i-273 for i in T_TESS_registry[10][1]], 'c', label='10 PVT - 1 TESS')
plt.plot([i-273 for i in T_TESS_registry[10][10]], 'y', label='10 PVT - 10 TESS')
plt.legend(loc='lower center', ncol = 5, prop={'size': 10})
plt.grid()
plt.xlim([0, end_day*24])
plt.xticks(arange(1, 365*24*4, step=24*31*4), ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])  # Set label locations. 
#plt.xlabel('Time [h]')
plt.ylabel('TESS temperature, T$_{TESS}$, [°C]')
plt.show()

#
f24 = plt.figure()
#data = [[i-273 for i in T_TESS_registry[0][2]], [i-273 for i in T_TESS_registry[0][4]], [i-273 for i in T_TESS_registry[0][6]], [i-273 for i in T_TESS_registry[0][8]], [i-273 for i in T_TESS_registry[0][10]]]
data = [[t-273 for t in i] for i in T_TESS_registry[0]]
#labels = ['1 m$^{3}$', '2 m$^{3}$', '3 m$^{3}$', '4 m$^{3}$', '5 m$^{3}$', '6 m$^{3}$', '7 m$^{3}$', '8 m$^{3}$', '9 m$^{3}$', '10 m$^{3}$']
labels = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10']
plt.xlabel('TESS volume, V, [m$^{3}$]')
plt.ylabel('TESS temperature, T$_{TESS}$, [°C]')
#labels = ['2 m$^{3}$', '4 m$^{3}$', '6 m$^{3}$', '8 m$^{3}$', '10 m$^{3}$']
plt.boxplot(data[1:], showfliers = False, labels = labels)     # , 'g', label=
#plt.legend(loc='lower center', ncol = 5)
plt.show()
#
#
#f25 = plt.figure()
#plt.plot([i*100 for i in SoC_BESS_registry[10][10]], 'y', label='10 PVT - 10 TESS')
#plt.plot([i*100 for i in SoC_BESS_registry[10][1]], 'c', label='10 PVT - 1 TESS')
#plt.plot([i*100 for i in SoC_BESS_registry[5][5]], 'r', label='5 PVT - 5 TESS')
#plt.plot([i*100 for i in SoC_BESS_registry[1][10]], 'b', label='1 PVT - 10 TESS')
#plt.plot([i*100 for i in SoC_BESS_registry[1][1]], 'g', label='1 PVT - 1 TESS')
#plt.legend(loc='lower center', ncol = 5)
#plt.show()
#
#
#f26 = plt.figure()
#data = [[i*100 for i in SoC_BESS_registry[1][1] if i > 0.2], [i*100 for i in SoC_BESS_registry[10][1] if i > 0.2], [i*100 for i in SoC_BESS_registry[5][5] if i > 0.2], [i*100 for i in SoC_BESS_registry[1][10] if i > 0.2], [i*100 for i in SoC_BESS_registry[10][10] if i > 0.2]]
#labels = ['1 PVT - 1 TESS', '10 PVT - 1 TESS', '5 PVT - 5 TESS', '1 PVT - 10 TESS', '10 PVT - 10 TESS']
#
#plt.boxplot(data, showfliers = False, labels = labels)     # , 'g', label=
#
##plt.legend(loc='lower center', ncol = 5)
#plt.show()
#
#
#

thermal_performance = zeros(shape=(len(n_modules_list[0]), len(TESS_vol_list[0])))

for n in range(len(n_modules_list[0])-1):
    for vol in range(len(TESS_vol_list[0])-1):
#        print(100*(Qdot_TESS_registry[n + 1,vol + 1] - 1000*vol*4200*(T_TESS_0 - T_sup)/(1000*3600))/(Qdot_SC_TESS_registry[n + 1,vol + 1] + Qdot_HP_TESS_registry[n + 1,vol + 1]))
        thermal_performance[n+1, vol+1] = 100*(Qdot_TESS_registry[n + 1,vol + 1] - 1000*vol*4200*(T_TESS_0 - T_sup)/(1000*3600))/(Qdot_SC_TESS_registry[n + 1,vol + 1] + Qdot_HP_TESS_registry[n + 1,vol + 1])

fig27, (ax27) = plt.subplots(1, 1, constrained_layout=True, sharey=True)
a27 = plt.imshow(thermal_performance, cmap='plasma', interpolation='bicubic', origin='lower', vmin = 50, vmax = 70)   # , vmin = -15, vmax = 35
cbar27 = fig27.colorbar(a27, ax=ax27)
#cbar27.set_ticks(arange(-15, 36, step=5), [round(i,0) for i in arange(-15, 36, step=5)])
cbar27.set_label('Thermal performance, [%]')
#ax27.set_title('Year 2021')
#plt.xticks(arange(1, int((end_day - start_day)*24*4), step=24*31*4), ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])  # Set label locations. 
ax27.set_xlabel('TESS volume, [m$^{3}$]')
ax27.set_xlim([1,10])
ax27.set_ylim([1,10])
#plt.yticks(arange(0, int(y_f/dy), step=10), [round(i,1) for i in arange(0, y_f, step=0.1)])
ax27.set_ylabel('Number of PVT modules, [-]')     



thermal_electric_performance = empty(shape=(len(n_modules_list), len(TESS_vol_list)), dtype='float')

for n in range(len(n_modules_list)-1):
    for vol in range(len(TESS_vol_list)-1):
        
        thermal_electric_performance[n+1, vol+1] = (Qdot_PVT_Network_registry[n + 1,vol + 1] + Qdot_SC_TESS_registry[n + 1,vol + 1])/P_PVT_registry[n + 1,vol + 1]

fig28, (ax28) = plt.subplots(1, 1, constrained_layout=True, sharey=True)
a28 = plt.imshow(thermal_electric_performance, cmap='plasma', interpolation='bicubic', origin='lower')   # , vmin = -15, vmax = 35
cbar28 = fig28.colorbar(a28, ax=ax28)
#cbar28.set_ticks(arange(-15, 36, step=5), [round(i,0) for i in arange(-15, 36, step=5)])
cbar28.set_label('Termal/electric PVT output ratio, [-]')
#ax28.set_title('Year 2021')
#plt.xticks(arange(1, int((end_day - start_day)*24*4), step=24*31*4), ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])  # Set label locations. 
ax28.set_xlabel('TESS volume, [m$^{3}$]')
ax28.set_xlim([1,10])
ax28.set_ylim([1,10])
#plt.yticks(arange(0, int(y_f/dy), step=10), [round(i,1) for i in arange(0, y_f, step=0.1)])
ax28.set_ylabel('Number of PVT modules, [-]') 


Net_grid = [i - j for i, j in zip(P_Grid_consumed_registry, P_Grid_returned_registry)]

fig29, (ax29) = plt.subplots(1, 1, constrained_layout=True, sharey=True)
a29 = plt.imshow(Net_grid, cmap='plasma', interpolation='bicubic', origin='lower')   # , vmin = -15, vmax = 35
cbar29 = fig29.colorbar(a29, ax=ax29)
#cbar29.set_ticks(arange(-15, 36, step=5), [round(i,0) for i in arange(-15, 36, step=5)])
cbar29.set_label('Net energy consumption, E$_{net}$, [kWh]')
#ax29.set_title('Year 2021')
#plt.xticks(arange(1, int((end_day - start_day)*24*4), step=24*31*4), ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])  # Set label locations. 
ax29.set_xlabel('TESS volume, [m$^{3}$]')
#ax29.set_xlim([1,10])
#ax29.set_ylim([1,10])
#plt.yticks(arange(0, int(y_f/dy), step=10), [round(i,1) for i in arange(0, y_f, step=0.1)])
ax29.set_ylabel('Number of PVT modules, [-]') 


electric_performance = 100*Qdot_TESS_registry/P_HP_registry

fig30, (ax30) = plt.subplots(1, 1, constrained_layout=True, sharey=True)
a30 = plt.imshow(electric_performance, cmap='plasma', interpolation='bicubic', origin='lower', vmin = 65, vmax = 90)   # , vmin = -15, vmax = 35
cbar30 = fig30.colorbar(a30, ax=ax30)
#cbar30.set_ticks(arange(-15, 36, step=5), [round(i,0) for i in arange(-15, 36, step=5)])
cbar30.set_label('Electric performance, [%]')
#ax30.set_title('Year 2021')
#plt.xticks(arange(1, int((end_day - start_day)*24*4), step=24*31*4), ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])  # Set label locations. 
ax30.set_xlabel('TESS volume, [m$^{3}$]')
ax30.set_xlim([1,10])
ax30.set_ylim([1,10])
#plt.yticks(arange(0, int(y_f/dy), step=10), [round(i,1) for i in arange(0, y_f, step=0.1)])
ax30.set_ylabel('Number of PVT modules, [-]') 


Net_CO2 = [(i - j)*CO2_factor for i, j in zip(P_Grid_consumed_registry, P_Grid_returned_registry)]

fig31, (ax31) = plt.subplots(1, 1, constrained_layout=True, sharey=True)
a31 = plt.imshow(Net_CO2, cmap='plasma', interpolation='bicubic', origin='lower')   # , vmin = -15, vmax = 35
cbar31 = fig31.colorbar(a31, ax=ax31)
#cbar31.set_ticks(arange(-15, 36, step=5), [round(i,0) for i in arange(-15, 36, step=5)])
cbar31.set_label('Total equivalent emissions, [kg$_{CO_{2},eq}$]')
#ax31.set_title('Year 2021')
#plt.xticks(arange(1, int((end_day - start_day)*24*4), step=24*31*4), ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])  # Set label locations. 
ax31.set_xlabel('TESS volume, [m$^{3}$]')
#ax31.set_xlim([1,10])
#ax31.set_ylim([1,10])
#plt.yticks(arange(0, int(y_f/dy), step=10), [round(i,1) for i in arange(0, y_f, step=0.1)])
ax31.set_ylabel('Number of PVT modules, [-]') 


#import matplotlib.pyplot as plt
#
#
f1 = plt.figure()
plt.plot([i-273 for i in T_amb], 'r', label='Ambient temperature, $T_{amb}$')
plt.plot([i-273 for i in T_set], 'g', label='Setpoint temperature inside the house, $T_{set}$')
#plt.plot([i-273 for i in T_amb[0:4*h]], 'r', label='Ambient temperature, $T_{amb}$')
plt.plot([i-273 for i in T_in[1:]], 'b', label='Temperature inside the house, $T_{in}$')
plt.legend(loc='lower center')
plt.grid()
plt.xlim([0, end_day*24])
plt.xticks(arange(1, 365*24*4, step=24*31*4), ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])  # Set label locations. 
#plt.xlabel('Time [h]')
plt.ylabel('Temperature [°C]')
#plt.title('Temperature inside the house')
plt.show()   

#
f2 = plt.figure()
plt.plot([i/1000 for i in Qdot_HP], 'g', label='$\dot{Q}_{HP}$')
plt.plot([i/1000 for i in Qdot_TESS], 'r', label='$\dot{Q}_{TESS}$')
plt.plot([i/1000 for i in Qdot_PVT_Network], 'b', label='$\dot{Q}_{PVT}$')
plt.legend(loc='upper center', ncol = 3)
plt.grid()
plt.xlim([0, end_day*24])
plt.xticks(arange(1, 365*24*4, step=24*31*4), ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])  # Set label locations. 
#plt.xlabel('Time [h]')
plt.ylabel('Thermal power [kW]')
#plt.title('Temperature inside the house')
plt.show()   
#
#f3 = plt.figure()
#plt.plot([i-273 for i in T_registry[4]], 'g', label='$T_{4}$')
#plt.plot([i-273 for i in T_registry[3]], 'r', label='$T_{3}$')
#plt.plot([i-273 for i in T_registry[2]], 'b', label='$T_{2}$')
#plt.plot([i-273 for i in T_registry[1]], 'c', label='$T_{1}$')
##plt.plot([i-273 for i in T_registry[0]], 'k', label='$T_{0}$')
#plt.legend(loc='lower center', ncol = 4)
#plt.grid()
#plt.xlim([0, end_day*24])
#plt.xticks(arange(1, 365*24*4, step=24*31*4), ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])  # Set label locations. 
##plt.xlabel('Time [h]')
#plt.ylabel('Temperature [°C]')
##plt.title('Temperature inside the house')
#plt.show()   
#
f4 = plt.figure()
plt.plot([-i/1000 for i in Qdot_Losses], linewidth=0.5)
plt.grid()
plt.xlim([0, end_day*24])
plt.ylim([0,3])
plt.xticks(arange(1, 365*24*4, step=24*31*4), ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])  # Set label locations. 
#plt.xlabel('Time [h]')
plt.ylabel('Thermal losses, $\dot{Q}_{L}$, [kW]')
#plt.title('Temperature inside the house')
plt.show()  
#
f5 = plt.figure()
plt.plot([i-273 for i in T_TESS])
plt.grid()
plt.xlim([0, end_day*24])
plt.xticks(arange(1, 365*24*4, step=24*31*4), ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])  # Set label locations. 
#plt.xlabel('Time [h]')
plt.ylabel('Temperature in the TESS, $T_{TESS}$, [°C]')
#plt.title('Temperature in the TESS')
plt.show()
#
##f6 = plt.figure()
##plt.plot([i/3600 for i in G[start_day*24*4:end_day*24*4]], 'c', label='$G$')
##plt.plot([i-273 for i in T_amb[start_day*24*4:end_day*24*4]], 'r', label='$T_{amb}$')
###plt.plot([i-273 for i in fluid_registry[0]], 'g', label='$T_{glass}$')
###plt.plot([i-273 for i in fluid_registry[1]], 'b', label='$T_{PV}$')
###plt.plot([i-273 for i in fluid_registry[2]], 'c', label='$T_{a}$')
##plt.plot([i for i in fluid_registry[4]], 'b', label='$T_{out}$')
##plt.plot([i for i in fluid_registry[3]], 'y', label='$T_{f}$')
##
##
##plt.grid()
##f6.legend(loc='lower center', ncol = 3)
##plt.show()
#

#T_start = 178
#T_end = 182
T_start = 324
T_end = 328
f6, ax1 = plt.subplots()

#ax1.set_xlabel('June 27 - June 28 - June 29 - June 30 - July 1')
ax1.set_xlabel('November 20 - November 21 - November 22 - November 23 - November 24')
#ax1.set_xticks(arange(1, 365*24*4, step=24*31*4), ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])  # Set label locations. 
#ax1.set_xticks(arange(1, 365*24*4, step=24*31*4))
#ax1.set_xticks([])
ax1.set_ylabel('Global irradiance, G, [W/m$^{2}$]') # , color='b'  # we already handled the x-label with ax1
ax1.plot([i/3600 for i in G[0][T_start*96:(T_end+1)*96]], 'teal', label='$G$')
ax1.set_ylim([0, 1000])
#ax1.set_yticks(arange(0, 1001, step=30))
#ax1.tick_params(axis='y', labelcolor='b')

ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

ax2.set_ylabel('Temperature, [°C]') # , color='r'
ax2.plot([i-273 for i in T_amb[0][T_start*96:(T_end+1)*96]], 'red', label='$T_{amb}$')
ax2.plot([i for i in fluid_registry[0][0][T_start*96:(T_end+1)*96]], 'olive', label='$T_{glass}$')
ax2.plot([i for i in fluid_registry[1][0][T_start*96:(T_end+1)*96]], 'blue', label='$T_{PV}$')
ax2.plot([i for i in fluid_registry[2][0][T_start*96:(T_end+1)*96]], 'navy', label='$T_{a}$')
ax2.plot([i for i in fluid_registry[4][0][T_start*96:(T_end+1)*96]], 'magenta', label='$T_{net}$')
ax2.plot([i for i in fluid_registry[3][0][T_start*96:(T_end+1)*96]], 'yellow', label='$T_{f}$')
ax2.plot(T_tank_PVT_registry[0][T_start*96:(T_end+1)*96], 'black', label='$T_{PVT}$') # T_{tank, PVT}$
#T_tank_PVT_registry
ax2.set_xticks(arange(1, (T_end-T_start+1)*96, step=6*4))
ax2.set_xlim([1,(T_end-T_start+1)*96])
#ax2.set_ylim([0, 1001])
#ax2.set_xticks(arange(0, int(t_final/2)+1, step=300))
#ax2.set_yticks(arange(0, 1001, step=100))
ax2.tick_params(axis='y', labelcolor='r')

f6.legend(loc='upper center', ncol = 8, fontsize="10")

f6.tight_layout()  # otherwise the right y-label is slightly clipped
ax1.grid()
plt.show()


#Imbalance = [Load - PV for Load, PV in zip(P_Load, P_PV_av[1:])]
##T_start = 180
##T_end = 186
#T_start = 0
#T_end = 6
#
#
#
#plt.figure()
#plt.plot(Imbalance[T_start*96:(T_end+1)*96], label='$P_{Grid}$')
#plt.xlim([1,(T_end-T_start+1)*96])
#plt.xticks(arange(1, (T_end-T_start+1)*96, step=4*24))  # Set label locations. 
##plt.xlabel('Jun 29 - Jun 30 - Jul 1 - Jul 2 - Jul 3 - Jul 4 - Jul 5 - Jul 6')
#plt.xlabel('Jan 1 - Jan 2 - Jan 3 - Jan 4 - Jan 5 - Jan 6 - Jan 7 - Jan 8')
#plt.ylabel('Grid power, $P_{Grid}$, [kW]')
##plt.legend(loc='lower center', ncol = 5)
#plt.grid()
#plt.show()



f7 = plt.figure()
plt.plot(P_Grid, 'y', label='$P_{Grid}$')
plt.plot(P_HP, 'r', label='$P_{HP}$')
plt.plot(P_Load, 'c', label='$P_{Load}$')
plt.plot(P_BESS, 'k', label='$P_{BESS}$')
plt.plot(P_PVT, 'b', label='$P_{PVT}$')
plt.xlim([0, end_day*24])
plt.xticks(arange(1, 365*24*4, step=24*31*4), ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])  # Set label locations. 
plt.ylabel('Electric power, $P$, [kW]')
plt.legend(loc='lower center', ncol = 5)
plt.grid()
plt.show()
#
#
#f8 = plt.figure()
#plt.plot([i*100 for i in SoC_BESS], label='$SoC_{BESS}$')
#plt.xlim([0, end_day*24])
#plt.xticks(arange(1, 365*24*4, step=24*31*4), ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])  # Set label locations. 
#plt.ylabel('State-of-Charge of the BESS, $SoC_{BESS}$, [%]')
##plt.legend(loc='lower center', ncol = 1)
#plt.grid()
#plt.show()
#
#f8 = plt.figure()
#plt.plot(Qdot_TESS_SD, label='$\dot{Q}_{TESS}^{SD}$')
#plt.xlim([0, end_day*24])
#plt.xticks(arange(1, 365*24*4, step=24*31*4), ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])  # Set label locations. 
#plt.ylabel('TESS self-discharge themal power, $\dot{Q}_{TESS}^{SD}$, [W]')
##plt.legend(loc='lower center', ncol = 1)
#plt.grid()
#plt.show()
#
