# -*- coding: utf-8 -*-
#############################   Code description   ############################


## This code defines a library with a 
## Created by: Joel Alpízar Castillo.
## TU Delft
## Version: 1.0

###############################################################################
###########################   Imported modules   ##############################
import time
from numpy import append, array, arange, cos, pi

###############################################################################
###########################   Functions definition  ###########################

###############################################################################
###########################   PVT model definition  ###########################

def T_sky(T_amb, RH = 70, N = 0):
    from math import exp
#    return 0.0552*(T_amb)**1.5 + 2.652*N
    return 0.62+0.056*(6.11*(RH/100)*exp(17.63*T_amb/(T_amb+243)))**0.5
#    return 0.7

def Wind_forced_convection(v):
    
    # Duffie and Beckman
#    return 2.8 + 3*v

    # McAdams correlation
    if v < 5:
        return 5.7 + 3.8*v
    else:   # 5 <= v < 10
        return 6.47 + v**0.78

def radiative_heat_transfer_coefficient(T_amb, T_s, epsilon):
    from scipy.constants import Stefan_Boltzmann as sigma    
    return epsilon*sigma*((T_s + 273.15)**2 + T_sky(T_amb)**2)*((T_s + 273.15) + T_sky(T_amb))

def Thermal_diffusivity(k, density, c):
    
    return k/(density*c)

def radiative_heat_transfer_coefficient_two_plates(T_p1, epsilon_1, T_p2, epsilon_2):
    from scipy.constants import Stefan_Boltzmann as sigma   
    
    return sigma*((T_p1 + 273.15)**2 + (T_p2 + 273.15)**2)*((T_p1 + 273.15) + (T_p2 + 273.15))*(epsilon_1**-1 + epsilon_2**-1 - 1)**-1
    
def conductive_heat_transfer_coefficient(L, k):
    
    return sum([L[i]/k[i] for i in range(len(k))])**-1

def Rayleigh(T_p1, T_p2, L_gap, beta = 3.4e-3, visc = 1.562e-5, k = 0.024, density = 1.184, c = 1005, g = 9.8):
    
#    alpha = k/(density*c)
    
    return (g*beta*(abs(T_p1 - T_p2))*L_gap**3)/(visc*Thermal_diffusivity(k, density, c))

def Nusselt(T_p1, T_p2, L_gap, tilt = 0):
    # https://watermark-silverchair-com.tudelft.idm.oclc.org/189_1.pdf?token=AQECAHi208BE49Ooan9kkhW_Ercy7Dm3ZL_9Cf3qfKAc485ysgAABG0wggRpBgkqhkiG9w0BBwagggRaMIIEVgIBADCCBE8GCSqGSIb3DQEHATAeBglghkgBZQMEAS4wEQQMmeIjkuilybvmOu6eAgEQgIIEIDb4_PENA7Dbd2jCBVnYUS4X2W1seHzJALTJjlFqGBx5NOs9rkDi9MnddnNFwaz0kftoRS7-16r-HX-jRJHnl_0ugd8v4E_LiwpyeS3MQQlipGe2JkW4SXOHHjRu8y2t-NWC9c1azO1WcGJn74lN1Ei7AYML6z7FwAsyr88yyeO8gUOr9wXzatVyBZtjT_lkdXPWhumDfkeFiPLr1O73iZWkhM2SGCISff4K5t09xsHVpROcPIXUPuj9X6rkhJZ7QenkX6XrpnRs91HEhrUELQrfhknDSL1y5wqBlO1WmP5x_mEyyKkFl_V5ycEvecaZKz_zRNpJGc7M8ajSkNiSJXW1aKnCXG3-vHSE6ErV7MzReqUO-_RcZedWuD11pHtgf-7gdy_rSo0BWWEfA_baqbZlYdeAvzenhhePnjgfA4Bema6lNm0fJvlVfbOr0nDablq1IenqMmbK6whUiWKgBdbzKMWLubjVWWTH3fB3jcRjhxEhjKgHoYcbfZkmbt8bXpgOl4YTZWphgD3goTsfgeN0-4hZxuzuZQIXQDY41pGjXRjUnN8XM3jP3jubPAz8ZPcRlPYIBfbq6FThKNMTirn80E77eEmeraKM57UqfMVVAMzRiUi5Vb-sQWWvI49xprdBk9L7jRL0Xw93MCMloKebfscnRBO7icrnyOhEJShNXp8NT-h4EiVbL-b069LYYdENkO-P14xAUQk3pU7eKkw8_JNeExhHFCRVJn_tVo8audf8hEQA0IYrjxPtxMWLl5QbPKI3f-A2_KlWXmTJDquaSSwTK-EjF4cm1V8_h_9GUMAZe9M9qQcuDFsV81Z8or8reRfpGIvPrJfeM9RqczCeqdVNnBfLECrc8ZCOWg5LeHTgB98KqXbXo7bKklFN1Iji9jZ8UuTDoI2eqSnnVPb9HcPkSXu8rDeNJTEXqxNMQGe9AHDxzJhAeycc5H8idMPVE1NFRTOrZEdlbG356bCW_kZKvwCPE3IrKFrE9v0lzJSO4V5f4uG98W1jqCJU3oHvlmGTnmdk-jQLZnvKLNZac7uwxIBHXWGCHPtqtGnFXe_f_WkJC0-W3FQjWwC3AN_mULEXnVVwMsIzClF58wwpUK-Tv0WSAzRQC94XhZJklXhbAoMew_RcaYHBOz0p5j8fxNWIWifLF-rBopwoddm2_sU-NXNFzT3tXBJZIeSixcgi_FRYVBY3N8eZ7uT8hfF-ckX54SfSb710WJ7ZKgKD7gPXNlSK8KEbymMU0ozWCaVMCK0Q0UW3GrgeGa7lEa6Z4fDa_-GQ13jcnVxMWYHt8M-fvBoQQqES9c4_p6BbSbgd8OK3wqvXCDGCTiqs3xyY5y7mj0HeBHlYHWMkGTlNFno1eAixBbP_K_fwqXj1YULF5Ff4h2o94FFXPjNagw
    # having 0° < tilt < 60°, and 0 < Ra < 105
    from math import sin, cos, radians
    
    return 1 + 1.44*(1 - 1708/(Rayleigh(T_p1, T_p2, L_gap)*cos(radians(tilt))))*(1 - 1708*(sin(radians(1.6*tilt))**1.6)/(Rayleigh(T_p1, T_p2, L_gap)*cos(radians(tilt))))*((Rayleigh(T_p1, T_p2, L_gap)*cos(radians(tilt))/5830)**(1/3)-1)
    
def Reynolds(m_dot, D_tube, density = 1000, visc = 1.562e-5):
    from math import pi
    
    return 4*m_dot/(density*pi*D_tube*visc)

def Prandlt(visc = 0.8927e-6, k = 0.6071, density = 1000, c = 4200):
    return visc/Thermal_diffusivity(k, density, c)

def fluid_conductive_heat_transfer(m_dot, k_f, D_tube):
    
    if m_dot == 0:
        return 2*k_f/D_tube
    elif Reynolds(m_dot, D_tube) < 2300:
        return 4.36*k_f/D_tube
    else: 
        return 0.023*k_f*(Reynolds(m_dot, D_tube)**0.8)*(Prandlt()**0.4)/D_tube
    

def air_gap_convection(T_p1, T_p2, L_gap, k):
    return Nusselt(T_p1, T_p2, L_gap)*k/L_gap
    

def PVT_model(T_amb, G, v, T_f_in, T_glass_0, T_PV_0, T_a_0, T_f_0, n_STC, N_tubes, D_tube, n_HE = 0.8, A_collector = 1.48, A_PV = 1.46, m_f_dot = 0, Len_collector = 1, c_f = 3800, dt = 1):
#    from math import pi
    

    A_glass = A_collector
#    A_gap = A_collector    
    A_PV = A_PV    
    A_a = A_collector
#    A_t_cross = D_tube*Len_collector*N_tubes
#    A_t_surf = 0.5*pi*D_tube*Len_collector*N_tubes
    A_t = 1.48 #A_t_cross + A_t_surf

    rho_glass = 2200
    rho_PV = 2330    
    rho_a = 2699
    rho_f = 1050

    L_glass = 0.004
    L_PV = 0.0004
    L_a = 0.001
    L_gap = 0.02
    L_ins = 0.04
    
    L_PV_glass = 0.003
    L_PV_EVA = 0.0005
    L_PV_tedlar = 0.0001

    k_PV_glass = 1.8
    k_PV_EVA = 0.35
    k_PV_tedlar = 0.2
    k_air = 0.024
    k_f = 0.6071
    k_ins = 0.035
    

    if m_f_dot == 0:
        m_f_dot = rho_f*2.77e-5    # kg/s 
    
    c_glass = 670
    c_pv = 900
    c_a = 800
#    c_f = 4200

    alpha_glass = 0.1
    alpha_PV = 0.9
    tau_glass = 0.93**1
    
    epsilon_glass = 0.9
    epsilon_PV = 0.96
    
    m_glass = rho_glass*A_glass*L_glass
    m_PV = rho_PV*A_PV*L_PV
    m_a = rho_a*A_a*L_a
    m_f = 0.65*rho_f/1000 # N_tubes*rho_f*(0.125*pi*D_tube**2)*Len_collector    

    h_glass_conv = Wind_forced_convection(v)
    h_glass_r = radiative_heat_transfer_coefficient(T_sky(T_amb), T_glass_0, epsilon_glass)    
    h_gap = (conductive_heat_transfer_coefficient([L_PV_glass, L_PV_EVA], [k_PV_glass, k_PV_EVA])**-1 + air_gap_convection(T_glass_0, T_PV_0, L_gap, k_air)**-1)**-1
    h_glassPV_r = radiative_heat_transfer_coefficient_two_plates(T_glass_0, epsilon_glass, T_PV_0, epsilon_PV)
    h_PVa_cond = conductive_heat_transfer_coefficient([L_PV_EVA, L_PV_tedlar], [k_PV_EVA, k_PV_tedlar])
    h_af = fluid_conductive_heat_transfer(m_f_dot, k_f, D_tube)
    h_a_cond = conductive_heat_transfer_coefficient([L_ins],[k_ins])
        

    T_glass = T_glass_0 + dt*A_glass/(m_glass*c_glass)*(h_glass_conv*(T_amb - T_glass_0) + h_glass_r*(T_sky(T_amb) - T_glass_0) + h_gap*(T_PV_0 - T_glass_0) + h_glassPV_r*(T_PV_0 - T_glass_0) + alpha_glass*G)
    T_PV = T_PV_0 + (A_PV*dt/(m_PV*c_pv))*(h_gap*(T_glass_0 - T_PV_0) + h_glassPV_r*(T_glass_0 - T_PV_0) + h_PVa_cond*(T_a_0 - T_PV_0) + alpha_PV*tau_glass*G*(1 - PV_efficiency(T_PV_0, n_STC)))
    T_a = T_a_0 + (dt/(m_a*c_a))*(A_a*h_PVa_cond*(T_PV_0 - T_a_0) + h_af*A_t*(T_f_0 - T_a_0) + h_a_cond*A_a*(T_amb - T_a_0))
    T_f = T_f_0 + (dt/(m_f*c_f))*(h_af*A_t*(T_a_0 - T_f_0) + 2*m_f_dot*c_f*(T_f_in - T_f_0)) 

#    h_glass_conv_registry.append(h_glass_conv)
#    h_glass_r_registry.append(h_glass_r)
#    h_gap_registry.append(h_gap)
#    h_glassPV_r_registry.append(h_glassPV_r)
#    h_PVa_cond_registry.append(h_PVa_cond)
#    h_af_registry.append(h_af)
#    h_a_cond_registry.append(h_a_cond)
#    glass_absorb_registry.append(alpha_glass*G)
#    PV_absorb_registry.append(alpha_PV*tau_glass*G*(1 - PV_efficiency(T_PV_0, n_STC)))
#    massflow_absorb_registry.append(2*m_f_dot*c_f)
    
#    h_glass_conv_registry.append(A_glass*h_glass_conv*(T_amb - T_glass_0))
#    h_glass_r_registry.append(A_glass*h_glass_r*(T_sky(T_amb) - T_glass_0))
#    h_gap_registry.append(A_glass*h_gap*(T_PV_0 - T_glass_0))
#    h_glassPV_r_registry.append(A_glass*h_glassPV_r*(T_PV_0 - T_glass_0))
#    h_PVa_cond_registry.append(A_a*h_PVa_cond*(T_a_0 - T_PV_0))
#    h_af_registry.append(h_af*A_t*(T_f_0 - T_a_0))
#    h_a_cond_registry.append(h_a_cond*A_a*(T_amb - T_a_0))
#    glass_absorb_registry.append(A_glass*alpha_glass*G)
#    PV_absorb_registry.append(A_glass*alpha_PV*tau_glass*G*(1 - PV_efficiency(T_PV_0, n_STC)))
#    massflow_absorb_registry.append(2*m_f_dot*c_f*(T_f_in - T_f_0))    
    
    Q_PVT_dot = n_HE*2*m_f_dot*c_f*(T_f - T_f_in)
    
    return [Q_PVT_dot, T_glass, T_PV, T_a, T_f]


def Qdot_PVT(PVT_active, TESS_charge, T_amb, G, v, T_in_Network, T_glass_0, T_PV_0, T_a_0, T_f_0, T_tank_PVT, T0_TESS, n_STC = 0.184, N_tubes = 0, D_tube = 0.009, n_HE = 0.8, A_collector = 1.48, A_PV = 1.46, m_f_dot = 0, Len_collector = 1, t_end = int(0.25*3600), n_modules = 1, recirculation = True):

#    T_tank_PVT_registry = []
#    Q_PVT_dot_registry = []
#    Qdot_PVT2Network_registry = []
#    T_network_registry = []
#    import matplotlib.pyplot as plt
    
    if recirculation:

        for step in range(t_end):
        
            [Q_PVT_dot, T_glass_0, T_PV_0, T_a_0, T_f_0] = PVT_model(T_amb, G, v, T_tank_PVT, T_glass_0, T_PV_0, T_a_0, T_f_0, n_STC, n_modules*N_tubes, D_tube, n_HE, n_modules*A_collector, n_modules*A_PV, m_f_dot, n_modules*Len_collector)                     
            [T_out_Network, T_tank_PVT, Qdot_PVT2Network, Qdot_PVT2TESS] = Update_PVT_tank(T_fluid_out(T_tank_PVT, T_f_0), T_in_Network, T_tank_PVT, T0_TESS, PVT_active, TESS_charge)         
            
            
#            T_tank_PVT_registry.append(T_tank_PVT)
#            Q_PVT_dot_registry.append(Q_PVT_dot)
#            Qdot_PVT2Network_registry.append(Qdot_PVT2Network)
#            T_network_registry.append(T_out_Network)
            
            

#        plt.figure(1)
#        plt.plot(T_tank_PVT_registry)
#        plt.title('T_tank_PVT_registry')
#        plt.figure(2)
#        plt.plot(Q_PVT_dot_registry)
#        plt.title('Q_PVT_dot_registry')
#        plt.figure(3)
#        plt.plot(Qdot_PVT2Network_registry)
#        plt.title('Qdot_PVT2Network_registry')
#        plt.figure(4)
#        plt.plot(T_network_registry)
#        plt.title('T_network')        
        
        return [Q_PVT_dot, [T_glass_0, T_PV_0, T_a_0, T_f_0], T_tank_PVT, T_out_Network, Qdot_PVT2Network, Qdot_PVT2TESS]            
    else:
        
        return False
    
        
#    return [Q_PVT_dot, T_glass_0, T_PV_0, T_a_0, T_f_0]

#    return Q_PVT_dot


def T_fluid_out(T_in, T_f):
    return (T_f - 0.5*T_in)/0.5

def PV_efficiency(T_PV, n_STC, beta_PV = -0.003, T_ref = 25):
    return n_STC*(1 - beta_PV*(T_PV - T_ref))

def PV_out (T_PV, G, A_PV = 60*0.156*0.156, n_STC = 0.184, kW = True):
    if kW:
        return PV_efficiency(T_PV, n_STC)*G*A_PV/1000
    
    else:
        return PV_efficiency(T_PV, n_STC)*G*A_PV

def Update_PVT_tank(T_in_PVT, T_in_Network, T0, T0_TESS, PVT_active = True, TESS_charge = False, mdot = 0.029085, m = 200, c = 4200, efficiency_TESS = 0.8, dt = 1):
    # Assuming mass conservation, using energy balance equatiion.
    
    # Thermal power from the PVT.
    [T_out, Qdot_PVT] = Heat_exchanger(T_in_PVT, T0, 0.029085, efficiency = 1, c1 = 3800, mode = "TESS")
    
    # Thermal power supplied to the thermal network.
    if PVT_active == True:
        [T_out_Network, T_in_tank, Qdot_Network, Qdot_tank] = Heat_exchanger(T_in_Network, T0, 0.22, 0.22, mode = "L2L")            
    else:
        T_out_Network = T_in_Network
        Qdot_Network = 0
        Qdot_tank = 0
        
    if TESS_charge == True:
        [T_out_TESS, Qdot_TESS] = Heat_exchanger(T0, T0_TESS, efficiency = efficiency_TESS, mode = "TESS")
        
    else:        
        Qdot_TESS = 0
        
    T_tank = T0 + (dt/(m*c))*(Qdot_PVT - Qdot_tank - Qdot_TESS/efficiency_TESS)
    
    return [T_out_Network, T_tank, Qdot_Network, Qdot_TESS] 
    
    
def Heat_exchanger(T1_in, T2_in, mdot1 = 0.22, mdot2 = 0.029085, c1 = 4200, c2 = 4200, efficiency = 0.8, mode = "L2L"):
# Assuming a small difference between T1_in and T2_in, a large contact surface and T2_in > T1_in.
    


    if mode == "L2L":
    #    T2_out = T1_in
    #    Qdot2 = mdot2*c2*(T2_out - T2_in)
    #    Qdot1 = -efficiency*Qdot2
    #    dT = Qdot1/(mdot1*c1)
    #    T1_out = T1_in + dT        
        Qdot1 = efficiency*mdot2*c2*(T2_in - T1_in)
        T1_out = T1_in + Qdot1/(mdot1*c1)
        
        return [T1_out, T1_in, Qdot1, Qdot1/efficiency] # T1_out, T2_out, Qdot1, Qdot2

    elif mode == "TESS":
    #    T1_out = T2_in
    #    Qdot2 = efficiency*mdot1*c1*(T1_in - T1_out)
        Qdot2 = efficiency*mdot1*c1*(T1_in - T2_in)
    
        return [T2_in, Qdot2]  # T1_out, Qdot2

    
#    else:
#        # Put something in here

    

###############################################################################
########################   Simple SC model definition  ########################

def Qdot_SolarCollector(active, G, A = 6, SC_eff = 0.45, dt = 1*3600):
    
    if active:
        return A*SC_eff*G/dt
    else:
        return 0


###############################################################################
#######################   Simple TESS model definition  #######################

def update_TESS(active, T_in, T_0, T_soil, Qdot_SC = 0, Qdot_HP = 0, Tmax = 95 + 273, Tmin = 50 + 273, mdot = 0.1, T_network = 50 + 273, efficiency = 0.8, m = 4000, c = 4200, dt=0.25*3600): # , Qdot_SD = 100

       
    Qdot_SD = Qdot_SD_TESS2Soil(T_soil, T_0 - 273, m)
    
    if active and T_0 >= Tmin:
        Qdot_TESS = mdot*c*(T_network - T_in)
    else:
        Qdot_TESS = 0

    
    if T_0 <= Tmin:         # TESS discharged, only charge
        T_new = T_0 + (Qdot_SC*efficiency + Qdot_HP*efficiency - Qdot_SD - Qdot_TESS/efficiency)*dt/(m*c)
#        if T_0 <= T_soil:   # Code for heat dissipation/absoprtion through the soil needed
#            T_new = T_0
#            
#        else:
#            T_new = T_0 + (Qdot_SC*efficiency + Qdot_HP*efficiency - Qdot_SD - Qdot_TESS/efficiency)*dt/(m*c)
            
    elif T_0 <= Tmax:       # TESS available for charge and discharge
        
        T_new = T_0 + (Qdot_SC*efficiency + Qdot_HP*efficiency - Qdot_SD - Qdot_TESS/efficiency)*dt/(m*c)
        
    else:                   # TESS fully charged, only discharge
        
        T_new = T_0 + (- Qdot_SD - Qdot_TESS/efficiency)*dt/(m*c)        
        
    return [T_new, Qdot_TESS, Qdot_SD]     

def Qdot_SD_TESS2Soil(T_soil, T_TESS, m = 4000, dy = 0.1, l_TESS = 3.8, w_TESS = 2.8, h_TESS = 1.6, t_TESS = 0.4, k_styro = 0.033, h_styro = 0): # t_TESS = 0.314
    
    m0 = 4000
    
    volume_factor = m/m0
    
    U = (t_TESS/k_styro)**(-1)      # 1/h_styro + 
    dA = 2*(l_TESS + w_TESS)*dy*(volume_factor**(0.5))
    
    Qdot_top = U*(l_TESS*w_TESS*volume_factor)*(T_TESS - T_soil[0])
    Qdot_sides = 0
    
    for T_y in T_soil[int(t_TESS/dy):int((h_TESS - t_TESS)/dy) + 1]:
        Qdot_sides += U*dA*(T_TESS - T_y)
    
    Qdot_bottom = U*(l_TESS*w_TESS*volume_factor)*(T_TESS - T_soil[-1])
    
    return (Qdot_top + Qdot_sides + Qdot_bottom)
    

###############################################################################
########################   Simple HP model definition  ########################
    
def HP_Power_0(active, P_in = 2.7*1000, COP = 4.1):
    if active:
        return [P_in/1000, P_in*COP]
    else:
        return [0,0]

###############################################################################
#########################   Nikos HP model definition  ########################
    
global COP_registry    
COP_registry = []
    
def HP_Power(active, T_amb, T_ret, m_dot = 0.22, HE_eff = 0.8, T_sup = 273 + 53, c = 4200):
    from math import e
    
    if active:
        
        COP = 7.90471*e**(-0.024*(T_ret - T_amb))
        COP_registry.append(COP)
#        print(COP)
#        Q_dot = HE_eff*m_dot*c*(T_sup - T_ret)
        
        return [m_dot*c*(T_sup - T_ret)/(1000*COP), HE_eff*m_dot*c*(T_sup - T_ret)]
    else:
        return [0,0]

###############################################################################
#######################   Simple BESS model definition  #######################


## Enphase IQ3 https://enphase.com/download/iq-battery-3-data-sheet
#
#def update_BESS(SoC_0, P_BESS, P_Load, P_PV = 0, SoCmax = 0.9, SoCmin = 0.2, P_BESS_max = 1.28, P_Grid_max = 0, Capacity_BESS = 3.36, charge_efficiency = 0.943, discharge_efficiency = 0.943, P_SD = 0, dt = 0.25):
#
#    
#    E_BESS_0 = Capacity_BESS*SoC_0
#    
#    if P_BESS > 0:
#        E_BESS = E_BESS_0 - P_BESS*dt/discharge_efficiency
#    elif P_BESS <= 0:
#        E_BESS = E_BESS_0 - P_BESS*dt*charge_efficiency
#    
#    E_BESS = E_BESS*(1-P_SD)
#    SoC_BESS = E_BESS /  Capacity_BESS
#    
#    return SoC_BESS
#    
#
#def BESS_perm_min(SoC, Capacity_BESS = 3.36, SoCmax = 0.9, P_BESS_max = 1.28, dt = 0.25):
#    from numpy import clip
#    
#    return clip(Capacity_BESS*(SoC - SoCmax)/dt, -P_BESS_max, P_BESS_max)
#    
#    
#def BESS_perm_max(SoC, Capacity_BESS = 3.36, SoCmin = 0.2, P_BESS_max = 1.28, dt = 0.25):
#    from numpy import clip
#    
#    return clip(Capacity_BESS*(SoC - SoCmin)/dt, -P_BESS_max, P_BESS_max)
        
###############################################################################
        
def BESS_perm_min(SoC, Capacity_BESS = 3.36, SoCmax = 0.9, P_BESS_max = 1.28, dt = 0.25):
    from numpy import clip
    
    return clip(Capacity_BESS*(SoC - SoCmax)/dt, -P_BESS_max, P_BESS_max)
    
    
def BESS_perm_max(SoC, Capacity_BESS = 3.36, SoCmin = 0.2, P_BESS_max = 1.28, dt = 0.25):
    from numpy import clip
    
    return clip(Capacity_BESS*(SoC - SoCmin)/dt, -P_BESS_max, P_BESS_max)

# Enphase IQ3 https://enphase.com/download/iq-battery-3-data-sheet
def update_BESS(SoC_0, P_Load, P_PV, SoCmax = 0.9, SoCmin = 0.2, P_BESS_max = 1.28, P_Grid_max = 0, Capacity_BESS = 3.36, charge_efficiency = 0.943, discharge_efficiency = 0.943, P_SD = 0, dt = 0.25):
    
    
    E_BESS_0 = Capacity_BESS*SoC_0
    
    if SoC_0 >= SoCmax:                  # Battery can only discharge

        if (P_Load-P_PV) <= P_Grid_max:        # No peakshaving needed
            P_BESS = 0
            E_BESS = E_BESS_0*(1-P_SD)
            SoC_BESS = E_BESS/Capacity_BESS
            P_Grid = P_Load - P_PV 
            
#            state = 1

        
        else:                                                       # Peakshaving needed
            if P_Load - P_PV - P_Grid_max  <= BESS_perm_max(SoC_0):    # Below the BESS max power
                P_BESS = P_Load - P_PV - P_Grid_max
                E_BESS = E_BESS_0 - P_BESS*dt/discharge_efficiency
                E_BESS = E_BESS*(1-P_SD)
                SoC_BESS = E_BESS/Capacity_BESS
                P_Grid = P_Grid_max
                
#                state = 2
                
            else:                                                       # Above the BESS max power
                P_BESS = BESS_perm_max(SoC_0)
                E_BESS = E_BESS_0 - P_BESS*dt/discharge_efficiency
                E_BESS = E_BESS*(1-P_SD)
                SoC_BESS = E_BESS/Capacity_BESS
                P_Grid = -BESS_perm_max(SoC_0) - P_PV + P_Load     
                
#                state = 3
    
    elif SoC_0 > SoCmin:                  # Battery can charge and discharge
        
        if (P_Load-P_PV) <= P_Grid_max:        # No peakshaving needed
            
            if P_Load >= P_PV:                    # PV below demanded load
                P_BESS = 0
                E_BESS = E_BESS_0*(1-P_SD)
                SoC_BESS = E_BESS/Capacity_BESS
                P_Grid = P_Load - P_PV     

#                state = 4

            else:                                                   # Surplus of PV power
                if P_PV - P_Load <= -BESS_perm_min(SoC_0):    # Below the BESS max power
                    P_BESS = P_Load - P_PV
                    E_BESS = E_BESS_0 - P_BESS*dt*charge_efficiency
                    E_BESS = E_BESS*(1-P_SD)
                    SoC_BESS = E_BESS/Capacity_BESS
                    P_Grid = 0
                    
#                    state = 5

                    
                else:                                                       # Above the BESS max power
                    P_BESS = BESS_perm_min(SoC_0)
                    E_BESS = E_BESS_0 - P_BESS*dt*charge_efficiency
                    E_BESS = E_BESS*(1-P_SD)
                    SoC_BESS = E_BESS /  Capacity_BESS
                    P_Grid = -BESS_perm_min(SoC_0) - P_PV + P_Load
                    
#                    state = 6

                    
        else:                                                       # Peakshaving needed
            if P_Load - P_PV - P_Grid_max <= BESS_perm_max(SoC_0):    # Below the BESS max power
                P_BESS = P_Load - P_PV - P_Grid_max
                E_BESS = E_BESS_0 - P_BESS*dt/discharge_efficiency
                E_BESS = E_BESS*(1-P_SD)
                SoC_BESS = E_BESS /  Capacity_BESS
                P_Grid = P_Grid_max
                
#                state = 7

                
            else:                                                       # Above the BESS max power
                P_BESS = BESS_perm_max(SoC_0)
                E_BESS = E_BESS_0 - P_BESS*dt/discharge_efficiency
                E_BESS = E_BESS*(1-P_SD)
                SoC_BESS = E_BESS /  Capacity_BESS
                P_Grid = -BESS_perm_max(SoC_0) - P_PV + P_Load  
                
#                state = 8

    
    else: # self.BESS.SoC[1,i] <= SoCmin:                 # Battery can only charge     
        
        if P_Load >= P_PV:                        # PV below demanded load
            P_BESS = 0
            E_BESS = E_BESS_0*(1-P_SD)
            SoC_BESS = E_BESS/Capacity_BESS
            P_Grid = P_Load - P_PV 
            
#            state = 9
             
        else:                                                       # Surplus of PV power
            
            if P_PV - P_Load <= -BESS_perm_min(SoC_0):    # Below the BESS max power
                P_BESS = P_Load - P_PV
                E_BESS = E_BESS_0 - P_BESS*dt*charge_efficiency
                E_BESS = E_BESS*(1-P_SD)
                SoC_BESS = E_BESS /  Capacity_BESS
                P_Grid = 0
                
#                state = 10

                
                
            else:                                                       # Above the BESS max power
                P_BESS = BESS_perm_min(SoC_0)
                E_BESS = E_BESS_0 - P_BESS*dt*charge_efficiency
                E_BESS = E_BESS*(1-P_SD)
                SoC_BESS = E_BESS /  Capacity_BESS
                P_Grid = -BESS_perm_min(SoC_0) - P_PV + P_Load      
                
#                state = 11
                
                
    return [SoC_BESS, P_Grid, P_BESS]    

###############################################################################
##########################   House model definition  ##########################

def heat_loss(U, Area, T_in, T_out):
    return U * Area * (T_out - T_in)
    
def new_house_Temperature(T_0, Qdot_Losses, mc_T, Qdot_HP = 0, Qdot_TESS = 0, Qdot_SC = 0, dt=0.25*3600):
    T_new = T_0 + dt*(Qdot_Losses + Qdot_HP + Qdot_TESS + Qdot_SC)/(mc_T)
    
    return T_new

def R_convective(m_dot, k_f = 0.6071, D_tube = 13e-3, L_tube = 1):
    
    return 1/(pi*D_tube*L_tube*fluid_conductive_heat_transfer(m_dot, k_f, D_tube))

def R_conductive(D_in, D_out, k_tube, L_tube = 1):
    from math import log
    
    return log(D_out/D_in)/(2*pi*L_tube*k_tube)

def Celcius2Farenheit(T_C):
    return T_C*9/5 + 32

def Infiltration_Airflow(T_in, T_out, U = 1, A_exposed = 119.6, C_stack = 0.015, C_wind = 0.0012):
    from math import sqrt
    
    U *= 2.23694    # m/s to mph
    A_exposed *= 10.7639    # m3 to ft3

    dT = Celcius2Farenheit(T_in) - Celcius2Farenheit(T_out)          # °C to °F
    
    A_unit = 0.01   # 
    
    if dT > 0:    
        return A_exposed*A_unit*sqrt(C_stack*abs(dT) + C_wind*U**2)*0.47194745/1000

    else:
        return 0

def Ventilation_Airflow(A_conditioned = 120, N_bedrooms = 1):
    A_conditioned *= 10.7639    # m3 to ft3
    return (0.03*A_conditioned + 7.5*(1+N_bedrooms))*0.47194745/1000

def Infiltration_Ventilation_Heat_loss(T_in, T_out, U = 1, A_exposed = 119.6, C_stack = 0.015, C_wind = 0.0012, A_conditioned = 120, N_bedrooms = 1, c_air = 1012, density_air = 1.293):
    
    return -c_air*density_air*(Infiltration_Airflow(T_in, T_out, U, A_exposed, C_stack, C_wind) + Ventilation_Airflow(A_conditioned, N_bedrooms))*(T_in - T_out)
    

def pipe_losses(T_f, T_out, m_dot, diameters = [13, 15, 40], k_f = 0.6071, k_tubes = [398, 0.038], L_tube = 1):
    
    R_total = R_convective(m_dot, k_f, diameters[0], L_tube = 1)
    
    
    for i in range(len(k_tubes)):
        
        R_total += R_conductive(diameters[i], diameters[i+1], k_tubes[i], L_tube)
    
    return (T_f - T_out) / R_total

def House_Thermal_Losses(T_0_in, T_amb, U = 1, A_exposed = 119.6, C_stack = 0.015, C_wind = 0.0012, A_conditioned = 120, N_bedrooms = 1, c_air = 1012, density_air = 1.293):

    #############################   Paremeters ################################
    
    # Convective heat transfer coefficients [W/m^2K]
    h_air_wall   = 0.9#24       # Indoor air -> walls, scaled to R-value of a C-label house
    h_wall_atm   = 0.9#34       # Walls -> atmosphere, scaled to R-value of a C-label house
    h_air_window = 25            # Indoor air -> windows
    h_window_atm = 32            # Windows -> atmosphere
    h_air_roof   = 12            # Indoor air -> roof
    h_roof_atm   = 38            # Roof -> atmosphere
    
    ## House
    
    # Air
    c_air        = 1005.4        # Specific heat of air at 273 K [J/kgK]
    airDensity   = 1.025         # Densiity of air at 293 K [kg/m^3]
    kAir         = 0.0257        # Thermal conductivity of air at 293 K [W/mK]
    
    
    # Windows (glass)
    n1_window     = 3            # Number of windows in room 1
    n2_window     = 2            # Number of windows in room 2
    n3_window     = 2            # Number of windows in room 3
    n4_window     = 1            # Number of windows in room 4
    htWindows     = 1            # Height of windows [m]
    widWindows    = 1            # Width of windows [m]
    windows_area  = (n1_window + n2_window + n3_window + n4_window) * htWindows * widWindows
    LWindow       = 0.004        # Thickness of a single window pane [m]
    LCavity       = 0.014        # Thickness of the cavity between the double glass window [m]  
    windowDensity = 2500         # Density of glass [kg/m^3]
    c_window      = 840          # Specific heat of glass [J/kgK]
    kWindow       = 0.8          # Thermal conductivity of glass [W/mK]
    U_windows = ((1/h_air_window) + (LWindow/kWindow) + (LCavity/kAir) + (LWindow/kWindow) + (1/h_window_atm))**-1
    m_windows = windowDensity * windows_area * LWindow
    
    # Walls (concrete)
    lenHouse    = 15             # House length [m]
    widHouse    = 8              # House width [m]
    htHouse     = 2.6            # House height [m]
    LWall       = 0.25           # Wall thickness [m]
    wallDensity = 2400           # Density [kg/m^3]
    c_wall      = 750            # Specific heat [J/kgK]
    kWall       = 0.14           # Thermal conductivity [W/mK]
    walls_area = 2*(lenHouse + widHouse) * htHouse - windows_area
    U_wall = ((1/h_air_wall) + (LWall /kWall) + (1/h_wall_atm))**-1
    m_walls = wallDensity * walls_area * LWall
    
    # Roof (glass fiber)
    pitRoof     = 40/180/pi      # Roof pitch (40 deg)
    LRoof       = 0.2            # Roof thickness [m]
    roofDensity = 2440           # Density of glass fiber [kg/m^3]
    c_roof      = 835            # Specific heat of glass fiber [J/kgK]
    kRoof       = 0.04           # Thermal conductivity of glass fiber [W/mK]
    roof_Area = 2 * (widHouse/(2*cos(pitRoof))*lenHouse)
    U_roof = ((1/h_air_roof) + (LRoof/kRoof) + (1/h_roof_atm))**-1
    m_roof = roofDensity * roof_Area * LRoof
    
    
    m_air = airDensity * lenHouse * widHouse * htHouse
    
    mc_T = m_air*c_air + m_roof*c_roof + m_windows*c_window + m_walls*c_wall
     
    
    ############################   Calculations ###############################    
    ################### Thermal carrier ######################
    
    ## Thermal losses
    # Roof
    Qdot_roof = heat_loss(U_roof, roof_Area, T_0_in, T_amb)
    
    # Windows
    Qdot_windows = heat_loss(U_windows, windows_area, T_0_in, T_amb)
    
    # Walls
    Qdot_wall = heat_loss(U_wall, walls_area, T_0_in, T_amb)
    
    # Infiltration and ventilation
    Qdot_iv = Infiltration_Ventilation_Heat_loss(T_0_in - 273, T_amb - 273, U, (walls_area + windows_area), C_stack, C_wind, roof_Area, N_bedrooms, c_air, airDensity)    
    
    
    return [Qdot_roof + Qdot_windows + Qdot_wall + Qdot_iv, mc_T]


#global enable_HP_2_TESS
#enable_HP_2_TESS = [True]

def Thermal_Electrical_model(Thermal_Components, T_in_0, T4, T_TESS_0, T_tank_PVT, T_PVT_layers, TESS_min_op_T, TESS_max_op_T, SoC_0_BESS, T_amb, T_set, G, Tsoil, P_Load, t_final, n_modules = 0, P_BESS_max = 1.28, P_Grid_max = 0, Capacity_BESS = 3.36, PVT_2_TESS = True, HP_2_TESS = True, Network_losses = False, m = 4000, T_sup = 50 + 273, T_TESS_max = 273 + 95, m_dot = 0.22, c_f = 4200, A_PVT = 1.48, dt = 60*15): 
    
    
    TESS_active = True
                
    T_in = [T_in_0, T_in_0]
    T_TESS = [T_TESS_0]
    
    
    fluid_registry = [[], [], [], [], [T4 - 273]]
    T_registry = [[], [], [], [], []]
    T_tank_PVT_registry = []
    Qdot_PVT_registry = []
    Qdot_SC_TESS = []
    Qdot_PVT_Network = []
    P_PVT = []
    A_PVT = A_PVT*n_modules 
    
    Qdot_TESS = []  
    Qdot_TESS_SD = []
    enable_HP_2_TESS = True
    
    P_HP = []
    Qdot_HP = []
    Qdot_HP_TESS = []
    
    Qdot_Losses = []
    
    SoC_BESS = [SoC_0_BESS]  
    P_BESS = []
    
    P_Grid = []  
    

    for step in range(t_final): 
    
    
    ########## Thermal    
        TESS_HP_control = False
    
        [Qdot_Losses_i, mc_T] = House_Thermal_Losses(T_in[-1], T_amb[step])     
        Qdot_Losses.append(Qdot_Losses_i)  
        Q_dot_network = -(mc_T*(T_in[-1] - T_in[-2])/dt - Qdot_Losses_i) #  
        T0 = T4 + Q_dot_network/(m_dot*c_f)    
        T_registry[0].append(T0)   
        
        if Network_losses:
            Q_dot_network += -pipe_losses(T0, T_in[-1], m_dot)
        else:
            Q_dot_network += 0
        
        
        T1 = T4 + Q_dot_network/(m_dot*c_f)
        T_registry[1].append(T1)
            
        # PVT
        if Thermal_Components[0]:    
            if T_in[-1] < T_set[step] and T1 < T_sup and T_tank_PVT > T1 - 273:
                PVT_active = True
                
                
            else:
                PVT_active = False
        
#        PVT_2_TESS = True, HP_2_TESS = True
        
            if T_tank_PVT > T_TESS[-1] - 273 + 2 and Thermal_Components[1] == True and PVT_2_TESS:            
                PVT_TESS_charge = True
        
            else:
                PVT_TESS_charge = False
        
                
            [Q_dot_SC, T_PVT_layers, T_tank_PVT, T2, Qdot_PVT2Network, Qdot_PVT2TESS] = Qdot_PVT(PVT_active, PVT_TESS_charge, int(T_amb[step] - 273), G[step]/3600, 1, T1-273, T_PVT_layers[0], T_PVT_layers[1], T_PVT_layers[2], T_PVT_layers[3], T_tank_PVT, T_TESS[-1]-273, n_HE = 1, m_f_dot = 0.029085, t_end = int(15*60), n_modules = n_modules)
            T2 += 273     
            
            Qdot_SC_TESS.append(Qdot_PVT2TESS)     
            Qdot_PVT_registry.append(Q_dot_SC)   
            Qdot_PVT_Network.append(Qdot_PVT2Network)
            T_tank_PVT_registry.append(T_tank_PVT)
    
            fluid_registry[0].append(T_PVT_layers[0])
            fluid_registry[1].append(T_PVT_layers[1])
            fluid_registry[2].append(T_PVT_layers[2])
            fluid_registry[3].append(T_PVT_layers[3])
            fluid_registry[4].append(T2 - 273)
            
            P_PVT.append(PV_out (T_PVT_layers[1], G[step]/3600, A_PVT))
    
        else:
            T2 = T1
            Qdot_SC_TESS.append(0)     
            Qdot_PVT_registry.append(0)   
            Qdot_PVT_Network.append(0)
            T_tank_PVT_registry.append(0) 
    
            fluid_registry[0].append(0)
            fluid_registry[1].append(0)
            fluid_registry[2].append(0)
            fluid_registry[3].append(0)
            fluid_registry[4].append(0)      
            
            P_PVT.append(0)
    
        T_registry[2].append(T2)
        
        # TESS
        if Thermal_Components[1]:

        # Control the minimum TESS SoC during operation            
#            if T_TESS[-1] < 90 + 273:
            if T_TESS[-1] < TESS_min_op_T:                
                enable_HP_2_TESS = True
            elif T_TESS[-1] > TESS_max_op_T:
                enable_HP_2_TESS = False
            
            
            if T_in[-1] < T_set[step] and T2 < T_sup and T_TESS[-1] > T2:  
                TESS_active = True
#                print('TESS - Case 1')                
            
            elif T_in[-1] > T_set[step] and T_TESS[-1] < T_TESS_max - 5:
                TESS_active = False
#                print('TESS - Case 2')  
                
                if HP_2_TESS and enable_HP_2_TESS:
#                    print('TESS - Case 3')
                    
                    TESS_HP_control = True
                    HP_state = HP_Power(HP_2_TESS, T_amb[step], T_TESS[-1], m_dot, T_sup = T_TESS_max)                    
                    P_HP.append(HP_state[0])
                    Qdot_HP.append(0)
                    Qdot_HP_TESS.append(HP_state[1])
                
                else:
#                    print('TESS - Case 4')
                    TESS_HP_control = True
                    P_HP.append(0)
                    Qdot_HP.append(0)        
                    Qdot_HP_TESS.append(0)                    
            
            else:
#                print('TESS - Case 5')  
                TESS_active = False
#                TESS_HP_control = True
#                P_HP.append(0)
#                Qdot_HP.append(0)
#                Qdot_HP_TESS.append(0)                
                
                
            TESS_state = update_TESS(TESS_active, T2, T_TESS[-1], T_soil = Tsoil[step], Qdot_SC = Qdot_SC_TESS[-1], Qdot_HP = Qdot_HP_TESS[-1], mdot = m_dot, m = m)
            T_TESS.append(TESS_state[0])
            Qdot_TESS.append(TESS_state[1])
            Qdot_TESS_SD.append(TESS_state[2])
            T3 = T4 + (Q_dot_network + Qdot_PVT_Network[-1] + Qdot_TESS[-1])/(m_dot*c_f) 
        
        else:
#            print('TESS - Case 6')  
            T_TESS.append(0)
            Qdot_TESS.append(0)
            T3 = T2
        
        T_registry[3].append(T3)       
        
        # HP  
        if Thermal_Components[2]:  
            if T_in[-1] < T_set[step] and T3 < T_sup: 
#                print('HP - Case 1')
                HP_active = True
                
                HP_state = HP_Power(HP_active, T_amb[step], T3, m_dot)
                P_HP.append(HP_state[0])
                Qdot_HP.append(HP_state[1])   
                Qdot_HP_TESS.append(0)
                
                
#                T4 += (Q_dot_network + Qdot_PVT_Network[-1] + Qdot_TESS[-1] + Qdot_HP[-1])/(m_dot*c_f)  

            elif not(Thermal_Components[1]) or not(TESS_HP_control):
#                print('HP - Case 2')                
                P_HP.append(0)
                Qdot_HP.append(0)
                Qdot_HP_TESS.append(0)
                
#                T4 = T3
                                 
            
        else:
#            print('HP - Case 3')
            P_HP.append(0)
            Qdot_HP.append(0)
            Qdot_HP_TESS.append(0)
#            T4 = T3
    
        
        T4 += (Q_dot_network + Qdot_PVT_Network[-1] + Qdot_TESS[-1] + Qdot_HP[-1])/(m_dot*c_f)  
        T_registry[4].append(T4)   
    

        
        T_in.append(new_house_Temperature(T_in[-1], Qdot_Losses_i, mc_T, Qdot_HP = Qdot_HP[-1], Qdot_TESS = Qdot_TESS[-1], Qdot_SC = Qdot_PVT_Network[-1]))
     
    ########## Electric
    
        # PVT above
        
        
        
        # BESS
        [SoC_BESS_state, P_Grid_state, P_BESS_state] = update_BESS(SoC_BESS[-1], P_Load[step] + P_HP[-1], P_PVT[-1], P_BESS_max = P_BESS_max, P_Grid_max = P_Grid_max, Capacity_BESS = Capacity_BESS)
        P_BESS.append(P_BESS_state)    
        SoC_BESS.append(SoC_BESS_state)
        P_Grid.append(P_Grid_state)      
        
    
    return [T_in, Qdot_PVT_Network, Qdot_SC_TESS, Qdot_TESS, Qdot_HP, Qdot_HP_TESS, T_registry, Qdot_Losses, T_TESS, Qdot_TESS_SD, fluid_registry, T_tank_PVT_registry, P_Grid, P_PVT, P_HP, P_BESS, SoC_BESS]
    
#def Thermal_Electrical_model(T_amb, T_0_in, T_set, T_soil, P_PV_av, T4_0, P_Load, G, T_0_TESS, SoC_0_BESS, P_BESS, T_sup = 50, external_control = False, SC_active = False, TESS_active = False, HP_active = False, PV_curtail = 1, m_dot = 0.22, c_f = 4200, dt = 0.25):
#    
#    [Qdot_Losses, mc_T] = House_Thermal_Losses(T_0_in, T_amb)
#
#    if external_control:
#        
#        print('put some code here for the external control')
#
#
#    else:   # Rule-Based control
#        
#        # PVT
#        T1 = T4_0 + Qdot_Losses/(m_dot*c_f)
#        if T1 < T_sup:
#            SC_active = True
#
#        Qdot_PVT = Qdot_SolarCollector(SC_active, G)
#
#
#        # TESS     
#        T2 = T4_0 + (Qdot_Losses + Qdot_PVT)/(m_dot*c_f)         
#        if T2 < T_sup:
#            TESS_active = True
#            
#        Qdot_SC_TESS = Qdot_SolarCollector(not SC_active, G)
#        TESS_state = update_TESS(TESS_active, T2, T_0_TESS, T_soil = T_soil, Qdot_SC = Qdot_SC_TESS, dt = dt)
#        T_TESS = TESS_state[0]
#        Qdot_TESS =  TESS_state[1]
#
#            
#        # HP   
#        T3 = T4_0 + (Qdot_Losses + Qdot_PVT + Qdot_TESS)/(m_dot*c_f)  
#            
#        if T3 < T_sup:
#            HP_active = True
#
#        HP_state = HP_Power(HP_active, T_amb, T3)
#        P_HP = HP_state[0]
#        Qdot_HP = HP_state[1]   
#        
#        T4 = T4_0 + (Qdot_Losses + Qdot_PVT + Qdot_TESS + Qdot_HP)/(m_dot*c_f)        
#    
#    # Thermal demand
#    T_in = new_house_Temperature(T_0_in, Qdot_Losses, mc_T, Qdot_HP = Qdot_HP, Qdot_TESS = Qdot_TESS, Qdot_SC = Qdot_PVT, dt = dt)  
#
#
#
#    ################### Electric carrier ######################
#
#    # PV
#    P_PV = P_PV_av*PV_curtail
#
#    # BESS
##   [SoC_BESS_state, P_Grid_state, P_BESS_state] = update_BESS(SoC_BESS, P_Load + P_HP/1000, P_PV)
#    SoC_BESS = update_BESS(SoC_0_BESS, P_BESS, P_Load + P_HP)
#    
#    
#    # Grid balance
#    
#    P_Grid = P_Load + P_HP - P_PV - P_BESS
#    
#    return [T_in, T4]
##    return [Qdot_PVT, Qdot_SC_TESS, Qdot_TESS, Qdot_HP, Qdot_Losses, T_TESS, T_in, P_HP, P_PV, P_BESS, P_Grid, SoC_BESS, T4]