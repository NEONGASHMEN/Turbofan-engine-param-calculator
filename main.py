import math
import matplotlib as mpl
from matplotlib import pyplot as plt
import numpy as np
from functions import *
from turbofan_param_calculator import *

#All in SI units

gamma_a = 1.4
Cp_a = 1000.0

#Inlet
Ma = 0.78
Ta = 218.8
Pa = 23842.0
P_drop_2 = 0.98
T_drop_2 = 1.0
mdot_2 = 173.0
Pt_2, Tt_2, mdot_2 = get_inlet_params(Ma, Ta, Pa, P_drop_2, T_drop_2, mdot_2, gamma_a, Cp_a)
print("Inlet Pt, Tt & mdot:")
print(Pt_2, Tt_2, mdot_2)
print("------------------------------------------\n")

#Fan
PR_fan = 1.4
eta_fan = 0.9
Pt_21, Tt_21, mdot_21, W_fan = get_compressor_output(Tt_2, Pt_2, PR_fan, eta_fan, mdot_2, gamma_a, Cp_a)
print("Fan Output - Pt, Tt, mdot & W:")
print(Pt_21, Tt_21, mdot_21, W_fan)
print("------------------------------------------\n")

#Bypass
BPR = 12
mdot_21, mdot_13 = get_bypass_flow_rates(BPR, mdot_21)
print("Bypass flow - mdot_core & mdot_bypass:")
print(mdot_21, mdot_13)
print("------------------------------------------\n")

#LPC
PR_LPC = 1.7
eta_LPC = 0.92
Pt_25, Tt_25, mdot_25, W_LPC = get_compressor_output(Tt_21, Pt_21, PR_LPC, eta_LPC, mdot_21, gamma_a, Cp_a)
print("LPC Output - Pt, Tt, mdot & W:")
print(Pt_25, Tt_25, mdot_25, W_LPC)
print("------------------------------------------\n")

#HPC
PR_HPC = 12.5
eta_HPC = 0.92
Pt_3, Tt_3, mdot_3, W_HPC = get_compressor_output(Tt_25, Pt_25, PR_HPC, eta_HPC, mdot_25, gamma_a, Cp_a)
print("HPC Output - Pt, Tt, mdot & W:")
print(Pt_3, Tt_3, mdot_3, W_HPC)
print("------------------------------------------\n")

#CC
LHV = 43e6
eta_cc = 0.995
Cp_g = 1150
gamma_g = 1.33
P_drop_cc = 0.96
Pt_4, Tt_4, mdot_4, Q, mdot_fuel = get_cc_output(Tt_3, Pt_3, LHV, eta_cc, mdot_3, P_drop_cc, 1400.0, gamma_g, Cp_g)
print("CC Output - Pt, Tt, mdot, Q & fuel_flow:")
print(Pt_4, Tt_4, mdot_4, Q, mdot_fuel)
print("------------------------------------------\n")

#HPT
eta_m_HPshaft = 0.99
eta_HPT = 0.9
Pt_45, Tt_45, mdot_45, W_HPT = get_turbine_output(Tt_4, Pt_4, eta_HPT, eta_m_HPshaft, W_HPC, mdot_4, gamma_g, Cp_g)
print("HPT Output - Pt, Tt, mdot & W:")
print(Pt_45, Tt_45, mdot_45, W_HPT)
print("------------------------------------------\n")

#LPT
eta_m_LPshaft = 0.99
eta_LPT = 0.9
Pt_5, Tt_5, mdot_5, W_LPT = get_turbine_output(Tt_45, Pt_45, eta_LPT, eta_m_LPshaft, W_LPC+W_fan, mdot_45, gamma_g, Cp_g)
print("LPT Output - Pt, Tt, mdot & W:")
print(Pt_5, Tt_5, mdot_5, W_LPT)
print("------------------------------------------\n")

#Core Nozzle
eta_n_core = 0.98
R_g = 287.0
P_8, T_8, M_8, mdot_8 = get_nozzle_output(Tt_5, Pt_5, eta_n_core, mdot_5, Pa, gamma_g, Cp_g, R_g)
print("Nozzle Output - P, T, M & mdot:")
print(P_8, T_8, M_8, mdot_8)
print("------------------------------------------\n")

#Bypass Nozzle
eta_n_bypass = 0.98
Tt_16 = Tt_21
Pt_16 = Pt_21
R_a = R_g
P_16, T_16, M_16, mdot_16 = get_nozzle_output(Tt_16, Pt_16, eta_n_bypass, mdot_13, Pa, gamma_a, Cp_a, R_a)
print("Bypass Nozzle Output - P, T, M & mdot:")
print(P_16, T_16, M_16, mdot_16)
print("------------------------------------------\n")

#Engine performance
A_core = get_nozzle_area(M_8*math.sqrt(gamma_g*R_g*T_8), mdot_8, P_8, T_8, R_g)
F_core = get_thrust(Ma, Ta, Pa, R_a, gamma_a, M_8, T_8, P_8, R_g, gamma_g, mdot_8, A_core)
A_bypass = get_nozzle_area(M_16*math.sqrt(gamma_a*R_a*T_16), mdot_16, P_16, T_16, R_a)
F_bypass = get_thrust(Ma, Ta, Pa, R_a, gamma_a, M_16, T_16, P_16, R_a, gamma_a, mdot_16, A_bypass)
SFC = get_sfc(F_core + F_bypass, mdot_fuel)
print("Core thrust: ", F_core)
print("Bypass thrust: ", F_bypass)
print("Total thrust: ", F_core + F_bypass)
print("Specific Fuel Consumption: ", SFC)
print("\n")

 
#TS Diagram
plt.style.use('dark_background')
TS_fig = plt.figure()
TS_fig.canvas.manager.set_window_title("TS Diagram")
plt.xlabel("S (J/K)")
plt.ylabel("Tt (K)")

S1 = 0
Tta = Ta*get_isentropic_vals(["M",Ma], "TR", gamma_a)
Pta = Pa*get_isentropic_vals(["M",Ma], "PR", gamma_a)
S2 = S1 + get_dS(Tta, Tt_2, Pta, Pt_2, R_a, Cp_a)
S21 = S2 + get_dS(Tt_2, Tt_21, Pt_2, Pt_21, R_a, Cp_a)
S25 = S21 + get_dS(Tt_21, Tt_25, Pt_21, Pt_25, R_a, Cp_a)
S3 = S25 + get_dS(Tt_25, Tt_3, Pt_25, Pt_3, R_a, Cp_a)
S4 = S3 + get_dS(Tt_3, Tt_4, Pt_3, Pt_4, R_g, Cp_g)
S45 = S4 + get_dS(Tt_4, Tt_45, Pt_4, Pt_45, R_g, Cp_g)
S5 = S45 + get_dS(Tt_45, Tt_5, Pt_45, Pt_5, R_g, Cp_g)
Tt_8 = T_8*get_isentropic_vals(["M",M_8], "TR", gamma_g)
Pt_8 = P_8*get_isentropic_vals(["M",M_8], "PR", gamma_g)
S8 = S5 + get_dS(Tt_5, Tt_8, Pt_5, Pt_8, R_g, Cp_g)
S = np.array([S1,S2,S21,S25,S3,S4,S45,S5])
T = np.array([Tta,Tt_2,Tt_21,Tt_25,Tt_3,Tt_4,Tt_45,Tt_5])

plt.scatter(S, T, c='w')
state_names = ["1","2","21","25","3","4","45","5,7"]
for i in range(len(S)):
    plt.annotate(state_names[i],[S[i],T[i]])
plt.plot(S, T, c='w')

overshoot = 0.0
single_line_range = 200.0
for i in range(4):
    CP_line_y = np.linspace(T[i]-overshoot,T[i]+single_line_range,100)
    CP_line_x = np.zeros(len(CP_line_y))
    CP_line_x[0] = S[i]
    for j in range(len(CP_line_y)-1):
        CP_line_x[j+1] = CP_line_x[j] + get_dS(CP_line_y[j], CP_line_y[j+1], 1, 1, R_a, Cp_a)
    plt.plot(CP_line_x, CP_line_y, c='#606661')

single_line_range = 900.0
CP_line_y = np.linspace(T[4]-overshoot,T[4]+single_line_range,100)
CP_line_x = np.zeros(len(CP_line_y))
CP_line_x[0] = S[4]
for j in range(len(CP_line_y)-1):
    CP_line_x[j+1] = CP_line_x[j] + get_dS(CP_line_y[j], CP_line_y[j+1], 1, 1, R_a, Cp_a)
plt.plot(CP_line_x, CP_line_y, c='#606661')

single_line_range = 700.0
CP_line_y = np.linspace(T[5]+overshoot,T[5]-single_line_range,100)
CP_line_x = np.zeros(len(CP_line_y))
CP_line_x[0] = S[5]
for j in range(len(CP_line_y)-1):
    CP_line_x[j+1] = CP_line_x[j] + get_dS(CP_line_y[j], CP_line_y[j+1], 1, 1, R_a, Cp_a)
plt.plot(CP_line_x, CP_line_y, c='#606661')


single_line_range = 300.0
for i in range(6,8):
    CP_line_y = np.linspace(T[i]+overshoot,T[i]-single_line_range,100)
    CP_line_x = np.zeros(len(CP_line_y))
    CP_line_x[0] = S[i]
    for j in range(len(CP_line_y)-1):
        CP_line_x[j+1] = CP_line_x[j] + get_dS(CP_line_y[j], CP_line_y[j+1], 1, 1, R_a, Cp_a)
    plt.plot(CP_line_x, CP_line_y, c='#606661')

#BPR & FPR Variation
F_list = []
SFC_list = []
BPR = [10,11,12]
FPR = [1.4,1.5]

const = ["gamma_air", 1.4, "gamma_mixture", 1.33, "Cp_air", 1000.0, "Cp_mixture", 1150.0, "R_air", 287.0, "R_mixture", 287.0]
inlet = ["Ma", 0.78, "Pa", 23842.0, "Ta", 218.8, "P_drop", 0.98, "T_drop", 1.0, "mdot", 173.0]
fan = ["PR", 1.4, "eta", 0.9, "BPR", 12]
LPC = ["PR", 1.7, "eta", 0.92]
HPC = ["PR", 12.5, "eta", 0.92]
CC = ["LHV", 43e6, "eta", 0.995, "P_drop", 0.96, "CC_exit_temp", 1400.0]
HPT = ["eta", 0.9, "eta_shaft", 0.99]
LPT = ["eta", 0.9, "eta_shaft", 0.99]
nozzle_c = ["eta", 0.98]
nozzle_b = ["eta", 0.98]

print(80*"#")
for i in FPR:
    for j in BPR:
        fan[5] = j
        fan[1] = i
        F, SFC, eta_therm, eta_prop = turbofan_param_calculator(const, inlet, fan, LPC, HPC, CC, HPT, LPT, nozzle_c, nozzle_b)
        print("Thrust for FPR = ", i, " & BPR = ", j, " is; ", F)
        print("SFC for FPR = ", i, " & BPR = ", j, " is; ", SFC)
print(80*"#")


#HPC PR & Tcc variation 
HPC_PR = [12, 13, 14]
Tcc = [1400.0, 1500.0]
fan[5] = 12
fan[1] = 1.4
print(80*"#")
for i in Tcc:
    for j in HPC_PR:
        CC[7] = i
        HPC[1] = j
        F, SFC, eta_therm, eta_prop = turbofan_param_calculator(const, inlet, fan, LPC, HPC, CC, HPT, LPT, nozzle_c, nozzle_b)
        print("Thrust for Tcc = ", i, " & HPC_PR = ", j, " is; ", F)
        print("SFC for Tcc = ", i, " & HPC_PR = ", j, " is; ", SFC)
print(80*"#")

#Thrust, SFC, Thermal & Propulsive efficiency vs FPR, BPR, OPR, TIT
CC[7] = 1400.0
HPC[1] = 12.5
x_axis = np.linspace(1.4,1.5,51)
y_axis_F = np.zeros(len(x_axis))
y_axis_SFC = np.zeros(len(x_axis))
y_axis_eta_therm = np.zeros(len(x_axis))
y_axis_eta_prop = np.zeros(len(x_axis))
for i in range(len(x_axis)):
    fan[1] = x_axis[i]
    y_axis_F[i], y_axis_SFC[i], y_axis_eta_therm[i], y_axis_eta_prop[i] = turbofan_param_calculator(const, inlet, fan, LPC, HPC, CC, HPT, LPT, nozzle_c, nozzle_b)
    #print(x_axis[i],"--------")
FPR_fig = plt.figure()
FPR_fig.canvas.manager.set_window_title("Performance wrt Fan Pressure Ratio")
FPR_fig.set_figwidth(11)
FPR_fig.set_figheight(6)
plt.subplot(1,2,1)
plt.title("Thrust v/s FPR")
plt.xlabel("FPR")
plt.ylabel("Thrust (N)")
plt.plot(x_axis, y_axis_F, c='w')
plt.subplot(1,2,2)
plt.title("SFC v/s FPR")
plt.xlabel("FPR")
plt.ylabel("SFC (kg/Ns)")
plt.plot(x_axis, y_axis_SFC, c='w')
FPR_eta_fig = plt.figure()
FPR_eta_fig.canvas.manager.set_window_title("Efficiency wrt Fan Pressure Ratio")
FPR_eta_fig.set_figwidth(11)
FPR_eta_fig.set_figheight(6)
plt.subplot(1,2,1)
plt.title("Thermal Efficiency v/s FPR")
plt.xlabel("FPR")
plt.ylabel("Thermal Efficiency")
plt.plot(x_axis, y_axis_eta_therm, c='w')
plt.subplot(1,2,2)
plt.title("Propulsive Efficiency v/s FPR")
plt.xlabel("FPR")
plt.ylabel("Propulsive Efficiency")
plt.plot(x_axis, y_axis_eta_prop, c='w')


fan[1] = 1.4
x_axis = np.linspace(10,12,51)
y_axis_F = np.zeros(len(x_axis))
y_axis_SFC = np.zeros(len(x_axis))
y_axis_eta_therm = np.zeros(len(x_axis))
y_axis_eta_prop = np.zeros(len(x_axis))
for i in range(len(x_axis)):
    fan[5] = x_axis[i]
    y_axis_F[i], y_axis_SFC[i], y_axis_eta_therm[i], y_axis_eta_prop[i] = turbofan_param_calculator(const, inlet, fan, LPC, HPC, CC, HPT, LPT, nozzle_c, nozzle_b)
BPR_fig = plt.figure()
BPR_fig.canvas.manager.set_window_title("Performance wrt Bypass Ratio")
BPR_fig.set_figwidth(11)
BPR_fig.set_figheight(6)
plt.subplot(1,2,1)
plt.title("Thrust v/s BPR")
plt.xlabel("BPR")
plt.ylabel("Thrust (N)")
plt.plot(x_axis, y_axis_F, c='w')
plt.subplot(1,2,2)
plt.title("SFC v/s BPR")
plt.xlabel("BPR")
plt.ylabel("SFC (kg/Ns)")
plt.plot(x_axis, y_axis_SFC, c='w')
BPR_eta_fig = plt.figure()
BPR_eta_fig.canvas.manager.set_window_title("Efficiency wrt Bypass Ratio")
BPR_eta_fig.set_figwidth(11)
BPR_eta_fig.set_figheight(6)
plt.subplot(1,2,1)
plt.title("Thermal Efficiency v/s BPR")
plt.xlabel("BPR")
plt.ylabel("Thermal Efficiency")
plt.plot(x_axis, y_axis_eta_therm, c='w')
plt.subplot(1,2,2)
plt.title("Propulsive Efficiency v/s BPR")
plt.xlabel("BPR")
plt.ylabel("Propulsive Efficiency")
plt.plot(x_axis, y_axis_eta_prop, c='w')


fan[5] = 12
x_axis = np.linspace(1400.0,1500.0,51)
y_axis_F = np.zeros(len(x_axis))
y_axis_SFC = np.zeros(len(x_axis))
y_axis_eta_therm = np.zeros(len(x_axis))
y_axis_eta_prop = np.zeros(len(x_axis))
for i in range(len(x_axis)):
    CC[7] = x_axis[i]
    y_axis_F[i], y_axis_SFC[i], y_axis_eta_therm[i], y_axis_eta_prop[i] = turbofan_param_calculator(const, inlet, fan, LPC, HPC, CC, HPT, LPT, nozzle_c, nozzle_b)
CC_fig = plt.figure()
CC_fig.canvas.manager.set_window_title("Performance wrt Turbine Inlet Temperature")
CC_fig.set_figwidth(11)
CC_fig.set_figheight(6)
plt.subplot(1,2,1)
plt.title("Thrust v/s TIT")
plt.xlabel("TIT (K)")
plt.ylabel("Thrust (N)")
plt.plot(x_axis, y_axis_F, c='w')
plt.subplot(1,2,2)
plt.title("SFC v/s TIT")
plt.xlabel("TIT (K)")
plt.ylabel("SFC (kg/Ns)")
plt.plot(x_axis, y_axis_SFC, c='w')
CC_eta_fig = plt.figure()
CC_eta_fig.canvas.manager.set_window_title("Efficiency wrt Turbine Inlet Temperature")
CC_eta_fig.set_figwidth(11)
CC_eta_fig.set_figheight(6)
plt.subplot(1,2,1)
plt.title("Thermal Efficiency v/s TIT")
plt.xlabel("TIT (K)")
plt.ylabel("Thermal Efficiency")
plt.plot(x_axis, y_axis_eta_therm, c='w')
plt.subplot(1,2,2)
plt.title("Propulsive Efficiency v/s TIT")
plt.xlabel("TIT (K)")
plt.ylabel("Propulsive Efficiency")
plt.plot(x_axis, y_axis_eta_prop, c='w')


CC[7] = 1400.0
x_axis = np.linspace(12,14,51)
y_axis_F = np.zeros(len(x_axis))
y_axis_SFC = np.zeros(len(x_axis))
y_axis_eta_therm = np.zeros(len(x_axis))
y_axis_eta_prop = np.zeros(len(x_axis))
for i in range(len(x_axis)):
    HPC[1] = x_axis[i]
    y_axis_F[i], y_axis_SFC[i], y_axis_eta_therm[i], y_axis_eta_prop[i] = turbofan_param_calculator(const, inlet, fan, LPC, HPC, CC, HPT, LPT, nozzle_c, nozzle_b)
HPC_fig = plt.figure()
HPC_fig.canvas.manager.set_window_title("Performance wrt PR in HPC")
HPC_fig.set_figwidth(11)
HPC_fig.set_figheight(6)
plt.subplot(1,2,1)
plt.title("Thrust v/s PR in HPC")
plt.xlabel("PR")
plt.ylabel("Thrust (N)")
plt.plot(x_axis, y_axis_F, c='w')
plt.subplot(1,2,2)
plt.title("SFC v/s PR in HPC")
plt.xlabel("PR")
plt.ylabel("SFC (kg/Ns)")
plt.plot(x_axis, y_axis_SFC, c='w')
HPC_eta_fig = plt.figure()
HPC_eta_fig.canvas.manager.set_window_title("Efficiency wrt PR in HPC")
HPC_eta_fig.set_figwidth(11)
HPC_eta_fig.set_figheight(6)
plt.subplot(1,2,1)
plt.title("Thermal Efficiency v/s PR in HPC")
plt.xlabel("PR")
plt.ylabel("Thermal Efficiency")
plt.plot(x_axis, y_axis_eta_therm, c='w')
plt.subplot(1,2,2)
plt.title("Propulsive Efficiency v/s PR in HPC")
plt.xlabel("PR")
plt.ylabel("Propulsive Efficiency")
plt.plot(x_axis, y_axis_eta_prop, c='w')


plt.show()
