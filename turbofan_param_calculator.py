import math
from functions import *

#every argument is an array of the form ["name", value] pairs
#constants require "gamma_air", "gamma_mixture", "Cp_air", "Cp_mixture", "R_air" & "R_mixture"
#inlet require "Ma", "Ta", "Pa", "P_drop", "T_drop" & "mdot"
#fan require "PR", "eta" and "BPR"
#LPC require "PR" & "eta"
#HPC require "PR" & "eta"
#CC require "LHV", "eta", "P_drop" & "CC_exit_temp"
#HPT require "eta" & "eta_shaft"
#LPT require "eta" & "eta_shaft"
#core_nozzle require "eta"
#bypass_nozzle require "eta"


def turbofan_param_calculator(constants, inlet, fan, LPC, HPC, CC, HPT, LPT, core_nozzle, bypass_nozzle):


    gamma_a = constants[constants.index("gamma_air")+1]
    gamma_g = constants[constants.index("gamma_mixture")+1]
    Cp_a = constants[constants.index("Cp_air")+1]
    Cp_g = constants[constants.index("Cp_mixture")+1]
    R_a = constants[constants.index("R_air")+1]
    R_g = constants[constants.index("R_mixture")+1]

    #-------------------------Inlet------------------------------------#
    Ma = inlet[inlet.index("Ma")+1]
    Ta = inlet[inlet.index("Ta")+1]
    Pa = inlet[inlet.index("Pa")+1]
    P_drop_2 = inlet[inlet.index("P_drop")+1]
    T_drop_2 = inlet[inlet.index("T_drop")+1]
    mdot_2 = inlet[inlet.index("mdot")+1]
    Pt_2, Tt_2, mdot_2 = get_inlet_params(Ma, Ta, Pa, P_drop_2, T_drop_2, mdot_2, gamma_a, Cp_a)
    #print(get_inlet_params(0.78, 218.8, 23482.0, 0.98, 1.0, 173.0, 1.4, 1000.0))
    #print(Pt_2, Tt_2)

    #-------------------------Fan--------------------------------------#
    PR_fan = fan[fan.index("PR")+1]
    eta_fan = fan[fan.index("eta")+1]
    Pt_21, Tt_21, mdot_21, W_fan = get_compressor_output(Tt_2, Pt_2, PR_fan, eta_fan, mdot_2, gamma_a, Cp_a)
    #print(Tt_2, Pt_2, PR_fan, eta_fan, mdot_2, gamma_a, Cp_a)
    #print(get_compressor_output(Tt_2, Pt_2, PR_fan, eta_fan, mdot_2, gamma_a, Cp_a))
    #print(PR_fan, eta_fan)
    #print(Pt_21, Tt_21)

    #-----------------------Bypass-------------------------------------#
    BPR = fan[fan.index("BPR")+1]
    mdot_21, mdot_13 = get_bypass_flow_rates(BPR, mdot_21)

    #-------------------------LPC--------------------------------------#
    PR_LPC = LPC[LPC.index("PR")+1]
    eta_LPC = LPC[LPC.index("eta")+1]
    Pt_25, Tt_25, mdot_25, W_LPC = get_compressor_output(Tt_21, Pt_21, PR_LPC, eta_LPC, mdot_21, gamma_a, Cp_a)
    #print(Pt_25, Tt_25, mdot_25)

    #-------------------------HPC--------------------------------------#
    PR_HPC = HPC[HPC.index("PR")+1]
    eta_HPC = HPC[HPC.index("eta")+1]
    Pt_3, Tt_3, mdot_3, W_HPC = get_compressor_output(Tt_25, Pt_25, PR_HPC, eta_HPC, mdot_25, gamma_a, Cp_a)
    #print(PR_HPC, eta_HPC)
    #print(Pt_3, Tt_3, mdot_3)

    #--------------------------CC--------------------------------------#
    LHV = CC[CC.index("LHV")+1]
    eta_cc = CC[CC.index("eta")+1]
    P_drop_cc = CC[CC.index("P_drop")+1]
    T_4 = CC[CC.index("CC_exit_temp")+1]
    Pt_4, Tt_4, mdot_4, Q, mdot_fuel = get_cc_output(Tt_3, Pt_3, LHV, eta_cc, mdot_3, P_drop_cc, T_4, gamma_g, Cp_g)

    #-------------------------HPT--------------------------------------#
    eta_m_HPshaft = HPT[HPT.index("eta_shaft")+1]
    eta_HPT = HPT[HPT.index("eta")+1]
    Pt_45, Tt_45, mdot_45, W_HPT = get_turbine_output(Tt_4, Pt_4, eta_HPT, eta_m_HPshaft, W_HPC, mdot_4, gamma_g, Cp_g)

    #-------------------------LPT--------------------------------------#
    eta_m_LPshaft = LPT[LPT.index("eta_shaft")+1]
    eta_LPT = LPT[LPT.index("eta")+1]
    Pt_5, Tt_5, mdot_5, W_LPT = get_turbine_output(Tt_45, Pt_45, eta_LPT, eta_m_LPshaft, W_LPC+W_fan, mdot_45, gamma_g, Cp_g)

    #---------------------Core--Nozzle---------------------------------#
    eta_n_core = core_nozzle[core_nozzle.index("eta")+1]
    P_8, T_8, M_8, mdot_8 = get_nozzle_output(Tt_5, Pt_5, eta_n_core, mdot_5, Pa, gamma_g, Cp_g, R_g)
    A_core = get_nozzle_area(M_8*math.sqrt(gamma_g*R_g*T_8), mdot_8, P_8, T_8, R_g)
    F_core = get_thrust(Ma, Ta, Pa, R_a, gamma_a, M_8, T_8, P_8, R_g, gamma_g, mdot_8, A_core)
    #print(P_8, T_8, M_8, mdot_8)

    #-------------------Bypass--Nozzle---------------------------------#
    Tt_16 = Tt_21
    Pt_16 = Pt_21
    eta_n_bypass = bypass_nozzle[bypass_nozzle.index("eta")+1]
    P_16, T_16, M_16, mdot_16 = get_nozzle_output(Tt_16, Pt_16, eta_n_bypass, mdot_13, Pa, gamma_a, Cp_a, R_a)
    A_bypass = get_nozzle_area(M_16*math.sqrt(gamma_a*R_a*T_16), mdot_16, P_16, T_16, R_a)
    F_bypass = get_thrust(Ma, Ta, Pa, R_a, gamma_a, M_16, T_16, P_16, R_a, gamma_a, mdot_16, A_bypass)

    #--------------------Engine--Data---------------------------------#
    F = F_core + F_bypass
    SFC = get_sfc(F, mdot_fuel)
    eta_therm = get_thermal_efficiency(mdot_8, mdot_16, mdot_fuel, M_8*math.sqrt(gamma_g*R_g*T_8), M_16*math.sqrt(gamma_a*R_a*T_16), Ma*math.sqrt(gamma_a*R_a*Ta), LHV)
    eta_therm = get_thermal_efficiency(mdot_8, mdot_16, mdot_fuel, F_core, F_bypass, Ma*math.sqrt(gamma_a*R_a*Ta), LHV)
    eta_prop = get_propulsive_efficiency(mdot_8, mdot_16, F_core, F_bypass, Ma*math.sqrt(gamma_a*R_a*Ta))

    return F, SFC, eta_therm, eta_prop

















