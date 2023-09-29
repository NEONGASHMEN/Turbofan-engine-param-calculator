import math

glb = 1.0

def get_inlet_params(Ma, Ta, Pa, P_drop_inlet, T_drop_inlet, mdot,  gamma, Cp):
    Tta = Ta*get_isentropic_vals(["M",Ma], "TR", gamma)
    Pta = Pa*get_isentropic_vals(["M",Ma], "PR", gamma)
    return Pta*P_drop_inlet, Tta*T_drop_inlet, mdot

def get_compressor_output(Tt_i, Pt_i, PR, eta, mdot, gamma, Cp):
    Pt_o = Pt_i*PR
    Tt_o_ideal = Tt_i*get_isentropic_vals(["PR",PR], "TR", gamma)
    Tt_o = Tt_i + (1/eta)*(Tt_o_ideal - Tt_i)
    W = mdot*Cp*(Tt_o - Tt_i)
    return Pt_o, Tt_o, mdot, W

def get_bypass_flow_rates(BPR, mdot_total):
    mdot_c = mdot_total/(1 + BPR)
    mdot_b = mdot_total - mdot_c
    return mdot_c, mdot_b
    
def get_cc_output(Tt_i, Pt_i, LHV, eta, mdot_i, P_drop, Tt_o, gamma, Cp):
    Pt_o = P_drop*Pt_i
    mdot_fuel = mdot_i*Cp*(Tt_o - Tt_i)/(LHV*eta)
    W = mdot_fuel*LHV*eta
    return Pt_o, Tt_o, mdot_fuel+mdot_i, W, mdot_fuel

def get_turbine_output(Tt_i, Pt_i, eta_t, eta_m, W, mdot, gamma, Cp):
    Tt_o = Tt_i - W/(eta_m*mdot*Cp)
    Tt_o_ideal = Tt_i - (Tt_i - Tt_o)/eta_t
    Pt_o = Pt_i*get_isentropic_vals(["TR",Tt_o_ideal/Tt_i], "PR", gamma)
    return Pt_o, Tt_o, mdot, W/eta_m

def get_nozzle_output(Tt_i, Pt_i, eta, mdot, P_a, gamma, Cp, R):
    Pt_o = Pt_i
    Tt_o = Tt_i
    P_o = Pt_o*((1-(1/eta)*(gamma-1)/(gamma+1))**(gamma/(gamma-1)))
    if Pt_o/P_o <= Pt_o/P_a:
        M = 1.0
        TR = get_isentropic_vals(["M",1], "TR", gamma)
        T_o = Tt_o/TR
        return P_o, T_o, M, mdot
    else:
        PR = Pt_o/P_a
        TR = get_isentropic_vals(["PR",PR], "TR", gamma)
        T_o = Tt_o/TR
        V = math.sqrt(2*Cp*eta*Tt_o*(1-(PR**((1-gamma)/gamma))))
        M = V/math.sqrt(gamma*R*T_o)
        #print(M)
        #print(mdot,"mdot")
        return P_a, T_o, M, mdot
    #Not sure about unchoked cndtn

def get_nozzle_area(V, mdot, P, T, R):
    A = mdot/(V*P/(R*T))
    return A

def get_thrust(Mi, Ta, Pa, R_a, gamma_a, Mo, To, Po, R_o, gamma_o, mdot_o, A):
    Vi = Mi*math.sqrt(gamma_a*R_a*Ta)
    Vo = Mo*math.sqrt(gamma_o*R_o*To)
    F = mdot_o*(Vo - Vi) + A*(Po - Pa)
    return F

def get_sfc(F, mdot):
    return mdot/F

def get_isentropic_vals(arg1, arg2, gamma):
    if arg1[0] == "TR" and arg2 == "M":
        TR = arg1[1]
        return math.sqrt((2/(gamma-1))*(TR-1))
    elif arg1[0] == "PR"and arg2 == "M":
        PR = arg1[1]
        TR = PR**((gamma-1)/gamma)
        return math.sqrt((2/(gamma-1))*(TR-1))
    elif arg1[0] == "TR" and arg2 == "PR":
        TR = arg1[1]
        return TR**(gamma/(gamma-1))
    elif arg1[0] == "M" and arg2 == "PR":
        M = arg1[1]
        TR = 1 + ((gamma-1)/2)*M**2
        return TR**(gamma/(gamma-1))
    elif arg1[0] == "PR" and arg2 == "TR":
        PR = arg1[1]
        return PR**((gamma-1)/gamma)
    elif arg1[0] == "M" and arg2 == "TR":
        M = arg1[1]
        return 1 + ((gamma-1)/2)*M**2
    else:
        print("Invalid Input!")
        return "Error"
 
def get_isentropic_velocity(M, T, R, gamma):
    V = M*math.sqrt(gamma*R*T)
    return V

def get_dS(Tt1, Tt2, Pt1, Pt2, R, Cp):
    dS = Cp*math.log(Tt2/Tt1) - R*math.log(Pt2/Pt1)
    return dS

def get_thermal_efficiency(mdot_core, mdot_bypass, mdot_fuel, F_core, F_bypass, V_o, LHV):
    V_j_core_eff = F_core/mdot_core + V_o
    V_j_bypass_eff = F_bypass/mdot_bypass + V_o
    KE = 0.5*mdot_core*(V_j_core_eff**2 - V_o**2) + 0.5*mdot_bypass*(V_j_bypass_eff**2 - V_o**2)
    eta = KE/(mdot_fuel*LHV)
    return eta

def get_propulsive_efficiency(mdot_core, mdot_bypass, F_core, F_bypass, V_o):
    V_j_core_eff = F_core/mdot_core + V_o
    V_j_bypass_eff = F_bypass/mdot_bypass + V_o
    thrust_pow = mdot_core*(V_j_core_eff - V_o)*V_o + mdot_bypass*(V_j_bypass_eff - V_o)*V_o
    KE = 0.5*mdot_core*(V_j_core_eff**2 - V_o**2) + 0.5*mdot_bypass*(V_j_bypass_eff**2 - V_o**2)
    eta = thrust_pow/KE
    return eta

