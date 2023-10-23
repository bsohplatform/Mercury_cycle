from dataclasses import dataclass
import CoolProp.CoolProp as CP
from HX_module import HX
from COMP_module import COMP
from Hg_property import HG_Props
from math import log
import matplotlib.pyplot as plt


@dataclass
class Props:
    T_sat: float = 0.0 
    T: float = 0.0
    P: float = 0.0
    h: float = 0.0
    s: float = 0.0
    x: float = 0.0
    Q: float = 0.0
    m: float = 0.0
    fluid: str = ''
    Cp_gas: float = 0.0
    Cp_liq: float = 0.0
    Cp: float = 0.0

@dataclass
class inputs:
    comp_eff: float = 0.0
    mech_eff: float = 1.0
    dT_evap = 0.0
    dP_evap = 0.0
    dT_cond = 0.0
    dP_cond = 0.0
    DSH = 0.0
    evap_out_x = 0.0
    DSC = 0.0
    cond_out_x = 0.0
    N_element = 30

@dataclass
class outputs:
    COP: float = 0.0
    W_comp: float = 0.0
    Q_cond: float = 0.0
    Q_evap: float = 0.0



class Mercury_cycle:
    def __init__(self):
        print('---------------------------계산시작---------------------------------')
        
    def preprocess(self, primary, secondary):
       secondary['evap_in'].Cp = CP.PropsSI("C","T",secondary['evap_in'].T,"P",secondary['evap_in'].P,secondary['evap_in'].fluid)
       
       secondary['cond_in'].Cp = CP.PropsSI("C","T",secondary['cond_in'].T,"P",secondary['cond_in'].P,secondary['cond_in'].fluid)
       
       secondary['cond_in'].Q = secondary['cond_in'].Cp*(secondary['cond_out'].T - secondary['cond_in'].T)*secondary['cond_in'].m
       primary['cond_in'].Q = -secondary['cond_in'].Q
       return(secondary)
       
    def solver(self, primary, secondary, inputs, outputs):
        evap_a = 1
        
        evap_P_lb = 1.65e-7
        evap_P_ub = HG_Props.TP_sat(secondary['evap_in'].T - inputs.dT_evap)
        
        while evap_a:
            primary['evap_out'].P = 0.5*(evap_P_lb+evap_P_ub)
            primary['evap_in'].P = primary['evap_out'].P/(1-inputs.dP_evap)
            
            primary['evap_in'].T_sat = HG_Props.PT_sat(primary['evap_in'].P)
            primary['evap_out'].T_sat = HG_Props.PT_sat(primary['evap_out'].P)
            
            primary['evap_out'].Cp_gas =  HG_Props.Cp_gas(primary['evap_out'].T_sat)
            primary['evap_out'].T = primary['evap_out'].T_sat + inputs.DSH
            if inputs.DSH == 0:
                primary['evap_out'].h = (HG_Props.H_gas(primary['evap_out'].T_sat)-HG_Props.H_liq(primary['evap_out'].T_sat))*inputs.evap_out_x+HG_Props.H_liq(primary['evap_out'].T_sat)
                primary['evap_out'].s = (HG_Props.S_gas(primary['evap_out'].T_sat)-HG_Props.S_liq(primary['evap_out'].T_sat))*inputs.evap_out_x+HG_Props.S_liq(primary['evap_out'].T_sat)
            else:
                primary['evap_out'].h = primary['evap_out'].Cp_gas*(primary['evap_out'].T - primary['evap_out'].T_sat) + HG_Props.H_gas(primary['evap_out'].T_sat)
                primary['evap_out'].s = primary['evap_out'].Cp_gas*log(primary['evap_out'].T/primary['evap_out'].T_sat) + HG_Props.S_gas(primary['evap_out'].T_sat)
            
            cond_a = 1
            cond_P_lb = HG_Props.TP_sat(secondary['cond_in'].T + inputs.dT_cond)
            cond_P_ub = 10000.0
            
            while cond_a:
                primary['cond_in'].P = 0.5*(cond_P_lb + cond_P_ub)
                primary['cond_out'].P = primary['cond_in'].P*(1-inputs.dP_cond)
                
                primary['cond_out'].T_sat = HG_Props.PT_sat(primary['cond_out'].P)
                primary['cond_in'].T_sat = HG_Props.PT_sat(primary['cond_in'].P)
                primary['cond_in'].Cp_gas = HG_Props.Cp_gas(primary['cond_in'].T_sat)
                primary['cond_out'].Cp_liq = HG_Props.Cp_liq(primary['cond_out'].T_sat)
                primary['cond_out'].T = primary['cond_out'].T_sat - inputs.DSC
                if inputs.DSC == 0:
                    primary['cond_out'].h = (HG_Props.H_gas(primary['cond_out'].T_sat)-HG_Props.H_liq(primary['cond_out'].T_sat))*inputs.cond_out_x+HG_Props.H_liq(primary['cond_out'].T_sat)
                else:
                    primary['cond_out'].h = HG_Props.H_liq(primary['cond_out'].T_sat) - primary['cond_out'].Cp_liq*(HG_Props.H_gas(primary['cond_out'].T_sat)-primary['cond_out'].T)
                
                comp = COMP(primary['evap_out'], primary['cond_in'], inputs.comp_eff)
                [primary['cond_in'].T, primary['cond_in'].h] = comp()
                
                primary['cond_in'].m = secondary['cond_in'].Q/(primary['cond_in'].h - primary['cond_out'].h)
                primary['evap_in'].m = primary['cond_in'].m
                
                cond = HX(primary['cond_in'], primary['cond_out'], secondary['cond_in'], secondary['cond_out'])
                cond.PHE(inputs.N_element)
                
                cond_err = (inputs.dT_cond - cond.T_pp)/inputs.dT_cond
                
                if cond.T_rvs == 1:
                    cond_P_lb = primary['cond_in'].P
                else:
                    if cond.T_rvs < 0:
                        cond_P_ub = primary['cond_in'].P
                    else:
                        cond_P_lb = primary['cond_in'].P
                
                if abs(cond_err) < 1.0e-4:
                    cond_a = 0
                elif (cond_P_ub - cond_P_lb)/(0.5*(cond_P_ub + cond_P_lb)) < 1.0e-4:
                    cond_a = 0
            
            primary['evap_in'].h = primary['cond_out'].h
            if HG_Props.H_liq(primary['evap_in'].T_sat) < primary['evap_in'].h < HG_Props.H_gas(primary['evap_in'].T_sat):
                primary['evap_in'].T = primary['evap_in'].T_sat
            else:
                primary['evap_in'].Cp_liq = HG_Props.Cp_liq(primary['evap_in'].T_sat)
                primary['evap_in'].T = (HG_Props.H_liq(primary['evap_in'].T_sat)-primary['evap_in'].h)/primary['evap_in'].m/primary['evap_in'].Cp_liq+primary['evap_in'].T_sat
                    
            primary['evap_in'].Q = primary['evap_in'].m*(primary['evap_out'].h - primary['evap_in'].h)
            secondary['evap_in'].Q = -primary['evap_in'].Q
            secondary['evap_in'].m = secondary['evap_in'].Q/(secondary['evap_out'].T - secondary['evap_in'].T)/secondary['evap_in'].Cp
            
            evap = HX(primary['evap_in'], primary['evap_out'], secondary['evap_in'], secondary['evap_out'])
            evap.PHE(inputs.N_element)
            
            evap_err = (inputs.dT_evap - evap.T_pp)/inputs.dT_evap
            
            if evap.T_rvs == 1:
                evap_P_lb = primary['evap_out'].P
            else:
                if evap.T_rvs < 0:
                    evap_P_lb = primary['evap_out'].P
                else:
                    evap_P_ub = primary['evap_out'].P
            
            if abs(evap_err) < 1.0e-4:
                evap_a = 0
            elif (evap_P_ub - evap_P_lb)/(0.5*(evap_P_ub + evap_P_lb)) < 1.0e-4:
                evap_a = 0
        
        outputs.w_comp = comp.W_s*primary['evap_in'].m
        outputs.Q_cond = abs(secondary['cond_in'].Q)
        outputs.Q_evap = abs(secondary['evap_in'].Q)
        outputs.COP = outputs.Q_cond/outputs.w_comp
        
        return(primary, secondary, outputs)
    
    def Plot_diamgram(primary):
        T_dome = range(-39,800)
        s_liq_dome = [HG_Props.S_liq(t+273.15) for t in T_dome]
        s_gas_dome = [HG_Props.S_gas(t+273.15) for t in T_dome]
        P_dome = [HG_Props.TP_sat(t+273.15) for t in T_dome]
        h_liq_dome = [HG_Props.H_liq(t+273.15) for t in T_dome]
        h_gas_dome = [HG_Props.H_gas(t+273.15) for t in T_dome]
                
        incond_hg_s_gas = HG_Props.S_gas(primary['cond_in'].T_sat)
        incond_hg_s_liq = HG_Props.S_liq(primary['cond_in'].T_sat)
        incond_hg_h_gas = HG_Props.H_gas(primary['cond_in'].T_sat)
        incond_hg_h_liq = HG_Props.H_liq(primary['cond_in'].T_sat)
        outcond_hg_s_gas = HG_Props.S_gas(primary['cond_out'].T_sat)
        outcond_hg_s_liq = HG_Props.S_liq(primary['cond_out'].T_sat)
        outcond_hg_h_gas = HG_Props.H_gas(primary['cond_out'].T_sat)
        outcond_hg_h_liq = HG_Props.H_liq(primary['cond_out'].T_sat)

        inevap_hg_s_gas = HG_Props.S_gas(primary['evap_in'].T_sat)
        inevap_hg_s_liq = HG_Props.S_liq(primary['evap_in'].T_sat)
        inevap_hg_h_gas = HG_Props.H_gas(primary['evap_in'].T_sat)
        inevap_hg_h_liq = HG_Props.H_liq(primary['evap_in'].T_sat)
        outevap_hg_s_gas = HG_Props.S_gas(primary['evap_out'].T_sat)
        outevap_hg_s_liq = HG_Props.S_liq(primary['evap_out'].T_sat)
        outevap_hg_h_gas = HG_Props.H_gas(primary['evap_out'].T_sat)
        outevap_hg_h_liq = HG_Props.H_liq(primary['evap_out'].T_sat)
        
        s_vec = [primary['evap_out'].s, primary['cond_in'].s]
        T_vec = [primary['evap_out'].T, primary['cond_in'].T]
        
        if primary['cond_in'].h > incond_hg_h_gas: # 과열증기
            primary['cond_in'].Cp_gas = HG_Props.Cp_gas(primary['cond_in'].T)
            primary['cond_in'].s = primary['cond_in'].Cp_gas*log(primary['cond_in'].T/primary['cond_in'].T_sat)+incond_hg_s_gas
        else:
            primary['cond_in'].s = (incond_hg_s_gas - incond_hg_s_liq)*(primary['cond_in'].h - incond_hg_h_liq)/(incond_hg_h_gas - incond_hg_h_liq)+incond_hg_s_liq
        
        if primary['cond_out'].h < outcond_hg_h_liq:
            primary['cond_out'].Cp_liq = HG_Props.Cp_liq(primary['cond_out'].T)
            primary['cond_out'].s = outcond_hg_s_liq - primary['cond_out'].Cp_liq*log(primary['cond_out'].T_sat/primary['cond_out'].T)
            s_vec.append(outcond_hg_s_liq)
            T_vec.append(primary['cond_out'].T_sat)
        else:
            primary['cond_out'].s = (outcond_hg_s_gas - outcond_hg_s_liq)*(primary['cond_out'].h - outcond_hg_h_liq)/(outcond_hg_h_gas - outcond_hg_h_liq)+outcond_hg_s_liq
        
        s_vec.append(primary['cond_out'].s)
        T_vec.append(primary['cond_out'].T)
            
        primary['evap_in'].s = (inevap_hg_s_gas - inevap_hg_s_liq)*(primary['evap_in'].h - inevap_hg_h_liq)/(inevap_hg_h_gas - inevap_hg_h_liq) + inevap_hg_s_liq
        s_vec.append(primary['evap_in'].s)
        T_vec.append(primary['evap_in'].T)
        
        if primary['evap_out'].h > outevap_hg_h_gas:
            s_vec.append(outevap_hg_s_liq)
            T_vec.append(primary['evap_out'].T_sat)
        
        s_vec.append(primary['evap_out'].s)
        T_vec.append(primary['evap_out'].T)
        
        h_vec = [primary['evap_out'].h, primary['cond_in'].h, primary['cond_out'].h, primary['evap_in'].h, primary['evap_out'].h]
        P_vec = [primary['evap_out'].P, primary['cond_in'].P, primary['cond_out'].P, primary['evap_in'].h, primary['evap_out'].P]
        
        fig_ph, ax_ph = plt.subplots()
        ax_ph.plot([i/1.0e3 for i in h_liq_dome], [i for i in P_dome], [i/1.0e3 for i in h_gas_dome[::-1]], [i for i in P_dome[::-1]])
        ax_ph.plot([i/1.0e3 for i in h_vec], [i for i in P_vec], 'b--')
        ax_ph.set_xlabel('Enthalpy [kJ/kg]',fontsize = 15)
        ax_ph.set_ylabel('Pressure [kPa]',fontsize = 15)
        ax_ph.set_title('Pressure-Enthalpy Diagram')
        ax_ph.tick_params(axis = 'x', labelsize = 13)
        ax_ph.tick_params(axis = 'y', labelsize = 13)
        
        fig_ts, ax_ts = plt.subplots()
        ax_ts.plot([i/1.0e3 for i in s_liq_dome], [i for i in T_dome], [i/1.0e3 for i in s_gas_dome[::-1]], [i for i in T_dome[::-1]])
        ax_ts.plot([i/1.0e3 for i in s_vec], [i-273.15 for i in T_vec], 'b--')
        ax_ts.set_xlabel('Entropy [kJ/kg-K]',fontsize = 15)
        ax_ts.set_ylabel('Temperature [℃]',fontsize = 15)
        ax_ts.set_title('Temperature-Entropy Diagram')
        ax_ts.tick_params(axis = 'x', labelsize = 13)
        ax_ts.tick_params(axis = 'y', labelsize = 13)
        
        fig_ph.save_fig('.\Ph_diagram.png',dpi=300)
        fig_ts.save_fig('.\Ts_diagram.png',dpi=300)


if __name__ == '__main__':

    cond_in = Props()
    cond_out = Props()
    evap_in = Props()
    evap_out = Props()

    cond_in_hg = Props()
    cond_out_hg = Props()
    evap_in_hg = Props()
    evap_out_hg = Props()

    inputs = inputs()
    outputs = outputs()


    cond_in.T = 600 + 273.15
    cond_out.T = 620 + 273.15
    cond_in.m = 1.0
    cond_in.P = 1.013e5
    cond_out.P = 1.013e5
    cond_in.fluid = 'Air'

    evap_in.T = 400 + 273.15
    evap_in.T = 380 + 273.15
    evap_in.P = 1.013e5
    evap_in.fluid = 'Air'

    inputs.dT_evap = 5.0
    inputs.dT_cond = 5.0
    inputs.comp_eff = 0.4
    inputs.mech_eff = 0.5
    inputs.DSH = 0.0
    inputs.DSC = 0.0
    inputs.evap_out_x = 0.4
    inputs.cond_in_x = 0.2
    
    primary = {"cond_in":cond_in_hg, "cond_out":cond_out_hg, "evap_in":evap_in_hg, "evap_out":evap_out_hg}
    secondary = {"cond_in":cond_in, "cond_out":cond_out, "evap_in":evap_in, "evap_out":evap_out}
    
    HG = Mercury_cycle()
    secondary = HG.preprocess(primary, secondary)
    (primary, secondary, outputs) = HG.solver(primary, secondary, inputs, outputs)
    HG.Plot_diamgram(primary)
