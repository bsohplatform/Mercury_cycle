import numpy as np
import CoolProp.CoolProp as CP
from Hg_property import HG_Props

class HX:
    def __init__(self, primary_in, primary_out, secondary_in, secondary_out):
        self.primary_in = primary_in
        self.primary_out = primary_out
        self.secondary_in = secondary_in
        self.secondary_out = secondary_out
        
    def PHE(self, N_element: int):
        h_primary = np.zeros(shape=(N_element+1))
        T_primary = np.zeros(shape=(N_element+1))
        P_primary = np.zeros(shape=(N_element+1))
        T_secondary = np.zeros(shape=(N_element+1))
        h_secondary = np.zeros(shape=(N_element+1))
        P_secondary = np.zeros(shape=(N_element+1))
        dT = np.zeros(shape=(N_element+1))
        
        P_primary[0]= self.primary_in.P
        T_primary[0]= self.primary_in.T
        h_primary[0]= self.primary_in.h
        h_secondary[0]= self.secondary_out.h
        T_secondary[0]= self.secondary_out.T
        P_secondary[0]= self.secondary_out.P
        
        dT[0] = T_primary[0] - T_secondary[0]
        
        if (dT[0] < 0.0 and self.primary_in.Q < 0.0) or (dT[0] > 0.0 and self.primary_in.Q > 0.0):
            self.T_rvs = 1
        else:
            self.T_rvs = 0
    
        for i in range(N_element):
            if self.T_rvs == 1:
                break
            
            h_primary[i+1] = h_primary[i] + self.primary_in.Q/N_element/self.primary_in.m
            P_primary[i+1] = P_primary[i] - (self.primary_in.P - self.primary_out.P)/N_element
            T_sat_now = HG_Props.PT_sat(P_primary[i+1])
            h_gas = HG_Props.H_gas(T_sat_now)
            h_liq = HG_Props.H_liq(T_sat_now)
            
            if h_gas < h_primary[i+1]: # 과열 기체
                T_primary[i+1] = h_primary[i+1]/self.primary_in.Cp_gas
            elif h_liq > h_primary[i+1]: # 과냉 액체
                T_primary[i+1] = h_primary[i+1]/self.primary_out.Cp_liq
            else:
                T_primary[i+1] = T_sat_now
            
            P_secondary[i+1] = P_secondary[i] + (self.secondary_in.P - self.secondary_out.P)/N_element
            h_secondary[i+1] = h_secondary[i] + self.primary_in.Q/N_element/self.secondary_in.m
            T_secondary[i+1] = CP.PropsSI("T","H",h_secondary[i+1],"P",P_secondary[i+1],self.secondary_in.fluid)
    
            dT[i+1] = T_primary[i+1] - T_secondary[i+1]
            if (dT[i+1] < 0.0 and self.primary_in.Q < 0.0) or (dT[i+1] > 0.0 and self.primary_in.Q > 0.0):
                self.T_rvs = 1
            else:
                self.T_rvs = 0
      
        self.T_pp = min(abs(dT))
        
        self.primary_out.T = T_primary[N_element]
        self.primary_out.h = h_primary[N_element]
        
        return(T_secondary)