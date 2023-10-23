from math import log
from Hg_property import HG_Props


class COMP:
    def __init__(self, evap_out, cond_in, comp_eff):
        self.evap_out = evap_out
        self.cond_in = cond_in
        self.comp_eff = comp_eff
        
    def __call__(self):
        cond_in_T_sat = HG_Props.PT_sat(self.cond_in.P)
        cond_in_s_gas = HG_Props.S_gas(cond_in_T_sat)
    
        if cond_in_s_gas < self.evap_out.s:
            cond_in_T_ideal = cond_in_T_sat
            while 1:
                cond_in_s_ideal = self.cond_in.Cp_gas*log(cond_in_T_ideal/self.evap_out.T)+self.evap_out.s
                dsdT = self.cond_in.Cp_gas/cond_in_T_ideal
                cond_in_T1_ideal = cond_in_T_ideal - (cond_in_s_ideal - self.evap_out.s)/dsdT;
                if abs(cond_in_T1_ideal - cond_in_T_ideal) < 1.0e-4:
                    break
                else:
                    cond_in_T_ideal = cond_in_T1_ideal
            
            cond_in_h_ideal = self.cond_in.Cp_gas*(cond_in_T_ideal - cond_in_T_sat)+HG_Props.H_gas(cond_in_T_sat)
            self.cond_in.h = (cond_in_h_ideal - self.evap_out.h)/self.comp_eff + self.evap_out.h
            self.cond_in.T = (self.cond_in.h - self.evap_out.h)/self.cond_in.Cp_gas + self.evap_out.T
        else:
            cond_in_s_liq = HG_Props.S_liq(cond_in_T_sat)
            x = 0.5
            while 1:
                cond_in_s_ideal = (cond_in_s_gas - cond_in_s_liq)*x + cond_in_s_liq
                dsdx = cond_in_s_gas - cond_in_s_liq
                x1 = x - (cond_in_s_ideal - self.evap_out.s)/dsdx
                if abs(x1-x) < 1.0e-4:
                    break
                else:
                    x = x1
            
            cond_in_h_ideal = (HG_Props.H_gas(cond_in_T_sat)-HG_Props.H_liq(cond_in_T_sat))*x + HG_Props.H_liq(cond_in_T_sat)
            self.cond_in.h = (cond_in_h_ideal - self.evap_out.h)/self.comp_eff + self.evap_out.h
            self.cond_in.T = cond_in_T_sat
        
        self.W_s = (self.cond_in.h - self.evap_out.h)
        return(self.cond_in.T, self.cond_in.h)