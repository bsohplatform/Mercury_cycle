from math import log
from Hg_property import HG_Props


class COMP:
    def __init__(self, evap_out, cond_in, comp_eff):
        self.evap_out = evap_out
        self.cond_in = cond_in
        self.comp_eff = comp_eff
        
    def __call__(self):
        cond_in_s_gas = HG_Props.S_gas(self.cond_in.T_sat)
    
        if cond_in_s_gas < self.evap_out.s:
            mid_Tsat_lb = self.evap_out.T_sat
            mid_Tsat_ub = self.cond_in.T_sat
            while 1:
                mid_Tsat = 0.5*(mid_Tsat_lb + mid_Tsat_ub)
                mid_s_gas = HG_Props.S_gas(mid_Tsat)
                s_err = (self.evap_out.s - mid_s_gas)/self.evap_out.s
                if s_err < 0.0:
                    mid_Tsat_lb = mid_Tsat
                else:
                    mid_Tsat_ub = mid_Tsat
                if abs(s_err) <= 1.0e-6:
                    break
                    
            
            mid_Cp = HG_Props.Cp_gas(mid_Tsat)
            gamma = mid_Cp/(mid_Cp-8314.0/200.592)
            cond_in_T_ideal = self.evap_out.T*(self.cond_in.P/self.evap_out.P)**((gamma-1)/gamma)
            cond_in_h_ideal = mid_Cp*(cond_in_T_ideal - self.cond_in.T_sat)+HG_Props.H_gas(mid_Tsat)
            self.cond_in.h = (cond_in_h_ideal - self.evap_out.h)/self.comp_eff + self.evap_out.h
            self.cond_in.T = (self.cond_in.h - self.evap_out.h)/self.cond_in.Cp_gas + self.evap_out.T
        else:
            cond_in_s_liq = HG_Props.S_liq(self.cond_in.T_sat)
            x = 0.5
            while 1:
                cond_in_s_ideal = (cond_in_s_gas - cond_in_s_liq)*x + cond_in_s_liq
                dsdx = cond_in_s_gas - cond_in_s_liq
                x1 = x - (cond_in_s_ideal - self.evap_out.s)/dsdx
                if abs(x1-x) < 1.0e-4:
                    break
                else:
                    x = x1
            
            cond_in_h_ideal = (HG_Props.H_gas(self.cond_in.T_sat)-HG_Props.H_liq(self.cond_in.T_sat))*x + HG_Props.H_liq(self.cond_in.T_sat)
            self.cond_in.h = (cond_in_h_ideal - self.evap_out.h)/self.comp_eff + self.evap_out.h
            self.cond_in.T = self.cond_in.T_sat
        
        self.W_s = (self.cond_in.h - self.evap_out.h)
        return(self.cond_in.T, self.cond_in.h)