from math import exp, log10

class HG_Props:
    @staticmethod
    def TP_sat(T):
        Tcrt = 1764
        Pcrt = 167.0e3 #kPa
        a1 = -4.57618368
        a2 = -1.40726277
        a3 = 2.36263541
        a4 = -31.0889985
        a5 = 58.0183959
        a6 = -27.6304546
        tau = 1-T/Tcrt
        P = Pcrt*exp((Tcrt/T)*(a1*tau+a2*tau**1.89+a3*tau**2+a4*tau**8+a5*tau**8.5+a6*tau**9))
        
        return(P)
    
    @staticmethod
    def PT_sat(P):
        T_lb = 0.0
        T_ub = 1764
        
        while 1:
            T = 0.5*(T_lb+T_ub)
            P1  = HG_Props.TP_sat(T)
            
            P_err = P-P1
            if P_err < 0:
                T_ub = T
            else:
                T_lb = T
                
            if abs(P_err) < 1.0e-6:
                break
            
        return(T)
    
    @staticmethod
    def TdPdT_sat(T):
        Tcrt = 1764
        a1 = -4.57618368
        a2 = -1.40726277
        a3 = 2.36263541
        a4 = -31.0889985
        a5 = 58.0183959
        a6 = -27.6304546
        tau = 1-T/Tcrt
        dtaudT = -1/Tcrt
        P = HG_Props.TP_sat(T)
        dPdT = P*((-Tcrt/T**2)*(a1*tau+a2*tau**1.89+a3*tau**2+a4*tau**8+a5*tau**8.5+a6*tau**9)+(Tcrt/T)*(a1+a2*tau**0.89+a3*tau+a4*tau**7+a5*tau**7.5+a6*tau**8)*dtaudT)
        
        return(dPdT)
    
    @staticmethod
    def H_liq(T):
        [theta, null] = HG_Props.correct_T(T)
        hl = 7.25939*theta-1.36651e-3*theta**2+8.0906e-7*theta**3-1636.13
        hl = hl*4184.0/200.592
        
        return(hl)
    
    @staticmethod
    def H_gas(T):
        [theta, null] = HG_Props.correct_T(T)
        B = HG_Props.second_virial(T)
        dB = HG_Props.dsecond_virial(T)
        T0 = 629.8818
        [theta0, dthetadT0] = HG_Props.correct_T(T0)
        dB0 = HG_Props.dsecond_virial(T0)
        B0 = HG_Props.second_virial(T0)
        hg = 13648.676+4.96797*theta+2022.0*exp(-7136.5/theta)*(B-theta*dB)+0.2503*B0+15.25*dB0+22.53*(T0-theta0)*(dthetadT0+1)+1.419e4*(dthetadT0-1)+4.97*((T-theta)-(T0-theta0))
        hg = hg*4184.0/200.592
        
        return(hg)
        
    @staticmethod
    def Cp_liq(T):
        [theta,dthetadT] = HG_Props.correct_T(T)
        Cpl = 7.25939-2.73302e-3*theta+2.42718e-6*theta**2-2.294e8/theta**2*exp(-7136.5/theta)+6.55*(dthetadT-1)+44585/theta*exp(-7136.5/theta)
        Cpl = Cpl*4184.0/200.592
                
        return(Cpl)
    
    @staticmethod
    def Cp_gas(T):
        [theta, null] = HG_Props.correct_T(T)
        P = HG_Props.TP_sat(T)*7.50062
        B = HG_Props.second_virial(T)
        dB = HG_Props.dsecond_virial(T)
        ddB = -1310/theta**3*(56.4-B)-655/theta**2*dB
        Cpg = 4.96797-3.186e-5*P*theta*ddB
        Cpg = Cpg*4184.0/200.592
        
        return(Cpg)
    
    @staticmethod
    def S_liq(T):
        [theta, null] = HG_Props.correct_T(T)
        T0 = 629.8818
        [theta0, dthetadT0] = HG_Props.correct_T(T0)
        dB0 = HG_Props.dsecond_virial(T0)
        B0 = HG_Props.second_virial(T0)
        int_term = HG_Props.integral_term(theta)
        sl = 16.71536*log10(theta)-2.73302e-3*theta+1.21359e-6*theta**2-4.511*(7136.5/theta+1)*exp(-7136.5/theta)-22.559734-0.02422*dB0-4.359e-4*B0+6.55*int_term+11.44*log10(T0/theta0)-00.003577*(T0-theta0)*dthetadT0-22.53*(dthetadT0-1)
        sl = sl*4184.0/200.592
        
        return(sl)
    
    @staticmethod
    def S_gas(T):
        [theta, null] = HG_Props.correct_T(T)
        P = HG_Props.TP_sat(T)*7.50062
        dB = HG_Props.dsecond_virial(T)
        sg = 11.439185*log10(theta)-4.575674*log10(P)+26.6702+11.44*log10(T/theta)-3.186e-5*P*dB
        sg = sg*4184.0/200.592
        
        return(sg)
    
    @staticmethod
    def correct_T(T):
        theta = T
        
        while 1:
            f = 0.6381+(1-4.809e-3)*theta+1.1096e-5*theta**2-7.481e-9*theta**3 - T
            dfdtheta = 1-4.809e-3 + 2*1.1096e-5*theta-3*7.481e-9*theta**2
            theta1 = theta - f/dfdtheta
            
            if abs(theta1-theta) < 1.0e-4:
                break
            else:
                theta  = theta1
        
        dthetadT = 1/dfdtheta
        
        return(theta, dthetadT)
    
    @staticmethod
    def second_virial(T):
        [theta, null] = HG_Props.correct_T(T)
        B = 56.4-43.82*exp(655/theta)
        
        return(B)
    
    @staticmethod
    def dsecond_virial(T):
        [theta, null] = HG_Props.correct_T(T)
        B = HG_Props.second_virial(T)
        dB = 655/theta**2*(56.4-B)
        
        return(dB)
    
    @staticmethod
    def integral_term(theta):
        int_term = 0.01107*log10(theta)+0.6381/theta-1.1096e-5*theta+3.7405e-9*theta**2-0.0264958
        
        return(int_term)

if __name__ == '__main__':
    print(HG_Props.S_gas(1350.))