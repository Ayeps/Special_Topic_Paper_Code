function [Wbar,Vbar,rhobar,state_feasible] = ZF_SDR_initial_point(Pmax,H,F,G,Q,sigma_s,sigma_FU,sigma_PE,P_FU,P_PE,J,K,L,P,Nt,Gamma1,Gamma2,Gamma3,xi1,xi2)


            [Wbar,Vbar,rhobar,state_feasible] = ZF_SDR_minPower(3,Pmax,H,F,G,Q,sigma_s,sigma_FU,sigma_PE,P_FU,P_PE,J,K,L,P,Nt,Gamma1,Gamma2,Gamma3,xi1,xi2);
            
            if state_feasible == 0
                [Wbar,Vbar,rhobar,state_feasible] = ZF_SDR_minPower(2,Pmax,H,F,G,Q,sigma_s,sigma_FU,sigma_PE,P_FU,P_PE,J,K,L,P,Nt,Gamma1,Gamma2,Gamma3,xi1,xi2);
            end
            
            if state_feasible == 0
                [Wbar,Vbar,rhobar,state_feasible] = ZF_SDR_minPower(1,Pmax,H,F,G,Q,sigma_s,sigma_FU,sigma_PE,P_FU,P_PE,J,K,L,P,Nt,Gamma1,Gamma2,Gamma3,xi1,xi2);
            end
            
            if state_feasible == 0
                [Wbar,Vbar,rhobar,state_feasible] = ZF_SDR_minPower(0.2,Pmax,H,F,G,Q,sigma_s,sigma_FU,sigma_PE,P_FU,P_PE,J,K,L,P,Nt,Gamma1,Gamma2,Gamma3,xi1,xi2);
            end
            
            if state_feasible == 0
                [Wbar,Vbar,rhobar,state_feasible] = ZF_SDR_minPower(0.05,Pmax,H,F,G,Q,sigma_s,sigma_FU,sigma_PE,P_FU,P_PE,J,K,L,P,Nt,Gamma1,Gamma2,Gamma3,xi1,xi2);
            end
            if state_feasible == 0
                [Wbar,Vbar,rhobar,state_feasible] = ZF_SDR_minPower(10^(-3),Pmax,H,F,G,Q,sigma_s,sigma_FU,sigma_PE,P_FU,P_PE,J,K,L,P,Nt,Gamma1,Gamma2,Gamma3,xi1,xi2);
            end
            if state_feasible == 0
                [Wbar,Vbar,rhobar,state_feasible] = ZF_SDR_minPower(10^(-5),Pmax,H,F,G,Q,sigma_s,sigma_FU,sigma_PE,P_FU,P_PE,J,K,L,P,Nt,Gamma1,Gamma2,Gamma3,xi1,xi2);
            end
end


            