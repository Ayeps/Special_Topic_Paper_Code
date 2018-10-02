function [err_cvx_value_ZF,err_cvx_nan_ZF,state_feasible_ZF_SDR,state_eta_nan_ZF_SDR,secure_EE_ZF_SDR] = ZF_SCA_SDR(Omega,h,f,g,q,Nt,J,K,L,P,zeta,xi1,xi2,a1tilde,a2tilde,Gamma1,Gamma2inWatt,Gamma3,Pmax,P_RF,sigma_s,sigma_FU,sigma_PE,P_FU,P_PE)

% function [err_cvx_value_ZF,err_cvx_nan_ZF,state_feasible_ZF_SDR,state_eta_nan_ZF_SDR,secure_EE_ZF_SDR,w_rank_SEE,v_rank_SEE] = ZF_SCA_SDR(Omega,h,f,g,q,Nt,J,K,L,P,zeta,xi1,xi2,a1tilde,a2tilde,Gamma1,Gamma2inWatt,Gamma3,Pmax,P_RF,sigma_s,sigma_FU,sigma_PE,P_FU,P_PE)

state_eta_nan_ZF_SDR =0;
secure_EE_ZF_SDR = zeros(numel(Gamma2inWatt),1)
err_cvx_value_ZF=0 ;
err_cvx_nan_ZF = 0;

% w_rank_SEE = zeros(numel(Gamma2inWatt),L);
% v_rank_SEE = zeros(numel(Gamma2inWatt),L);


H =zeros(Nt,Nt,L,L);
F = zeros(Nt,Nt,L,L,J);
G = zeros(Nt,Nt,L,L,K);
Q=zeros(Nt,Nt,L,P);

%%%%%%%%%%%%   Generate channel matrix first
for ii=1:L
    for jj=1:L
        H(:,:,ii,jj) = h(:,ii,jj)*h(:,ii,jj)';
        for jjj= 1:J
            F(:,:,ii,jj,jjj) = f(:,ii,jj,jjj)*f(:,ii,jj,jjj)';
        end
        for kk = 1:K
            G(:,:,ii,jj,kk) = g(:,ii,jj,kk)*g(:,ii,jj,kk)';
        end
    end
    
    for qq= 1:P
        Q(:,:,ii,qq) = q(:,ii,qq)*q(:,ii,qq)';
    end
end


for Gamma2_counter = 1 :numel(Gamma2inWatt)
    Gamma2 = (Gamma2inWatt(Gamma2_counter));
    
    % user power minimization to get an intiail point
    [Wbar,Vbar,rhobar,state_feasible_ZF_SDR] = ZF_SDR_initial_point(Pmax,H,F,G,Q,sigma_s,sigma_FU,sigma_PE,P_FU,P_PE,J,K,L,P,Nt,Gamma1,Gamma2,Gamma3,xi1,xi2);
    
    %After power minimization, we need to check whether it works or not
    if state_feasible_ZF_SDR == 0
        break
    end
    
    %After deriving Wbar and Vbar, now we set initial point of slack variables 
    [phibar,bbar,lambdabar,thetabar] = ZF_SDR_slk_variable(J,K,L,H,F,G,Wbar,Vbar,rhobar,P_RF,sigma_s,sigma_FU,P_FU,zeta,a1tilde,a2tilde,Omega);
    
    
    %To run CVX
    state_of_convergence =0;
    previous_optimal_value=0;
    
    while state_of_convergence ~= 2
        
        cvx_begin

        variable eta
        variable lambda
        variable theta
        variable rho(L,1)
        variable W1(Nt,Nt) Hermitian
        variable W2(Nt,Nt) Hermitian
        variable V1(Nt,Nt) Hermitian
        variable V2(Nt,Nt) Hermitian
        variable R_FU(L,1)
        variable phi(L,1)
        variable a(L,1)
        variable b(L,1)
        variable c(L,1)
        maximize eta
        subject to
        W1==hermitian_semidefinite(Nt);
        W2==hermitian_semidefinite(Nt);
        V1==hermitian_semidefinite(Nt);
        V2==hermitian_semidefinite(Nt);
        eta <= (exp(lambdabar - thetabar))*(1+lambda - lambdabar - theta + thetabar);
        exp(lambda) <= (a1tilde*(R_FU(1) ) + a2tilde*( R_FU(2)  ));
        (exp(thetabar))*(1+theta - thetabar) >= P_RF + zeta*(real(trace(W1)+trace(W2)+trace(V1)+trace(V2)));
        (exp(a(1)))*1000 <= Omega*(real( trace(W1*H(:,:,1,1))) )*1000;
        (exp(a(2)))*1000 <= Omega*(real( trace(W2*H(:,:,2,2))) )*1000;
        ((exp(bbar(1)))*(1+b(1)-bbar(1)))*1000 >= Omega*( real(trace( H(:,:,2,1)*W2) + trace(H(:,:,1,1)*V1)  + trace(H(:,:,2,1)*V2) ) )*1000;
        ((exp(bbar(2)))*(1+b(2)-bbar(2)))*1000 >= Omega*( real(trace( H(:,:,1,2)*W1) + trace(H(:,:,1,2)*V1)  + trace(H(:,:,2,2)*V2) ) )*1000;
        for jjjjj=1:J
            real(trace( W1*F(:,:,1,1,jjjjj) ))  == 0;
            real(trace( W2*F(:,:,2,2,jjjjj) ))  == 0;
            ( real( trace(W2*F(:,:,2,1,jjjjj)) + trace(V1*F(:,:,1,1,jjjjj)) + trace(V2*F(:,:,2,1,jjjjj))) )*1000000 >= (Gamma2/xi2 - P_PE(1,jjjjj) -sigma_PE(1,jjjjj)^2)*1000000;
            ( real( trace(W1*F(:,:,1,2,jjjjj)) + trace(V1*F(:,:,1,2,jjjjj)) + trace(V2*F(:,:,2,2,jjjjj))) )*1000000 >= (Gamma2/xi2 - P_PE(2,jjjjj) -sigma_PE(2,jjjjj)^2)*1000000;
        end
        for kkk=1:K
            real(trace( W1*G(:,:,1,1,kkk) ))  == 0;
            real(trace( W2*G(:,:,2,2,kkk) ))  == 0;
        end
        ( real(trace(W1) + trace(V1)) )*1000000 <= Pmax(1)*1000000;
        ( real(trace(W2) + trace(V2)) )*1000000 <= Pmax(2)*1000000;
        for ii=1:L
            R_FU(ii) <= log2(1 + exp(phibar(ii)) ) +(exp(phibar(ii)))*(phi(ii) - phibar(ii))/((1+exp(phibar(ii)))*log(2) );
            (exp(phi(ii) -a(ii)))*(exp(b(ii)) + Omega*P_FU(ii) + Omega*sigma_FU(ii)^2 + Omega*(sigma_s^2)*exp(-c(ii)) ) <= 1 ;
            exp(c(ii)) <= rho(ii);
            ( real( trace(W1*H(:,:,1,ii)) + trace(W2*H(:,:,2,ii))+ trace(V1*H(:,:,1,ii)) + trace(V2*H(:,:,2,ii)) ) )*1000000  >= (Gamma1/xi1*(pow_p((1-rho(ii)),-1) ) - P_FU(ii) - sigma_FU(ii)^2)*1000000 ;
            
            for pp=1:P
                (real(  trace(W1*Q(:,:,1,pp)) + trace(W2*Q(:,:,2,pp))+ trace(V1*Q(:,:,1,pp)) + trace(V2*Q(:,:,2,pp)) ))*(1000000) <= Gamma3*(1000000)
            end
        end
        cvx_end
        if (isnan(cvx_optval)==1) || (cvx_optval < previous_optimal_value)
            state_eta_nan_ZF_SDR =1
            
            % The optimal value should be higher than the previous one
            % If it is lower, then the value must be wrong
            % This is to record how many times that happens.
            if (cvx_optval < previous_optimal_value)
                err_cvx_value_ZF = err_cvx_value_ZF + 1
            end
            
            % CVX might be fail / nan
            % This is to record how many times that happens.
            if (isnan(cvx_optval)==1)
                err_cvx_nan_ZF = err_cvx_nan_ZF + 1
            end
            break
        else
            % If it's not nan and the optimal value is higher
            % Then record the optimal value
            current_optimal_value = cvx_optval
        end
        
        % Convergence Criteria, < 1% and 2 times
        if ( (current_optimal_value - previous_optimal_value)/previous_optimal_value <= 0.01 ) && (current_optimal_value > 0)
            state_of_convergence = state_of_convergence + 1
            
        end
        
        if state_of_convergence == 2
            break
        else
            %  push the slk variables higher
            phibar(1) = log( ( real(trace(W1*H(:,:,1,1))) )/( real( trace(W2*H(:,:,2,1) + V1*H(:,:,1,1) + V2*H(:,:,2,1) ) ) + P_FU(1) + sigma_FU(1)^2 + (sigma_s^2)/(rho(1)) ) );
            phibar(2) = log( ( real(trace(W2*H(:,:,2,2))) )/( real( trace(W1*H(:,:,1,2) + V1*H(:,:,1,2) + V2*H(:,:,2,2) ) ) + P_FU(2) + sigma_FU(2)^2 + (sigma_s^2)/(rho(2)) ) );
            bbar(1) = log(Omega*( real(trace( H(:,:,2,1)*W2) + trace(H(:,:,1,1)*V1)  + trace(H(:,:,2,1)*V2) ) ));
            bbar(2) = log(Omega*( real(trace( H(:,:,1,2)*W1) + trace(H(:,:,1,2)*V1)  + trace(H(:,:,2,2)*V2) ) ));
            thetabar = log( P_RF + zeta*(real(trace(W1)+trace(W2)+trace(V1)+trace(V2))) );
            lambdabar= lambda;

            previous_optimal_value = current_optimal_value;
%             phibar=phi;
%             bbar=b;
%             thetabar= theta;        
        end
    end % To this point, it means that we have a convergent solution
    if state_eta_nan_ZF_SDR ==1
        break
    else
        secure_EE_ZF_SDR(Gamma2_counter) = current_optimal_value;
%         w_rank_SEE(Gamma2_counter,:) = [ rank(W1/(real(trace(W1))),0.01)  ; rank(W2/(real(trace(W2))),0.01)];
%         v_rank_SEE(Gamma2_counter,:) = [ rank(V1/(real(trace(V1))),0.01)  ; rank(W2/(real(trace(V2))),0.01)];
    end
    
end 


end