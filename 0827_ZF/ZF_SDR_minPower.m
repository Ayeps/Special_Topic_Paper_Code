function [Wbar,Vbar,rhobar,state_feasible] = ZF_SDR_minPower(FU_SINR_threshold,Pmax,H,F,G,Q,sigma_s,sigma_FU,sigma_PE,P_FU,P_PE,J,K,L,P,Nt,Gamma1,Gamma2,Gamma3,xi1,xi2)
%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%

Wbar = zeros(Nt,Nt,L);
Vbar = zeros(Nt,Nt,L);
rhobar = zeros(L,1);

                    cvx_begin
                    variable W1bar(Nt,Nt) Hermitian
                    variable W2bar(Nt,Nt) Hermitian
                    variable V1bar(Nt,Nt) Hermitian
                    variable V2bar(Nt,Nt) Hermitian
                    variable rhobar(L,1) 
                    minimize real(trace(W1bar) + trace(V1bar) + trace(W2bar) + trace(V2bar))
                    subject to
                    W1bar==hermitian_semidefinite(Nt);
                    W2bar==hermitian_semidefinite(Nt);
                    V1bar==hermitian_semidefinite(Nt);
                    V2bar==hermitian_semidefinite(Nt);
                    ((real(trace(W1bar*H(:,:,1,1)))))*10000 >= (FU_SINR_threshold*( real(trace( H(:,:,2,1)*W2bar) + trace(H(:,:,1,1)*V1bar)  + trace(H(:,:,2,1)*V2bar) )  + P_FU(1) + sigma_FU(1)^2 + ((sigma_s^2)*pow_p( rhobar(1) , -1  )  ) ))*10000;
                    ((real(trace(W2bar*H(:,:,2,2)))))*10000 >=  (FU_SINR_threshold*(real(trace( H(:,:,1,2)*W1bar) + trace(H(:,:,1,2)*V1bar)  + trace(H(:,:,2,2)*V2bar) )  + P_FU(2) + sigma_FU(2)^2 + ((sigma_s^2)*pow_p( rhobar(2) , -1  )  ) ))*10000;
                    for jjj=1:J
                        real(trace(W1bar*F(:,:,1,1,jjj)))  == 0 ;
                        real(trace(W2bar*F(:,:,2,2,jjj)))  == 0 ;
                        ( real( trace(W2bar*F(:,:,2,1,jjj)) + trace(V1bar*F(:,:,1,1,jjj)) + trace(V2bar*F(:,:,2,1,jjj))) )*1000000 >= (Gamma2/xi2 - P_PE(1,jjj) -sigma_PE(1,jjj)^2)*1000000;  
                        ( real( trace(W1bar*F(:,:,1,2,jjj)) + trace(V1bar*F(:,:,1,2,jjj)) + trace(V2bar*F(:,:,2,2,jjj))) )*1000000 >= (Gamma2/xi2 - P_PE(2,jjj) -sigma_PE(2,jjj)^2)*1000000;  
                    end
                    ( real(trace(W1bar) + trace(V1bar)) )*1000000 <= Pmax(1)*1000000;
                    ( real(trace(W2bar) + trace(V2bar)) )*1000000 <= Pmax(2)*1000000;
                    
                    for ii=1:L           
                        ( real( trace(W1bar*H(:,:,1,ii)) + trace(W2bar*H(:,:,2,ii))+ trace(V1bar*H(:,:,1,ii)) + trace(V2bar*H(:,:,2,ii)) ) )*1000000  >= ( Gamma1/xi1*( pow_p( (1-rhobar(ii)),-1 ) ) - P_FU(ii) - sigma_FU(ii)^2 )*1000000 ;
                        rhobar(ii) >= 0.00001;
                        rhobar(ii) <= 1 ;
                    end
                    for pp=1:P
                            (real( trace(W1bar*Q(:,:,1,pp)) + trace(W2bar*Q(:,:,2,pp))+ trace(V1bar*Q(:,:,1,pp)) + trace(V2bar*Q(:,:,2,pp)) ))*(1000000) <= Gamma3*(1000000)
                    end
                    for kk=1:K
                             real(trace(W1bar*G(:,:,1,1,kk)))  == 0 ;
                             real(trace(W2bar*G(:,:,2,2,kk)))  == 0 ;
                    end
                    cvx_end
                    
                    if (isnan(cvx_optval)==1) || (cvx_optval >= 100)
                        state_feasible =0
                       
                    else
                       state_feasible = 1
                       Wbar(:,:,1) = W1bar;
                       Wbar(:,:,2) = W2bar;
                       Vbar(:,:,1) = V1bar;
                       Vbar(:,:,2) = V2bar;
                    end
                    
                    
end






