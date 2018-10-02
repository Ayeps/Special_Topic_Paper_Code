function [phibar,bbar,lambdabar,thetabar] = ZF_SDR_slk_variable(J,K,L,H,F,G,Wbar,Vbar,rhobar,P_RF,sigma_s,sigma_FU,P_FU,zeta,a1tilde,a2tilde,Omega)

quad_wh_bar = ones(L,L); % to calcaulare w'*h first
quad_vh_bar = ones(L,L); % to calcaulare v'*h first

for ii=1:L
    for jj=1:L
        quad_wh_bar(ii,jj) = trace(Wbar(:,:,ii)*H(:,:,ii,jj));  
        quad_vh_bar(ii,jj) = trace(Vbar(:,:,ii)*H(:,:,ii,jj));  
    end
end
quad_wf_bar = ones(L,L,J); % to calcaulare w'*h first
quad_vf_bar = ones(L,L,J); % to calcaulare v'*h first
for ii=1:L
    for jj=1:L
       for jjjj=1:J
          quad_wf_bar(ii,jj,jjjj) = trace(Wbar(:,:,ii)*F(:,:,ii,jj,jjjj)); 
          quad_vf_bar(ii,jj,jjjj) = trace(Vbar(:,:,ii)*F(:,:,ii,jj,jjjj));
       end
    end
end
quad_wg_bar = ones(L,L,K); % to calcaulare w'*g first
quad_vg_bar = ones(L,L,K); % to calcaulare v'*g first
for ii=1:L
    for jj=1:L
       for kk=1:K
          quad_wg_bar(ii,jj,kk) = trace(Wbar(:,:,ii)*G(:,:,ii,jj,kk));
          quad_vg_bar(ii,jj,kk) = trace(Vbar(:,:,ii)*G(:,:,ii,jj,kk));
       end
    end
end
phibar = ones(L,1); 
FUratebar = ones(L,1);
bbar = ones(L,1);
for ii=1:L
   for jj=1:L
       phibar(ii) = real(log((rhobar(ii)*quad_wh_bar(ii,ii))/(rhobar(ii)*(quad_wh_bar([1:(ii-1),(ii+1):L],ii)+sum(quad_vh_bar(:,ii))+P_FU(ii) +sigma_FU(ii)^2)+sigma_s^2) ) );
       FUratebar(ii) = real(log2(1 +exp(phibar(ii))));
       bbar(ii) = real(log(Omega*(quad_wh_bar([1:(ii-1),(ii+1):L],ii) +sum(quad_vh_bar(:,ii)))));
   end
end

beamforming_power = ones(L,1); % the power of wbar
AN_power = ones(L,1); % the power of vbar
for ii=1:L
    beamforming_power(ii) = real(trace(Wbar(:,:,ii))); 
    AN_power(ii) = real(trace(Vbar(:,:,ii)));
end

lambdabar = real(log(a1tilde*(FUratebar(1)  )+a2tilde*(FUratebar(2)  )  ));
thetabar =real(log(P_RF +zeta*( sum(beamforming_power) + sum(AN_power)) ));

end