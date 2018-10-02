function [h,f] = given_FU(Nt,L,J,channel_of_users)


h=zeros(Nt,L,L); % channel of femto users
f= zeros(Nt,L,L,J);
%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%


for cellcounter=1:L 
    h(:,:,cellcounter) = channel_of_users(:,:,cellcounter , 1);
    f(:,:,cellcounter,:) = channel_of_users(:,:,cellcounter,  (2: J+1)   );
end

end

