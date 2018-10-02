function [channel_of_users,g,q] = all_users_channel(Nt,L,J,K,P,pathloss,distance_of_users,distance_of_EVEs)

number_of_total_users = J+1;

channel_of_users = zeros(Nt,L,L,number_of_total_users);
g=zeros(Nt,L,L,K);  % channel of eavesdroppers
q=zeros(Nt,L,P);
%%%%%%%%%%%%%%%%% The following is to define the vectors and the matrixs of
%%%%%%%%%%%%%%%%% channels
%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%
% This for loop is to randomly generate the distances for all users & EVEs

%The following is the channels for all the users, including one FU & J PEs. But we haven't chosen which one is gonne be the FU)
for ii= 1:L
    for jj=1:L
        for uu = 1:number_of_total_users
        if ii==jj
            channel_of_users(:,ii,jj,uu) = (randn(Nt,1) + (-1)^(0.5)*randn(Nt,1))/sqrt(2)*( (distance_of_users(ii,jj,uu))^(pathloss));
        else
            channel_of_users(:,ii,jj,uu) = (randn(Nt,1) + (-1)^(0.5)*randn(Nt,1))/sqrt(2)*( (distance_of_users(ii,jj,uu))^(pathloss));
        end
        end
    end
end

for ii=1:L
    for jj=1:L
       for kk = 1:K
            if ii==jj
                g(:,ii,jj,kk) = (randn(Nt,1) + (-1)^(0.5)*randn(Nt,1))/sqrt(2)*((distance_of_EVEs(ii,jj,kk))^(pathloss));
            else
                g(:,ii,jj,kk) = (randn(Nt,1) + (-1)^(0.5)*randn(Nt,1))/sqrt(2)*((distance_of_EVEs(ii,jj,kk))^(pathloss));
            end
       end
    end
end


for ii=1:L
    for jj=1:P
        q(:,ii,jj) = (randn(Nt,1) + (-1)^(0.5)*randn(Nt,1))/sqrt(2)*((randi([60 100]))^(pathloss));
    end
    
end
end

