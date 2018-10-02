function [h,f,g,q] = channel(Nt,L,J,K,P)
pathloss = -1.5;
h=zeros(Nt,L,L); % channel of femto users
% H =zeros(Nt,Nt,L,L);

for ii = 1:L
    for jj = 1:L
        if jj == ii
            h(:,ii,jj) = (randn(Nt,1) + (-1)^(0.5)*randn(Nt,1))/sqrt(2)*(10^(pathloss));
        else
            h(:,ii,jj) = (randn(Nt,1) + (-1)^(0.5)*randn(Nt,1))/sqrt(2)*(80^(pathloss));
        end
%         H(:,:,ii,jj) =  h(:,ii,jj)*h(:,ii,jj)';
    end
end
f = zeros(Nt,L,L,J); % channel of potential eavesdroppers
% F = zeros(Nt,Nt,L,L,J);
for ii = 1 :L
    for jj = 1:L 
        for kk = 1:J
            if ii==jj
                f(:,ii,jj,kk) = (randn(Nt,1) + (-1)^(0.5)*randn(Nt,1))/sqrt(2)*(10^(pathloss));
            else
                f(:,ii,jj,kk) = (randn(Nt,1) + (-1)^(0.5)*randn(Nt,1))/sqrt(2)*(80^(pathloss));
            end 
%             F(:,:,ii,jj,kk) = f(:,ii,jj,kk)*f(:,ii,jj,kk)';
        end
    end
end
g=zeros(Nt,L,L,K);  % channel of eavesdroppers
% G = zeros(Nt,Nt,L,L,K);
for ii=1:L
    for jj=1:L
       for kk = 1:K
            if ii==jj
                g(:,ii,jj,kk) = (randn(Nt,1) + (-1)^(0.5)*randn(Nt,1))/sqrt(2)*(10^(pathloss));
            else
                g(:,ii,jj,kk) = (randn(Nt,1) + (-1)^(0.5)*randn(Nt,1))/sqrt(2)*(80^(pathloss));
            end
%             G(:,:,ii,jj,kk) = g(:,ii,jj,kk)*g(:,ii,jj,kk)';
       end
    end
end
q=zeros(Nt,L,P);
% Q=zeros(Nt,Nt,L,P);
for ii=1:L
    for jj=1:P
        q(:,ii,jj) = (randn(Nt,1) + (-1)^(0.5)*randn(Nt,1))/sqrt(2)*(100^(pathloss));
%         Q(:,:,ii,jj) = q(:,ii,jj)*q(:,ii,jj)';
    end
end

end