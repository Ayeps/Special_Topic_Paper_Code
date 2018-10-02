function [h,f] = choose_FU_ch_gain(Nt,L,J,channel_of_users)


h=zeros(Nt,L,L); % channel of femto users
f= zeros(Nt,L,L,J);
%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%
channel_gain = zeros(L,J+1);
best_user = zeros(L,1);

%This for loop is to find which user has the smallest g'*h
for cellcounter=1:L  % to do the math in the cellcounter-th femto cell
    for usercounter = 1:(J+1) % to calaulate the channel inner product of the usercounter-th user and all K EVEs
        channel_gain(cellcounter,usercounter) = norm(channel_of_users(:,cellcounter,cellcounter,usercounter))
    end
    [~ , best_user(cellcounter)] = max( channel_gain(cellcounter,:) );
    h(:,:,cellcounter) = channel_of_users(:,:,cellcounter , best_user(cellcounter));
    f(:,:,cellcounter,:) = channel_of_users(:,:,cellcounter, [1:(best_user(cellcounter)-1) , (best_user(cellcounter)+1) : J+1 ]  );
end

end

