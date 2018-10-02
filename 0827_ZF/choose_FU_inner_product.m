function [h,f] = choose_FU_inner_product(Nt,L,J,K,g,channel_of_users)


h=zeros(Nt,L,L); % channel of femto users
f= zeros(Nt,L,L,J);
%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%
inner_product = zeros(L,J+1,K);
best_user = zeros(L,1);
max_inn_product = zeros(L,J);

%This for loop is to find which user has the smallest g'*h
for cellcounter=1:L  % to do the math in the cellcounter-th femto cell
    for usercounter = 1:(J+1) % to calaulate the channel inner product of the usercounter-th user and all K EVEs
        for EVEcounter = 1:K
           inner_product(cellcounter,usercounter,EVEcounter) = abs( channel_of_users(:,cellcounter,cellcounter,usercounter)'*g(:,cellcounter,cellcounter,EVEcounter) );
        end
        [ max_inn_product(cellcounter,usercounter) , ~] = max(  inner_product(cellcounter,usercounter,:)  );
    end
    [~ , best_user(cellcounter)] = min( max_inn_product(cellcounter,:) );
    h(:,:,cellcounter) = channel_of_users(:,:,cellcounter , best_user(cellcounter));
    f(:,:,cellcounter,:) = channel_of_users(:,:,cellcounter, [1:(best_user(cellcounter)-1) , (best_user(cellcounter)+1) : J+1 ]  );
end

end

