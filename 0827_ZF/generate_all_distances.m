function [distance_of_users,distance_of_EVEs] = generate_all_distances(femtocell_radius,shortest_dist_another_cell,longest_dist_another_cell,J,K,L)

number_of_total_users = J+1;

distance_of_users = zeros(L,L,number_of_total_users);
distance_of_EVEs = zeros(L,L,K);

%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%
% This for loop is to randomly generate the distances for all users & EVEs
for ii=1:L
    for jj= 1:L
        for kk=1:number_of_total_users
            if ii==jj
                distance_of_users(ii,jj,kk) = unidrnd(femtocell_radius);
            else
                distance_of_users(ii,jj,kk) = randi([shortest_dist_another_cell longest_dist_another_cell]);
            end
        end
        
        for kkk=1:K
            if ii==jj
                distance_of_EVEs(ii,jj,kkk) = unidrnd(femtocell_radius);
            else
                distance_of_EVEs(ii,jj,kkk) = randi([shortest_dist_another_cell longest_dist_another_cell]);
            end
        end
    end    
end


end

