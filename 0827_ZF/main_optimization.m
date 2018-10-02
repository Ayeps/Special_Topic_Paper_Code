clc;
clear
tic
format long;
%cvx_precision high
cvx_solver mosek
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
distance_realization_times = 1;
channel_realization_times = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
femtocell_radius = 10;
shortest_dist_another_cell = 40;
longest_dist_another_cell = 80;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pathloss = -1;
Omega=1; 

%Input data
Nt = 5;% number of BS antenna
J = 1;%number of potential eavesdroppers
K = 1;%number of eavesdroppers
L = 2 ;% number of famto cells
P = 3 ;%number of macro users
zeta = 5;%amplifier coefficient
xi1 = 0.5;%EH efficiency of FU
xi2 = 0.5;%EH efficiency of PE
a1tilde = 0.5;%weight parameter femto cell #1
a2tilde = 0.5;%weight parameter femto cell #2
Gamma1 = 10^(-3 + -30/10); %FU EH threshold -30dBm
Gamma3 = 10^(-3 + -20/10); %threshold of interference to Macro user  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PE EH constraint
Gamma2indBm = [-20,-30];%-10,-12,-14,-16,-18,-20,-22,-24,-26,-28,-30
Gamma2inWatt = 10.^(-3 + Gamma2indBm/10);

%total power dissipation of every single FBS
Pmax = 10^(-3 + 30/10)*ones(L,1); 

%total power dissipation of all FBSs
P_RF = 1; 

%AWGN noise
sigma_s =10^(-4); %noise power = 1e-6 = -50dBm
sigma_FU = (10^(-4))*ones(L,1); %noise power = 1e-8 = -50dBm
sigma_PE = (10^(-4))*ones(L,J); %noise power = 1e-8 = -50dBm
sigma_EVE = (10^(-4))*ones(L,K);%noise power = 1e-8 = -50dBm

%Interference from macro BS constraint
P_FU = (10^(-6))*ones(L,1); %interference power from MBS
P_PE = (10^(-6))*ones(L,J);
P_EVE = (10^(-6))*ones(L,K);

%%  For ZF
actual_realization_times_ZF_SDR = zeros(distance_realization_times,1);
secure_energy_efficiency_ZF_SDR = zeros(numel(Gamma2indBm),channel_realization_times,distance_realization_times);
sum_avgSEE_over_dist_ZF_SDR = zeros(numel(Gamma2indBm),1);
error_at_initial_ZF_SDR = 0;
error_at_cvx_ZF_SDR = 0;
error_at_cvx_value_ZF = 0 ;
error_at_cvx_nan_ZF = 0;



%% Begin Simulation

%Randomly Generate distances first.
%Then, given distances, simulate over different channel realizations
for distance_counter = 1:distance_realization_times
    
    [distance_of_users,distance_of_EVEs] = generate_all_distances(femtocell_radius,shortest_dist_another_cell,longest_dist_another_cell,J,K,L);
    
    sum_SEE_given_distance_ZF_SDR = zeros(numel(Gamma2indBm),1);
    
    for channel_counter = 1:channel_realization_times
        
        
        
        % This is to generate all USERs & EVEs' channel
        [channel_of_users,g,q] = all_users_channel(Nt,L,J,K,P,pathloss,distance_of_users,distance_of_EVEs);
        
        [h,f] = given_FU(Nt,L,J,channel_of_users);
        

        %% this is Zero-Forcing
        %%%%%%%%%
        %%%%%%%%%
        disp(">> bgein ZF_SCA_SDR")
        [err_cvx_value_ZF,err_cvx_nan_ZF,state_feasible_ZF_SDR,state_eta_nan_ZF_SDR,secure_EE_ZF_SDR] = ZF_SCA_SDR(Omega,h,f,g,q,Nt,J,K,L,P,zeta,xi1,xi2,a1tilde,a2tilde,Gamma1,Gamma2inWatt,Gamma3,Pmax,P_RF,sigma_s,sigma_FU,sigma_PE,P_FU,P_PE) 
        disp(">> done ZF_SCA_SDR")
        %%%%%%%%%
        %%%%%%%%%
        
        if state_feasible_ZF_SDR ==0
            error_at_initial_ZF_SDR = error_at_initial_ZF_SDR + 1;
        end
        if state_eta_nan_ZF_SDR ==1
            error_at_cvx_ZF_SDR = error_at_cvx_ZF_SDR + 1;
            error_at_cvx_value_ZF  = error_at_cvx_value_ZF + err_cvx_value_ZF;
            error_at_cvx_nan_ZF  = error_at_cvx_nan_ZF + err_cvx_nan_ZF;
            
        end
        if (state_eta_nan_ZF_SDR ~= 1) && (state_feasible_ZF_SDR ~= 0)
            secure_energy_efficiency_ZF_SDR(:,channel_counter,distance_counter) = secure_EE_ZF_SDR;
            sum_SEE_given_distance_ZF_SDR = sum_SEE_given_distance_ZF_SDR + secure_energy_efficiency_ZF_SDR(:,channel_counter,distance_counter);
            actual_realization_times_ZF_SDR(distance_counter) = actual_realization_times_ZF_SDR(distance_counter) + 1;
            
        end
        

    end %this corresponees to channel realization times (channel_counter)
    avg_SEE_given_distance_ZF_SDR = sum_SEE_given_distance_ZF_SDR/(max(1,actual_realization_times_ZF_SDR(distance_counter)))
    sum_avgSEE_over_dist_ZF_SDR = sum_avgSEE_over_dist_ZF_SDR + avg_SEE_given_distance_ZF_SDR
   
end

toc

avg_SEE_PEEH_ZF_SDR = sum_avgSEE_over_dist_ZF_SDR/distance_realization_times

%save('PEEH0823test.mat');




