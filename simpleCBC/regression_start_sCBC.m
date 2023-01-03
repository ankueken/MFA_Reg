%% Flux estimation from 13C labeling data by constrained regression
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% non-compartmented model of CBC %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc;
addpath('..')

%% load data (Treves et al.)
[Data_mean_abs, Data_std_abs, Data_mean_norm, ~,AbsoluteLevel,~,Data_std_abs_original] = load_experimental_data;

%% ceate model and convert to model with irreversible rxns only
network_simpleCBC
model_rev=model;
model=convertToIrreversible(model);

%% Build matrix S combining MID data 
algea_c=1:5; %chlamy
algea_s=6:10; %sorokiniana
algea_o=11:15; %ohadii
algea_o_eil=16:20; %ohadii EI

%chlamy
[S_t0_a1,S_t5_a1,S_t10_a1,S_t20_a1] = build_differential_matrix_sCBC(Data_mean_norm(algea_c,:),model);
%sorokiniana
[S_t0_a2,S_t5_a2,S_t10_a2,S_t20_a2] = build_differential_matrix_sCBC(Data_mean_norm(algea_s,:),model);
%ohadii
[S_t0_a3,S_t5_a3,S_t10_a3,S_t20_a3] = build_differential_matrix_sCBC(Data_mean_norm(algea_o,:),model);
%ohadii EIL
[S_t0_a4,S_t5_a4,S_t10_a4,S_t20_a4] = build_differential_matrix_sCBC(Data_mean_norm(algea_o_eil,:),model);

%% index of used metabolites MIDs in MID data table
rubp_MIDs_used = find(contains(Data_mean_abs.Properties.VariableNames,'RuBP'));
fbp_MIDs_used = find(contains(Data_mean_abs.Properties.VariableNames,'FBP'));
f6p_MIDs_used = find(contains(Data_mean_abs.Properties.VariableNames,'f6p'));
g6p_MIDs_used = find(contains(Data_mean_abs.Properties.VariableNames,'G6P'));
g1p_MIDs_used = find(contains(Data_mean_abs.Properties.VariableNames,'g1p'));
adpg_MIDs_used = find(contains(Data_mean_abs.Properties.VariableNames,'ADPG'));

%% vector including measured starch synthesis rates
sta = [19 51 85 74];
d_sta = sta*0.3;

%% run regression
% CHLAMY
[v_a1, residuals_a1, Chi_a1, df_a1, EXITFLAG_a1, H_a1] = regression_sCBC(S_t0_a1,S_t5_a1,S_t10_a1,S_t20_a1, Data_mean_abs(1:5,:),...
    Data_std_abs(1,:), Data_std_abs(2,:), Data_std_abs(3,:), Data_std_abs(4,:),Data_std_abs(5,:),...
    rubp_MIDs_used, f6p_MIDs_used, g6p_MIDs_used, fbp_MIDs_used, g1p_MIDs_used, adpg_MIDs_used, model,...
    sta(1),d_sta(1));

disp('expected interval for chi-square:')
disp([chi2inv(0.025,df_a1) chi2inv(0.975,df_a1)])
disp('')
if Chi_a1<chi2inv(0.975,df_a1)
    disp('Chlamy: fit accepted')
else
    disp('Chlamy: fit rejected')
end

% SOROKINIANA
[v_a2, residuals_a2, Chi_a2, df_a2, EXITFLAG_a2, H_a2] = regression_sCBC(S_t0_a2,S_t5_a2,S_t10_a2,S_t20_a2, Data_mean_abs(6:10,:), ...
    Data_std_abs(6,:),Data_std_abs(7,:), Data_std_abs(8,:), Data_std_abs(9,:),Data_std_abs(10,:),...
    rubp_MIDs_used, f6p_MIDs_used, g6p_MIDs_used, fbp_MIDs_used, g1p_MIDs_used, adpg_MIDs_used, model,...
    sta(2),d_sta(2));

disp('')
if Chi_a2<chi2inv(0.975,df_a2)
    disp('Sorokiniana: fit accepted')
else
    disp('Sorokiniana: fit rejected')
end

% OHADII
[v_a3, residuals_a3, Chi_a3, df_a3, EXITFLAG_a3, H_a3] = regression_sCBC(S_t0_a3,S_t5_a3,S_t10_a3,S_t20_a3, Data_mean_abs(11:15,:), ...
    Data_std_abs(11,:), Data_std_abs(12,:), Data_std_abs(13,:), Data_std_abs(14,:),Data_std_abs(15,:),...
    rubp_MIDs_used, f6p_MIDs_used, g6p_MIDs_used, fbp_MIDs_used, g1p_MIDs_used, adpg_MIDs_used, model,...
    sta(3),d_sta(3));

disp('')
if Chi_a3<chi2inv(0.975,df_a3)
    disp('Ohadii: fit accepted')
else
    disp('Ohadii: fit rejected')
end

% OHADII EIL
[v_a4, residuals_a4, Chi_a4, df_a4, EXITFLAG_a4, H_a4] = regression_sCBC(S_t0_a4,S_t5_a4,S_t10_a4,S_t20_a4, Data_mean_abs(16:20,:), ...
    Data_std_abs(16,:), Data_std_abs(17,:), Data_std_abs(18,:), Data_std_abs(19,:),Data_std_abs(20,:),...
    rubp_MIDs_used, f6p_MIDs_used, g6p_MIDs_used, fbp_MIDs_used, g1p_MIDs_used, adpg_MIDs_used, model,...
    sta(4),d_sta(4));
    
disp('')
if Chi_a4<chi2inv(0.975,df_a4)
    disp('Ohadii eil: fit accepted')
else
    disp('Ohadii eil: fit rejected')
end

bootstrap_non_parametric_sCBC

v_a1=get_net_flux(v_a1,model,model_rev);
v_a2=get_net_flux(v_a2,model,model_rev);
