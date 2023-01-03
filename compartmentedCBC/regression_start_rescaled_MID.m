clear
close all
addpath('..')

% load experimental data
% MID absolute and relative (mean, std)
% pool size measurements (mean, std)
[Data_mean_abs, Data_std_abs, Data_mean_norm,~,AbsoluteLevel] = load_experimental_data(true);

algea_c=1:5; %chlamy
algea_s=6:10; %sorokiniana
algea_o=11:15; %ohadii
algea_o_eil=16:20; %ohadii eil

network_compCBC % load network structure
model_rev=model; % split model into irreversible but also keep the original model in variable model_rev
model=convertToIrreversible(model);

%% vectors including measured starch and sucrose synthesis rates
sta = [19 51 85 74]; % mean 
d_sta = sta*0.3; % std
suc = [1e-7 17 35 61];
d_suc = suc.*0.3;

%% specify idx of used MIDs
rubp_MIDs_used = find(contains(Data_mean_abs.Properties.VariableNames,'RuBP'));
fbp_MIDs_used = find(contains(Data_mean_abs.Properties.VariableNames,'FBP'));
f6p_MIDs_used = find(contains(Data_mean_abs.Properties.VariableNames,'f6p'));
g6p_MIDs_used = find(contains(Data_mean_abs.Properties.VariableNames,'G6P'));
g1p_MIDs_used = find(contains(Data_mean_abs.Properties.VariableNames,'g1p'));
adpg_MIDs_used = find(contains(Data_mean_abs.Properties.VariableNames,'ADPG'));
udpg_MIDs_used = find(contains(Data_mean_abs.Properties.VariableNames,'UDPG'));
model.lb(1)=1e-3;
for sample_comp=1 % sample 50 parameterization of c and calculate fit, take the best fit and do bootstrap

%% create matrix that combines the MID data
%chlamy
[S_t0_a1,S_t5_a1,S_t10_a1,S_t20_a1,comp_a1(:,sample_comp)] = build_differential_matrix_comp(Data_mean_norm(algea_c,:),model,[0.38 0.36 0.18 0 0.001]);
%sorokiniana
[S_t0_a2,S_t5_a2,S_t10_a2,S_t20_a2,comp_a2(:,sample_comp)] = build_differential_matrix_comp(Data_mean_norm(algea_s,:),model,[0 0.003 0 0 0]);
%ohadii
[S_t0_a3,S_t5_a3,S_t10_a3,S_t20_a3,comp_a3(:,sample_comp)] = build_differential_matrix_comp(Data_mean_norm(algea_o,:),model,[2.4 0.65 9999 0.002 1.68]);
%ohadii EIL
[S_t0_a4,S_t5_a4,S_t10_a4,S_t20_a4,comp_a4(:,sample_comp)] = build_differential_matrix_comp(Data_mean_norm(algea_o_eil,:),model);

%% run regression
% CHLAMY
[v_a1(:,sample_comp), residuals_a1(:,sample_comp), Chi_a1(sample_comp), df_a1, EXITFLAG_a1, H_a1] = regression_T40(S_t0_a1,S_t5_a1,S_t10_a1,S_t20_a1, Data_mean_abs(1:5,:),...
    Data_std_abs(1,:), Data_std_abs(2,:), Data_std_abs(3,:), Data_std_abs(4,:),Data_std_abs(5,:),...
    rubp_MIDs_used, f6p_MIDs_used, g6p_MIDs_used, fbp_MIDs_used, g1p_MIDs_used, adpg_MIDs_used, udpg_MIDs_used, model,...
    sta(1),d_sta(1),suc(1),d_suc(1));

% SOROKINIANA
[v_a2(:,sample_comp), residuals_a2(:,sample_comp), Chi_a2(sample_comp), df_a2, EXITFLAG_a2, H_a2] = regression_T40(S_t0_a2,S_t5_a2,S_t10_a2,S_t20_a2, Data_mean_abs(6:10,:), ...
    Data_std_abs(6,:),Data_std_abs(7,:), Data_std_abs(8,:), Data_std_abs(9,:),Data_std_abs(10,:),...
    rubp_MIDs_used, f6p_MIDs_used, g6p_MIDs_used, fbp_MIDs_used, g1p_MIDs_used, adpg_MIDs_used, udpg_MIDs_used, model,...
    sta(2),d_sta(2),suc(2),d_suc(2));

% OHADII
% [v_a3_t40(:,sample_comp), residuals_a3_t40(:,sample_comp), Chi_a3_t40(sample_comp), df_a3_t40, EXITFLAG_a3_t40, H_a3_t40] = regression_T40(S_t0_a3,S_t5_a3,S_t10_a3,S_t20_a3, Data_mean_abs(11:15,:), ...
%     Data_std_abs(11,:), Data_std_abs(12,:), Data_std_abs(13,:), Data_std_abs(14,:),Data_std_abs(15,:),...
%     rubp_MIDs_used, f6p_MIDs_used, g6p_MIDs_used, fbp_MIDs_used, g1p_MIDs_used, adpg_MIDs_used, udpg_MIDs_used, model,...
%     sta(3),d_sta(3),suc(3),d_suc(3));

[v_a3(:,sample_comp), residuals_a3(:,sample_comp), Chi_a3(sample_comp), df_a3, EXITFLAG_a3, H_a3] = regression_T20(S_t0_a3,S_t5_a3,S_t10_a3,S_t20_a3, Data_mean_abs(11:15,:), ...
    Data_std_abs(11,:), Data_std_abs(12,:), Data_std_abs(13,:), Data_std_abs(14,:),Data_std_abs(15,:),...
    rubp_MIDs_used, f6p_MIDs_used, g6p_MIDs_used, fbp_MIDs_used, g1p_MIDs_used, adpg_MIDs_used, udpg_MIDs_used, model,...
    sta(3),d_sta(3),suc(3),d_suc(3));

% OHADII EIL
[v_a4(:,sample_comp), residuals_a4(:,sample_comp), Chi_a4(sample_comp), df_a4, EXITFLAG_a4, H_a4] = regression_T40(S_t0_a4,S_t5_a4,S_t10_a4,S_t20_a4, Data_mean_abs(16:20,:), ...
    Data_std_abs(16,:), Data_std_abs(17,:), Data_std_abs(18,:), Data_std_abs(19,:),Data_std_abs(20,:),...
    rubp_MIDs_used, f6p_MIDs_used, g6p_MIDs_used, fbp_MIDs_used, g1p_MIDs_used, adpg_MIDs_used, udpg_MIDs_used, model,...
    sta(4),d_sta(4),suc(4),d_suc(4));

end

%% find best c parameterization and do bootstrap for chlamy, sorokiniana, ohadii LL with tp 0-20sec
[~,order_chi_a1] = sort(Chi_a1);
[~,order_chi_a2] = sort(Chi_a2);
[~,order_chi_a3] = sort(Chi_a3);
[~,order_chi_a4] = sort(Chi_a4);

v_a1=v_a1(:,order_chi_a1(1)); residuals_a1=residuals_a1(:,order_chi_a1(1)); Chi_a1=Chi_a1(order_chi_a1(1)); comp_a1=comp_a1(:,order_chi_a1(1));
v_a2=v_a2(:,order_chi_a2(1)); residuals_a2=residuals_a2(:,order_chi_a2(1)); Chi_a2=Chi_a2(order_chi_a2(1)); comp_a2=comp_a2(:,order_chi_a2(1));
v_a3=v_a3(:,order_chi_a3(1)); residuals_a3=residuals_a3(:,order_chi_a3(1)); Chi_a3=Chi_a3(order_chi_a3(1)); comp_a3=comp_a3(:,order_chi_a3(1));
v_a4=v_a4(:,order_chi_a4(1)); residuals_a4=residuals_a4(:,order_chi_a4(1)); Chi_a4=Chi_a4(order_chi_a4(1)); comp_a4=comp_a4(:,order_chi_a4(1));

disp('fit chlamy:');...
disp(Chi_a1<chi2inv(0.975,df_a1))
disp('fit soro:');...
disp(Chi_a2<chi2inv(0.975,df_a2))
disp('fit ohadii LL 0-20sec:');...
disp(Chi_a3<chi2inv(0.975,df_a3))
disp('fit ohadii EIL:');...
disp(Chi_a4<chi2inv(0.975,df_a4))

% since no chi-square was significant for ohadii under eil we commented out
% that part in the bootstrap
% bootstrap_non_parametric

