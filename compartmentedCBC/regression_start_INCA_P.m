clear
close all
addpath('..')

% load experimental data
% MID absolute and relative (mean, std)
% pool size measurements (mean, std)
[~, ~, Data_mean_norm, Data_std_norm, ~] = load_experimental_data;

%% specify idx of used MIDs
rubp_MIDs_used = find(contains(Data_mean_norm.Properties.VariableNames,'RuBP'));
fbp_MIDs_used = find(contains(Data_mean_norm.Properties.VariableNames,'FBP'));
f6p_MIDs_used = find(contains(Data_mean_norm.Properties.VariableNames,'f6p'));
g6p_MIDs_used = find(contains(Data_mean_norm.Properties.VariableNames,'G6P'));
g1p_MIDs_used = find(contains(Data_mean_norm.Properties.VariableNames,'g1p'));
adpg_MIDs_used = find(contains(Data_mean_norm.Properties.VariableNames,'ADPG'));
udpg_MIDs_used = find(contains(Data_mean_norm.Properties.VariableNames,'UDPG'));
MIDs_used_all=[rubp_MIDs_used, f6p_MIDs_used, g6p_MIDs_used, fbp_MIDs_used, g1p_MIDs_used, adpg_MIDs_used, udpg_MIDs_used];

INCA_P = load('../INCA_compCBC/AbsoluteLevelINCA.mat');
AbsoluteLevel = INCA_P.AbsoluteLevel_mean;

R_rubp=repmat(AbsoluteLevel.RuBP,1,5)';
R_f6p=repmat(AbsoluteLevel.F6P,1,5)';
R_g6p=repmat(AbsoluteLevel.G6P,1,5)';
R_fbp=repmat(AbsoluteLevel.FBP,1,5)';
R_g1p=repmat(AbsoluteLevel.G1P,1,5)';
R_adpg=repmat(AbsoluteLevel.ADPG_Levels_nmol_gDW__,1,5)';
R_udpg=repmat(AbsoluteLevel.UDPG,1,5)';

R=[repmat(R_rubp(:),1,length(rubp_MIDs_used)) repmat(R_f6p(:),1,length(f6p_MIDs_used)) ...
    repmat(R_g6p(:),1,length(g6p_MIDs_used)) repmat(R_fbp(:),1,length(fbp_MIDs_used)) ...
    repmat(R_g1p(:),1,length(g1p_MIDs_used)) repmat(R_adpg(:),1,length(adpg_MIDs_used)) repmat(R_udpg(:),1,length(udpg_MIDs_used))];

Data_mean_abs = Data_mean_norm;
Data_mean_abs{:,MIDs_used_all} = R.*Data_mean_norm{:,MIDs_used_all};

Data_std_abs = Data_std_norm;
Data_std_abs{:,MIDs_used_all} = R.*Data_std_norm{:,MIDs_used_all};

algea_c=1:5; %chlamy
algea_s=6:10; %sorokiniana
algea_o=11:15; %ohadii

network_compCBC % load network structure
model_rev=model; % split model into irreversible but also keep the original model in variable model_rev
model=convertToIrreversible(model);
model.lb(1) = 1e-3; % TCO2 should be active

for sample_comp=1:50 % sample 50 parameterization of c and calculate fit, take the best fit and do bootstrap

%% create matrix that combines the MID data
%chlamy
[S_t0_a1,S_t5_a1,S_t10_a1,S_t20_a1,comp_a1(:,sample_comp)] = build_differential_matrix_comp(Data_mean_norm(algea_c,:),model);
%sorokiniana
[S_t0_a2,S_t5_a2,S_t10_a2,S_t20_a2,comp_a2(:,sample_comp)] = build_differential_matrix_comp(Data_mean_norm(algea_s,:),model);
%ohadii
[S_t0_a3,S_t5_a3,S_t10_a3,S_t20_a3,comp_a3(:,sample_comp)] = build_differential_matrix_comp(Data_mean_norm(algea_o,:),model);

%% vectors including measured starch and sucrose synthesis rates
sta = [19 51 85]; % mean 
d_sta = sta*0.3; % std
suc = [1e-7 17 35];
d_suc = suc.*0.3;

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
[v_a3(:,sample_comp), residuals_a3(:,sample_comp), Chi_a3(sample_comp), df_a3, EXITFLAG_a3, H_a3] = regression_T20(S_t0_a3,S_t5_a3,S_t10_a3,S_t20_a3, Data_mean_abs(11:15,:), ...
    Data_std_abs(11,:), Data_std_abs(12,:), Data_std_abs(13,:), Data_std_abs(14,:),Data_std_abs(15,:),...
    rubp_MIDs_used, f6p_MIDs_used, g6p_MIDs_used, fbp_MIDs_used, g1p_MIDs_used, adpg_MIDs_used, udpg_MIDs_used, model,...
    sta(3),d_sta(3),suc(3),d_suc(3));

end

%% find best c parameterization and do bootstrap for chlamy, sorokiniana, ohadii LL with tp 0-20sec
[~,order_chi_a1] = sort(Chi_a1);
[~,order_chi_a2] = sort(Chi_a2);
[~,order_chi_a3] = sort(Chi_a3);

v_a1=v_a1(:,order_chi_a1(1)); residuals_a1=residuals_a1(:,order_chi_a1(1)); Chi_a1=Chi_a1(order_chi_a1(1)); comp_a1=comp_a1(:,order_chi_a1(1));
v_a2=v_a2(:,order_chi_a2(1)); residuals_a2=residuals_a2(:,order_chi_a2(1)); Chi_a2=Chi_a2(order_chi_a2(1)); comp_a2=comp_a2(:,order_chi_a2(1));
v_a3=v_a3(:,order_chi_a3(1)); residuals_a3=residuals_a3(:,order_chi_a3(1)); Chi_a3=Chi_a3(order_chi_a3(1)); comp_a3=comp_a3(:,order_chi_a3(1));

disp('fit chlamy:');...
disp(Chi_a1<chi2inv(0.975,df_a1))
disp('fit soro:');...
disp(Chi_a2<chi2inv(0.975,df_a2))
disp('fit ohadii LL 0-20sec:');...
disp(Chi_a3<chi2inv(0.975,df_a3))

bootstrap_non_parametric

