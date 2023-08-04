clear
close all
addpath('..')

% load experimental data
% MID absolute and relative (mean, std)
% pool size measurements (mean, std)
 [Data_mean_abs, Data_std_abs, Data_mean_norm, Data_std_norm,AbsoluteLevel,AbsoluteLevel_std,Data_std_abs_original] = load_experimental_data; % run with original time series
%[Data_mean_abs, Data_std_abs, Data_mean_norm, Data_std_norm, AbsoluteLevel,AbsoluteLevel_std,Data_std_abs_original] = load_experimental_data(false,true); % run with dense time series obtained from data interpolation

algea_c=find(strcmp(Data_mean_abs.Alga,'Chlamy'));
algea_s=find(strcmp(Data_mean_abs.Alga,'Sorokiniana'));
algea_o=find(strcmp(Data_mean_abs.Alga,'Ohadii') & strcmp(Data_mean_abs.Light,'100 uE'));
algea_o_eil=find(strcmp(Data_mean_abs.Alga,'Ohadii') & strcmp(Data_mean_abs.Light,'3000 uE'));

network_compCBC % load network structure
model_rev=model;
model=convertToIrreversible(model);
model.lb(1)=1e-3; % TCO2 should be active

%% specify idx of used MIDs
rubp_MIDs_used = find(contains(Data_mean_abs.Properties.VariableNames,'RuBP'));
fbp_MIDs_used = find(contains(Data_mean_abs.Properties.VariableNames,'FBP'));
f6p_MIDs_used = find(contains(Data_mean_abs.Properties.VariableNames,'f6p'));
g6p_MIDs_used = find(contains(Data_mean_abs.Properties.VariableNames,'G6P'));
g1p_MIDs_used = find(contains(Data_mean_abs.Properties.VariableNames,'g1p'));
adpg_MIDs_used = find(contains(Data_mean_abs.Properties.VariableNames,'ADPG'));
udpg_MIDs_used = find(contains(Data_mean_abs.Properties.VariableNames,'UDPG'));

%% vectors including measured starch and sucrose synthesis rates
sta = [19 51 85 74]; % mean
d_sta = sta*0.3; % std
suc = [1e-5 17 35 61];
d_suc = suc.*0.3;

MIDs_used_all=[rubp_MIDs_used, f6p_MIDs_used, g6p_MIDs_used, fbp_MIDs_used, g1p_MIDs_used, adpg_MIDs_used, udpg_MIDs_used];

PMEAN_a1 = [eye(length(rubp_MIDs_used))*AbsoluteLevel.RuBP(1)             zeros(length(rubp_MIDs_used),length(f6p_MIDs_used))           zeros(length(rubp_MIDs_used),length(g6p_MIDs_used))   zeros(length(rubp_MIDs_used),length(fbp_MIDs_used)) zeros(length(rubp_MIDs_used),length(g1p_MIDs_used)) zeros(length(rubp_MIDs_used),length(adpg_MIDs_used)) zeros(length(rubp_MIDs_used),length(udpg_MIDs_used));
    zeros(length(f6p_MIDs_used),length(rubp_MIDs_used))   eye(length(f6p_MIDs_used))*AbsoluteLevel.F6P(1)                          zeros(length(f6p_MIDs_used),length(g6p_MIDs_used))    zeros(length(f6p_MIDs_used),length(fbp_MIDs_used))  zeros(length(f6p_MIDs_used),length(g1p_MIDs_used))  zeros(length(f6p_MIDs_used),length(adpg_MIDs_used))  zeros(length(f6p_MIDs_used),length(udpg_MIDs_used)) ;
    zeros(length(g6p_MIDs_used),length(rubp_MIDs_used))   zeros(length(g6p_MIDs_used),length(f6p_MIDs_used))            eye(length(g6p_MIDs_used))*AbsoluteLevel.G6P(1)                  zeros(length(g6p_MIDs_used),length(fbp_MIDs_used))  zeros(length(g6p_MIDs_used),length(g1p_MIDs_used))  zeros(length(g6p_MIDs_used),length(adpg_MIDs_used))  zeros(length(g6p_MIDs_used),length(udpg_MIDs_used)) ;
    zeros(length(fbp_MIDs_used),length(rubp_MIDs_used))   zeros(length(fbp_MIDs_used),length(f6p_MIDs_used))            zeros(length(fbp_MIDs_used),length(g6p_MIDs_used))    eye(length(fbp_MIDs_used))*AbsoluteLevel.FBP(1)                zeros(length(fbp_MIDs_used),length(g1p_MIDs_used))  zeros(length(fbp_MIDs_used),length(adpg_MIDs_used))  zeros(length(fbp_MIDs_used),length(udpg_MIDs_used)) ;
    zeros(length(g1p_MIDs_used),length(rubp_MIDs_used))   zeros(length(g1p_MIDs_used),length(f6p_MIDs_used))            zeros(length(g1p_MIDs_used),length(g6p_MIDs_used))    zeros(length(g1p_MIDs_used),length(fbp_MIDs_used))  eye(length(g1p_MIDs_used))*AbsoluteLevel.G1P(1)                zeros(length(g1p_MIDs_used),length(adpg_MIDs_used))  zeros(length(g1p_MIDs_used),length(udpg_MIDs_used)) ;
    zeros(length(adpg_MIDs_used),length(rubp_MIDs_used))  zeros(length(adpg_MIDs_used),length(f6p_MIDs_used))           zeros(length(adpg_MIDs_used),length(g6p_MIDs_used))   zeros(length(adpg_MIDs_used),length(fbp_MIDs_used)) zeros(length(adpg_MIDs_used),length(g1p_MIDs_used)) eye(length(adpg_MIDs_used))*AbsoluteLevel.ADPG_Levels_nmol_gDW__(1)              zeros(length(adpg_MIDs_used),length(udpg_MIDs_used));
    zeros(length(udpg_MIDs_used),length(rubp_MIDs_used))  zeros(length(udpg_MIDs_used),length(f6p_MIDs_used))           zeros(length(udpg_MIDs_used),length(g6p_MIDs_used))   zeros(length(udpg_MIDs_used),length(fbp_MIDs_used)) zeros(length(udpg_MIDs_used),length(g1p_MIDs_used)) zeros(length(udpg_MIDs_used),length(adpg_MIDs_used)) eye(length(udpg_MIDs_used))*AbsoluteLevel.UDPG(1)             ];

PVar_a1 = [eye(length(rubp_MIDs_used))*AbsoluteLevel_std.RuBP(1)             zeros(length(rubp_MIDs_used),length(f6p_MIDs_used))           zeros(length(rubp_MIDs_used),length(g6p_MIDs_used))   zeros(length(rubp_MIDs_used),length(fbp_MIDs_used)) zeros(length(rubp_MIDs_used),length(g1p_MIDs_used)) zeros(length(rubp_MIDs_used),length(adpg_MIDs_used)) zeros(length(rubp_MIDs_used),length(udpg_MIDs_used));
    zeros(length(f6p_MIDs_used),length(rubp_MIDs_used))   eye(length(f6p_MIDs_used))*AbsoluteLevel_std.F6P(1)                          zeros(length(f6p_MIDs_used),length(g6p_MIDs_used))    zeros(length(f6p_MIDs_used),length(fbp_MIDs_used))  zeros(length(f6p_MIDs_used),length(g1p_MIDs_used))  zeros(length(f6p_MIDs_used),length(adpg_MIDs_used))  zeros(length(f6p_MIDs_used),length(udpg_MIDs_used)) ;
    zeros(length(g6p_MIDs_used),length(rubp_MIDs_used))   zeros(length(g6p_MIDs_used),length(f6p_MIDs_used))            eye(length(g6p_MIDs_used))*AbsoluteLevel_std.G6P(1)                  zeros(length(g6p_MIDs_used),length(fbp_MIDs_used))  zeros(length(g6p_MIDs_used),length(g1p_MIDs_used))  zeros(length(g6p_MIDs_used),length(adpg_MIDs_used))  zeros(length(g6p_MIDs_used),length(udpg_MIDs_used)) ;
    zeros(length(fbp_MIDs_used),length(rubp_MIDs_used))   zeros(length(fbp_MIDs_used),length(f6p_MIDs_used))            zeros(length(fbp_MIDs_used),length(g6p_MIDs_used))    eye(length(fbp_MIDs_used))*AbsoluteLevel_std.FBP(1)                zeros(length(fbp_MIDs_used),length(g1p_MIDs_used))  zeros(length(fbp_MIDs_used),length(adpg_MIDs_used))  zeros(length(fbp_MIDs_used),length(udpg_MIDs_used)) ;
    zeros(length(g1p_MIDs_used),length(rubp_MIDs_used))   zeros(length(g1p_MIDs_used),length(f6p_MIDs_used))            zeros(length(g1p_MIDs_used),length(g6p_MIDs_used))    zeros(length(g1p_MIDs_used),length(fbp_MIDs_used))  eye(length(g1p_MIDs_used))*AbsoluteLevel_std.G1P(1)                zeros(length(g1p_MIDs_used),length(adpg_MIDs_used))  zeros(length(g1p_MIDs_used),length(udpg_MIDs_used)) ;
    zeros(length(adpg_MIDs_used),length(rubp_MIDs_used))  zeros(length(adpg_MIDs_used),length(f6p_MIDs_used))           zeros(length(adpg_MIDs_used),length(g6p_MIDs_used))   zeros(length(adpg_MIDs_used),length(fbp_MIDs_used)) zeros(length(adpg_MIDs_used),length(g1p_MIDs_used)) eye(length(adpg_MIDs_used))*AbsoluteLevel_std.ADPG_Levels_nmol_gDW__(1)              zeros(length(adpg_MIDs_used),length(udpg_MIDs_used));
    zeros(length(udpg_MIDs_used),length(rubp_MIDs_used))  zeros(length(udpg_MIDs_used),length(f6p_MIDs_used))           zeros(length(udpg_MIDs_used),length(g6p_MIDs_used))   zeros(length(udpg_MIDs_used),length(fbp_MIDs_used)) zeros(length(udpg_MIDs_used),length(g1p_MIDs_used)) zeros(length(udpg_MIDs_used),length(adpg_MIDs_used)) eye(length(udpg_MIDs_used))*AbsoluteLevel_std.UDPG(1)             ];

PMEAN_a2 = [eye(length(rubp_MIDs_used))*AbsoluteLevel.RuBP(2)             zeros(length(rubp_MIDs_used),length(f6p_MIDs_used))           zeros(length(rubp_MIDs_used),length(g6p_MIDs_used))   zeros(length(rubp_MIDs_used),length(fbp_MIDs_used)) zeros(length(rubp_MIDs_used),length(g1p_MIDs_used)) zeros(length(rubp_MIDs_used),length(adpg_MIDs_used)) zeros(length(rubp_MIDs_used),length(udpg_MIDs_used));
    zeros(length(f6p_MIDs_used),length(rubp_MIDs_used))   eye(length(f6p_MIDs_used))*AbsoluteLevel.F6P(2)                          zeros(length(f6p_MIDs_used),length(g6p_MIDs_used))    zeros(length(f6p_MIDs_used),length(fbp_MIDs_used))  zeros(length(f6p_MIDs_used),length(g1p_MIDs_used))  zeros(length(f6p_MIDs_used),length(adpg_MIDs_used))  zeros(length(f6p_MIDs_used),length(udpg_MIDs_used)) ;
    zeros(length(g6p_MIDs_used),length(rubp_MIDs_used))   zeros(length(g6p_MIDs_used),length(f6p_MIDs_used))            eye(length(g6p_MIDs_used))*AbsoluteLevel.G6P(2)                  zeros(length(g6p_MIDs_used),length(fbp_MIDs_used))  zeros(length(g6p_MIDs_used),length(g1p_MIDs_used))  zeros(length(g6p_MIDs_used),length(adpg_MIDs_used))  zeros(length(g6p_MIDs_used),length(udpg_MIDs_used)) ;
    zeros(length(fbp_MIDs_used),length(rubp_MIDs_used))   zeros(length(fbp_MIDs_used),length(f6p_MIDs_used))            zeros(length(fbp_MIDs_used),length(g6p_MIDs_used))    eye(length(fbp_MIDs_used))*AbsoluteLevel.FBP(2)                zeros(length(fbp_MIDs_used),length(g1p_MIDs_used))  zeros(length(fbp_MIDs_used),length(adpg_MIDs_used))  zeros(length(fbp_MIDs_used),length(udpg_MIDs_used)) ;
    zeros(length(g1p_MIDs_used),length(rubp_MIDs_used))   zeros(length(g1p_MIDs_used),length(f6p_MIDs_used))            zeros(length(g1p_MIDs_used),length(g6p_MIDs_used))    zeros(length(g1p_MIDs_used),length(fbp_MIDs_used))  eye(length(g1p_MIDs_used))*AbsoluteLevel.G1P(2)                zeros(length(g1p_MIDs_used),length(adpg_MIDs_used))  zeros(length(g1p_MIDs_used),length(udpg_MIDs_used)) ;
    zeros(length(adpg_MIDs_used),length(rubp_MIDs_used))  zeros(length(adpg_MIDs_used),length(f6p_MIDs_used))           zeros(length(adpg_MIDs_used),length(g6p_MIDs_used))   zeros(length(adpg_MIDs_used),length(fbp_MIDs_used)) zeros(length(adpg_MIDs_used),length(g1p_MIDs_used)) eye(length(adpg_MIDs_used))*AbsoluteLevel.ADPG_Levels_nmol_gDW__(2)              zeros(length(adpg_MIDs_used),length(udpg_MIDs_used));
    zeros(length(udpg_MIDs_used),length(rubp_MIDs_used))  zeros(length(udpg_MIDs_used),length(f6p_MIDs_used))           zeros(length(udpg_MIDs_used),length(g6p_MIDs_used))   zeros(length(udpg_MIDs_used),length(fbp_MIDs_used)) zeros(length(udpg_MIDs_used),length(g1p_MIDs_used)) zeros(length(udpg_MIDs_used),length(adpg_MIDs_used)) eye(length(udpg_MIDs_used))*AbsoluteLevel.UDPG(2)             ];

PVar_a2 = [eye(length(rubp_MIDs_used))*AbsoluteLevel_std.RuBP(2)             zeros(length(rubp_MIDs_used),length(f6p_MIDs_used))           zeros(length(rubp_MIDs_used),length(g6p_MIDs_used))   zeros(length(rubp_MIDs_used),length(fbp_MIDs_used)) zeros(length(rubp_MIDs_used),length(g1p_MIDs_used)) zeros(length(rubp_MIDs_used),length(adpg_MIDs_used)) zeros(length(rubp_MIDs_used),length(udpg_MIDs_used));
    zeros(length(f6p_MIDs_used),length(rubp_MIDs_used))   eye(length(f6p_MIDs_used))*AbsoluteLevel_std.F6P(2)                          zeros(length(f6p_MIDs_used),length(g6p_MIDs_used))    zeros(length(f6p_MIDs_used),length(fbp_MIDs_used))  zeros(length(f6p_MIDs_used),length(g1p_MIDs_used))  zeros(length(f6p_MIDs_used),length(adpg_MIDs_used))  zeros(length(f6p_MIDs_used),length(udpg_MIDs_used)) ;
    zeros(length(g6p_MIDs_used),length(rubp_MIDs_used))   zeros(length(g6p_MIDs_used),length(f6p_MIDs_used))            eye(length(g6p_MIDs_used))*AbsoluteLevel_std.G6P(2)                  zeros(length(g6p_MIDs_used),length(fbp_MIDs_used))  zeros(length(g6p_MIDs_used),length(g1p_MIDs_used))  zeros(length(g6p_MIDs_used),length(adpg_MIDs_used))  zeros(length(g6p_MIDs_used),length(udpg_MIDs_used)) ;
    zeros(length(fbp_MIDs_used),length(rubp_MIDs_used))   zeros(length(fbp_MIDs_used),length(f6p_MIDs_used))            zeros(length(fbp_MIDs_used),length(g6p_MIDs_used))    eye(length(fbp_MIDs_used))*AbsoluteLevel_std.FBP(2)                zeros(length(fbp_MIDs_used),length(g1p_MIDs_used))  zeros(length(fbp_MIDs_used),length(adpg_MIDs_used))  zeros(length(fbp_MIDs_used),length(udpg_MIDs_used)) ;
    zeros(length(g1p_MIDs_used),length(rubp_MIDs_used))   zeros(length(g1p_MIDs_used),length(f6p_MIDs_used))            zeros(length(g1p_MIDs_used),length(g6p_MIDs_used))    zeros(length(g1p_MIDs_used),length(fbp_MIDs_used))  eye(length(g1p_MIDs_used))*AbsoluteLevel_std.G1P(2)                zeros(length(g1p_MIDs_used),length(adpg_MIDs_used))  zeros(length(g1p_MIDs_used),length(udpg_MIDs_used)) ;
    zeros(length(adpg_MIDs_used),length(rubp_MIDs_used))  zeros(length(adpg_MIDs_used),length(f6p_MIDs_used))           zeros(length(adpg_MIDs_used),length(g6p_MIDs_used))   zeros(length(adpg_MIDs_used),length(fbp_MIDs_used)) zeros(length(adpg_MIDs_used),length(g1p_MIDs_used)) eye(length(adpg_MIDs_used))*AbsoluteLevel_std.ADPG_Levels_nmol_gDW__(2)              zeros(length(adpg_MIDs_used),length(udpg_MIDs_used));
    zeros(length(udpg_MIDs_used),length(rubp_MIDs_used))  zeros(length(udpg_MIDs_used),length(f6p_MIDs_used))           zeros(length(udpg_MIDs_used),length(g6p_MIDs_used))   zeros(length(udpg_MIDs_used),length(fbp_MIDs_used)) zeros(length(udpg_MIDs_used),length(g1p_MIDs_used)) zeros(length(udpg_MIDs_used),length(adpg_MIDs_used)) eye(length(udpg_MIDs_used))*AbsoluteLevel_std.UDPG(2)             ];

PMEAN_a3 = [eye(length(rubp_MIDs_used))*AbsoluteLevel.RuBP(3)             zeros(length(rubp_MIDs_used),length(f6p_MIDs_used))           zeros(length(rubp_MIDs_used),length(g6p_MIDs_used))   zeros(length(rubp_MIDs_used),length(fbp_MIDs_used)) zeros(length(rubp_MIDs_used),length(g1p_MIDs_used)) zeros(length(rubp_MIDs_used),length(adpg_MIDs_used)) zeros(length(rubp_MIDs_used),length(udpg_MIDs_used));
    zeros(length(f6p_MIDs_used),length(rubp_MIDs_used))   eye(length(f6p_MIDs_used))*AbsoluteLevel.F6P(3)                          zeros(length(f6p_MIDs_used),length(g6p_MIDs_used))    zeros(length(f6p_MIDs_used),length(fbp_MIDs_used))  zeros(length(f6p_MIDs_used),length(g1p_MIDs_used))  zeros(length(f6p_MIDs_used),length(adpg_MIDs_used))  zeros(length(f6p_MIDs_used),length(udpg_MIDs_used)) ;
    zeros(length(g6p_MIDs_used),length(rubp_MIDs_used))   zeros(length(g6p_MIDs_used),length(f6p_MIDs_used))            eye(length(g6p_MIDs_used))*AbsoluteLevel.G6P(3)                  zeros(length(g6p_MIDs_used),length(fbp_MIDs_used))  zeros(length(g6p_MIDs_used),length(g1p_MIDs_used))  zeros(length(g6p_MIDs_used),length(adpg_MIDs_used))  zeros(length(g6p_MIDs_used),length(udpg_MIDs_used)) ;
    zeros(length(fbp_MIDs_used),length(rubp_MIDs_used))   zeros(length(fbp_MIDs_used),length(f6p_MIDs_used))            zeros(length(fbp_MIDs_used),length(g6p_MIDs_used))    eye(length(fbp_MIDs_used))*AbsoluteLevel.FBP(3)                zeros(length(fbp_MIDs_used),length(g1p_MIDs_used))  zeros(length(fbp_MIDs_used),length(adpg_MIDs_used))  zeros(length(fbp_MIDs_used),length(udpg_MIDs_used)) ;
    zeros(length(g1p_MIDs_used),length(rubp_MIDs_used))   zeros(length(g1p_MIDs_used),length(f6p_MIDs_used))            zeros(length(g1p_MIDs_used),length(g6p_MIDs_used))    zeros(length(g1p_MIDs_used),length(fbp_MIDs_used))  eye(length(g1p_MIDs_used))*AbsoluteLevel.G1P(3)                zeros(length(g1p_MIDs_used),length(adpg_MIDs_used))  zeros(length(g1p_MIDs_used),length(udpg_MIDs_used)) ;
    zeros(length(adpg_MIDs_used),length(rubp_MIDs_used))  zeros(length(adpg_MIDs_used),length(f6p_MIDs_used))           zeros(length(adpg_MIDs_used),length(g6p_MIDs_used))   zeros(length(adpg_MIDs_used),length(fbp_MIDs_used)) zeros(length(adpg_MIDs_used),length(g1p_MIDs_used)) eye(length(adpg_MIDs_used))*AbsoluteLevel.ADPG_Levels_nmol_gDW__(3)              zeros(length(adpg_MIDs_used),length(udpg_MIDs_used));
    zeros(length(udpg_MIDs_used),length(rubp_MIDs_used))  zeros(length(udpg_MIDs_used),length(f6p_MIDs_used))           zeros(length(udpg_MIDs_used),length(g6p_MIDs_used))   zeros(length(udpg_MIDs_used),length(fbp_MIDs_used)) zeros(length(udpg_MIDs_used),length(g1p_MIDs_used)) zeros(length(udpg_MIDs_used),length(adpg_MIDs_used)) eye(length(udpg_MIDs_used))*AbsoluteLevel.UDPG(3)             ];

PVar_a3 = [eye(length(rubp_MIDs_used))*AbsoluteLevel_std.RuBP(3)             zeros(length(rubp_MIDs_used),length(f6p_MIDs_used))           zeros(length(rubp_MIDs_used),length(g6p_MIDs_used))   zeros(length(rubp_MIDs_used),length(fbp_MIDs_used)) zeros(length(rubp_MIDs_used),length(g1p_MIDs_used)) zeros(length(rubp_MIDs_used),length(adpg_MIDs_used)) zeros(length(rubp_MIDs_used),length(udpg_MIDs_used));
    zeros(length(f6p_MIDs_used),length(rubp_MIDs_used))   eye(length(f6p_MIDs_used))*AbsoluteLevel_std.F6P(3)                          zeros(length(f6p_MIDs_used),length(g6p_MIDs_used))    zeros(length(f6p_MIDs_used),length(fbp_MIDs_used))  zeros(length(f6p_MIDs_used),length(g1p_MIDs_used))  zeros(length(f6p_MIDs_used),length(adpg_MIDs_used))  zeros(length(f6p_MIDs_used),length(udpg_MIDs_used)) ;
    zeros(length(g6p_MIDs_used),length(rubp_MIDs_used))   zeros(length(g6p_MIDs_used),length(f6p_MIDs_used))            eye(length(g6p_MIDs_used))*AbsoluteLevel_std.G6P(3)                  zeros(length(g6p_MIDs_used),length(fbp_MIDs_used))  zeros(length(g6p_MIDs_used),length(g1p_MIDs_used))  zeros(length(g6p_MIDs_used),length(adpg_MIDs_used))  zeros(length(g6p_MIDs_used),length(udpg_MIDs_used)) ;
    zeros(length(fbp_MIDs_used),length(rubp_MIDs_used))   zeros(length(fbp_MIDs_used),length(f6p_MIDs_used))            zeros(length(fbp_MIDs_used),length(g6p_MIDs_used))    eye(length(fbp_MIDs_used))*AbsoluteLevel_std.FBP(3)                zeros(length(fbp_MIDs_used),length(g1p_MIDs_used))  zeros(length(fbp_MIDs_used),length(adpg_MIDs_used))  zeros(length(fbp_MIDs_used),length(udpg_MIDs_used)) ;
    zeros(length(g1p_MIDs_used),length(rubp_MIDs_used))   zeros(length(g1p_MIDs_used),length(f6p_MIDs_used))            zeros(length(g1p_MIDs_used),length(g6p_MIDs_used))    zeros(length(g1p_MIDs_used),length(fbp_MIDs_used))  eye(length(g1p_MIDs_used))*AbsoluteLevel_std.G1P(3)                zeros(length(g1p_MIDs_used),length(adpg_MIDs_used))  zeros(length(g1p_MIDs_used),length(udpg_MIDs_used)) ;
    zeros(length(adpg_MIDs_used),length(rubp_MIDs_used))  zeros(length(adpg_MIDs_used),length(f6p_MIDs_used))           zeros(length(adpg_MIDs_used),length(g6p_MIDs_used))   zeros(length(adpg_MIDs_used),length(fbp_MIDs_used)) zeros(length(adpg_MIDs_used),length(g1p_MIDs_used)) eye(length(adpg_MIDs_used))*AbsoluteLevel_std.ADPG_Levels_nmol_gDW__(3)              zeros(length(adpg_MIDs_used),length(udpg_MIDs_used));
    zeros(length(udpg_MIDs_used),length(rubp_MIDs_used))  zeros(length(udpg_MIDs_used),length(f6p_MIDs_used))           zeros(length(udpg_MIDs_used),length(g6p_MIDs_used))   zeros(length(udpg_MIDs_used),length(fbp_MIDs_used)) zeros(length(udpg_MIDs_used),length(g1p_MIDs_used)) zeros(length(udpg_MIDs_used),length(adpg_MIDs_used)) eye(length(udpg_MIDs_used))*AbsoluteLevel_std.UDPG(3)             ];

PMEAN_a4 = [eye(length(rubp_MIDs_used))*AbsoluteLevel.RuBP(4)             zeros(length(rubp_MIDs_used),length(f6p_MIDs_used))           zeros(length(rubp_MIDs_used),length(g6p_MIDs_used))   zeros(length(rubp_MIDs_used),length(fbp_MIDs_used)) zeros(length(rubp_MIDs_used),length(g1p_MIDs_used)) zeros(length(rubp_MIDs_used),length(adpg_MIDs_used)) zeros(length(rubp_MIDs_used),length(udpg_MIDs_used));
    zeros(length(f6p_MIDs_used),length(rubp_MIDs_used))   eye(length(f6p_MIDs_used))*AbsoluteLevel.F6P(4)                          zeros(length(f6p_MIDs_used),length(g6p_MIDs_used))    zeros(length(f6p_MIDs_used),length(fbp_MIDs_used))  zeros(length(f6p_MIDs_used),length(g1p_MIDs_used))  zeros(length(f6p_MIDs_used),length(adpg_MIDs_used))  zeros(length(f6p_MIDs_used),length(udpg_MIDs_used)) ;
    zeros(length(g6p_MIDs_used),length(rubp_MIDs_used))   zeros(length(g6p_MIDs_used),length(f6p_MIDs_used))            eye(length(g6p_MIDs_used))*AbsoluteLevel.G6P(4)                  zeros(length(g6p_MIDs_used),length(fbp_MIDs_used))  zeros(length(g6p_MIDs_used),length(g1p_MIDs_used))  zeros(length(g6p_MIDs_used),length(adpg_MIDs_used))  zeros(length(g6p_MIDs_used),length(udpg_MIDs_used)) ;
    zeros(length(fbp_MIDs_used),length(rubp_MIDs_used))   zeros(length(fbp_MIDs_used),length(f6p_MIDs_used))            zeros(length(fbp_MIDs_used),length(g6p_MIDs_used))    eye(length(fbp_MIDs_used))*AbsoluteLevel.FBP(4)                zeros(length(fbp_MIDs_used),length(g1p_MIDs_used))  zeros(length(fbp_MIDs_used),length(adpg_MIDs_used))  zeros(length(fbp_MIDs_used),length(udpg_MIDs_used)) ;
    zeros(length(g1p_MIDs_used),length(rubp_MIDs_used))   zeros(length(g1p_MIDs_used),length(f6p_MIDs_used))            zeros(length(g1p_MIDs_used),length(g6p_MIDs_used))    zeros(length(g1p_MIDs_used),length(fbp_MIDs_used))  eye(length(g1p_MIDs_used))*AbsoluteLevel.G1P(4)                zeros(length(g1p_MIDs_used),length(adpg_MIDs_used))  zeros(length(g1p_MIDs_used),length(udpg_MIDs_used)) ;
    zeros(length(adpg_MIDs_used),length(rubp_MIDs_used))  zeros(length(adpg_MIDs_used),length(f6p_MIDs_used))           zeros(length(adpg_MIDs_used),length(g6p_MIDs_used))   zeros(length(adpg_MIDs_used),length(fbp_MIDs_used)) zeros(length(adpg_MIDs_used),length(g1p_MIDs_used)) eye(length(adpg_MIDs_used))*AbsoluteLevel.ADPG_Levels_nmol_gDW__(4)              zeros(length(adpg_MIDs_used),length(udpg_MIDs_used));
    zeros(length(udpg_MIDs_used),length(rubp_MIDs_used))  zeros(length(udpg_MIDs_used),length(f6p_MIDs_used))           zeros(length(udpg_MIDs_used),length(g6p_MIDs_used))   zeros(length(udpg_MIDs_used),length(fbp_MIDs_used)) zeros(length(udpg_MIDs_used),length(g1p_MIDs_used)) zeros(length(udpg_MIDs_used),length(adpg_MIDs_used)) eye(length(udpg_MIDs_used))*AbsoluteLevel.UDPG(4)             ];

PVar_a4 = [eye(length(rubp_MIDs_used))*AbsoluteLevel_std.RuBP(4)             zeros(length(rubp_MIDs_used),length(f6p_MIDs_used))           zeros(length(rubp_MIDs_used),length(g6p_MIDs_used))   zeros(length(rubp_MIDs_used),length(fbp_MIDs_used)) zeros(length(rubp_MIDs_used),length(g1p_MIDs_used)) zeros(length(rubp_MIDs_used),length(adpg_MIDs_used)) zeros(length(rubp_MIDs_used),length(udpg_MIDs_used));
    zeros(length(f6p_MIDs_used),length(rubp_MIDs_used))   eye(length(f6p_MIDs_used))*AbsoluteLevel_std.F6P(4)                          zeros(length(f6p_MIDs_used),length(g6p_MIDs_used))    zeros(length(f6p_MIDs_used),length(fbp_MIDs_used))  zeros(length(f6p_MIDs_used),length(g1p_MIDs_used))  zeros(length(f6p_MIDs_used),length(adpg_MIDs_used))  zeros(length(f6p_MIDs_used),length(udpg_MIDs_used)) ;
    zeros(length(g6p_MIDs_used),length(rubp_MIDs_used))   zeros(length(g6p_MIDs_used),length(f6p_MIDs_used))            eye(length(g6p_MIDs_used))*AbsoluteLevel_std.G6P(4)                  zeros(length(g6p_MIDs_used),length(fbp_MIDs_used))  zeros(length(g6p_MIDs_used),length(g1p_MIDs_used))  zeros(length(g6p_MIDs_used),length(adpg_MIDs_used))  zeros(length(g6p_MIDs_used),length(udpg_MIDs_used)) ;
    zeros(length(fbp_MIDs_used),length(rubp_MIDs_used))   zeros(length(fbp_MIDs_used),length(f6p_MIDs_used))            zeros(length(fbp_MIDs_used),length(g6p_MIDs_used))    eye(length(fbp_MIDs_used))*AbsoluteLevel_std.FBP(4)                zeros(length(fbp_MIDs_used),length(g1p_MIDs_used))  zeros(length(fbp_MIDs_used),length(adpg_MIDs_used))  zeros(length(fbp_MIDs_used),length(udpg_MIDs_used)) ;
    zeros(length(g1p_MIDs_used),length(rubp_MIDs_used))   zeros(length(g1p_MIDs_used),length(f6p_MIDs_used))            zeros(length(g1p_MIDs_used),length(g6p_MIDs_used))    zeros(length(g1p_MIDs_used),length(fbp_MIDs_used))  eye(length(g1p_MIDs_used))*AbsoluteLevel_std.G1P(4)                zeros(length(g1p_MIDs_used),length(adpg_MIDs_used))  zeros(length(g1p_MIDs_used),length(udpg_MIDs_used)) ;
    zeros(length(adpg_MIDs_used),length(rubp_MIDs_used))  zeros(length(adpg_MIDs_used),length(f6p_MIDs_used))           zeros(length(adpg_MIDs_used),length(g6p_MIDs_used))   zeros(length(adpg_MIDs_used),length(fbp_MIDs_used)) zeros(length(adpg_MIDs_used),length(g1p_MIDs_used)) eye(length(adpg_MIDs_used))*AbsoluteLevel_std.ADPG_Levels_nmol_gDW__(4)              zeros(length(adpg_MIDs_used),length(udpg_MIDs_used));
    zeros(length(udpg_MIDs_used),length(rubp_MIDs_used))  zeros(length(udpg_MIDs_used),length(f6p_MIDs_used))           zeros(length(udpg_MIDs_used),length(g6p_MIDs_used))   zeros(length(udpg_MIDs_used),length(fbp_MIDs_used)) zeros(length(udpg_MIDs_used),length(g1p_MIDs_used)) zeros(length(udpg_MIDs_used),length(adpg_MIDs_used)) eye(length(udpg_MIDs_used))*AbsoluteLevel_std.UDPG(4)             ];


Data_std_Px=Data_std_norm;
for i=1:5
Data_std_Px{i,MIDs_used_all} = (PVar_a1*Data_std_norm{i,MIDs_used_all}')' + ...
    (PVar_a1*(Data_mean_norm{i,MIDs_used_all})')' + ...
    Data_std_norm{i,MIDs_used_all}*(PMEAN_a1);...% a1 tpi
end
for i=6:10
Data_std_Px{i,MIDs_used_all} = (PVar_a2*Data_std_norm{i,MIDs_used_all}')' + ...
    (PVar_a2*(Data_mean_norm{i,MIDs_used_all})')' + ...
    Data_std_norm{i,MIDs_used_all}*(PMEAN_a2);...% a1 tpi
end
for i=11:15
Data_std_Px{i,MIDs_used_all} = (PVar_a3*Data_std_norm{i,MIDs_used_all}')' + ...
    (PVar_a3*(Data_mean_norm{i,MIDs_used_all})')' + ...
    Data_std_norm{i,MIDs_used_all}*(PMEAN_a3);...% a1 tpi
end
for i=16:20
Data_std_Px{i,MIDs_used_all} = (PVar_a4*Data_std_norm{i,MIDs_used_all}')' + ...
    (PVar_a4*(Data_mean_norm{i,MIDs_used_all})')' + ...
    Data_std_norm{i,MIDs_used_all}*(PMEAN_a4);...% a1 tpi
end

for sample_comp=1  % parameterization of c fixed to estimates from INCA
    disp(sample_comp)
%% create matrix that combines the MID data
%chlamy
[S_a1,comp_a1(:,sample_comp)] = build_differential_matrix_comp(Data_mean_norm(algea_c,:),model,[0.38 0.36 0.18 0 0.001]);
%sorokiniana
[S_a2,comp_a2(:,sample_comp)] = build_differential_matrix_comp(Data_mean_norm(algea_s,:),model,[0 0.003 0 0 0]);
%ohadii
[S_a3,comp_a3(:,sample_comp)] = build_differential_matrix_comp(Data_mean_norm(algea_o,:),model,[2.4 0.65 9999 0.002 1.68]);
%ohadii EIL
[S_a4,comp_a4(:,sample_comp)] = build_differential_matrix_comp(Data_mean_norm(algea_o_eil,:),model);

%% run regression
% CHLAMY
[v_a1(:,sample_comp), residuals_a1(:,sample_comp), Chi_a1(sample_comp), df_a1, EXITFLAG_a1, H_a1, P_a1(:,sample_comp)] = QP_regression_var_P(S_a1, Data_mean_norm(algea_c,:),...
    Data_std_Px(algea_c,:), rubp_MIDs_used, f6p_MIDs_used, g6p_MIDs_used, fbp_MIDs_used, g1p_MIDs_used, adpg_MIDs_used, udpg_MIDs_used, model,...
    sta(1),d_sta(1),suc(1),d_suc(1),...
    AbsoluteLevel.RuBP(1), AbsoluteLevel.F6P(1), AbsoluteLevel.G6P(1), AbsoluteLevel.FBP(1), AbsoluteLevel.G1P(1), AbsoluteLevel.ADPG_Levels_nmol_gDW__(1), AbsoluteLevel.UDPG(1),...
    AbsoluteLevel_std.RuBP(1), AbsoluteLevel_std.F6P(1), AbsoluteLevel_std.G6P(1), AbsoluteLevel_std.FBP(1), AbsoluteLevel_std.G1P(1), AbsoluteLevel_std.ADPG_Levels_nmol_gDW__(1), AbsoluteLevel_std.UDPG(1));

% SOROKINIANA
[v_a2(:,sample_comp), residuals_a2(:,sample_comp), Chi_a2(sample_comp), df_a2, EXITFLAG_a2, H_a2, P_a2(:,sample_comp)] = QP_regression_var_P(S_a2, Data_mean_norm(algea_s,:), ...
    Data_std_Px(algea_s,:), rubp_MIDs_used, f6p_MIDs_used, g6p_MIDs_used, fbp_MIDs_used, g1p_MIDs_used, adpg_MIDs_used, udpg_MIDs_used, model,...
    sta(2),d_sta(2),suc(2),d_suc(2),...
    AbsoluteLevel.RuBP(2), AbsoluteLevel.F6P(2), AbsoluteLevel.G6P(2), AbsoluteLevel.FBP(2), AbsoluteLevel.G1P(2), AbsoluteLevel.ADPG_Levels_nmol_gDW__(2), AbsoluteLevel.UDPG(2),...
    AbsoluteLevel_std.RuBP(2), AbsoluteLevel_std.F6P(2), AbsoluteLevel_std.G6P(2), AbsoluteLevel_std.FBP(2), AbsoluteLevel_std.G1P(2), AbsoluteLevel_std.ADPG_Levels_nmol_gDW__(2), AbsoluteLevel_std.UDPG(2));

% OHADII LL
[v_a3(:,sample_comp), residuals_a3(:,sample_comp), Chi_a3(sample_comp), df_a3, EXITFLAG_a3, H_a3, P_a3(:,sample_comp)] = QP_regression_var_P(S_a2, Data_mean_norm(algea_o,:), ...
    Data_std_Px(algea_o,:), rubp_MIDs_used, f6p_MIDs_used, g6p_MIDs_used, fbp_MIDs_used, g1p_MIDs_used, adpg_MIDs_used, udpg_MIDs_used, model,...
    sta(3),d_sta(3),suc(3),d_suc(3),...
    AbsoluteLevel.RuBP(3), AbsoluteLevel.F6P(3), AbsoluteLevel.G6P(3), AbsoluteLevel.FBP(3), AbsoluteLevel.G1P(3), AbsoluteLevel.ADPG_Levels_nmol_gDW__(3), AbsoluteLevel.UDPG(3),...
    AbsoluteLevel_std.RuBP(3), AbsoluteLevel_std.F6P(3), AbsoluteLevel_std.G6P(3), AbsoluteLevel_std.FBP(3), AbsoluteLevel_std.G1P(3), AbsoluteLevel_std.ADPG_Levels_nmol_gDW__(3), AbsoluteLevel_std.UDPG(3),length(algea_o));

% OHADII EIL
[v_a4(:,sample_comp), residuals_a4(:,sample_comp), Chi_a4(sample_comp), df_a4, EXITFLAG_a4, H_a4, P_a4(:,sample_comp)] = QP_regression_var_P(S_a2, Data_mean_norm(algea_o_eil,:), ...
    Data_std_Px(algea_o_eil,:), rubp_MIDs_used, f6p_MIDs_used, g6p_MIDs_used, fbp_MIDs_used, g1p_MIDs_used, adpg_MIDs_used, udpg_MIDs_used, model,...
    sta(4),d_sta(4),suc(4),d_suc(4),...
    AbsoluteLevel.RuBP(4), AbsoluteLevel.F6P(4), AbsoluteLevel.G6P(4), AbsoluteLevel.FBP(4), AbsoluteLevel.G1P(4), AbsoluteLevel.ADPG_Levels_nmol_gDW__(4), AbsoluteLevel.UDPG(4),...
    AbsoluteLevel_std.RuBP(4), AbsoluteLevel_std.F6P(4), AbsoluteLevel_std.G6P(4), AbsoluteLevel_std.FBP(4), AbsoluteLevel_std.G1P(4), AbsoluteLevel_std.ADPG_Levels_nmol_gDW__(4), AbsoluteLevel_std.UDPG(4));
end

%% find best c parameterization and do bootstrap
[~,order_chi_a1] = sort(Chi_a1);
[~,order_chi_a2] = sort(Chi_a2);
[~,order_chi_a3] = sort(Chi_a3);
[~,order_chi_a4] = sort(Chi_a4);

v_a1=v_a1(:,order_chi_a1(1)); residuals_a1=residuals_a1(:,order_chi_a1(1)); Chi_a1=Chi_a1(order_chi_a1(1)); comp_a1=comp_a1(:,order_chi_a1(1)); P_a1=P_a1(:,order_chi_a1(1));
v_a2=v_a2(:,order_chi_a2(1)); residuals_a2=residuals_a2(:,order_chi_a2(1)); Chi_a2=Chi_a2(order_chi_a2(1)); comp_a2=comp_a2(:,order_chi_a2(1)); P_a2=P_a2(:,order_chi_a2(1));
v_a3=v_a3(:,order_chi_a3(1)); residuals_a3=residuals_a3(:,order_chi_a3(1)); Chi_a3=Chi_a3(order_chi_a3(1)); comp_a3=comp_a3(:,order_chi_a3(1)); P_a3=P_a3(:,order_chi_a3(1));
v_a4=v_a4(:,order_chi_a4(1)); residuals_a4=residuals_a4(:,order_chi_a4(1)); Chi_a4=Chi_a4(order_chi_a4(1)); comp_a4=comp_a4(:,order_chi_a4(1)); P_a4=P_a4(:,order_chi_a4(1));

disp('fit chlamy:');...
disp(Chi_a1<chi2inv(0.975,df_a1))
disp('fit soro:');...
disp(Chi_a2<chi2inv(0.975,df_a2))
disp('fit ohadii LL 0-20sec:');...
disp(Chi_a3<chi2inv(0.975,df_a3))
disp('fit ohadii EIL:');...
disp(Chi_a4<chi2inv(0.975,df_a4))

AbsoluteLevel.RuBP./[v_a1(2);v_a2(2);v_a3(2);v_a4(2)]
table(model.rxns, v_a1, v_a2, v_a3, v_a4);

v_a1_net = get_net_flux(v_a1,model,model_rev);
v_a2_net = get_net_flux(v_a2,model,model_rev);
v_a3_net = get_net_flux(v_a3,model,model_rev);

% bootstrap_non_parametric_var_P


