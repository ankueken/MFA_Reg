%% Bootstrap results to get CI

nboot=1000;

chiBoot_a1 = Inf(1, nboot);
chiBoot_a2 = Inf(1, nboot);
v_boot_a1  = nan(length(v_a1), nboot);
v_boot_a2  = nan(length(v_a2), nboot);

MIDs_used_all=[rubp_MIDs_used, f6p_MIDs_used, g6p_MIDs_used, fbp_MIDs_used, g1p_MIDs_used, adpg_MIDs_used];

A1=repmat([repmat(AbsoluteLevel.RuBP(1),length(rubp_MIDs_used),1);...
    repmat(AbsoluteLevel.F6P(1),length(f6p_MIDs_used),1);...
    repmat(AbsoluteLevel.G6P(1),length(g6p_MIDs_used),1);...
    repmat(AbsoluteLevel.FBP(1),length(fbp_MIDs_used),1);...
    repmat(AbsoluteLevel.G1P(1),length(g1p_MIDs_used),1);...
    repmat(AbsoluteLevel.ADPG_Levels_nmol_gDW__(1),length(adpg_MIDs_used),1)],1,4);

A2=repmat([repmat(AbsoluteLevel.RuBP(2),length(rubp_MIDs_used),1);...
    repmat(AbsoluteLevel.F6P(2),length(f6p_MIDs_used),1);...
    repmat(AbsoluteLevel.G6P(2),length(g6p_MIDs_used),1);...
    repmat(AbsoluteLevel.FBP(2),length(fbp_MIDs_used),1);...
    repmat(AbsoluteLevel.G1P(2),length(g1p_MIDs_used),1);...
    repmat(AbsoluteLevel.ADPG_Levels_nmol_gDW__(2),length(adpg_MIDs_used),1)],1,4);

for i=1:nboot
    %% Resample data by bootstrap
    % Resample the residuals with replacement to generate a new set of residuals.
    
    disp('Bootstrap:'); ...
        disp(i)
    residuals_a1_norm = reshape(residuals_a1,length(MIDs_used_all),4)./A1;
    boot_residuals_a1 = residuals_a1_norm(randi(numel(residuals_a1_norm),size(residuals_a1_norm)));
    
    residuals_a2_norm = reshape(residuals_a2,length(MIDs_used_all),4)./A2;
    boot_residuals_a2 = residuals_a2_norm(randi(numel(residuals_a2_norm),size(residuals_a2_norm)));
    
    %% Add the resampled residuals to the predicted response to generate a bootstrap data set
    Data_mean_boot = Data_mean_abs;
    Data_mean_boot{2:5,MIDs_used_all} = Data_mean_abs{2:5,MIDs_used_all}+Data_mean_abs{2:5,MIDs_used_all}.*boot_residuals_a1';
    Data_mean_boot{7:10,MIDs_used_all} = Data_mean_abs{7:10,MIDs_used_all}+Data_mean_abs{7:10,MIDs_used_all}.*boot_residuals_a2';
    
    Names=unique(cellfun(@(x) x(1:end-2),Data_mean_boot.Properties.VariableNames(4:end),'UniformOutput',false));
    Data_norm_boot=Data_mean_boot;
    for m=1:length(Names)
        MIDs = find(contains(Data_mean_boot.Properties.VariableNames,Names{m}));
        Data_norm_boot{:,MIDs}=Data_mean_boot{:,MIDs}./sum(Data_mean_boot{:,MIDs},2);
    end
    
    %% Treat the bootstrap dataset as an independent replicate experiment,
    %  and fit it to the model to calculate new estimates of model parameters.
    
    %chlamy
    [S_t0_a1,S_t5_a1,S_t10_a1,S_t20_a1] = build_differential_matrix_sCBC(Data_norm_boot(algea_c,:),model);
    
    [v_a1_b, ~, Chi_a1_b, ~, EXITFLAG_a1] = regression_sCBC(S_t0_a1,S_t5_a1,S_t10_a1,S_t20_a1, Data_mean_boot(1:5,:),...
        Data_std_abs(1,:), Data_std_abs(2,:), Data_std_abs(3,:), Data_std_abs(4,:),Data_std_abs(5,:),...
        rubp_MIDs_used, f6p_MIDs_used, g6p_MIDs_used, fbp_MIDs_used, g1p_MIDs_used, adpg_MIDs_used, model,...
        sta(1),d_sta(1));
    
    if EXITFLAG_a1==1
        v_boot_a1(:,i) = v_a1_b;
        chiBoot_a1(:,i) = Chi_a1_b;
        clear v_a1_b Chi_a1_b
    end
    
    %sorokiniana
    [S_t0_a2,S_t5_a2,S_t10_a2,S_t20_a2] = build_differential_matrix_sCBC(Data_norm_boot(algea_s,:),model);
    
    [v_a2_b, ~, Chi_a2_b, df_a2, EXITFLAG_a2] = regression_sCBC(S_t0_a2,S_t5_a2,S_t10_a2,S_t20_a2, Data_mean_boot(6:10,:), ...
        Data_std_abs(6,:),Data_std_abs(7,:), Data_std_abs(8,:), Data_std_abs(9,:),Data_std_abs(10,:),...
        rubp_MIDs_used, f6p_MIDs_used, g6p_MIDs_used, fbp_MIDs_used, g1p_MIDs_used, adpg_MIDs_used, model,...
        sta(2),d_sta(2));
    
    if EXITFLAG_a2==1
        v_boot_a2(:,i) = v_a2_b;
        chiBoot_a2(:,i) = Chi_a2_b;
        clear v_a2_b Chi_a2_b
    end
end

%% To calculate the 95% bootstrap CIs, compute the 97.5th and the 2.5th
% % percentile values of each parameter from the bootstrap distributions.

v_boot_a1=get_net_flux_v_boot(v_boot_a1,17,model,model_rev);
v_boot_a2=get_net_flux_v_boot(v_boot_a2,17,model,model_rev);

vbootCI_a1 = prctile(v_boot_a1',[2.5 97.5])';
vbootCI_a2 = prctile(v_boot_a2',[2.5 97.5])';

chibootCI_a1 = prctile(chiBoot_a1',[2.5 97.5]);
chibootCI_a2 = prctile(chiBoot_a2',[2.5 97.5]);

chi2cdf(chibootCI_a1,df_a1)
chi2cdf(chibootCI_a2,df_a2)
