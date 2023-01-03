%% Bootstrap results to get CI

nboot=1000; % number of boot samples

chiBoot_a1 = Inf(1, nboot);
chiBoot_a2 = Inf(1, nboot);
chiBoot_a3 = Inf(1, nboot);
% chiBoot_a4 = Inf(1, nboot);

v_boot_a1  = nan(length(v_a1), nboot);
v_boot_a2  = nan(length(v_a2), nboot);
v_boot_a3  = nan(length(v_a3), nboot);
% v_boot_a4  = nan(length(v_a4), nboot);

MIDs_used_all=[rubp_MIDs_used, f6p_MIDs_used, g6p_MIDs_used, fbp_MIDs_used, g1p_MIDs_used, adpg_MIDs_used, udpg_MIDs_used];

A1=repmat([repmat(AbsoluteLevel.RuBP(1),length(rubp_MIDs_used),1);...
    repmat(AbsoluteLevel.F6P(1),length(f6p_MIDs_used),1);...
    repmat(AbsoluteLevel.G6P(1),length(g6p_MIDs_used),1);...
    repmat(AbsoluteLevel.FBP(1),length(fbp_MIDs_used),1);...
    repmat(AbsoluteLevel.G1P(1),length(g1p_MIDs_used),1);...
    repmat(AbsoluteLevel.ADPG_Levels_nmol_gDW__(1),length(adpg_MIDs_used),1);...
    repmat(AbsoluteLevel.UDPG(1),length(udpg_MIDs_used),1)],1,4);

A2=repmat([repmat(AbsoluteLevel.RuBP(2),length(rubp_MIDs_used),1);...
    repmat(AbsoluteLevel.F6P(2),length(f6p_MIDs_used),1);...
    repmat(AbsoluteLevel.G6P(2),length(g6p_MIDs_used),1);...
    repmat(AbsoluteLevel.FBP(2),length(fbp_MIDs_used),1);...
    repmat(AbsoluteLevel.G1P(2),length(g1p_MIDs_used),1);...
    repmat(AbsoluteLevel.ADPG_Levels_nmol_gDW__(2),length(adpg_MIDs_used),1);...
    repmat(AbsoluteLevel.UDPG(2),length(udpg_MIDs_used),1)],1,4);

A3=repmat([repmat(AbsoluteLevel.RuBP(3),length(rubp_MIDs_used),1);...
    repmat(AbsoluteLevel.F6P(3),length(f6p_MIDs_used),1);...
    repmat(AbsoluteLevel.G6P(3),length(g6p_MIDs_used),1);...
    repmat(AbsoluteLevel.FBP(3),length(fbp_MIDs_used),1);...
    repmat(AbsoluteLevel.G1P(3),length(g1p_MIDs_used),1);...
    repmat(AbsoluteLevel.ADPG_Levels_nmol_gDW__(3),length(adpg_MIDs_used),1);...
    repmat(AbsoluteLevel.UDPG(3),length(udpg_MIDs_used),1)],1,3);

% A4=repmat([repmat(AbsoluteLevel.RuBP(4),length(rubp_MIDs_used),1);...
%     repmat(AbsoluteLevel.F6P(4),length(f6p_MIDs_used),1);...
%     repmat(AbsoluteLevel.G6P(4),length(g6p_MIDs_used),1);...
%     repmat(AbsoluteLevel.FBP(4),length(fbp_MIDs_used),1);...
%     repmat(AbsoluteLevel.G1P(4),length(g1p_MIDs_used),1);...
%     repmat(AbsoluteLevel.ADPG_Levels_nmol_gDW__(4),length(adpg_MIDs_used),1);...
%     repmat(AbsoluteLevel.UDPG(4),length(udpg_MIDs_used),1)],1,3);

for i=1:nboot
%% Resample data by bootstrap
% Resample the residuals with replacement to generate a new set of residuals. 
% 
residuals_a1_norm = reshape(residuals_a1,length(MIDs_used_all),4)./A1;
boot_residuals_a1 = residuals_a1_norm(randi(numel(residuals_a1_norm),size(residuals_a1_norm)));

residuals_a2_norm = reshape(residuals_a2,length(MIDs_used_all),4)./A2;
boot_residuals_a2 = residuals_a2_norm(randi(numel(residuals_a2_norm),size(residuals_a2_norm)));

residuals_a3_norm = reshape(residuals_a3,length(MIDs_used_all),3)./A3;
boot_residuals_a3 = residuals_a3_norm(randi(numel(residuals_a3_norm),size(residuals_a3_norm)));

% residuals_a4_norm = reshape(residuals_a4,length(MIDs_used_all),3)./A4;
% boot_residuals_a4 = residuals_a4_norm(randi(numel(residuals_a4_norm),size(residuals_a4_norm)));

disp('Bootstrap:'); ...
    disp(i)
    
    %% Add the resampled residuals to the predicted response to generate a bootstrap data set    
    Data_mean_boot = Data_mean_abs;
    Data_mean_boot{2:5,MIDs_used_all} = Data_mean_abs{2:5,MIDs_used_all}+Data_mean_abs{2:5,MIDs_used_all}.*boot_residuals_a1';
    Data_mean_boot{7:10,MIDs_used_all} = Data_mean_abs{7:10,MIDs_used_all}+Data_mean_abs{7:10,MIDs_used_all}.*boot_residuals_a2';    
    Data_mean_boot{12:14,MIDs_used_all} = Data_mean_abs{12:14,MIDs_used_all}+Data_mean_abs{12:14,MIDs_used_all}.*boot_residuals_a3';    
%     Data_mean_boot{17:19,MIDs_used_all} = Data_mean_abs{17:19,MIDs_used_all}+Data_mean_abs{17:19,MIDs_used_all}.*boot_residuals_a4';    

    Names=unique(cellfun(@(x) x(1:end-2),Data_mean_boot.Properties.VariableNames(4:end),'UniformOutput',false));
    Data_norm_boot=Data_mean_boot;
    for m=1:length(Names)
        MIDs = find(contains(Data_mean_boot.Properties.VariableNames,Names{m}));
        Data_norm_boot{:,MIDs}=Data_mean_boot{:,MIDs}./sum(Data_mean_boot{:,MIDs},2);
    end

    %% Treat the bootstrap dataset as an independent replicate experiment,
    %  and fit it to the model to calculate new estimates of model parameters.
    
    %chlamy
    for j=1:10
    [S_t0_a1,S_t5_a1,S_t10_a1,S_t20_a1,comp_boot_a1(:,i)] = build_differential_matrix_comp(Data_norm_boot(algea_c,:),model,comp_a1);
    
        [v_a1_b, ~, Chi_a1_b, ~, EXITFLAG_a1] = regression_T40(S_t0_a1,S_t5_a1,S_t10_a1,S_t20_a1, Data_mean_boot(1:5,:),...
    Data_std_abs(1,:), Data_std_abs(2,:), Data_std_abs(3,:), Data_std_abs(4,:),Data_std_abs(5,:),...
    rubp_MIDs_used, f6p_MIDs_used, g6p_MIDs_used, fbp_MIDs_used, g1p_MIDs_used, adpg_MIDs_used, udpg_MIDs_used, model,...
    sta(1),d_sta(1),suc(1),d_suc(1));

        if EXITFLAG_a1==1 && Chi_a1_b<=chi2inv(0.975,df_a1)
            v_boot_a1(:,i) = v_a1_b;
            chiBoot_a1(:,i) = Chi_a1_b;
            break
        elseif EXITFLAG_a1==1 && Chi_a1_b<chiBoot_a1(:,i)
            v_boot_a1(:,i) = v_a1_b;
            chiBoot_a1(:,i) = Chi_a1_b;
        end
        clear v_a1_b Chi_a1_b
    end
    
    %sorokiniana
    for j=1:10
    [S_t0_a2,S_t5_a2,S_t10_a2,S_t20_a2,comp_boot_a2(:,i)] = build_differential_matrix_comp(Data_norm_boot(algea_s,:),model,comp_a2);
        
        [v_a2_b, ~, Chi_a2_b, df_a2, EXITFLAG_a2] = regression_T40(S_t0_a2,S_t5_a2,S_t10_a2,S_t20_a2, Data_mean_boot(6:10,:), ...
    Data_std_abs(6,:),Data_std_abs(7,:), Data_std_abs(8,:), Data_std_abs(9,:),Data_std_abs(10,:),...
    rubp_MIDs_used, f6p_MIDs_used, g6p_MIDs_used, fbp_MIDs_used, g1p_MIDs_used, adpg_MIDs_used, udpg_MIDs_used, model,...
    sta(2),d_sta(2),suc(2),d_suc(2));

        if EXITFLAG_a2==1 && Chi_a2_b<=chi2inv(0.975,df_a2)
            v_boot_a2(:,i) = v_a2_b;
            chiBoot_a2(:,i) = Chi_a2_b;
            break
        elseif EXITFLAG_a2==1 && Chi_a2_b<chiBoot_a2(:,i)
            v_boot_a2(:,i) = v_a2_b;
            chiBoot_a2(:,i) = Chi_a2_b;
        end
        clear v_a2_b Chi_a2_b
    end
    
    %ohadii LL 0-20sec
    for j=1:10
    [S_t0_a3,S_t5_a3,S_t10_a3,S_t20_a3,comp_boot_a3(:,i)] = build_differential_matrix_comp(Data_norm_boot(algea_o,:),model,comp_a3);
     
        [v_a3_b, ~, Chi_a3_b, ~, EXITFLAG_a3] = regression_T20(S_t0_a3,S_t5_a3,S_t10_a3,S_t20_a3, Data_mean_boot(11:15,:), ...
    Data_std_abs(11,:), Data_std_abs(12,:), Data_std_abs(13,:), Data_std_abs(14,:),Data_std_abs(15,:),...
    rubp_MIDs_used, f6p_MIDs_used, g6p_MIDs_used, fbp_MIDs_used, g1p_MIDs_used, adpg_MIDs_used, udpg_MIDs_used, model,...
    sta(3),d_sta(3),suc(3),d_suc(3));

        if EXITFLAG_a3==1 && Chi_a3_b<=chi2inv(0.975,df_a3)
            v_boot_a3(:,i) = v_a3_b;
            chiBoot_a3(:,i) = Chi_a3_b;
            break
        elseif EXITFLAG_a3==1 && Chi_a3_b<chiBoot_a3(:,i)
            v_boot_a3(:,i) = v_a3_b;
            chiBoot_a3(:,i) = Chi_a3_b;
        end
        clear v_a3_b Chi_a3_b
    end
    
%     %ohadii eil
%     for j=1:10   
%     [S_t0_a4,S_t5_a4,S_t10_a4,S_t20_a4,comp_boot_a4(:,i)] = build_differential_matrix_comp(Data_norm_boot(algea_o_eil,:),model,comp_a4);
%        
%         [v_a4_b, ~, Chi_a4_b, ~, EXITFLAG_a4] = regression_T40(S_t0_a4,S_t5_a4,S_t10_a4,S_t20_a4, Data_mean_boot(16:20,:), ...
%     Data_std_abs(16,:), Data_std_abs(17,:), Data_std_abs(18,:), Data_std_abs(19,:),Data_std_abs(20,:),...
%     rubp_MIDs_used, f6p_MIDs_used, g6p_MIDs_used, fbp_MIDs_used, g1p_MIDs_used, adpg_MIDs_used, udpg_MIDs_used, model,...
%     sta(4),d_sta(4),suc(4),d_suc(4));
%         if EXITFLAG_a4==1 && Chi_a4_b<=chi2inv(0.975,df_a4)
%             v_boot_a4(:,i) = v_a4_b;
%             chiBoot_a4(:,i) = Chi_a4_b;
%             break
%         elseif EXITFLAG_a4==1 && Chi_a4_b<chiBoot_a4(:,i)
%             v_boot_a4(:,i) = v_a4_b;
%             chiBoot_a4(:,i) = Chi_a4_b;
%         end
%         clear v_a4_b Chi_a4_b
%     end

end

%% To calculate the 95% bootstrap CIs, compute the 97.5th and the 2.5th
%  percentile values of each parameter from the bootstrap distributions.

chiBoot_a3(isinf(chiBoot_a3))=[];
v_boot_a3=v_boot_a3(:,~isinf(chiBoot_a3));

% we want to report net flux in the final results for easier
% comparison across algae 
v_boot_a1 = get_net_flux_v_boot(v_boot_a1,34,model,model_rev);
v_boot_a2 = get_net_flux_v_boot(v_boot_a2,34,model,model_rev);
v_boot_a3 = get_net_flux_v_boot(v_boot_a3,34,model,model_rev);

vbootCI_a1 = prctile(v_boot_a1',[2.5 97.5])';
vbootCI_a2 = prctile(v_boot_a2',[2.5 97.5])';
vbootCI_a3 = prctile(v_boot_a3',[2.5 97.5])';

chibootCI_a1 = prctile(chiBoot_a1',[2.5 97.5]);
chibootCI_a2 = prctile(chiBoot_a2',[2.5 97.5]);
chibootCI_a3 = prctile(chiBoot_a3',[2.5 97.5]);

chi2cdf(chibootCI_a1,df_a1)
sum(chiBoot_a1<=chi2inv(0.975,df_a1))/1000
chi2cdf(chibootCI_a2,df_a2)
sum(chiBoot_a2<=chi2inv(0.975,df_a2))/1000
chi2cdf(chibootCI_a3,df_a3)
sum(chiBoot_a3<=chi2inv(0.975,df_a3))/1000

v_a1 = get_net_flux(v_a1,model,model_rev);
v_a2 = get_net_flux(v_a2,model,model_rev);
v_a3 = get_net_flux(v_a3,model,model_rev);

save(strcat('Results_non-parametric_MID_rescaled_',datestr(now,'yyyymmdd_HH_MMPM'),'.mat'))