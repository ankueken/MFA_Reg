%% build differential matrix S combining MID data
function [S_t0,S_t5,S_t10,S_t20] = build_differential_matrix_sCBC(Data_mean,model)

RuBP_MID_norm = Data_mean{:,contains(Data_mean.Properties.VariableNames,'RuBP')};
PP_MID_norm = Data_mean{:,contains(Data_mean.Properties.VariableNames,'PP')};
FBP_MID_norm = Data_mean{:,contains(Data_mean.Properties.VariableNames,'FBP')};
F6P_MID_norm = Data_mean{:,contains(Data_mean.Properties.VariableNames,'f6p')};
G6P_MID_norm = Data_mean{:,contains(Data_mean.Properties.VariableNames,'G6P')};
G1P_MID_norm = Data_mean{:,contains(Data_mean.Properties.VariableNames,'g1p')};
T3P_MID_norm = Data_mean{:,contains(Data_mean.Properties.VariableNames,'DHAP')};
ADPG_MID_norm = Data_mean{:,contains(Data_mean.Properties.VariableNames,'ADPG')};

for tp=1:4 % for each time point
    
    S_RuBP=zeros(6,size(model.S,2)); % (#MIDs x #fluxes)
    S_F6P=zeros(7,size(model.S,2));
    S_G6P=zeros(7,size(model.S,2));
    S_G1P=zeros(7,size(model.S,2));
    S_FBP=zeros(7,size(model.S,2));
    S_ADPG=zeros(7,size(model.S,2));
    %%
    % *RuBP*
    %
    % $$\frac{dRuBP}{dt}=v_{10} -v_1$$
    %
    % $$\frac{{dRuBP}_{m+x} }{dt}=\frac{1}{P_{RuBP} }{\left(v_{10} ?{PP}_{m+x} -v_1
    % ?{RuBP}_{m+x} \right)}$$
    
    S_RuBP(:,find(strcmp(model.rxns,'RUBISCO_CO2')))=-RuBP_MID_norm(tp,:); % RuBP + CO2 -> 2* 3PGA
    S_RuBP(:,find(strcmp(model.rxns,'RUBISCO_O2')))=-RuBP_MID_norm(tp,:); % RuBP -> 3PGA + 2PG
    S_RuBP(:,find(strcmp(model.rxns,'PRK')))=PP_MID_norm(tp,:); % Ru5P -> RuBP
    %%
    % *F6P*
    %
    % $$\frac{dF6P}{dt}=v_4 +\overleftarrow{{\mathit{\mathbf{v}}}_5 } +\overleftarrow{v_{11}
    % } -\overrightarrow{v_5 } -\overrightarrow{v_{11} }$$
    
    S_F6P(:,find(strcmp(model.rxns,'TK3')))=-F6P_MID_norm(tp,:); %v5
    S_F6P(:,find(strcmp(model.rxns,'PGIp')))=-F6P_MID_norm(tp,:); %v11f
    %%
    % $$\frac{{dF6P}_{m+0} }{dt}=\frac{1}{P_{F6P} }{\left(v_4 ?{FBP}_{m+0} +\overleftarrow{v_{11}
    % } ?{G6P}_{m+0} -\left(v_5 +\overrightarrow{v_{11} } \right)?{F6P}_{m+0} \right)}$$
    
    S_F6P(:,find(strcmp(model.rxns,'PFPp')))=FBP_MID_norm(tp,:); % FBP -> F6P
    S_F6P(:,find(strcmp(model.rxns,'PGIp_rev')))=G6P_MID_norm(tp,:); %v11b
    %%
    % *G6P*
    %
    % $$\frac{dF6P}{dt}=v_4 +\overleftarrow{{\mathit{\mathbf{v}}}_5 } +\overleftarrow{v_{11}
    % } -\overrightarrow{v_5 } -\overrightarrow{v_{11} }$$
    
    S_G6P(:,find(strcmp(model.rxns,'PGIp')))=F6P_MID_norm(tp,:); %plastid
    S_G6P(:,find(strcmp(model.rxns,'PGMp_rev')))=G1P_MID_norm(tp,:);
    %%
    % $$\frac{{dF6P}_{m+0} }{dt}=\frac{1}{P_{F6P} }{\left(v_4 ?{FBP}_{m+0} +\overleftarrow{v_{11}
    % } ?{G6P}_{m+0} -\left(v_5 +\overrightarrow{v_{11} } \right)?{F6P}_{m+0} \right)}$$
    
    S_G6P(:,find(strcmp(model.rxns,'PGMp')))=-G6P_MID_norm(tp,:); %plastid
    S_G6P(:,find(strcmp(model.rxns,'PGIp_rev')))=-G6P_MID_norm(tp,:); %
    %%
    % *FBP*
    %
    % $$\frac{dFBP}{dt}=\overrightarrow{v_3 } -v_4 -\overleftarrow{v_3 }$$
    
    S_FBP(:,find(strcmp(model.rxns,'PFPp')))=-FBP_MID_norm(tp,:);
    S_FBP(:,find(strcmp(model.rxns,'FBAp_rev')))=-FBP_MID_norm(tp,:);
    %%
    % $$\frac{{dFBP}_{m+0} }{dt}=\frac{1}{P_{FBP} }{\left(\overrightarrow{v_3 }
    % ?{T3P}_{m+0} {T3P}_{m+0} -\left(v_4 +-\overleftarrow{v_3 } \right)?{FBP}_{m+0}
    % \right)}$$
    
    S_FBP(1,find(strcmp(model.rxns,'FBAp')))=(T3P_MID_norm(tp,1)*T3P_MID_norm(tp,1));
    
    %%
    % $$\frac{{dFBP}_{m+1} }{dt}=\frac{1}{P_{FBP} }{\left(\overrightarrow{v_3 }
    % ?\left({T3P}_{m+1} {T3P}_{m+0} \right)-\left(v_4 +-\overleftarrow{v_3 } \right)?{FBP}_{m+1}
    % \right)}$$
    
    S_FBP(2,find(strcmp(model.rxns,'FBAp')))=(T3P_MID_norm(tp,2))*(T3P_MID_norm(tp,1));
    %%
    % $$\frac{{dFBP}_{m+2} }{dt}=\frac{1}{P_{FBP} }{\left(\overrightarrow{v_3 }
    % ?\left({T3P}_{m+2} {T3P}_{m+0} +{T3P}_{m+1} {T3P}_{m+1} \right)-\left(v_4 +-\overleftarrow{v_3
    % } \right)?{FBP}_{m+2} \right)}$$
    
    S_FBP(3,find(strcmp(model.rxns,'FBAp')))=((T3P_MID_norm(tp,3))*(T3P_MID_norm(tp,1)))+((T3P_MID_norm(tp,2))*(T3P_MID_norm(tp,2)));
    %%
    % $$\frac{{dFBP}_{m+3} }{dt}=\frac{1}{P_{FBP} }{\left(\overrightarrow{v_3 }
    % ?\left({T3P}_{m+3} {T3P}_{m+0} +{T3P}_{m+2} {T3P}_{m+1} \right)-\left(v_4 +-\overleftarrow{v_3
    % } \right)?{FBP}_{m+3} \right)}$$
    
    S_FBP(4,find(strcmp(model.rxns,'FBAp')))=((T3P_MID_norm(tp,4))*(T3P_MID_norm(tp,1)))+((T3P_MID_norm(tp,3))*(T3P_MID_norm(tp,2)));
    %%
    % $$\frac{{dFBP}_{m+4} }{dt}=\frac{1}{P_{FBP} }{\left(\overrightarrow{v_3 }
    % ?\left({T3P}_{m+3} {T3P}_{m+1} +{T3P}_{m+2} {T3P}_{m+2} \right)-\left(v_4 +-\overleftarrow{v_3
    % } \right)?{FBP}_{m+4} \right)}$$
    
    S_FBP(5,find(strcmp(model.rxns,'FBAp')))=((T3P_MID_norm(tp,4))*(T3P_MID_norm(tp,2)))+((T3P_MID_norm(tp,3))*(T3P_MID_norm(tp,3)));
    %%
    % $$\frac{{dFBP}_{m+5} }{dt}=\frac{1}{P_{FBP} }{\left(\overrightarrow{v_3 }
    % ?\left({T3P}_{m+3} {T3P}_{m+2} \right)-\left(v_4 +-\overleftarrow{v_3 } \right)?{FBP}_{m+5}
    % \right)}$$
    
    S_FBP(6,find(strcmp(model.rxns,'FBAp')))=(T3P_MID_norm(tp,4))*(T3P_MID_norm(tp,3));
    %%
    % $$\frac{{dFBP}_{m+6} }{dt}=\frac{1}{P_{FBP} }{\left(\overrightarrow{v_3 }
    % ?\left({T3P}_{m+3} T3P_{m+3} \right)-\left(v_4 +-\overleftarrow{v_3 } \right)?{FBP}_{m+6}
    % \right)}$$
    
    S_FBP(7,find(strcmp(model.rxns,'FBAp')))=(T3P_MID_norm(tp,4))*(T3P_MID_norm(tp,4));
    
    %%
    % *ADPG*
    %
    S_ADPG(:,find(strcmp(model.rxns,'SS')))=-ADPG_MID_norm(tp,:); % ADPG ->
    S_ADPG(:,find(strcmp(model.rxns,'AGP')))=G1P_MID_norm(tp,:); % G1P -> ADPG
    %
    %%
    % *G1P*
    S_G1P(:,find(strcmp(model.rxns,'PGMp')))=G6P_MID_norm(tp,:); 

    S_G1P(:,find(strcmp(model.rxns,'AGP')))=-G1P_MID_norm(tp,:); 
    S_G1P(:,find(strcmp(model.rxns,'PGMp_rev')))=-G1P_MID_norm(tp,:); 
    
    % combine matrices of metabolites
    if tp==1
        S_t0 = [S_RuBP; S_F6P; S_G6P; S_FBP; S_G1P; S_ADPG];
    elseif tp==2
        S_t5 = [S_RuBP; S_F6P; S_G6P; S_FBP; S_G1P; S_ADPG];
    elseif tp==3
        S_t10 = [S_RuBP; S_F6P; S_G6P; S_FBP; S_G1P; S_ADPG];
    elseif tp==4
        S_t20 = [S_RuBP; S_F6P; S_G6P; S_FBP; S_G1P; S_ADPG];
    end
end
end