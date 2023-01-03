%% build matrix S combining MID data for compartmented CBC model
function [S_t0,S_t5,S_t10,S_t20,comp] = build_differential_matrix_comp(Data_mean,model,comp)
% INPUT:
% Data_mean: MID data
% model: network structure
% comp (optional): set of compartmentation parameter 

RuBP_MID_norm = Data_mean{:,contains(Data_mean.Properties.VariableNames,'RuBP')};
ADPG_MID_norm = Data_mean{:,contains(Data_mean.Properties.VariableNames,'ADPG')};
PP_MID_norm = Data_mean{:,contains(Data_mean.Properties.VariableNames,'PP')};
FBP_MID_norm = Data_mean{:,contains(Data_mean.Properties.VariableNames,'FBP')};
F6P_MID_norm = Data_mean{:,contains(Data_mean.Properties.VariableNames,'f6p')};
G6P_MID_norm = Data_mean{:,contains(Data_mean.Properties.VariableNames,'G6P')};
G1P_MID_norm = Data_mean{:,contains(Data_mean.Properties.VariableNames,'g1p')};
T3P_MID_norm = Data_mean{:,contains(Data_mean.Properties.VariableNames,'DHAP')};
AKG_MID_norm = Data_mean{:,contains(Data_mean.Properties.VariableNames,'x2_OG')};
UDPG_MID_norm = Data_mean{:,contains(Data_mean.Properties.VariableNames,'UDPG')};
GLUT_MID_norm = Data_mean{:,contains(Data_mean.Properties.VariableNames,'glut')};

% sample compartmentation parameter
comp_levels=table(ones(6,1),'VariableNames',{'cytosol'});
comp_levels.cytosol(2)=0.1 + (0.9-0.1).*rand(1,1);
comp_levels.cytosol(3)=0.1 + (0.9-0.1).*rand(1,1);
comp_levels.cytosol(4)=0.1 + (0.9-0.1).*rand(1,1);
comp_levels.cytosol(5)=0.1 + (0.9-0.1).*rand(1,1);
comp_levels.cytosol(6)=0.1 + (0.9-0.1).*rand(1,1);

if ~exist('comp')
    comp=comp_levels.cytosol(2:6);
else % if a set of compartmentation parameter is given only allow small deviation to find feasible solution
    comp_levels.cytosol(2)=(comp(1)*0.975) + ((comp(1)*1.025)-(comp(1)*0.975)).*rand(1,1);
    comp_levels.cytosol(3)=(comp(2)*0.975) + ((comp(2)*1.025)-(comp(2)*0.975)).*rand(1,1);
    comp_levels.cytosol(4)=(comp(3)*0.975) + ((comp(3)*1.025)-(comp(3)*0.975)).*rand(1,1);
    comp_levels.cytosol(5)=(comp(4)*0.975) + ((comp(4)*1.025)-(comp(4)*0.975)).*rand(1,1);
    comp_levels.cytosol(6)=(comp(5)*0.975) + ((comp(5)*1.025)-(comp(5)*0.975)).*rand(1,1);

end

N=model.S;

for tp=1:4 % for each time point build matrix S
    
    S_RuBP=zeros(6,size(N,2)); % (#MIDs x #fluxes)
    S_ADPG=zeros(7,size(N,2));
    S_G1P=zeros(7,size(N,2));
    S_G6P=zeros(7,size(N,2));
    S_FBP=zeros(7,size(N,2));
    S_F6P=zeros(7,size(N,2));
    S_UDPG=zeros(7,size(N,2));
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
    %
    %%
    % *ADPG*
    %
    S_ADPG(:,find(strcmp(model.rxns,'SS')))=-ADPG_MID_norm(tp,:); % ADPG ->
    S_ADPG(:,find(strcmp(model.rxns,'AGP')))=G1P_MID_norm(tp,:); % G1P -> ADPG
    %
    %%
    % *G1P*
    %
    S_G1P(:,find(strcmp(model.rxns,'PGMp')))=G6P_MID_norm(tp,:); %plastid
    S_G1P(:,find(strcmp(model.rxns,'PGMc')))=comp_levels.cytosol(6)*G6P_MID_norm(tp,:); %cytosol
%     S_G1P(:,find(strcmp(model.rxns,'GPU_rev')))=UDPG_MID_norm(tp,:); %
    
    S_G1P(:,find(strcmp(model.rxns,'AGP')))=-G1P_MID_norm(tp,:); %plastid
    S_G1P(:,find(strcmp(model.rxns,'PGMp_rev')))=-G1P_MID_norm(tp,:); %
    S_G1P(:,find(strcmp(model.rxns,'GPU')))=-comp_levels.cytosol(5)*G1P_MID_norm(tp,:); %cytosol
    S_G1P(:,find(strcmp(model.rxns,'PGMc_rev')))=-comp_levels.cytosol(5)*G1P_MID_norm(tp,:); %
    
    %%
    % *G6P*
    %
    S_G6P(:,find(strcmp(model.rxns,'PGIp')))=F6P_MID_norm(tp,:); %plastid
    S_G6P(:,find(strcmp(model.rxns,'PGMp_rev')))=G1P_MID_norm(tp,:);
    S_G6P(:,find(strcmp(model.rxns,'PGIc')))=comp_levels.cytosol(3)*F6P_MID_norm(tp,:);  %cytosol
    S_G6P(:,find(strcmp(model.rxns,'PGMc_rev')))=comp_levels.cytosol(5)*G1P_MID_norm(tp,:);
    
    S_G6P(:,find(strcmp(model.rxns,'PGMp')))=-G6P_MID_norm(tp,:); %plastid
    S_G6P(:,find(strcmp(model.rxns,'PGIp_rev')))=-G6P_MID_norm(tp,:); %
    S_G6P(:,find(strcmp(model.rxns,'PGMc')))=-comp_levels.cytosol(6)*G6P_MID_norm(tp,:); %cytosol
    S_G6P(:,find(strcmp(model.rxns,'PGIc_rev')))=-comp_levels.cytosol(6)*G6P_MID_norm(tp,:); %
    
    %%
    % *UDPG*
    %
%     S_UDPG(:,find(strcmp(model.rxns,'GPU_rev')))=-UDPG_MID_norm(tp,:); %
    S_UDPG(:,find(strcmp(model.rxns,'SPS')))=-UDPG_MID_norm(tp,:); %
    S_UDPG(:,find(strcmp(model.rxns,'GPU')))=comp_levels.cytosol(5)*G1P_MID_norm(tp,:); % G1P -> UDPG
    
    %%
    % *FBP*
    %
    % $$\frac{dFBP}{dt}=\overrightarrow{v_3 } -v_4 -\overleftarrow{v_3 }$$
    
    S_FBP(:,find(strcmp(model.rxns,'PFPp')))=-FBP_MID_norm(tp,:);
    S_FBP(:,find(strcmp(model.rxns,'FBAp_rev')))=-FBP_MID_norm(tp,:);
    S_FBP(:,find(strcmp(model.rxns,'PFPc')))=-comp_levels.cytosol(4)*FBP_MID_norm(tp,:);
    S_FBP(:,find(strcmp(model.rxns,'FBAc_rev')))=-comp_levels.cytosol(4)*FBP_MID_norm(tp,:);
    %%
    % $$\frac{{dFBP}_{m+0} }{dt}=\frac{1}{P_{FBP} }{\left(\overrightarrow{v_3 }
    % ?{T3P}_{m+0} {T3P}_{m+0} -\left(v_4 +-\overleftarrow{v_3 } \right)?{FBP}_{m+0}
    % \right)}$$
    
    S_FBP(1,find(strcmp(model.rxns,'FBAp')))=(T3P_MID_norm(tp,1)*T3P_MID_norm(tp,1));
    S_FBP(1,find(strcmp(model.rxns,'FBAc')))=(comp_levels.cytosol(2)*T3P_MID_norm(tp,1))*(comp_levels.cytosol(2)*T3P_MID_norm(tp,1));
    %%
    % $$\frac{{dFBP}_{m+1} }{dt}=\frac{1}{P_{FBP} }{\left(\overrightarrow{v_3 }
    % ?\left({T3P}_{m+1} {T3P}_{m+0} \right)-\left(v_4 +-\overleftarrow{v_3 } \right)?{FBP}_{m+1}
    % \right)}$$
    
    S_FBP(2,find(strcmp(model.rxns,'FBAp')))=(T3P_MID_norm(tp,2))*(T3P_MID_norm(tp,1));
    S_FBP(2,find(strcmp(model.rxns,'FBAc')))=(comp_levels.cytosol(2)*T3P_MID_norm(tp,2))*(comp_levels.cytosol(2)*T3P_MID_norm(tp,1));
    %%
    % $$\frac{{dFBP}_{m+2} }{dt}=\frac{1}{P_{FBP} }{\left(\overrightarrow{v_3 }
    % ?\left({T3P}_{m+2} {T3P}_{m+0} +{T3P}_{m+1} {T3P}_{m+1} \right)-\left(v_4 +-\overleftarrow{v_3
    % } \right)?{FBP}_{m+2} \right)}$$
    
    S_FBP(3,find(strcmp(model.rxns,'FBAp')))=((T3P_MID_norm(tp,3))*(T3P_MID_norm(tp,1)))+((T3P_MID_norm(tp,2))*(T3P_MID_norm(tp,2)));
    S_FBP(3,find(strcmp(model.rxns,'FBAc')))=((comp_levels.cytosol(2)*T3P_MID_norm(tp,3))*(comp_levels.cytosol(2)*T3P_MID_norm(tp,1)))+((comp_levels.cytosol(2)*T3P_MID_norm(tp,2))*(comp_levels.cytosol(2)*T3P_MID_norm(tp,2)));
    %%
    % $$\frac{{dFBP}_{m+3} }{dt}=\frac{1}{P_{FBP} }{\left(\overrightarrow{v_3 }
    % ?\left({T3P}_{m+3} {T3P}_{m+0} +{T3P}_{m+2} {T3P}_{m+1} \right)-\left(v_4 +-\overleftarrow{v_3
    % } \right)?{FBP}_{m+3} \right)}$$
    
    S_FBP(4,find(strcmp(model.rxns,'FBAp')))=((T3P_MID_norm(tp,4))*(T3P_MID_norm(tp,1)))+((T3P_MID_norm(tp,3))*(T3P_MID_norm(tp,2)));
    S_FBP(4,find(strcmp(model.rxns,'FBAc')))=((comp_levels.cytosol(2)*T3P_MID_norm(tp,4))*(comp_levels.cytosol(2)*T3P_MID_norm(tp,1)))+((comp_levels.cytosol(2)*T3P_MID_norm(tp,3))*(comp_levels.cytosol(2)*T3P_MID_norm(tp,2)));
    %%
    % $$\frac{{dFBP}_{m+4} }{dt}=\frac{1}{P_{FBP} }{\left(\overrightarrow{v_3 }
    % ?\left({T3P}_{m+3} {T3P}_{m+1} +{T3P}_{m+2} {T3P}_{m+2} \right)-\left(v_4 +-\overleftarrow{v_3
    % } \right)?{FBP}_{m+4} \right)}$$
    
    S_FBP(5,find(strcmp(model.rxns,'FBAp')))=((T3P_MID_norm(tp,4))*(T3P_MID_norm(tp,2)))+((T3P_MID_norm(tp,3))*(T3P_MID_norm(tp,3)));
    S_FBP(5,find(strcmp(model.rxns,'FBAc')))=((comp_levels.cytosol(2)*T3P_MID_norm(tp,4))*(comp_levels.cytosol(2)*T3P_MID_norm(tp,2)))+((comp_levels.cytosol(2)*T3P_MID_norm(tp,3))*(comp_levels.cytosol(2)*T3P_MID_norm(tp,3)));
    %%
    % $$\frac{{dFBP}_{m+5} }{dt}=\frac{1}{P_{FBP} }{\left(\overrightarrow{v_3 }
    % ?\left({T3P}_{m+3} {T3P}_{m+2} \right)-\left(v_4 +-\overleftarrow{v_3 } \right)?{FBP}_{m+5}
    % \right)}$$
    
    S_FBP(6,find(strcmp(model.rxns,'FBAp')))=(T3P_MID_norm(tp,4))*(T3P_MID_norm(tp,3));
    S_FBP(6,find(strcmp(model.rxns,'FBAc')))=(comp_levels.cytosol(2)*T3P_MID_norm(tp,4))*(comp_levels.cytosol(2)*T3P_MID_norm(tp,3));
    %%
    % $$\frac{{dFBP}_{m+6} }{dt}=\frac{1}{P_{FBP} }{\left(\overrightarrow{v_3 }
    % ?\left({T3P}_{m+3} T3P_{m+3} \right)-\left(v_4 +-\overleftarrow{v_3 } \right)?{FBP}_{m+6}
    % \right)}$$
    
    S_FBP(7,find(strcmp(model.rxns,'FBAp')))=(T3P_MID_norm(tp,4))*(T3P_MID_norm(tp,4));
    S_FBP(7,find(strcmp(model.rxns,'FBAc')))=(comp_levels.cytosol(2)*T3P_MID_norm(tp,4))*(comp_levels.cytosol(2)*T3P_MID_norm(tp,4));
    
    % *F6P*
    %
    % $$\frac{dF6P}{dt}=v_4 +\overleftarrow{{\mathit{\mathbf{v}}}_5 } +\overleftarrow{v_{11}
    % } -\overrightarrow{v_5 } -\overrightarrow{v_{11} }$$
    
    S_F6P(:,find(strcmp(model.rxns,'TK3')))=-F6P_MID_norm(tp,:); %plastid
    S_F6P(:,find(strcmp(model.rxns,'PGIp')))=-F6P_MID_norm(tp,:); %plastid
    S_F6P(:,find(strcmp(model.rxns,'PGIc')))=-comp_levels.cytosol(3)*F6P_MID_norm(tp,:); %cytosol
    S_F6P(:,find(strcmp(model.rxns,'SPS')))=-comp_levels.cytosol(3)*F6P_MID_norm(tp,:); %cytosol
   
    %%
    % $$\frac{{dF6P}_{m+0} }{dt}=\frac{1}{P_{F6P} }{\left(v_4 ?{FBP}_{m+0} +\overleftarrow{v_{11}
    % } ?{G6P}_{m+0} -\left(v_5 +\overrightarrow{v_{11} } \right)?{F6P}_{m+0} \right)}$$
    
    S_F6P(:,find(strcmp(model.rxns,'PFPp')))=FBP_MID_norm(tp,:); % FBP -> F6P
    S_F6P(:,find(strcmp(model.rxns,'PFPc')))=comp_levels.cytosol(4)*FBP_MID_norm(tp,:); % FBP -> F6P
    S_F6P(:,find(strcmp(model.rxns,'PGIp_rev')))=G6P_MID_norm(tp,:); %v11b
    S_F6P(:,find(strcmp(model.rxns,'PGIc_rev')))=comp_levels.cytosol(6)*G6P_MID_norm(tp,:); %v11b
    
    % combine submatrices created for individual metabolites
    if tp==1
        S_t0 = [S_RuBP; S_F6P; S_G6P; S_FBP; S_G1P; S_ADPG; S_UDPG];
    elseif tp==2
        S_t5 = [S_RuBP; S_F6P; S_G6P; S_FBP; S_G1P; S_ADPG; S_UDPG];
    elseif tp==3
        S_t10 = [S_RuBP; S_F6P; S_G6P; S_FBP; S_G1P; S_ADPG; S_UDPG];
    elseif tp==4
        S_t20 = [S_RuBP; S_F6P; S_G6P; S_FBP; S_G1P; S_ADPG; S_UDPG];
    end
end
end