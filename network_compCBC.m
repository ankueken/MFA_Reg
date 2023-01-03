%% Build the compartmented model

model=createModel;

% CBC
model = addReaction(model,'T_CO2','metaboliteList',{'CO2'},'stoichCoeffList',[1], 'reversible',false);
model = addReaction(model,'RUBISCO_CO2','metaboliteList',{'RUBP_p' 'CO2' '3PGA_p'},'stoichCoeffList',[-1 -1 2], 'reversible',false);
model = addReaction(model,'RUBISCO_O2','metaboliteList',{'RUBP_p' '3PGA_p' '2PG_p'},'stoichCoeffList',[-1 1 1], 'reversible',false);
model = addReaction(model,'Gapdhp','metaboliteList',{'3PGA_p' 'T3P_p'},'stoichCoeffList',[-1 1], 'reversible',false);
model = addReaction(model,'FBAp','metaboliteList',{'T3P_p' 'FBP_p'},'stoichCoeffList',[-2 1], 'reversible',false); %rev
model = addReaction(model,'PFPp','metaboliteList',{'FBP_p' 'F6P_p'},'stoichCoeffList',[-1 1], 'reversible',false); 
model = addReaction(model,'TK3','metaboliteList',{'F6P_p' 'E4P_p' 'EC2'},'stoichCoeffList',[-1 1 1], 'reversible',false); 
model = addReaction(model,'TK3_rev','metaboliteList',{'E4P_p' 'EC2' 'F6P_p' },'stoichCoeffList',[-1 -1 1], 'reversible',false); 
model = addReaction(model,'ALD','metaboliteList',{'T3P_p' 'E4P_p' 'SBP'},'stoichCoeffList',[-1 -1 1], 'reversible',false);
model = addReaction(model,'SBPase','metaboliteList',{'SBP' 'S7P_p'},'stoichCoeffList',[-1 1], 'reversible',false);
model = addReaction(model,'TK1','metaboliteList',{'T3P_p' 'EC2' 'PP_p'},'stoichCoeffList',[-1 -1 1], 'reversible',true);
model = addReaction(model,'TK2','metaboliteList',{'S7P_p' 'R5P_p' 'EC2'},'stoichCoeffList',[-1 1 1], 'reversible',true);
model = addReaction(model,'PPI','metaboliteList',{'R5P_p' 'PP_p'},'stoichCoeffList',[-1 1], 'reversible',true);
model = addReaction(model,'PRK','metaboliteList',{'PP_p' 'RUBP_p'},'stoichCoeffList',[-1 1], 'reversible',false);

% respiration
model = addReaction(model,'PGP','metaboliteList',{'2PG_p' 'GLY_p'},'stoichCoeffList',[-1 1], 'reversible',false);
model = addReaction(model,'GDC','metaboliteList',{'GLY_p' 'SER_p' 'CO2'},'stoichCoeffList',[-2 1 1], 'reversible',false);
model = addReaction(model,'SGA1','metaboliteList',{'SER_p' 'GA_p'},'stoichCoeffList',[-1 1], 'reversible',false);
model = addReaction(model,'GK','metaboliteList',{'GA_p' '3PGA_p'},'stoichCoeffList',[-1 1], 'reversible',true);

% starch biosynthesis
model = addReaction(model,'PGIp','metaboliteList',{'F6P_p' 'G6P_p'},'stoichCoeffList',[-1 1], 'reversible',false); %rev
model = addReaction(model,'PGMp','metaboliteList',{'G6P_p' 'G1P_p'},'stoichCoeffList',[-1 1], 'reversible',false); %rev
model = addReaction(model,'AGP','metaboliteList',{'G1P_p' 'ADPG'},'stoichCoeffList',[-1 1], 'reversible',false);
model = addReaction(model,'SS','metaboliteList',{'ADPG'},'stoichCoeffList',-1, 'reversible',false);

% sucrose synthesis cytosol
model = addReaction(model,'FBAc','metaboliteList',{'T3P_c' 'FBP_c'},'stoichCoeffList',[-2 1], 'reversible',false); %rev
model = addReaction(model,'PFPc','metaboliteList',{'FBP_c' 'F6P_c'},'stoichCoeffList',[-1 1], 'reversible',false); 
model = addReaction(model,'PGIc','metaboliteList',{'F6P_c' 'G6P_c'},'stoichCoeffList',[-1 1], 'reversible',false); %rev
model = addReaction(model,'PGMc','metaboliteList',{'G6P_c' 'G1P_c'},'stoichCoeffList',[-1 1], 'reversible',false); %rev
model = addReaction(model,'GPU','metaboliteList',{'G1P_c' 'UDPG'},'stoichCoeffList',[-1 1], 'reversible',false); %rev
model = addReaction(model,'SPS','metaboliteList',{'F6P_c' 'UDPG'},'stoichCoeffList',[-1 -1], 'reversible',false);

% TCA cycle
model = addReaction(model,'PGAMc','metaboliteList',{'3PGA_c' 'PEP_c'},'stoichCoeffList',[-1 1], 'reversible',false);
model = addReaction(model,'PKc','metaboliteList',{'PEP_c' },'stoichCoeffList',[-1 ], 'reversible',false);

% transport
model = addReaction(model,'T_3PGA','metaboliteList',{'3PGA_p' '3PGA_c'},'stoichCoeffList',[-1 1], 'reversible',true);
model = addReaction(model,'T_TP','metaboliteList',{'T3P_p' 'T3P_c'},'stoichCoeffList',[-1 1], 'reversible',true);
model = addReaction(model,'T_PEP','metaboliteList',{'PEP_p' 'PEP_c'},'stoichCoeffList',[-1 1], 'reversible',true);


% FA
model = addReaction(model,'PGAMp','metaboliteList',{'3PGA_p' 'PEP_p'},'stoichCoeffList',[-1 1], 'reversible',false);
model = addReaction(model,'PKp','metaboliteList',{'PEP_p' },'stoichCoeffList',[-1], 'reversible',false);

%% reversible reactions
model = addReaction(model,'FBAp_rev','metaboliteList',{'FBP_p' 'T3P_p'},'stoichCoeffList',[-1 2], 'reversible',false);

model = addReaction(model,'PGIp_rev','metaboliteList',{'G6P_p' 'F6P_p'},'stoichCoeffList',[-1 1], 'reversible',false);
model = addReaction(model,'PGMp_rev','metaboliteList',{'G1P_p' 'G6P_p'},'stoichCoeffList',[-1 1], 'reversible',false);
% 
model = addReaction(model,'FBAc_rev','metaboliteList',{'FBP_c' 'T3P_c'},'stoichCoeffList',[-1 2], 'reversible',false);
model = addReaction(model,'PGIc_rev','metaboliteList',{'G6P_c' 'F6P_c'},'stoichCoeffList',[-1 1], 'reversible',false);
model = addReaction(model,'PGMc_rev','metaboliteList',{'G1P_c' 'G6P_c'},'stoichCoeffList',[-1 1], 'reversible',false);

