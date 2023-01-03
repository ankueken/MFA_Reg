%% create simple CBC model struct

model=createModel;

% CBC
model = addReaction(model,'T_CO2','metaboliteList',{'CO2'},'stoichCoeffList',[1], 'reversible',false);
model = addReaction(model,'RUBISCO_CO2','metaboliteList',{'RUBP_p' 'CO2' '3PGA_p'},'stoichCoeffList',[-1 -1 2], 'reversible',false);
model = addReaction(model,'RUBISCO_O2','metaboliteList',{'RUBP_p' '3PGA_p'},'stoichCoeffList',[-1 1], 'reversible',false);
model = addReaction(model,'GAPDHp','metaboliteList',{'3PGA_p' 'T3P_p'},'stoichCoeffList',[-1 1], 'reversible',false);
model = addReaction(model,'FBAp','metaboliteList',{'T3P_p' 'FBP_p'},'stoichCoeffList',[-2 1], 'reversible',false); %rev
model = addReaction(model,'PFPp','metaboliteList',{'FBP_p' 'F6P_p'},'stoichCoeffList',[-1 1], 'reversible',false); 
model = addReaction(model,'TK3','metaboliteList',{'F6P_p' 'E4P_p' 'EC2'},'stoichCoeffList',[-1 1 1], 'reversible',false); % !
model = addReaction(model,'ALD','metaboliteList',{'T3P_p' 'E4P_p' 'SBP'},'stoichCoeffList',[-1 -1 1], 'reversible',false);
model = addReaction(model,'SBPase','metaboliteList',{'SBP' 'S7P_p'},'stoichCoeffList',[-1 1], 'reversible',false);
model = addReaction(model,'TK1','metaboliteList',{'T3P_p' 'EC2' 'PP_p'},'stoichCoeffList',[-1 -1 1], 'reversible',true);
model = addReaction(model,'TK2','metaboliteList',{'S7P_p' 'R5P_p' 'EC2'},'stoichCoeffList',[-1 1 1], 'reversible',true);
model = addReaction(model,'PPI','metaboliteList',{'R5P_p' 'PP_p'},'stoichCoeffList',[-1 1], 'reversible',true);
model = addReaction(model,'PRK','metaboliteList',{'PP_p' 'RUBP_p'},'stoichCoeffList',[-1 1], 'reversible',false);

model = addReaction(model,'PGIp','metaboliteList',{'F6P_p' 'G6P_p'},'stoichCoeffList',[-1 1], 'reversible',false); %rev
model = addReaction(model,'PGMp','metaboliteList',{'G6P_p' 'G1P_p'},'stoichCoeffList',[-1 1], 'reversible',false); %rev
model = addReaction(model,'AGP','metaboliteList',{'G1P_p' 'ADPG'},'stoichCoeffList',[-1 1], 'reversible',false);
model = addReaction(model,'SS','metaboliteList',{'ADPG'},'stoichCoeffList',-1, 'reversible',false);

model = addReaction(model,'FBAp_rev','metaboliteList',{'FBP_p' 'TP_p'},'stoichCoeffList',[-1 2], 'reversible',false);
model = addReaction(model,'PGIp_rev','metaboliteList',{'G6P_p' 'F6P_p'},'stoichCoeffList',[-1 1], 'reversible',false);
model = addReaction(model,'PGMp_rev','metaboliteList',{'G1P_p' 'G6P_p'},'stoichCoeffList',[-1 1], 'reversible',false);
