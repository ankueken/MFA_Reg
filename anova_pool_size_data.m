% perform anova to test if measured metabolite data reflect steady state
% for those metabolites whose MID enter the regression

AbsoluteLevel=readtable('MixerAbsoluteLevels.csv');

SST(1,1)=anova1(AbsoluteLevel.RuBP(1:15),AbsoluteLevel.LabelingTime_sec_(1:15),'off'); 
SST(1,2)=anova1(AbsoluteLevel.RuBP(16:30),AbsoluteLevel.LabelingTime_sec_(16:30),'off');
SST(1,3)=anova1(AbsoluteLevel.RuBP(31:45),AbsoluteLevel.LabelingTime_sec_(31:45),'off');
SST(1,4)=anova1(AbsoluteLevel.RuBP(31:42),AbsoluteLevel.LabelingTime_sec_(31:42),'off');
SST(1,5)=anova1(AbsoluteLevel.RuBP(46:60),AbsoluteLevel.LabelingTime_sec_(46:60),'off');

SST(2,1)=anova1(AbsoluteLevel.DHAP(1:15),AbsoluteLevel.LabelingTime_sec_(1:15),'off'); 
SST(2,2)=anova1(AbsoluteLevel.DHAP(16:30),AbsoluteLevel.LabelingTime_sec_(16:30),'off');
SST(2,3)=anova1(AbsoluteLevel.DHAP(31:45),AbsoluteLevel.LabelingTime_sec_(31:45),'off'); 
SST(2,4)=anova1(AbsoluteLevel.DHAP(31:42),AbsoluteLevel.LabelingTime_sec_(31:42),'off'); 
SST(2,5)=anova1(AbsoluteLevel.DHAP(46:60),AbsoluteLevel.LabelingTime_sec_(46:60),'off');

SST(3,1)=anova1(AbsoluteLevel.FBP(1:15),AbsoluteLevel.LabelingTime_sec_(1:15),'off'); 
SST(3,2)=anova1(AbsoluteLevel.FBP(16:30),AbsoluteLevel.LabelingTime_sec_(16:30),'off');
SST(3,3)=anova1(AbsoluteLevel.FBP(31:45),AbsoluteLevel.LabelingTime_sec_(31:45),'off'); 
SST(3,4)=anova1(AbsoluteLevel.FBP(31:42),AbsoluteLevel.LabelingTime_sec_(31:42),'off'); 
SST(3,5)=anova1(AbsoluteLevel.FBP(46:60),AbsoluteLevel.LabelingTime_sec_(46:60),'off');

SST(4,1)=anova1(AbsoluteLevel.F6P(1:15),AbsoluteLevel.LabelingTime_sec_(1:15),'off');
SST(4,2)=anova1(AbsoluteLevel.F6P(16:30),AbsoluteLevel.LabelingTime_sec_(16:30),'off');
SST(4,3)=anova1(AbsoluteLevel.F6P(31:45),AbsoluteLevel.LabelingTime_sec_(31:45),'off'); 
SST(4,4)=anova1(AbsoluteLevel.F6P(31:42),AbsoluteLevel.LabelingTime_sec_(31:42),'off'); 
SST(4,5)=anova1(AbsoluteLevel.F6P(46:60),AbsoluteLevel.LabelingTime_sec_(46:60),'off');

SST(5,1)=anova1(AbsoluteLevel.G6P(1:15),AbsoluteLevel.LabelingTime_sec_(1:15),'off');
SST(5,2)=anova1(AbsoluteLevel.G6P(16:30),AbsoluteLevel.LabelingTime_sec_(16:30),'off');
SST(5,3)=anova1(AbsoluteLevel.G6P(31:45),AbsoluteLevel.LabelingTime_sec_(31:45),'off'); 
SST(5,4)=anova1(AbsoluteLevel.G6P(34:42),AbsoluteLevel.LabelingTime_sec_(34:42),'off'); 
SST(5,5)=anova1(AbsoluteLevel.G6P(46:60),AbsoluteLevel.LabelingTime_sec_(46:60),'off');

SST(6,1)=anova1(AbsoluteLevel.PP(1:15),AbsoluteLevel.LabelingTime_sec_(1:15),'off'); 
SST(6,2)=anova1(AbsoluteLevel.PP(16:30),AbsoluteLevel.LabelingTime_sec_(16:30),'off');
SST(6,3)=anova1(AbsoluteLevel.PP(31:45),AbsoluteLevel.LabelingTime_sec_(31:45),'off'); 
SST(6,4)=anova1(AbsoluteLevel.PP(31:42),AbsoluteLevel.LabelingTime_sec_(31:42),'off'); 
SST(6,5)=anova1(AbsoluteLevel.PP(46:60),AbsoluteLevel.LabelingTime_sec_(46:60),'off');

SST(7,1)=anova1(AbsoluteLevel.G1P(1:15),AbsoluteLevel.LabelingTime_sec_(1:15),'off');
SST(7,2)=anova1(AbsoluteLevel.G1P(16:30),AbsoluteLevel.LabelingTime_sec_(16:30),'off');
SST(7,3)=anova1(AbsoluteLevel.G1P(31:45),AbsoluteLevel.LabelingTime_sec_(31:45),'off'); 
SST(7,4)=anova1(AbsoluteLevel.G1P(34:42),AbsoluteLevel.LabelingTime_sec_(34:42),'off'); 
SST(7,5)=anova1(AbsoluteLevel.G1P(46:60),AbsoluteLevel.LabelingTime_sec_(46:60),'off');

SST(8,1)=anova1(AbsoluteLevel.ADPG_Levels_nmol_gDW__(1:15),AbsoluteLevel.LabelingTime_sec_(1:15),'off');
SST(8,2)=anova1(AbsoluteLevel.ADPG_Levels_nmol_gDW__(16:30),AbsoluteLevel.LabelingTime_sec_(16:30),'off');
SST(8,3)=anova1(AbsoluteLevel.ADPG_Levels_nmol_gDW__(34:45),AbsoluteLevel.LabelingTime_sec_(34:45),'off'); 
SST(8,4)=anova1(AbsoluteLevel.ADPG_Levels_nmol_gDW__(34:42),AbsoluteLevel.LabelingTime_sec_(34:42),'off'); 
SST(8,5)=anova1(AbsoluteLevel.ADPG_Levels_nmol_gDW__(46:60),AbsoluteLevel.LabelingTime_sec_(46:60),'off');

SST(9,1)=anova1(AbsoluteLevel.UDPG(1:15),AbsoluteLevel.LabelingTime_sec_(1:15),'off');
SST(9,2)=anova1(AbsoluteLevel.UDPG(16:30),AbsoluteLevel.LabelingTime_sec_(16:30),'off');
SST(9,3)=anova1(AbsoluteLevel.UDPG(34:45),AbsoluteLevel.LabelingTime_sec_(34:45),'off'); 
SST(9,4)=anova1(AbsoluteLevel.UDPG(34:42),AbsoluteLevel.LabelingTime_sec_(34:42),'off'); 
SST(9,5)=anova1(AbsoluteLevel.UDPG(46:60),AbsoluteLevel.LabelingTime_sec_(46:60),'off');

Result_ANOVA_p_val=array2table(SST,'RowNames',{'RuBP'; 'DHAP'; 'FBP'; 'F6P'; 'G6P'; 'PP'; 'G1P'; 'ADPG'; 'UDPG'},...
    'VariableNames',{'C.reinhatdtii' 'C.sorokiniana' 'C.ohadii-LL-5tp' 'C.ohadii-LL-4tp' 'C.ohadii-EIL'});