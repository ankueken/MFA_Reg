function [Data_mean_abs, Data_std_abs_max, Data_mean_norm, Data_std_norm_max, AbsoluteLevel_mean, AbsoluteLevel_std, Data_std_abs_original] = load_experimental_data(scale)

% load experimental data MIDs [nmol/sample] 
% -> use fractions*pool size to get values in nmol/gDW
% calculate mean and std over samples
% for std = 0, max std over time 
%
% pool size measurements [nmol/gDW]
% calculate mean and std
%
% the function can take one optional input argument
%   if scale is set to 'true' the MID data will be scaled according to the
%   last time point

warning off
Data_MID=readtable('mixerIsotopomers.csv');
AbsoluteLevel=readtable('MixerAbsoluteLevels.csv');

%% normalize MID data, claculate MID data in nmol/gDW
Names=unique(cellfun(@(x) x(1:end-2),Data_MID.Properties.VariableNames(4:end),'UniformOutput',false));
Data_norm=Data_MID;
Data_nmol_gDW=Data_MID;

for i=1:length(Names)
    MIDs = find(contains(Data_MID.Properties.VariableNames,Names{i}));
    pool_size_idx = find(contains(AbsoluteLevel.Properties.VariableNames,Names{i},'IgnoreCase',true));
    
    Data_norm{:,MIDs}=Data_MID{:,MIDs}./sum(Data_MID{:,MIDs},2);
        
    if exist('scale') && scale==true
        relative_amount_a1 = Data_norm{13:15,MIDs(1)};
        relative_amount_a2 = Data_norm{28:30,MIDs(1)};
        relative_amount_a3 = Data_norm{43:45,MIDs(1)};
        relative_amount_a4 = Data_norm{58:60,MIDs(1)};
        Data_norm{1:15,MIDs(1)}=Data_norm{1:15,MIDs(1)}-repmat(Data_norm{13:15,MIDs(1)},5,1); % rescale chlamy
        Data_norm{16:30,MIDs(1)}=Data_norm{16:30,MIDs(1)}-repmat(Data_norm{28:30,MIDs(1)},5,1); % rescale soro
        Data_norm{31:45,MIDs(1)}=Data_norm{31:45,MIDs(1)}-repmat(Data_norm{43:45,MIDs(1)},5,1); % rescale ohadii
        Data_norm{46:60,MIDs(1)}=Data_norm{46:60,MIDs(1)}-repmat(Data_norm{58:60,MIDs(1)},5,1); % rescale ohadii eil
    end
           
    if ~isempty(pool_size_idx)
        Data_nmol_gDW{:,MIDs}=Data_norm{:,MIDs}.*AbsoluteLevel{:,pool_size_idx};
    else
        Data_nmol_gDW{:,MIDs}=nan(size((Data_MID{:,MIDs}./sum(Data_MID{:,MIDs},2))));
    end
end

%% calculate mean and std form MID replicates

repl_start = 1:3:size(Data_norm,1);
repl_end = 3:3:size(Data_norm,1);
Data_mean_norm=table();Data_mean_abs=table();Data_std_abs=table();
for i=1:length(repl_start)
    Data_mean_norm(i,1:size(Data_norm,2)) = Data_norm(repl_start(i),:);
    Data_mean_norm(i,4:size(Data_norm,2)) = array2table(mean(Data_norm{repl_start(i):repl_end(i),4:end}));
    Data_mean_abs(i,1:size(Data_nmol_gDW,2)) = Data_nmol_gDW(repl_start(i),:);
    Data_mean_abs(i,4:size(Data_nmol_gDW,2)) = array2table(mean(Data_nmol_gDW{repl_start(i):repl_end(i),4:end}));
    Data_std_abs(i,1:size(Data_nmol_gDW,2)) = Data_nmol_gDW(repl_start(i),:);
    Data_std_abs(i,4:size(Data_nmol_gDW,2)) = array2table(std(Data_nmol_gDW{repl_start(i):repl_end(i),4:end}));
    Data_std_norm(i,1:size(Data_norm,2)) = Data_norm(repl_start(i),:);
    Data_std_norm(i,4:size(Data_norm,2)) = array2table(std(Data_norm{repl_start(i):repl_end(i),4:end}));

end
Data_mean_norm.Properties.VariableNames=Data_norm.Properties.VariableNames;
Data_mean_abs.Properties.VariableNames=Data_nmol_gDW.Properties.VariableNames;
Data_std_abs.Properties.VariableNames=Data_nmol_gDW.Properties.VariableNames;
Data_std_norm.Properties.VariableNames=Data_norm.Properties.VariableNames;

Data_std_abs_original = Data_std_abs;
%% if std zero set it to minimum observed for that metabolite
algae_idx{1}=1:5;
algae_idx{2}=6:10;
algae_idx{3}=11:15;
algae_idx{4}=16:20;

temp=Data_std_abs{:,4:size(Data_std_abs,2)};
for i=1:size(temp,2)
    for a=1:4
        z=find(temp(algae_idx{a},i)==0);
        if ~isempty(z)
            temp_a=temp(algae_idx{a},i);
            min_std=max(temp_a(temp_a>0));
            if isempty(min_std)
                temp_a=temp(algae_idx{a},i-1);
                min_std=max(temp_a(temp_a>0));
            end
            temp_a(z)=min_std;
            temp(algae_idx{a},i)=temp_a;
        end
    end
end
Data_std_abs(:,4:size(Data_std_abs,2)) = array2table(temp);

%% use maximum std observed for a MID over time, assumption: no change in std over time
repl_start = 1:5:size(Data_std_abs,1);
repl_end =5:5:size(Data_std_abs,1);
Data_std_abs_max=table();
Data_std_abs_max(:,1:3) = Data_std_abs(:,1:3);

for i=1:length(repl_start)
    Data_std_abs_max(repl_start(i):repl_end(i),4:size(Data_std_abs,2)) = repmat(array2table(max(Data_std_abs{repl_start(i):repl_end(i),4:end})),repl_end(1),1);
end
Data_std_abs_max.Properties.VariableNames=Data_std_abs.Properties.VariableNames;

%% mean pool size over time points and replicates

% chlamy
AbsoluteLevel_mean(1,1:3) = AbsoluteLevel(1,1:3);
AbsoluteLevel_mean(1,4:size(AbsoluteLevel,2)) = array2table(mean(AbsoluteLevel{1:15,4:end}));
% sorokiniana
AbsoluteLevel_mean(2,1:3) = AbsoluteLevel(16,1:3);
AbsoluteLevel_mean(2,4:size(AbsoluteLevel,2)) = array2table(mean(AbsoluteLevel{16:30,4:end}));
% ohadii
AbsoluteLevel_mean(3,1:3) = AbsoluteLevel(31,1:3);
AbsoluteLevel_mean(3,4:size(AbsoluteLevel,2)) = array2table(mean(AbsoluteLevel{31:42,4:end}));
% ohadii eil
AbsoluteLevel_mean(4,1:3) = AbsoluteLevel(46,1:3);
AbsoluteLevel_mean(4,4:size(AbsoluteLevel,2)) = array2table(mean(AbsoluteLevel{46:60,4:end}));

AbsoluteLevel_mean.Properties.VariableNames=AbsoluteLevel.Properties.VariableNames;

AbsoluteLevel_std(1,1:3) = AbsoluteLevel(1,1:3);
AbsoluteLevel_std(1,4:size(AbsoluteLevel,2)) = array2table(std(AbsoluteLevel{1:15,4:end}));
% sorokiniana
AbsoluteLevel_std(2,1:3) = AbsoluteLevel(16,1:3);
AbsoluteLevel_std(2,4:size(AbsoluteLevel,2)) = array2table(std(AbsoluteLevel{16:30,4:end}));
% ohadii
AbsoluteLevel_std(3,1:3) = AbsoluteLevel(31,1:3);
AbsoluteLevel_std(3,4:size(AbsoluteLevel,2)) = array2table(std(AbsoluteLevel{31:42,4:end}));
% ohadii eil
AbsoluteLevel_std(4,1:3) = AbsoluteLevel(46,1:3);
AbsoluteLevel_std(4,4:size(AbsoluteLevel,2)) = array2table(std(AbsoluteLevel{46:60,4:end}));

AbsoluteLevel_std.Properties.VariableNames=AbsoluteLevel.Properties.VariableNames;

Data_std_norm_max=Data_std_abs_max;
for i=1:length(Names)
    MIDs = find(contains(Data_std_abs_max.Properties.VariableNames,Names{i}));
    pool_size_idx = find(contains(AbsoluteLevel_mean.Properties.VariableNames,Names{i},'IgnoreCase',true));
    
    if ~isempty(pool_size_idx)
        A=repmat(AbsoluteLevel_mean{:,pool_size_idx}',5,1);
        Data_std_norm_max{:,MIDs}=Data_std_abs_max{:,MIDs}./A(:);
    else
        Data_std_norm_max{:,MIDs}=nan(size(Data_std_norm_max{:,MIDs}));
    end
end


end