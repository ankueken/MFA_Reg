function  [b, residuals, Chi, df, EXITFLAG, H] = QP_regression(S, x_meas_t,...
    d_meas, rubp_MIDs, f6p_MIDs, g6p_MIDs, fbp_MIDs, g1p_MIDs, adpg_MIDs, udpg_MIDs, model,sta_n,d_sta,suc_n,d_suc,tp_excluded)
%% regression for time points 0-40sec
%% INPUT:
% S - cell array with matrix combining MID data for time point i (obtained
%     from build_differential_matrix_comp)
% x_meas_t - MID data table (absolute data)
% d_meas - standard deviation table 
% [met_name]_MIDs - index of MIDs in data table
% model - model struct 
% sta_n - starch synthesis rate
% d_sta - std of starch synthesis rate
% suc_n - sucrose synthesis rate
% d_suc - std of sucrose synthesis rate
% tp_excluded (optional) - index of rows in MID/std data table to time points that should not be considered in fit  
%
%% OUTPUT:
% b - flux vector estimated
% residuals - difference MID slope
% Chi - chi-square statistic
% df - degrees of freedom
% EXITFLAG - exit status from solver
% H - MID and time point wise contribution to total chi-squared

%% reduce to used MIDs
% absolute
d_meas = d_meas{:,[rubp_MIDs, f6p_MIDs, g6p_MIDs, fbp_MIDs, g1p_MIDs, adpg_MIDs, udpg_MIDs]};
if exist('tp_excluded') 
    if ~isempty(tp_excluded)
        d_meas(tp_excluded,:) = 0;
    end
end
d_meas = d_meas(2:end,:);
d_meas = reshape(d_meas',[1,numel(d_meas)]);

d_tp = x_meas_t.LabelingTime_sec_(2:end)-x_meas_t.LabelingTime_sec_(1:end-1);

% absolute
x_meas_t=x_meas_t{:,[rubp_MIDs, f6p_MIDs, g6p_MIDs, fbp_MIDs, g1p_MIDs, adpg_MIDs, udpg_MIDs]};

num_MIDs=length([[rubp_MIDs, f6p_MIDs, g6p_MIDs, fbp_MIDs, g1p_MIDs, adpg_MIDs, udpg_MIDs]]);

% load underlying network
N=model.S;

idx_e=1:num_MIDs*4;
idx_b=(1:size(N,2))+idx_e(end);

LB_b=model.lb*1e3;
UB_b=model.ub*1e3;

LB_e=-1e6;UB_e=1e6;
%LB_e=-1e9;UB_e=1e9; % for INCA P

LB=[];
LB(idx_e)=ones(idx_e(end),1)*LB_e;% difference MID slope to measured
LB(idx_b)=LB_b;  % flux bounds
LB(end+1:end+2)=[LB_e;LB_e]; % difference estimated and measured v_exch

UB=[];
UB(idx_e)=ones(idx_e(end),1)*UB_e;
UB(idx_b)=UB_b;
UB(end+1:end+2)=[UB_e;UB_e];

%% steady state and MID constraints

% fit MIDs per time point
% x(t_i+dt_i) - x(t_i) <= e + S_ti*b

beq = [...
    reshape(x_meas_t(2:end,:)',[numel(x_meas_t(2:end,:)),1])- reshape(x_meas_t(1:end-1,:)',[numel(x_meas_t(1:end-1,:)),1]);
    zeros(size(N,1),1);
    sta_n;
    suc_n];

S_delta_t=[];
for temp = 1:length(S)
    S_delta_t = [S_delta_t; S{temp}*d_tp(temp)];
end
%           e                b
Aeq = [eye(num_MIDs*length(S)) S_delta_t;
       zeros(size(N,1),num_MIDs*length(S)) N];
   
Aeq(end+1,:)=zeros(1,size(Aeq,2));
Aeq(end,idx_b(find(strcmp(model.rxns,'SS')))) = 1;
Aeq(:,end+1) = zeros(size(Aeq,1),1);
Aeq(end,end) = -1;

Aeq(end+1,:)=zeros(1,size(Aeq,2));
Aeq(end,idx_b(find(strcmp(model.rxns,'SPS')))) = 1;
Aeq(:,end+1) = zeros(size(Aeq,1),1);
Aeq(end,end) = -1;

%% objective
% linear part
f = zeros(size(Aeq,2),1);

% quadratic part
% (multiply measurement standard deviation with P to transform from
% relative to absolute)
d_MIDs = 1./(2*((d_meas').^2));
d_MIDs(d_meas'==0) = 0;

%        e                                   b
Q_MID = [diag(d_MIDs)                  zeros(num_MIDs*length(S),size(N,2));  % e
         zeros(size(N,2),num_MIDs*length(S))   zeros(size(N,2))]; % b
     
Q_MID(end+1,:)=zeros(1,size(Q_MID,2));
Q_MID(:,end+1) = zeros(size(Q_MID,1),1);
Q_MID(end,end) = 1/(d_sta.^2);

Q_MID(end+1,:)=zeros(1,size(Q_MID,2));
Q_MID(:,end+1) = zeros(size(Q_MID,1),1);
Q_MID(end,end) = 1/(d_suc.^2);

%% solve problem
% using matlab quadprog
OPTIONS=optimset('quadprog');
OPTIONS.Display = 'off';
OPTIONS.MaxIter = 5000;
[X,FVAL,EXITFLAG,output,lambda] = quadprog(Q_MID,f,[],[],Aeq,beq,LB,UB,[],OPTIONS);

if EXITFLAG==1
    df = length(find(diag(Q_MID~=0)))-rank(full(N));
    residuals = X(idx_e(d_MIDs~=0));
    if exist('tp_excluded') 
        H=reshape((residuals.^2).*d_MIDs(d_MIDs~=0),length(d_MIDs(d_MIDs~=0))/(length(S)-length(tp_excluded)),(length(S)-length(tp_excluded)));
    else
        H=reshape((residuals.^2).*d_MIDs(d_MIDs~=0),length(d_MIDs(d_MIDs~=0))/length(S),length(S));
    end
    b=X(idx_b);
    Chi = FVAL*2;
else
    Chi=nan;df=nan;residuals=nan(sum(d_MIDs~=0),1);H=nan(sum(d_MIDs~=0));b=nan(length(idx_b),1);
end

end