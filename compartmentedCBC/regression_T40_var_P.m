function  [b, residuals, Chi, df, EXITFLAG, H, Pest] = regression_T40_var_P(S_t0, S_t5, S_t10, S_t20, x_meas_t,...
    d_meas_t0, d_meas_t5, d_meas_t10, d_meas_t20,d_meas_t40,...
    rubp_MIDs, f6p_MIDs, g6p_MIDs, fbp_MIDs, g1p_MIDs, adpg_MIDs, udpg_MIDs, model,sta_n,d_sta,suc_n,d_suc,...
    PRuBP, PF6P, PG6P, PFBP, PG1P, PADPG, PUDPG, d_PRuBP, d_PF6P, d_PG6P, d_PFBP, d_PG1P, d_PADPG, d_PUDPG)
%% regression with P as variable for time points 0-40sec
%% INPUT:
% S_ti - matrix combining MID data for time point i (obtained from
% x_meas_t - MID data table (relative data)
% build_differential_matrix_sCBC.m)
% d_meas_ti - standard deviation entering matrix Q_MID in objective
% [met_name]_MIDs - index of MIDs in data table
% model - model struct from network_simpleCBC.m
% sta_n - starch synthesis rate
% d_sta - std of starch synthesis rate
% P[met_name] - mean pool size measured
% d_[met_name] - std pool size measured
% suc_n - sucrose synthesis rate
% d_suc - std of sucrose synthesis rate
%
%% OUTPUT:
% b - flux vector estimated
% residuals - difference MID slope
% Chi - chi-square statistic
% df - degrees of freedom
% EXITFLAG - exit status from solver
% H - MID and time point wise contribution to total chi-squared
% Pest - estimated pool size

%% reduce to used MIDs
% absolute
d_meas_t0 = d_meas_t0{:,[rubp_MIDs, f6p_MIDs, g6p_MIDs, fbp_MIDs, g1p_MIDs, adpg_MIDs, udpg_MIDs]};
d_meas_t5 = d_meas_t5{:,[rubp_MIDs, f6p_MIDs, g6p_MIDs, fbp_MIDs, g1p_MIDs, adpg_MIDs, udpg_MIDs]};
d_meas_t10 = d_meas_t10{:,[rubp_MIDs, f6p_MIDs, g6p_MIDs, fbp_MIDs, g1p_MIDs, adpg_MIDs, udpg_MIDs]};
d_meas_t20 = d_meas_t20{:,[rubp_MIDs, f6p_MIDs, g6p_MIDs, fbp_MIDs, g1p_MIDs, adpg_MIDs, udpg_MIDs]};
d_meas_t40 = d_meas_t40{:,[rubp_MIDs, f6p_MIDs, g6p_MIDs, fbp_MIDs, g1p_MIDs, adpg_MIDs, udpg_MIDs]};

d_tp = x_meas_t.LabelingTime_sec_(2:end)-x_meas_t.LabelingTime_sec_(1:end-1);

num_MIDs=length([[rubp_MIDs, f6p_MIDs, g6p_MIDs, fbp_MIDs, g1p_MIDs, adpg_MIDs, udpg_MIDs]]);

num_P=7;

% load underlying network
% network_INCA
N=model.S;

idx_e=1:num_MIDs*4; % deviation slope MID from measured 
idx_b=(1:size(N,2))+idx_e(end); % flux 
idx_P=(1:num_P)+idx_b(end); % pool size 
idx_tau=(1:num_P)+idx_P(end); % deviation pool size to measurement

LB_b=model.lb*1e3;
UB_b=model.ub*1e3;

LB_e=-1e6;UB_e=1e6;

LB_P = zeros(num_P,1);
UB_P = ones(num_P,1)*1e6;

LB=[];
LB(idx_e)=ones(idx_e(end),1)*LB_e;
LB(idx_b)=LB_b;
LB(idx_P)=LB_P;
LB(idx_tau)=ones(num_P,1)*LB_e;
LB(end+1:end+2)=[LB_e;LB_e];

UB=[];
UB(idx_e)=ones(idx_e(end),1)*UB_e;
UB(idx_b)=UB_b;
UB(idx_P)=UB_P;
UB(idx_tau)=ones(num_P,1)*UB_e;
UB(end+1:end+2)=[UB_e;UB_e];

%% matrix of pool sizes (P)
% MIDs x MIDs
% rubp g1p g6p fbp udpg adpg glut
P = [eye(length(rubp_MIDs))*PRuBP               zeros(length(rubp_MIDs),length(f6p_MIDs))           zeros(length(rubp_MIDs),length(g6p_MIDs))   zeros(length(rubp_MIDs),length(fbp_MIDs)) zeros(length(rubp_MIDs),length(g1p_MIDs)) zeros(length(rubp_MIDs),length(adpg_MIDs)) zeros(length(rubp_MIDs),length(udpg_MIDs));
    zeros(length(f6p_MIDs),length(rubp_MIDs))   eye(length(f6p_MIDs))*PF6P                          zeros(length(f6p_MIDs),length(g6p_MIDs))    zeros(length(f6p_MIDs),length(fbp_MIDs))  zeros(length(f6p_MIDs),length(g1p_MIDs))  zeros(length(f6p_MIDs),length(adpg_MIDs))  zeros(length(f6p_MIDs),length(udpg_MIDs)) ;
    zeros(length(g6p_MIDs),length(rubp_MIDs))   zeros(length(g6p_MIDs),length(f6p_MIDs))            eye(length(g6p_MIDs))*PG6P                  zeros(length(g6p_MIDs),length(fbp_MIDs))  zeros(length(g6p_MIDs),length(g1p_MIDs))  zeros(length(g6p_MIDs),length(adpg_MIDs))  zeros(length(g6p_MIDs),length(udpg_MIDs)) ;
    zeros(length(fbp_MIDs),length(rubp_MIDs))   zeros(length(fbp_MIDs),length(f6p_MIDs))            zeros(length(fbp_MIDs),length(g6p_MIDs))    eye(length(fbp_MIDs))*PFBP                zeros(length(fbp_MIDs),length(g1p_MIDs))  zeros(length(fbp_MIDs),length(adpg_MIDs))  zeros(length(fbp_MIDs),length(udpg_MIDs)) ;
    zeros(length(g1p_MIDs),length(rubp_MIDs))   zeros(length(g1p_MIDs),length(f6p_MIDs))            zeros(length(g1p_MIDs),length(g6p_MIDs))    zeros(length(g1p_MIDs),length(fbp_MIDs))  eye(length(g1p_MIDs))*PG1P                zeros(length(g1p_MIDs),length(adpg_MIDs))  zeros(length(g1p_MIDs),length(udpg_MIDs)) ;
    zeros(length(adpg_MIDs),length(rubp_MIDs))  zeros(length(adpg_MIDs),length(f6p_MIDs))           zeros(length(adpg_MIDs),length(g6p_MIDs))   zeros(length(adpg_MIDs),length(fbp_MIDs)) zeros(length(adpg_MIDs),length(g1p_MIDs)) eye(length(adpg_MIDs))*PADPG               zeros(length(adpg_MIDs),length(udpg_MIDs));
    zeros(length(udpg_MIDs),length(rubp_MIDs))  zeros(length(udpg_MIDs),length(f6p_MIDs))           zeros(length(udpg_MIDs),length(g6p_MIDs))   zeros(length(udpg_MIDs),length(fbp_MIDs)) zeros(length(udpg_MIDs),length(g1p_MIDs)) zeros(length(udpg_MIDs),length(adpg_MIDs)) eye(length(udpg_MIDs))*PUDPG             ];
 
%% steady state and MID constraints 

% fit MIDs per time point
% - P*x(t_i+dt_i) + P*x(t_0) <= e - S_t0*v
% P*x(t_i+dt_i) - Px(t_0) <= e + S_t0*v

beq = [...
    zeros(size(N,1),1); ...
    zeros(num_MIDs*4,1);...
    PRuBP;PF6P;PG6P;PFBP;PG1P;PADPG;PUDPG;...
    sta_n;suc_n];

%    epsilon per MID and time point, v(t0), v(t5), [...], RuBP(t0), F6p(t0), [...]
Aeq = [zeros(size(N,1),num_MIDs*4) N zeros(size(N,1),num_P) zeros(size(N,1),num_P);
       eye(num_MIDs*4) [S_t0*d_tp(1);S_t5*d_tp(2);S_t10*d_tp(3);S_t20*d_tp(4)] [                            -x_meas_t{2,rubp_MIDs}'-x_meas_t{1,rubp_MIDs}' zeros(length(rubp_MIDs),num_P-1);...
                                                                                  zeros(length(f6p_MIDs),1) -x_meas_t{2,f6p_MIDs}'+x_meas_t{1,f6p_MIDs}' zeros(length(f6p_MIDs),num_P-2) ;
                                                                                  zeros(length(g6p_MIDs),2) -x_meas_t{2,g6p_MIDs}'+x_meas_t{1,g6p_MIDs}' zeros(length(g6p_MIDs),num_P-3) ;
                                                                                  zeros(length(fbp_MIDs),3) -x_meas_t{2,fbp_MIDs}'+x_meas_t{1,fbp_MIDs}' zeros(length(fbp_MIDs),num_P-4) ;
                                                                                  zeros(length(g1p_MIDs),4) -x_meas_t{2,g1p_MIDs}'+x_meas_t{1,g1p_MIDs}'  zeros(length(g1p_MIDs),num_P-5) ;
                                                                                  zeros(length(adpg_MIDs),5) -x_meas_t{2,adpg_MIDs}'+x_meas_t{1,adpg_MIDs}'  zeros(length(adpg_MIDs),num_P-6) ;
                                                                                  zeros(length(udpg_MIDs),6) -x_meas_t{2,udpg_MIDs}'+x_meas_t{1,udpg_MIDs}';
                                                                                                                                                                    
                                                                                                             -x_meas_t{3,rubp_MIDs}'-x_meas_t{2,rubp_MIDs}' zeros(length(rubp_MIDs),num_P-1);...
                                                                                  zeros(length(f6p_MIDs),1) -x_meas_t{3,f6p_MIDs}'+x_meas_t{2,f6p_MIDs}' zeros(length(f6p_MIDs),num_P-2) ;
                                                                                  zeros(length(g6p_MIDs),2) -x_meas_t{3,g6p_MIDs}'+x_meas_t{2,g6p_MIDs}' zeros(length(g6p_MIDs),num_P-3) ;
                                                                                  zeros(length(fbp_MIDs),3) -x_meas_t{3,fbp_MIDs}'+x_meas_t{2,fbp_MIDs}' zeros(length(fbp_MIDs),num_P-4) ;
                                                                                  zeros(length(g1p_MIDs),4) -x_meas_t{3,g1p_MIDs}'+x_meas_t{2,g1p_MIDs}'  zeros(length(g1p_MIDs),num_P-5) ;
                                                                                  zeros(length(adpg_MIDs),5) -x_meas_t{3,adpg_MIDs}'+x_meas_t{2,adpg_MIDs}'  zeros(length(adpg_MIDs),num_P-6) ;
                                                                                  zeros(length(udpg_MIDs),6) -x_meas_t{3,udpg_MIDs}'+x_meas_t{2,udpg_MIDs}';

                                                                                                             -x_meas_t{4,rubp_MIDs}'-x_meas_t{3,rubp_MIDs}' zeros(length(rubp_MIDs),num_P-1);...
                                                                                  zeros(length(f6p_MIDs),1) -x_meas_t{4,f6p_MIDs}'+x_meas_t{3,f6p_MIDs}' zeros(length(f6p_MIDs),num_P-2) ;
                                                                                  zeros(length(g6p_MIDs),2) -x_meas_t{4,g6p_MIDs}'+x_meas_t{3,g6p_MIDs}' zeros(length(g6p_MIDs),num_P-3) ;
                                                                                  zeros(length(fbp_MIDs),3) -x_meas_t{4,fbp_MIDs}'+x_meas_t{3,fbp_MIDs}' zeros(length(fbp_MIDs),num_P-4) ;
                                                                                  zeros(length(g1p_MIDs),4) -x_meas_t{4,g1p_MIDs}'+x_meas_t{3,g1p_MIDs}'  zeros(length(g1p_MIDs),num_P-5) ;
                                                                                  zeros(length(adpg_MIDs),5) -x_meas_t{4,adpg_MIDs}'+x_meas_t{3,adpg_MIDs}'  zeros(length(adpg_MIDs),num_P-6) ;
                                                                                  zeros(length(udpg_MIDs),6) -x_meas_t{4,udpg_MIDs}'+x_meas_t{3,udpg_MIDs}';

                                                                                                            -x_meas_t{5,rubp_MIDs}'-x_meas_t{4,rubp_MIDs}' zeros(length(rubp_MIDs),num_P-1);...
                                                                                  zeros(length(f6p_MIDs),1) -x_meas_t{5,f6p_MIDs}'+x_meas_t{4,f6p_MIDs}' zeros(length(f6p_MIDs),num_P-2) ;
                                                                                  zeros(length(g6p_MIDs),2) -x_meas_t{5,g6p_MIDs}'+x_meas_t{4,g6p_MIDs}' zeros(length(g6p_MIDs),num_P-3) ;
                                                                                  zeros(length(fbp_MIDs),3) -x_meas_t{5,fbp_MIDs}'+x_meas_t{4,fbp_MIDs}' zeros(length(fbp_MIDs),num_P-4) ;
                                                                                  zeros(length(g1p_MIDs),4) -x_meas_t{5,g1p_MIDs}'+x_meas_t{4,g1p_MIDs}'  zeros(length(g1p_MIDs),num_P-5) ;
                                                                                  zeros(length(adpg_MIDs),5) -x_meas_t{5,adpg_MIDs}'+x_meas_t{4,adpg_MIDs}'  zeros(length(adpg_MIDs),num_P-6) ;
                                                                                  zeros(length(udpg_MIDs),6) -x_meas_t{5,udpg_MIDs}'+x_meas_t{4,udpg_MIDs}'] zeros(num_MIDs*4,num_P); 
    zeros(num_P,num_MIDs*4)  zeros(num_P,size(N,2))  eye(num_P) -eye(num_P)];

Aeq(end+1,:)=zeros(1,size(Aeq,2));
Aeq(end,idx_b(find(strcmp(model.rxns,'SS')))) = 1;
Aeq(:,end+1) = zeros(size(Aeq,1),1);
Aeq(end,end) = 1;

Aeq(end+1,:)=zeros(1,size(Aeq,2));
Aeq(end,idx_b(find(strcmp(model.rxns,'SPS')))) = 1;
Aeq(:,end+1) = zeros(size(Aeq,1),1);
Aeq(end,end) = 1;

%% objective
% linear part
f = zeros(size(Aeq,2),1);

% quadratic part 
% (multiply measurement standard deviation with P to transform from
% relative to absolute)
d_MIDs = [1./(2*((d_meas_t5').^2)); 1./(2*((d_meas_t10').^2)); 1./(2*((d_meas_t20').^2)); 1./(2*((d_meas_t40').^2))];

d_P = zeros(size(1./[d_PRuBP^2;d_PF6P^2;d_PG6P^2;d_PFBP^2;d_PG1P^2;d_PADPG^2;d_PUDPG^2]));

%        e                                   v              P           tau
Q_MID = [diag(d_MIDs)                  zeros(num_MIDs*4,size(N,2))  zeros(num_MIDs*4,num_P)  zeros(num_MIDs*4,num_P);  % e
         zeros(size(N,2),num_MIDs*4)   zeros(size(N,2))             zeros(size(N,2),num_P)   zeros(size(N,2),num_P); % v
         zeros(num_P,num_MIDs*4)       zeros(num_P,size(N,2))       zeros(num_P)             zeros(num_P,num_P); % P
         zeros(num_P,num_MIDs*4)       zeros(num_P,size(N,2))       zeros(num_P,num_P)       diag(d_P)];  % tau

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
[X,FVAL,EXITFLAG] = quadprog(Q_MID,f,[],[],Aeq,beq,LB',UB',[],OPTIONS);

if EXITFLAG==1
    df = length(find(diag(Q_MID~=0)))-rank(full(N));
    residuals = X(idx_e(d_MIDs~=0));
    H=reshape((residuals.^2).*d_MIDs(d_MIDs~=0),length(d_MIDs(d_MIDs~=0))/4,4);
    b=X(idx_b);
    Pest=X(idx_P);
    Chi = FVAL*2;
else
    Chi=nan;df=nan;residuals=nan(sum(d_MIDs~=0),1);H=nan(sum(d_MIDs~=0));b=nan(length(idx_b),1);Pest=nan(length(idx_P),1);
end

end