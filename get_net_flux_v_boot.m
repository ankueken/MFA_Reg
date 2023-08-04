function v_a3_net_all = get_net_flux_v_boot(v_boot,n,model,model_rev)
%% calculate net flux for reversible reactions and give back reduced flux distribution
% v_a3 - flux distribution with from irreversible model
% n - number of rxns in final net flux model
% model - model with irreversible rxns only obtained from
% applying convertToIrreversible function
% model_rev - model before splitting

v_a3_net_all=nan(n,1000);

% convert those flaged as reversible in match field, lb<0 in original model
% created by network_INCA.m
for r=1:size(v_boot,2)
    v_a3=v_boot(:,r);
    v_a3_net=convertIrrevFluxDistribution(v_a3,model.match);

    % convert those reversible rxns that where splitted manually in
    % network_INCA.m since they enter matrix S modeling the MIDs
    idx=find(contains(model_rev.rxns,'_rev'));

    for i=1:length(idx)
        name=model_rev.rxns{idx(i)}(1:end-4);
        idx_fb=find(contains(model_rev.rxns,name));
        v_a3_net(idx_fb(1)) = v_a3_net(idx_fb(1))-v_a3_net(idx_fb(2));
        v_a3_net(idx_fb(2)) = nan;

    end

    rxnNames=model_rev.rxns(~isnan(v_a3_net));
    v_a3_net(isnan(v_a3_net))=[];

    if ~isempty(v_a3_net)
        v_a3_net_all(:,r)=v_a3_net;
    end
end
end