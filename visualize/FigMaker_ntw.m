function [] = FigMaker_ntw(tbl,ntwObj,parameters,options)


%% load data
names = fieldnames(ntwObj);
for i = 1:length(names)
    iT = getfield(ntwObj(i),names{i});
    eval([names{i}, '=iT;']);
end

%% make figure
if strcmp(unique(tbl.Phase),{'probe_test'})    
    FigMaker_ntw_prb(tbl_prb_daily,parameters)
else
    FigMaker_ntw_trn_day_raw(tbl_trn_daily,options,parameters)
    FigMaker_ntw_trn_day_sum(tbl_trn_daily,options,parameters)
    FigMaker_ntw_trn_trial_sum(tbl,tbl_trn_daily,options,parameters)
%     FigMaker_ntw_trn_trial_sum_fit(tbl,tbl_trn_daily,options,parameters)
end

if ~isempty(options.champs)&&options.view_network_graph
    FigMaker_ntw_LocalNetwork(tbl,parameters,options)
    FigMaker_ntw_sortedLocalNetwork(tbl,parameters,options)
end

if options.view_global_network
    FigMaker_ntw_GlobalNetwork(tbl,parameters,options)
    FigMaker_ntw_sortedGlobalNetwork(tbl,parameters, options)
end




