function [] = FigMaker_cnv(tbl,cnvObj,parameters,options)


%% load data
names = fieldnames(cnvObj);
for i = 1:length(names)
    iT = getfield(cnvObj(i),names{i});
    eval([names{i}, '=iT;']);
end

%% make figure
if ismember(unique(tbl.Phase), {'training'})
    FigMaker_cnv_trn_day_raw(tbl_trn_daily,options,parameters)
    FigMaker_cnv_trn_day_sum(tbl_trn_daily,options,parameters)
    FigMaker_cnv_trn_trial_sum(tbl,tbl_trn_daily,options,parameters)
    FigMaker_cnv_trn_trial_sum_fit(tbl,tbl_trn_daily,options,parameters)
else % probe test
    FigMaker_cnv_prb_raw(tbl_prb_daily,parameters,options)
    FigMaker_cnv_prb_sum(tbl_prb_daily,parameters,options)
end