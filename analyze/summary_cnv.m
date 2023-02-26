function [filename,sheetnames] = summary_cnv(parameters,options,tbl)

%%
fprintf('\nAnalyze conventional indices......DONE\n')
cnvObj = analyze_cnv(parameters, tbl, options.order);

if options.order{'Conventional','Output'}    
    
    if ismember(unique(tbl.Phase),{'training','retraining','reversal_training'})
        sheetnames = {...
            'training(daily_mean)', 'training(summarize)',...
            'training(ANOVA)',...
            'training(HSD(Dxgn))','training(HSD(Gxdn))'...
            'training(HSD(G))','training(HSD(D))'};
        b = 1;
    elseif strcmp(unique(tbl.Phase),{'probe_test'})
        sheetnames = {...
            'probetest(daily_mean)', 'probetest(summarize)',...
            'probetest(ANOVA)',...
            'probetest(HSD(Hxgn))','probetest(HSD(Gxhn))'...
            'probetest(HSD(G))','probetest(HSD(H))'};
        b = 4;
    else
        b = 1;
    end
    
    filename = fullfile(parameters.PrjDir, 'ConventionalScore.xls');
    names = fieldnames(cnvObj);
    empt = [];
    a = 0; c = 0;    
    for i = b:length(names)
        
        Ti = getfield(cnvObj(i), names{i});
        
        if isstruct(Ti)
            names_Tj = fieldnames(Ti);
            for j = 1:length(names_Tj)
                Tj = getfield(Ti, {j}, names_Tj{j});
                a = a + 1;
                if isempty(Tj)
                    empt = [empt, a];
                else
                    c = c + 1;
                    dat = [Tj.Properties.VariableNames; table2cell(Tj)];
                    xlswrite(filename, dat, c)
                end
            end
            
        else
            a = a + 1;
            if isempty(Ti)
                empt = [empt, a];
            else
                c = c + 1;
                dat = [Ti.Properties.VariableNames; table2cell(Ti)];
                xlswrite(filename, dat, c)
            end
        end
    end
    
    sheetnames(empt) = [];
    if max(empt)==3
        sheetnames(max(empt):end) = [];
    end
    
    
%     % excel の各シートの名前を変更する
%     % http://goo.gl/8S2QXW    
    e = actxserver('Excel.Application'); % # open Activex server
    ewb = e.Workbooks.Open(filename); % # open file (enter full path!)
    for i = 1:length(sheetnames)
        ewb.Worksheets.Item(i).Name = sheetnames{i}; % rename #th sheet
        ewb.Save % # save to the same file
    end
    ewb.Close(false)
    e.Quit
    
    
    fprintf('\nOutput......DONE\n')
else
    fprintf('\nOutput......SKIP\n')
end


if options.order{'Conventional','View'}
    FigMaker_cnv(tbl, cnvObj, parameters, options)
    fprintf('\nVisualize......DONE\n')
else
    fprintf('\nVisualize......SKIP\n')
end