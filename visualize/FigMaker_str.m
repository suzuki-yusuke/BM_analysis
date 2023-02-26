function [] = FigMaker_str(strObj,tbl,options,parameters)

% load data
names = fieldnames(strObj);
for i = 1:length(names)
    iT = getfield(strObj(i),names{i});
    eval([names{i}, '=iT;']);
end


%%
if all(ismember(unique(tbl.Phase), {'training'}))
    
    strs = unique(tbl_trn_sum.Feature,'stable');
    daynums = unique(tbl_trn_daily.Day);
    daynames = tbl_trn_sum.Properties.VariableNames;
    daynames = daynames(end-length(daynums)+1:end);
    
    groups = unique(tbl_trn_sum.Group);
    
    [~,~,DayIdx] = unique(daynums,'stable');
    
    % drow figure
    f = figure('Position', [10,10,1800,800]);
    
    c = parula(3);
    c = mat2cell(c,ones(size(c,1),1),3);

    
    space = 1;
    for i = 1:length(groups)
        
        dat = tbl_trn_sum{strcmp(tbl_trn_sum.Group,groups{i}), daynames};
        dat = reshape([dat{:}],size(dat));        
        dat = dat ./ repmat(sum(dat),size(dat,1),1);   
        

        subplot(1,length(groups),i)
        
        b = bar(dat', space, 'stacked');
        set(b, {'FaceColor'},c)
        
        switch options.maze_dia
            case 1
                step = 1;
            case 3
                step = 2;
        end
        
        set(gca, 'FontSize',26, 'LineWidth',.5,...
            'XLim',[.5, max(DayIdx)+.5],...
            'YLim', [0,1],...
            'XTick',DayIdx(1:step:end),...
            'XTickLabel', daynums(1:step:end))

%         b = area(dat',0);        
%         set(b, {'FaceColor'},c)
        set(gca, 'FontSize',38, 'LineWidth',.5,...
            'XLim',[0.5, max(DayIdx)+.5],...
            'YLim', [0,1],...
            'XTick',DayIdx(1:2:end),...
            'XTickLabel', daynums(1:2:end))
    
        legend(gca,strs, 'Location','northeastoutside', 'Color','w', 'FontSize',20)
        axis square
        xlabel('Day', 'interpreter','none', 'FontWeight','bold')
        ylabel('Proportion', 'interpreter','none', 'FontWeight','bold')
        
        title(groups{i}, 'FontSize',18, 'FontWeight','bold', 'Interpreter','none');
        hold off
        
    end
    
        
    set(gcf,'color','w');
    F = getframe(f); F = F.cdata;
    imwrite(F,[parameters.PrjDir,'\str_trn.png'])
    close(f)
    
end