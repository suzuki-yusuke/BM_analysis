function [] = FigMaker_cnv_trn_trial_sum(tbl,tbl_trn_daily,options,parameters)

%%

switch options.maze_dia
    case 1
        step = 1;
        dcm = 200;
    case 3
        step = 2;
        dcm = 500;
end

rownames_daily = tbl_trn_daily.Properties.VariableNames;
cnvs = rownames_daily(end-2:end);
[groups,~,groupIdx] = unique(tbl.Group,'stable');
[daynums,~,dayIdx] = unique(tbl.Day);
[trials,~,trialIdx] = unique(tbl.Trial);
[distIdx,~,~] = unique([dayIdx,trialIdx],'rows','stable');
T = length(daynums)*length(trials);

NumSbj = arrayfun(@(x)length(find(groupIdx==x)),unique(groupIdx),'uniform',0);
NumSbj = cell2mat(NumSbj)./(length(daynums)*length(trials));


% draw figure for training
ylabelnames = {'No. of errors (times)', 'Latency (s)', 'Travel distance (cm)'};
f = figure('Position', [10,10,2200,640]);
cmap = mat2cell(lines(length(groups)),ones(length(groups),1),3);

trialnums = 1:(length(trials)*length(daynums));
for i = 1:length(cnvs)
    
    MED = zeros(length(groups),T);
    Q = cell(length(groups),1);
    I = Q;
    for j = 1:length(groups)

        a = tbl{strcmp(tbl.Group,groups{j}),cnvs{i}};
        a = reshape(a,NumSbj(j),T);
        
        % smoothing
        if sum(isnan(a(1,:)))>0
            b = arrayfun(@(x) smooth(a(x,1:find(isnan(a(x,:)),1)-1))',1:size(a,1),'uniform',0);            
            b = cat(1,b{:});
            a = [b,nan(size(b,1),sum(isnan(a(1,:))))];
            
        else
            b = arrayfun(@(x) smooth(a(x,:))',1:size(a,1),'uniform',0);
            b = cat(1,b{:});
            a = b;
        end
        
        
        MED(j,:) = median(a,1,'omitnan');
        
%         q = quantile(a,[.25,.75],1);
        q = [median(a,1,'omitnan')-mad(a,1); median(a,1,'omitnan')+mad(a,1)];
        
        
        [~,ind] = find(isnan(q));
        ind = unique(ind);
        q(:,ind) = [];        
                
        Q(j) = {[q(2,:),fliplr(q(1,:))]};        
        I(j) = {ind'};
    
    end

    subplot(1, length(cnvs), i);
    
    h0 = cellfun(@(x,y) patch([trialnums(~ismember(trialnums,y)),...
        fliplr(trialnums(~ismember(trialnums,y)))],x,'r',...
        'FaceAlpha',.3,'EdgeColor','none'),Q,I);
    set(h0,{'FaceColor'},cmap)
    hold on
    h1 = plot(1:T,MED,'-','LineWidth',3);
    set(h1,{'Color'},cmap)    
    
    ylabel(ylabelnames{i}, 'Interpreter','latex')
    xtick = 1:step*length(trials):T;
    xticklabels = num2cell(distIdx(xtick));
    xticklabels = cellfun(@num2str,xticklabels, 'uniform',0);
    
    set(gca, 'FontSize', 18, 'LineWidth', 3,...
        'XLim',[.5, T+.5],...
        'XTick',xtick, 'XTicklabel',xticklabels);
    
    switch i
        case 1
            legend(h1, groups, 'FontSize',18, 'Interpreter','none'); legend('boxoff')
            xlabel('Day', 'Interpreter','latex');
        case 2
            yticklabels = get(gca,'YTickLabels');
            yticks = 0:60:(60*floor(max(cellfun(@(x) str2double(x),yticklabels))/60)+60);
            yticklabels = arrayfun(@(x) num2str(x),yticks,'uniform',0);
            set(gca,'YTickLabel',yticklabels,'YTick',yticks);
        case 3
            yticklabels = get(gca,'YTickLabels');
            yticks = cellfun(@(x) str2double(x), yticklabels);
            yticks = 0:dcm:(100*floor(yticks(end)/100)+100);
            yticklabels = arrayfun(@(x) num2str(x),yticks,'uniform',0);
            set(gca,'YTickLabel',yticklabels,'YTick',yticks,'YLim',yticks([1,end]));
    end
    
    axis square
    box off
    
    hold off
    
end

F = getframe(f); F = F.cdata;
imwrite(F,[parameters.PrjDir,'\cnv_trn_summ_trial.png'])
close(f)

end