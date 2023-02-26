function [] = FigMaker_cnv_trn_day_raw(tbl_trn_daily,options,parameters)

%%

switch options.maze_dia
    case 1
        step = 1;
        dcm = 300;
    case 3
        step = 2;
        dcm = 1000;
end

rownames_daily = tbl_trn_daily.Properties.VariableNames;
cnvs = rownames_daily(end-2:end);
[groups,~,groupIdx] = unique(tbl_trn_daily.Group,'stable');
[daynums,~,dayIdx] = unique(tbl_trn_daily.Day);
NumDays = length(daynums);
NumSbj = arrayfun(@(x)length(find(groupIdx==x)), unique(groupIdx), 'Uniform', false);
NumSbj = cell2mat(NumSbj)./NumDays;

% drow figure for training
ylabelnames = {'No. of errors (times)', 'Latency (s)', 'Travel distance (cm)'};
f = figure('Position', [10,10,1800,640]);

cmap = colormap(lines(length(groups)));

for i = 1:length(cnvs)
    
    rawdata = cell(length(groups),1);
    ind = cell(2,1);
    catcolors = cell(1,length(groups));
    for j = 1:length(groups)
        a = tbl_trn_daily{strcmp(tbl_trn_daily.Group,groups(j)), cnvs(i)};
        b = [groupIdx(strcmp(tbl_trn_daily.Group, groups(j))),...
            dayIdx(strcmp(tbl_trn_daily.Group, groups(j)))+min(daynums)-1];
%         tf = ~isnan(a)&all(b>0,2);
        tf = all(b>0,2);
        rawdata(j) = {a(tf)};
        ind(j) = {b(tf,:)};
        catcolors(j) = {cmap(j,:)};
    end
    ind = vertcat(ind{:});
    rawdata = vertcat(rawdata{:});
    k = find(isnan(rawdata),1);
      
    catmarkers = repmat({'.'},1,length(groups));
    
    
    subplot(1, length(cnvs), i);
    [~,pos] = plotSpread(rawdata, 'distributionIdx',ind(:,2),...
        'categoryIdx',ind(:,1), 'categoryColors',catcolors,...
        'categoryMarkers', catmarkers, 'showMM',0);
    
    if ~isempty(k)
        pos = [pos(1:k-1,:);nan(size(rawdata,1)-size(pos,1),2);pos(k:end,:)];
    end
    
    
    % plot raw value    
    X = mat2cell(pos(:,1), NumSbj.*NumDays, 1);
    X = cellfun(@(x,y) reshape(x,y,NumDays), X,num2cell(NumSbj), 'uniform',0);
    X = vertcat(X{:})';
    Y = mat2cell(pos(:,2), NumSbj.*NumDays, 1);
    Y = cellfun(@(x,y) reshape(x,y,NumDays), Y,num2cell(NumSbj), 'uniform',0);
    E = cellfun(@(x) median(x,1,'omitnan'), Y, 'uniform',0);
    E = vertcat(E{:});
    Y = vertcat(Y{:})';
    
    mcolor_raw = mat2cell(lines(length(groups)),ones(length(groups),1),3);
    mcolor_raw = cellfun(@(x,y) repmat(x,y,1), mcolor_raw,num2cell(NumSbj), 'uniform',0);
    mcolor_raw = vertcat(mcolor_raw{:});
    mcolor_raw = mat2cell(mcolor_raw,ones(sum(NumSbj),1),3);
    lcolor = cellfun(@(x) x+(1-x).*.8, mcolor_raw, 'uniform',0);
    
    hold on
    % plot raw value
    h1 = plot(X,Y,'-');
    set(h1, {'Color'},lcolor);

    % plot representative value
    h2 = plot(daynums,E', '-', 'LineWidth',3);
    mcolor_rep = mat2cell(lines(length(groups)).*.7,ones(length(groups),1),3);
    set(h2, {'Color'},mcolor_rep)    
    
    % plot raw value
    c = vertcat(catcolors{:});
    scatter(pos(:,1),pos(:,2), 13, c(ind(:,1),:), 'fill')
        
    
    ylabel(ylabelnames{i}, 'Interpreter', 'latex')
    xtick = daynums(1:step:end);
    xticklabels = num2cell(xtick); xticklabels = cellfun(@num2str,xticklabels, 'uniform',0);
    ymax = reshape(rawdata, length(unique(tbl_trn_daily.SN)), NumDays);
    ymax = max(quantile(ymax,.95));
    
    ax = gca;
    ax.LineWidth = 3;
    
    set(ax, 'FontSize', 18, 'LineWidth', 2,...
        'XLim',[min(daynums)-.5, max(daynums)+.5],...        
        'XTick',xtick, 'XTicklabel',xticklabels);
    
    switch i
        case 1
            legend(h1, groups, 'FontSize',30, 'Interpreter','none'); legend('boxoff')
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
imwrite(F,[parameters.PrjDir,'\cnv_trn_raw.png'])
close(f)

end