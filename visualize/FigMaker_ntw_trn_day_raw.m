function [] = FigMaker_ntw_trn_day_raw(tbl_trn_daily,options,parameters)


switch options.maze_dia
    case 1
        step = 1;
    case 3
        step = 2;
end


%%
logdisp = 0; % ëŒêîï\é¶


rownames_daily = tbl_trn_daily.Properties.VariableNames;
ntws = rownames_daily(4:end);
ntws = ntws([8,9,11,13,12,10,14,15,1,3,5,4,2,6,7]);


% drow figure for training
f2 = figure('Position', [50,50,1480,900]);
f1 = figure('Position', [50,50,1480,900]);
ynames = {'Order','Shortest path length','Degree',...
    'Clustering coefficient','Density', 'Betweenness centrality', 'Closeness centrality',...
    'No. of stops','Order','Shortest path length','Degree',...
    'Clustering coefficient','Density','Betweenness centrality', 'Closeness centrality'};
toolong = ynames([2,4,6,7,10,12,14,15]);
ynames = ynames([8,9,11,13,12,10,14,15,1,3,5,4,2,6,7]);

yunits = {'$\overline{o}$','$\overline{n}$','$\overline{m}$','$\overline{\rho}$',...
    '$\overline{C_{WS}}$','$\overline{l}$','$\overline{x}$','$\overline{cl}$',...
    '$\overline{n}$','$\overline{m}$','$\overline{\rho}$',...
    '$\overline{C_{WS}}$','$\overline{l}$','$\overline{x}$','$\overline{cl}$'};



[groups,~,groupIdx] = unique(tbl_trn_daily.Group,'stable');
[daynums,~,dayIdx] = unique(tbl_trn_daily.Day);
NumDays = length(daynums);

ii2 = 1; ii1 = 0;

cmap = colormap(lines(length(groups)));

for i = 1:length(ntws)
    
    
    if ~isempty(strfind(ntws{i},'1'))
        figure(f2)
        ii2 = ii2 + 1; ii = ii2;
    else
        figure(f1)
        ii1 = ii1 + 1; ii = ii1;
    end
    
    %
    rawdata = cell(length(groups),1);
    ind = cell(2,1);
    catcolors = cell(1,length(groups));
    for j = 1:length(groups)
        a = tbl_trn_daily{strcmp(tbl_trn_daily.Group,groups(j)), ntws(i)};
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
    
    
    % ëŒêîï\é¶
    if logdisp
        rawdata = log(rawdata);
        if all(abs(rawdata(:))==Inf)
            continue
        else
            rawdata(isinf(rawdata)) = min(rawdata(~isinf(rawdata(:))));
        end
    end        
    
    
    subplot(3,3,ii);
    
    [~,pos] = plotSpread(rawdata, 'distributionIdx',ind(:,2),...
        'categoryIdx',ind(:,1), 'categoryColors',catcolors,...
        'categoryMarkers',catmarkers, 'showMM',0);
    
    if ~isempty(k)
        pos = [pos(1:k-1,:);nan(size(rawdata,1)-size(pos,1),2);pos(k:end,:)];
    end
    
    
    % plot raw value    
    NumSbj = histcounts(ind(:,1),.5:max(ind(:,1))+.5)./NumDays;
    X = mat2cell(pos(:,1), NumSbj.*NumDays,1);
    X = cellfun(@(x,y) reshape(x,y,NumDays), X,num2cell(NumSbj)', 'uniform',0);
    X = vertcat(X{:})';
    Y = mat2cell(pos(:,2), NumSbj.*NumDays,1);
    Y = cellfun(@(x,y) reshape(x,y,NumDays), Y,num2cell(NumSbj)', 'uniform',0);
    E = cellfun(@(x) median(x,1,'omitnan'), Y, 'uniform',0);
    E = vertcat(E{:});
    Y = vertcat(Y{:})';
    
    mcolor_raw = mat2cell(lines(length(groups)),ones(length(groups),1),3);
    mcolor_raw = cellfun(@(x,y) repmat(x,y,1), mcolor_raw,num2cell(NumSbj)', 'uniform',0);
    mcolor_raw = vertcat(mcolor_raw{:});
    mcolor_raw = mat2cell(mcolor_raw,ones(sum(NumSbj),1),3);
    lcolor = cellfun(@(x) x+(1-x).*.8, mcolor_raw, 'uniform',0);
    
    h1 = plot(X,Y, '-');
    set(h1, {'Color'},lcolor)
    hold on
    
    
    % plot representative value
    h2 = plot(daynums,E', '-', 'LineWidth',3);
    mcolor_rep = mat2cell(lines(length(groups)).*.7,ones(length(groups),1),3);
    set(h2, {'Color'},mcolor_rep)

    % plot raw value
    c = vertcat(catcolors{:});
    scatter(pos(:,1),pos(:,2), 13, c(ind(:,1),:), 'fill')
    
    if ismember(ntws{i}, {'n_1','o_2'})
        legend(h2, groups, 'FontSize',18, 'Interpreter','none'); legend('boxoff')
        xlabel('Day', 'Interpreter','latex');
    else
        xlabel('');
    end
    

    
    if ismember(ynames(i), toolong)
        str = ynames{i};
        k = strfind(str, ' '); k = k(1);
        str = {str(1:k(1)),str(k(1)+1:end)};
        ylabel([str(1), [str{2},' (',yunits{i},')']], 'Interpreter','latex')
    else
        ylabel([ynames{i}, ' (', yunits{i}, ')'], 'Interpreter','latex')
    end
    
    if logdisp
        Y = rawdata;
    else
        Y = tbl_trn_daily{:, ntws{i}};
    end
    
    if iscell(Y)
        Y = cell2mat(Y);
    end
    
    if min(Y)==max(Y)
        ylim = [-1,1];
    else
        ylim = quantile(Y, [.01,.99]);
    end
    xlim = [min(daynums)-.5, max(daynums)+.5];
    
    xtick = daynums(1:step:end);
    xticklabels = num2cell(xtick); xticklabels = cellfun(@num2str,xticklabels, 'uniform',0);
    set(gca, 'FontSize', 18, 'LineWidth', 2,...
        'XLim',xlim, 'XTick',xtick, 'XTickLabel',xticklabels, 'YLim',ylim);
    axis square
    box off
    
    hold off
    
end

set(f2, 'Position',[100,20,1300,1000])
set(f1, 'Position',[100,20,1300,1000])
hold off
figure(f2)
figure(f1)

F = getframe(f1); F = F.cdata;
imwrite(F,[parameters.PrjDir,'\ntw-d_trn_raw.png'])
close(f1)

F = getframe(f2); F = F.cdata;
imwrite(F,[parameters.PrjDir,'\ntw-s_trn_raw.png'])
close(f2)