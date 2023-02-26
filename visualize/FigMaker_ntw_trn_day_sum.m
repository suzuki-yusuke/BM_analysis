function [] = FigMaker_ntw_trn_day_sum(tbl_trn_daily,options,parameters)


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
f2 = figure('Position', [50,50,900,1480]);
f1 = figure('Position', [50,50,900,1480]);
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
NumSbj = arrayfun(@(x)length(find(groupIdx==x)), unique(groupIdx), 'Uniform', false);
NumSbj = cell2mat(NumSbj)./NumDays;

ii2 = 1; ii1 = 0;

cmap = mat2cell(lines(length(groups)),ones(length(groups),1),3);

for i = 1:length(ntws)
    
    
    if ~isempty(strfind(ntws{i},'1'))
        figure(f2)
        ii2 = ii2 + 1; ii = ii2;
    else
        figure(f1)
        ii1 = ii1 + 1; ii = ii1;
    end
    
    rawdata = cell(length(groups),1);
    for j = 1:length(groups)
        a = tbl_trn_daily{strcmp(tbl_trn_daily.Group,groups(j)), ntws(i)};
        b = [groupIdx(strcmp(tbl_trn_daily.Group, groups(j))),...
            dayIdx(strcmp(tbl_trn_daily.Group, groups(j)))+min(daynums)-1];
%         tf = ~isnan(a)&all(b>0,2);
        tf = all(b>0,2);
        rawdata(j) = {a(tf)};
    end
    
    
    % ëŒêîï\é¶
    if logdisp
        rawdata = cellfun(@log,rawdata,'uniform',0);
        tf = cellfun(@(x) all(abs(x(:))==Inf),rawdata);
        if all(tf)
            continue
        else
            a = min(cellfun(@(x) min(x(~isinf(x))),rawdata));
            for j = 1:length(groups)
                b = rawdata{j};
                b(isinf(b)) = a;
                rawdata(j) = {b};
            end
        end
    end
    
    
    rawdata = cellfun(@(x,y) reshape(x,y,length(daynums)),rawdata,num2cell(NumSbj),'uniform',0);
    MED = cellfun(@(x) median(x,1,'omitnan'),rawdata,'uniform',0);
    
    Q = cellfun(@(x) quantile(x,[.25,.75],1),rawdata,'uniform',0);
    [~,ind] = cellfun(@(x) find(isnan(x)),Q,'uniform',0);
    ind = cellfun(@(x) unique(x),ind,'uniform',0);
    Q = cellfun(@(x) reshape(x(~isnan(x)),2,sum(~isnan(x(:)))/2),Q,'uniform',0);
    Q = cellfun(@(x) [x(2,:),fliplr(x(1,:))], Q, 'uniform',0);
%     MAD = cellfun(@(x) mad(x,1,1),rawdata,'uniform',0);
%     MAD = cellfun(@(x,y) [x+y,fliplr(x-y)], MED,MAD, 'uniform',0);
    MED = cat(1,MED{:});
    
    
    
    subplot(4,2,ii);
    h0 = cellfun(@(x,y) patch([daynums(~ismember(daynums,y));...
        flipud(daynums(~ismember(daynums,y)))],x,'r',...
        'FaceAlpha',.3,'EdgeColor','none'),Q,ind);
    set(h0,{'FaceColor'},cmap)
    
    hold on
    
    h1 = plot(daynums,MED,'-','LineWidth',4);
    set(h1,{'Color'},cmap)
    
    xlabel('Day', 'Interpreter','none','FontWeight','bold');
    

    
    if ismember(ynames(i), toolong)
        str = ynames{i};
        k = strfind(str, ' '); k = k(1);
        str = {str(1:k(1)),str(k(1)+1:end)};
        ylabel(str, 'Interpreter','none','FontWeight','bold')
    else
        ylabel(ynames{i}, 'Interpreter','none','FontWeight','bold')
    end
    
%     title(yunits{i},'Interpreter','latex');
        
    if ismember(ntws{i}, {'n_1','o_2'})
        legend(h1, groups, 'FontSize',18, 'Interpreter','none');
        legend('boxoff')
    end    
    
    if min(MED)==max(MED)
        ylim = [-1,1];
    else
        ylim = [min([Q{:}]).*1.05, max([Q{:}]).*1.05];
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
set(f2, 'color','w')
set(f1, 'color','w')
hold off
figure(f2)
figure(f1)

F = getframe(f1); F = F.cdata;
imwrite(F,[parameters.PrjDir,'\ntw-d_trn_summ.png'])
close(f1)

F = getframe(f2); F = F.cdata;
imwrite(F,[parameters.PrjDir,'\ntw-s_trn_summ.png'])
close(f2)