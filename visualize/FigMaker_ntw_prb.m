function [] = FigMaker_ntw_prb(tbl_prb_daily,parameters)


%%
f2 = figure('Position', [50,50,900,1480]);
f1 = figure('Position', [50,50,900,1480]);
ii2 = 1; ii1 = 0;

rownames_daily = tbl_prb_daily.Properties.VariableNames;
ntws = rownames_daily(3:end);
ntws = ntws([8,9,11,13,12,10,14,15,1,3,5,4,2,6,7]);
ynames = {'Order','Shortest path length','Degree',...
    'Clustering coefficient','Density', 'Betweenness centrality', 'Closeness centrality',...
    'No. of stops','Order','Shortest path length','Degree',...
    'Clustering coefficient','Density','Betweenness centrality', 'Closeness centrality'};
toolong = ynames([2,4,6,7,10,12,14,15]);
ynames = ynames([8,9,11,13,12,10,14,15,1,3,5,4,2,6,7]);

%     [yunits, ~] = strtok(ntws,'_');
yunits = {'$\overline{o}$','$\overline{n}$','$\overline{m}$','$\overline{\rho}$',...
    '$\overline{C_{WS}}$','$\overline{l}$','$\overline{x}$','$\overline{cl}$',...
    '$\overline{n}$','$\overline{m}$','$\overline{\rho}$',...
    '$\overline{C_{WS}}$','$\overline{l}$','$\overline{x}$','$\overline{cl}$'};

groups = unique(tbl_prb_daily.Group,'stable');
cmap = colormap(lines(length(groups)));
legnames = cell(length(groups)*2,1);
for i = 1:2:length(groups)*2
    legnames(i:i+1) = {['Individuals of ',groups{(i+1)/2}],...
        ['median of ',groups{(i+1)/2}]};
end


for i = 1:length(ntws)
        
    if ~isempty(strfind(ntws{i},'1'))
        figure(f2)
        ii2 = ii2 + 1; ii = ii2;
    else
        figure(f1)
        ii1 = ii1 + 1; ii = ii1;
    end
    
    
    data = tbl_prb_daily{:, ntws{i}};
    if iscell(data)
        data = cell2mat(data);
    end
    
    if all(data==0)
        continue
    end
    
    subplot(4, 2, ii);
    
hs = zeros(length(groups),1);
    for j = 1:length(groups)
        
        Y = data(strcmp(tbl_prb_daily.Group, groups{j}));              
        catIdx = ones(length(Y),1).*j;
        distIdx = catIdx;
                
        subplot(4, 2, ii);
        
        [~,pos] = plotSpread(Y, 'distributionIdx',distIdx,...
            'categoryIdx',catIdx,...
            'categoryColors',{cmap(j,:)},...
            'categoryMarkers', {'none','.'},...
            'showMM',0);
        
        scatter(pos(:,1),pos(:,2), 20, cmap(j,:), 'filled');
        hold on
        
        h = plot(j, median(Y), 's', 'MarkerSize',13, 'LineWidth',6, 'Color', cmap(j, :).*.7);
        hs(j) = h;
                                   
        
        if ismember(ynames(i), toolong)
            str = ynames{i};
            k = strfind(str, ' '); k = k(1);
            str = {str(1:k(1)),str(k(1)+1:end)};
            ylabel(str, 'Interpreter','none')
        else
            ylabel(ynames{i}, 'Interpreter','none')
        end
        
        ylim = quantile(data,[.01,.99]);
        if ylim(1)==ylim(2)
            ylim = [-1,1];
        end
        set(gca, 'FontSize', 18, 'LineWidth', 2,...
            'XLim', [0.5, length(groups)+.5],...
            'XTick', 1:length(groups),...
            'YLim', ylim)        
        box off
        axis square
        
        
    end
    
    set(gca, 'XTickLabel',groups, 'FontWeight','bold')
    xaxisproperties= get(gca, 'XAxis');
    xaxisproperties.TickLabelInterpreter = 'none'; % latex for x-axis
    
end

hold off
% set(hs(1),'Position',p0);
set(f2, 'color','w')
set(f1, 'color','w')

F = getframe(f1); F = F.cdata;
imwrite(F,[parameters.PrjDir,'\ntw-d_prb.png'])
close(f1)

F = getframe(f2); F = F.cdata;
imwrite(F,[parameters.PrjDir,'\ntw-s_prb.png'])
close(f2)