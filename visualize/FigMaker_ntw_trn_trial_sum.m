function [] = FigMaker_ntw_trn_trial_sum(tbl,tbl_trn_daily,options,parameters)

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

[groups,~,groupIdx] = unique(tbl.Group,'stable');
[daynums,~,dayIdx] = unique(tbl.Day);
[trials,~,trialIdx] = unique(tbl.Trial);
[distIdx,~,~] = unique([dayIdx,trialIdx],'rows','stable');
T = length(daynums)*length(trials);

NumSbj = arrayfun(@(x)length(find(groupIdx==x)),unique(groupIdx),'uniform',0);
NumSbj = cell2mat(NumSbj)./(length(daynums)*length(trials));


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


cmap = mat2cell(lines(length(groups)),ones(length(groups),1),3);
ii1 = 0;
ii2 = 0;
trialnums = 1:(length(trials)*length(daynums));
for i = 1:length(ntws)
    
    
    if ~isempty(strfind(ntws{i},'1'))
        figure(f2)
        ii2 = ii2 + 1; ii = ii2;
    else
        figure(f1)
        ii1 = ii1 + 1; ii = ii1;
    end
    
    MED = zeros(length(groups),T);
    Q = cell(length(groups),1);
    I = Q;
    
    X = tbl{:,ntws{i}};
    if logdisp % ëŒêîï\é¶
        X = log(X); X(isinf(X)) = min(X(~isinf(X)));
    end
    
    for j = 1:length(groups)
        a = X(strcmp(tbl.Group,groups{j}));
        a = reshape(a,NumSbj(j),T);        
        
        MED(j,:) = median(a,1);
        q = quantile(a,[.25,.75],1);
        
        [~,ind] = find(isnan(q));
        ind = unique(ind);
        q(:,ind) = [];
        
        Q(j) = {[q(2,:),fliplr(q(1,:))]};
        I(j) = {ind'};
    end
    
    
    subplot(3,3,ii);
    h0 = cellfun(@(x,y) patch([trialnums(~ismember(trialnums,y)),...
        fliplr(trialnums(~ismember(trialnums,y)))],x,'r',...
        'FaceAlpha',.3,'EdgeColor','none'),Q,I);
    set(h0,{'FaceColor'},cmap)
    hold on
    h1 = plot(1:T,MED,'-','LineWidth',3);
    set(h1,{'Color'},cmap)
        

    if ismember(ynames(i), toolong)
        str = ynames{i};
        k = strfind(str, ' '); k = k(1);
        str = {str(1:k(1)),str(k(1)+1:end)};
        ylabel([str(1), [str{2},' (',yunits{i},')']], 'Interpreter','latex')
    else
        ylabel([ynames{i}, ' (', yunits{i}, ')'], 'Interpreter','latex')
    end
  
    if ismember(ntws{i}, {'n_1','o_2'})
        legend(h1, groups, 'FontSize',18, 'Interpreter','none'); legend('boxoff')
        xlabel('Day', 'Interpreter','latex');
    else
        xlabel('');
    end    
     
    if min(MED)==max(MED)
        ylim = [-1,1];
    else
        ylim = [min([Q{:}]).*1.05, max([Q{:}]).*1.05];
    end
    
    xtick = 1:step*length(trials):T;
    xticklabels = num2cell(distIdx(xtick)); xticklabels = cellfun(@num2str,xticklabels, 'uniform',0);
    set(gca, 'FontSize', 18, 'LineWidth', 2,...
        'XLim',[.5, T+.5], 'XTick',xtick, 'XTickLabel',xticklabels, 'YLim',ylim);
    
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
imwrite(F,[parameters.PrjDir,'\ntw-d_trn_summ_trial.png'])
close(f1)

F = getframe(f2); F = F.cdata;
imwrite(F,[parameters.PrjDir,'\ntw-s_trn_summ_trial.png'])
close(f2)