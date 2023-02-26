function [] = FigMaker_cnv_trn_day_sum(tbl_trn_daily,options,parameters)

%%

switch options.maze_dia
    case 1
        step = 1;
        dcm = 200;
        font_size = 18;
    case 3
        step = 2;
        dcm = 500;
        font_size = 18;
end

% set_model = 'Geometric'; % 'Exponential', 'Generalized Pareto', 'Geometric'
set_model = "Exponential";
% set_model = 'Pareto'; % 'Exponential', 'Generalized Pareto', 'Geometric'
% set_model = 'Generalized Pareto'; % 'Exponential', 'Generalized Pareto', 'Geometric'

rownames_daily = tbl_trn_daily.Properties.VariableNames;
cnvs = rownames_daily(end-2:end);
[groups,~,groupIdx] = unique(tbl_trn_daily.Group,'stable');
[daynums,~,dayIdx] = unique(tbl_trn_daily.Day);
NumDays = length(daynums);
NumSbj = arrayfun(@(x)length(find(groupIdx==x)), unique(groupIdx), 'Uniform', false);
NumSbj = cell2mat(NumSbj)./NumDays;

% draw figure for training
ylabelnames = {'No. of errors (times)', 'Latency (s)', 'Travel distance (cm)'};
f0 = figure('Position', [10,10,2200,640]);
cmap = mat2cell(lines(length(groups)),ones(length(groups),1),3);

for i = 1:length(cnvs)    
    
    Q = cell(1,length(groups));
    I = Q;
    
    for j = 1:length(groups)
        a = tbl_trn_daily{strcmp(tbl_trn_daily.Group,groups(j)), cnvs(i)};
        b = [groupIdx(strcmp(tbl_trn_daily.Group, groups(j))),...
            dayIdx(strcmp(tbl_trn_daily.Group, groups(j)))+min(daynums)-1];
%         tf = ~isnan(a)&all(b>0,2);
        tf = all(b>0,2);
        X = reshape(a(tf),NumSbj(j),length(daynums));
        
        [~,ind] = find(isnan(X));
        ind = unique(ind);
        X(:,ind) = [];                
        
%         Q(j) = {quantile(X,[.25,.5,.75],1)};
        Q(j) = {[median(X,1)-mad(X,1);median(X,1);median(X,1)+mad(X,1)]};
%         SEM = std(X,[],1)./sqrt(size(X,1)); % Standard Error
%         ts = tinv([0.025  0.975],size(X,1)-1); % T-Score
%         Q(j) = {[mean(X,1)+ts(1).*SEM; mean(X,1); mean(X,1)+ts(2).*SEM]}; % 95%CI
        I(j) = {ind};
        
    end
        
    %%
    figure(f0)
    subplot(1, length(cnvs), i);
    h0 = cellfun(@(x,y) patch([daynums(~ismember(daynums,y));...
        flipud(daynums(~ismember(daynums,y)))],[x(3,:),fliplr(x(1,:))],'r',...
        'FaceAlpha',.3,'EdgeColor','none'),Q,I);
    set(h0,{'FaceColor'},cmap)
    
    hold on
    
    h1 = cellfun(@(x,y) plot(daynums(~ismember(daynums,y)),x(2,:),'-','LineWidth',6),Q,I);
    set(h1,{'Color'},cmap);
%     set(h1,{'MarkerFaceColor'},cmap);
    
    ylabel(ylabelnames{i}, 'Interpreter','none', 'FontWeight','bold', 'FontSize',font_size)
%     xtick = daynums(1:step:end);
    xtick = daynums;
    xticklabels = num2cell(xtick);
    xticklabels = cellfun(@num2str,xticklabels, 'uniform',0);
    xticklabels(2:2:end) = repmat({' '}, 1,floor(length(daynums)/2));
    
    xlabel('Day', 'Interpreter','none', 'FontWeight','bold', 'FontSize',font_size);
    
    set(gca, 'FontSize',font_size, 'LineWidth',2,...
        'XLim',[min(daynums)-.5, max(daynums)+.5],...
        'XTick',xtick, 'XTicklabel',xticklabels);
        
    legend(h1, groups, 'FontSize',18, 'Interpreter','none'); legend('boxoff')
    
    switch i
        case 1
            
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
            set(gca,'YTickLabel',yticklabels,'YTick',yticks,'YLim',yticks([1,end]));
    end
    
    axis square
    box off
    
    hold off
        
                        
    
end
set(gcf,'color','w');
F = getframe(f0); F = F.cdata;
imwrite(F,[parameters.PrjDir,'\cnv_trn_summ.png'])
close(f0)

end