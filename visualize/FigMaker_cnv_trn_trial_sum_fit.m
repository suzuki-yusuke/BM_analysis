function [] = FigMaker_cnv_trn_trial_sum_fit(tbl,tbl_trn_daily,options,parameters)

%%

switch options.maze_dia
    case 1
        step = 1;
        dcm = 200;
    case 3
        step = 2;
        dcm = 500;
end

% 'Exponential', 'Generalized Pareto', 'Geometric'
% set_model = "Linear";
% set_model = 'Geometric';
set_model = "Exponential";
% set_model = 'Pareto';
% set_model = 'Generalized Pareto';
% set_model = 'HalfNormal';

normalize_frag = 1;

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
f0 = figure('Position', [10,10,2200,640]);
f1 = figure('Position', [10,10,2200,640]);
f2 = figure('Position', [10,10,2200,840]);
f3 = figure('Position', [10,10,840,2200]);
f4 = figure('Position', [10,10,2200,840]);
f5 = figure('Position', [10,10,2200,840]);
cmap = mat2cell(lines(length(groups)),ones(length(groups),1),3);


trialnums = 1:T;

P = cell(length(groups),length(cnvs));
T_t = P;
D = P;
Y = cell(length(groups),length(cnvs),2);
h3 = zeros(length(cnvs),1);
h4 = zeros(length(cnvs),1);
h5 = zeros(length(cnvs),1);
for i = 1:length(cnvs)
    
    Q = cell(length(groups),3);
    
    for j = 1:length(groups)
        
        X = tbl{strcmp(tbl.Group,groups{j}),cnvs{i}};
        X = reshape(X,NumSbj(j),T);
        
        if (i>1)&&normalize_frag % normalize latency and travel distance b/w BM1 & BM3
            X = mat2cell(X,ones(size(X,1),1),size(X,2));
            X = cellfun(@(x) (x-min(x))./(max(x)-min(x)),X,'uniform',0);
            X = cat(1,X{:});
        end
        
        [~,ind] = find(isnan(X));
        ind = unique(ind);
        X(:,ind) = [];
        
        Q(j,1) = {median(X,1)};
        X = mat2cell(X,ones(size(X,1),1),size(X,2));
        
        switch set_model
            case 'Linear'
                X = cat(1,X{:});
                X = cumsum(X,2);
                Q(j,1) = {median(X,1)};
                
                X = mat2cell(X,ones(size(X,1),1),size(X,2));
                
                initial_learning_trials = 1:9;
                
                X0 = X;
                X = cellfun(@(x) x(initial_learning_trials),X,'uniform',0);
                
                E = [ones(size(X{1}));1:length(X{1})]';
                c = cellfun(@(x) E\x', X, 'uniform',0);
                P(j,i) = {cellfun(@(x) x(2),c)};
                y = cellfun(@(x) x(2).*trialnums+x(1),c,'uniform',0);
                D(j,i) = {cellfun(@(x,y) y(length(x))-x(end),X0,y)};
                Y(j,i,1) = {cat(1,y{:})};
                T_t(j,i) = {cellfun(@(x) find(cumsum(x./sum(x))>.8,1),y)};
                
                c = E\Q{j,1}(initial_learning_trials)';
                y = c(2).*trialnums+c(1);
                Y(j,i,2) = {y};
                
            case 'Geometric'
                modelFun = @(a,x) a(1)*(1-a(1)).^x+a(2);
                nlModel = cellfun(@(x) fitnlm(1:length(x),x./sum(x),modelFun,[.1,0]),X,'uniform',0);
                p = cellfun(@(x) x.Coefficients.Estimate, nlModel,'uniform',0);
                p = [p{:}]';
                P(j,i) = {p(:,1)};
                %                 modelFun = @(a,x) a(1)*(1-a(1)).^x;
                %                 nlModel = cellfun(@(x) fitnlm(1:length(x),x./sum(x),modelFun,.1),X,'uniform',0);
                %                 p = cellfun(@(x) x.Coefficients.Estimate, nlModel);
                %                 P(j,i) = {p};
                y = cellfun(@(x) predict(x,trialnums'),nlModel,'uniform',0);
                y = cellfun(@(x,y) x.*sum(y),y,X,'uniform',0);
                offset = 0;
                y = cat(2,y{:})' + offset;
                D(j,i) = {arrayfun(@(x) y(x,end)-y(x,1),1:NumSbj(j))'};
                %                 T_converge(j,i) = {arrayfun(@(x) icdf(set_model,.8,x),p)};
                T_t(j,i) = {cellfun(@(x) find(cumsum(x./sum(x))>.5,1),y)};
                
                %                 nlModel = fitnlm(1:length(z),Q{j,1}./sum(Q{j,1}),modelFun,[.1,0]);
                %                 yy = predict(nlModel,trialnums')'.*sum(a);
                %                 y = yy;
                
            case 'Exponential'
                p = cellfun(@(x) fit((1:length(x))',x','exp1'),X,'uniform',0);
                P(j,i) = {cellfun(@(x) -x.b, p)};
                y = cellfun(@(x) x.a.*exp(x.b.*(1:1000)),p,'uniform',0);
                T_t(j,i) = {cellfun(@(x) log(2)/(-x.b),p)}; % half-life
%                 T_t(j,i) = {cellfun(@(x) find(cumsum(x./sum(x))>.5,1),y)};
                y = cat(1,y{:}); y(:,T+1:end) = [];
                Y(j,i,1) = {y};
                D(j,i) = {arrayfun(@(x) y(x,end)-y(x,1),1:NumSbj(j))'};
                %                 D(j,i) = {arrayfun(@(x) mean(y(x,end-2:end))-mean(y(x,1:3)),1:NumSbj(j))'};
                
                p = fit((1:length(X{1}))',Q{j,1}','exp1');
%                 p = fit((1:18)',Q{j,1}(1:18)','exp1');
                y = p.a.*exp(p.b.*trialnums);
                p.a;
                p.b;
                Y(j,i,2) = {y};
                
            case 'Pareto'
                %                 modelFun = @(a,x) (a(1)/a(2))./((x./a(2)).^(a(1)+1));
                modelFun = @(a,x) (a(1)*a(2)^a(1))./(x.^(a(1)+1));
                nlModel = cellfun(@(x) fitnlm(1:length(x),x./sum(x),modelFun,[.01,.01]),X,'uniform',0);
                p = cellfun(@(x) x.Coefficients.Estimate, nlModel,'uniform',0);
                p = [p{:}]';
                P(j,i) = {p(:,1)};
                y = cellfun(@(x) predict(x,(1:T)'),nlModel,'uniform',0);
                T_t(j,i) = {cellfun(@(x) find(cumsum(x./sum(x))>.5,1),y)};
                
                offset = cellfun(@(x) quantile(x,.25),X,'uniform',0);
                y = cellfun(@(x,y,z) x.*sum(y)+z*0,y,X,offset,'uniform',0);
                
                y = cat(2,y{:})'; y(:,T+1:end) = [];
                Y(j,i,1) = {y};
                D(j,i) = {arrayfun(@(x) y(x,end)-y(x,1),1:NumSbj(j))'};
                
                nlModel = fitnlm(1:length(Q{j,1}),Q{j,1}./sum(Q{j,1}),modelFun,[.01,.01]);
                offset = quantile(Q{j,1},.25);
                y = predict(nlModel,trialnums')';
                y = y.*sum(Q{j,1}) + offset*0;
                Y(j,i,2) = {y};
            case 'Generalized Pareto'
                theta = 1;
                p = cellfun(@(x) [gpfit(x+1),theta],X,'uniform',0);
                p = cat(1,p{:});
                P(j,i) = {p(:,2)};
                y = cellfun(@(x) gppdf(1:T,x(1),x(2),x(3)),p,'uniform',0);
                T_t(j,i) = {cellfun(@(x) find(cumsum(x./sum(x))>.5,1),y)};
                offset = cellfun(@(x) quantile(x,.25),X,'uniform',0);
                y = cellfun(@(x,y,z) smooth(x.*sum(x+1))'+z.*0,y,X,offset,'uniform',0);
                y = cat(1,y{:}); y(:,T+1:end) = [];
                Y(j,i,1) = {y};
                
                
                p = [gpfit(Q{j,1}+1),theta];
                y = gppdf(1:T,p(1),p(2),p(3));
                offset = quantile(Q{j,1},.25);
                y = smooth(y.*sum(Q{j,1}+1))' + offset*0;
                Y(j,i,2) = {y};
            case 'HalfNormal'
                p = cellfun(@(x) fitdist(x','HalfNormal'),X,'uniform',0);
                P(j,i) = {cellfun(@(x) x.sigma,p)};
                y = cellfun(@(x) pdf(x,1:T),p,'uniform',0);
                y = cellfun(@(a,b) a.*sum(b),y,X,'uniform',0);
                y = cat(1,y{:}); y(:,T+1:end) = [];
                Y(j,i,1) = {y};
                D(j,i) = {arrayfun(@(x) y(x,end)-y(x,1),1:NumSbj(j))'};
                T_t(j,i) = {cellfun(@(x) icdf('HalfNormal',.5,x.mu,x.sigma),p)};
                
                p = fitdist(Q{j,1}','HalfNormal');
                y = pdf(p,1:T);
                y = y.*sum(Q{j,1});
                Y(j,i,2) = {y};
        end
        
        q = quantile(Y{j,i,1},[.25,.5,.75],1);
        Q(j,2) = {q(2,:)};
        Q(j,3) = {[q(3,:),fliplr(q(1,:))]};
        
        % Q(j,2) = {mean(y,1)};
        % Q(j,3) = {[mean(y,1)+std(y,1),fliplr(mean(y,1)-std(y,1))]};
        
    end
    
    
    %%
    figure(f0)
    subplot(1, length(cnvs), i);
    
    hold on
    cellfun(@(x,y) plot(trialnums,x,'-','Color',y),Y(:,i,1),cmap)
    ylabel(ylabelnames{i}, 'Interpreter','none')
    xtick = 1:step*length(trials):T;
    xticklabels = num2cell(distIdx(xtick));
    xticklabels = cellfun(@num2str,xticklabels, 'uniform',0);
    
    set(gca, 'FontSize', 18, 'LineWidth', 2,...
        'XLim',[.5, T+.5],...
        'XTick',xtick, 'XTicklabel',xticklabels);
    
    legend(gca, groups, 'FontSize',18, 'Interpreter','none'); legend('boxoff')
    xlabel('Day', 'Interpreter','none');
    
    switch i
        case 1
                        
        case 2
%             yticklabels = get(gca,'YTickLabels');
%             yticks = 0:60:(60*floor(max(cellfun(@(x) str2double(x),yticklabels))/60)+60);
%             yticklabels = arrayfun(@(x) num2str(x),yticks,'uniform',0);
            %             set(gca,'YTickLabel',yticklabels,'YTick',yticks);
        case 3
%             yticklabels = get(gca,'YTickLabels');
%             yticks = cellfun(@(x) str2double(x), yticklabels);
%             yticks = 0:dcm:(100*floor(yticks(end)/100)+100);
%             yticklabels = arrayfun(@(x) num2str(x),yticks,'uniform',0);
            %             set(gca,'YTickLabel',yticklabels,'YTick',yticks,'YLim',yticks([1,end]));
    end
    
    axis square
    grid on
    box off
    hold off
    
    
    
    %%
    figure(f1)
    subplot(1, length(cnvs), i);
    h = cellfun(@(x,y) patch([trialnums,...
        fliplr(trialnums)],x,'r',...
        'FaceAlpha',.3,'EdgeColor','none'),Q(:,3));
    set(h,{'FaceColor'},cmap)
    hold on
    h1 = plot(trialnums,cat(1,Q{:,2}),'-','LineWidth',3);
    set(h1,{'Color'},cmap)
    
    cellfun(@(x,y) scatter(1:length(x),x,80,y,'filled'),Q(:,1),cmap);
    
    ylabel(ylabelnames{i}, 'Interpreter','none', 'FontWeight','bold')
    xtick = 1:step*length(trials):T;
    xticklabels = num2cell(distIdx(xtick));
    xticklabels = cellfun(@num2str,xticklabels, 'uniform',0);
    
    set(gca, 'FontSize', 18, 'LineWidth', 2,...
        'XLim',[.5, T+.5],...
        'XTick',xtick, 'XTicklabel',xticklabels);
    
    switch i
        case 1
            legend(h1, groups, 'FontSize',18, 'Interpreter','none'); legend('boxoff')
            xlabel('Day', 'Interpreter','none', 'FontWeight','bold');
        case 2
            yticklabels = get(gca,'YTickLabels');
            yticks = 0:60:(60*floor(max(cellfun(@(x) str2double(x),yticklabels))/60)+60);
            yticklabels = arrayfun(@(x) num2str(x),yticks,'uniform',0);
            %             set(gca,'YTickLabel',yticklabels,'YTick',yticks);
        case 3
            yticklabels = get(gca,'YTickLabels');
            yticks = cellfun(@(x) str2double(x), yticklabels);
            yticks = 0:dcm:(100*floor(yticks(end)/100)+100);
            yticklabels = arrayfun(@(x) num2str(x),yticks,'uniform',0);
            %             set(gca,'YTickLabel',yticklabels,'YTick',yticks,'YLim',yticks([1,end]));
    end
    
    axis square
    grid on
    box off
    hold off
    
    
    %%
    figure(f2)
    subplot(1, length(cnvs), i);
    hold on
       
    
    cellfun(@(x,y) scatter(1:length(x),x,80,y,'filled'),Q(:,1),cmap);    
    ylabel(ylabelnames{i}, 'Interpreter','none', 'FontWeight','bold')    
        
%     set(gca,'YGrid','on','XGrid','on','XMinorGrid','on')    
        
    xlabel('Day', 'Interpreter','none', 'FontWeight','bold');
    
    yl = get(gca,'YLim');
    
    switch i
        case 1            
            ytick = 0:2:max(yl)+1;
        case 2
            if normalize_frag
                yl = [0,1];
                ytick = 0:0.2:1;
                set(gca,'YLim',[0,1])
                ylabel('Normalized latency', 'Interpreter','none')
            else
                yl = [0,120];
                ytick = 0:30:120;
                ylabel('Latency (s)', 'Interpreter','none')
            end            
                        
%             yticklabels = get(gca,'YTickLabels');
%             yticks = 0:60:(60*floor(max(cellfun(@(x) str2double(x),yticklabels))/60)+60);
%             yticklabels = arrayfun(@(x) num2str(x),yticks,'uniform',0);
            %             set(gca,'YTickLabel',yticklabels,'YTick',yticks);
        case 3
            if normalize_frag
                yl = [0,1];
                ytick = 0:0.2:1;
                set(gca,'YLim',[0,1])
                ylabel('Normalized tarvel distance', 'Interpreter','none')
            else
                yl = [0,2500];
                ytick = 0:500:2500;
                ylabel('Tarvel distance (cm)', 'Interpreter','none')
            end
            
%             yticklabels = get(gca,'YTickLabels');
%             yticks = cellfun(@(x) str2double(x), yticklabels);
%             yticks = 0:dcm:(100*floor(yticks(end)/100)+100);
%             yticklabels = arrayfun(@(x) num2str(x),yticks,'uniform',0);
            %             set(gca,'YTickLabel',yticklabels,'YTick',yticks,'YLim',yticks([1,end]));
    end
    
    %     xtick = 1:step*length(trials):T;
    xtick = 1:length(trials):T;
    xticklabels = num2cell(distIdx(xtick));
    xticklabels = cellfun(@num2str,xticklabels, 'uniform',0);
    xticklabels = cellfun(@num2str,xticklabels, 'uniform',0);
    xticklabels(2:2:end) = repmat({' '},1,floor(length(xticklabels)/2));           
    
    [Xg,Yg] = meshgrid(ytick,1:T); Zg = zeros(size(Xg));
    mesh(Yg,Xg,Zg,'FaceColor','none','LineStyle',':',...
        'EdgeColor',[.7,.7,.7],'LineWidth',2);    
    
    cellfun(@(x,y) scatter(1:length(x),x,80,y,'filled'),Q(:,1),cmap);
    h1 = plot(trialnums,cat(1,Y{:,i,2}),'-','LineWidth',3);
    legend(h1, groups, 'FontSize',18, 'Interpreter','none'); legend('boxoff')    
    set(h1,{'Color'},cmap)    
    
    set(gca, 'XLim',[.5, T+.5],'XTick',xtick,'XTicklabel',xticklabels,...
        'YLim',[0,max(yl)],...
        'TickLength',[0.02,0.025])    
    set(gca, 'FontSize', 18, 'LineWidth', 2);
    set(gca ,'Layer', 'Top')
    set(gcf,'color','w');

    axis square        
    box off
    hold off
    
    
    
    
    %%
    figure(f3)
    h3(i) = subplot(length(cnvs),1,i);
    h = boxplot(cat(1,P{:,i}),groups(groupIdx(1:sum(NumSbj))),'Widths',0.7);
    set(h,{'LineWidth'},{4})
%     if i==1
%         p0 = get(gca, 'Position');
%     else
%         p1 = get(gca, 'Position');
%         set(gca,'Position',[p1(1),p1(2),p1(3),p1(4)])
%     end    
    ylabel('É¿', 'FontWeight','bold')
    set(gca,'FontSize',18, 'LineWidth',2)
%     title(ylabelnames{i})
    set(gcf,'color','w');
    axis square
    
    xtickangle(45)
    
    figure(f4)
    h4(i) = subplot(1,length(cnvs),i);
    h = boxplot(cat(1,T_t{:,i}),groups(groupIdx(1:sum(NumSbj))));
    set(h,{'LineWidth'},{4})
    if i==1
        p0 = get(gca, 'Position');
    else
        p1 = get(gca, 'Position');
        set(gca,'Position',[p1(1),p1(2),p1(3),p1(4)])
    end
    set(gca,'FontSize',18, 'LineWidth',2)
    ylabel('Cumulative trial', 'FontWeight','bold')
    title(ylabelnames{i})
    set(gcf,'color','w');
    axis square    
    
    figure(f5)
    h5(i) = subplot(1,length(cnvs),i);
    h = boxplot(cat(1,D{:,i}),groups(groupIdx(1:sum(NumSbj))));    
    set(h,{'LineWidth'},{4})
    if i==1
        p0 = get(gca, 'Position');
    else
        p1 = get(gca, 'Position');
        set(gca,'Position',[p1(1),p1(2),p1(3),p1(4)])
    end
    
    set(gca,'FontSize',18, 'LineWidth',2)
    ylabel('Difference', 'FontWeight','bold')
    title(ylabelnames{i})
    set(gcf,'color','w');
    axis square
    
    
end

F = getframe(f0); F = F.cdata;
imwrite(F,[parameters.PrjDir,'\cnv_trn_summ_trial_fit_raw.png'])
close(f0)

F = getframe(f1); F = F.cdata;
imwrite(F,[parameters.PrjDir,'\cnv_trn_summ_trial_fit.png'])
close(f1)

F = getframe(f2); F = F.cdata;
imwrite(F,[parameters.PrjDir,'\cnv_trn_summ_trial_fit_group.png'])
close(f2)

% set(h3(1),'Position',p0);
F = getframe(f3); F = F.cdata;
imwrite(F,[parameters.PrjDir,'\cnv_trn_summ_trial_fit_beta.png'])
close(f3)

set(h4(1),'Position',p0);
F = getframe(f4); F = F.cdata;
imwrite(F,[parameters.PrjDir,'\cnv_trn_summ_trial_fit_t.png'])
close(f4)

set(h5(1),'Position',p0);
F = getframe(f5); F = F.cdata;
imwrite(F,[parameters.PrjDir,'\cnv_trn_summ_trial_fit_d.png'])
close(f5)


a = table(groups(groupIdx(1:sum(NumSbj))));
a.Properties.VariableNames = {'Group'};
b = array2table(reshape(cat(1,P{:}),sum(NumSbj),length(cnvs)));
b.Properties.VariableNames = cnvs;
c = [a,b];
writetable(c,[parameters.PrjDir,'\cnv_trn_summ_trial_fit_beta.txt'],'delimiter','\t')

b = array2table(reshape(cat(1,T_t{:}),sum(NumSbj),length(cnvs)));
b.Properties.VariableNames = cnvs;
c = [a,b];
writetable(c,[parameters.PrjDir,'\cnv_trn_summ_trial_fit_t.txt'],'delimiter','\t')

b = array2table(reshape(cat(1,D{:}),sum(NumSbj),length(cnvs)));
b.Properties.VariableNames = cnvs;
c = [a,b];
writetable(c,[parameters.PrjDir,'\cnv_trn_summ_trial_fit_d.txt'],'delimiter','\t')

end