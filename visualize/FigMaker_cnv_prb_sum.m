function [] = FigMaker_cnv_prb_sum(tbl_prb_daily,parameters,options)

%%
NumHoles = length(unique(tbl_prb_daily.Holes));
[groups,~,ig] = unique(tbl_prb_daily.Group,'stable');


f = figure('Position', [10, 10, 900, 900]);


cmap1 = colormap(lines(length(groups)));

goal = 1;
RotIdx = circshift(1:NumHoles,5,2);


MED = cell(1,length(groups));
Q = MED;
xticknames = cellfun(@num2str, num2cell(-150:30:180), 'uniform',0);
xticknames(strcmp(xticknames,'0')) = {'Target'};
xticknames(~ismember(1:length(xticknames), 2:2:length(xticknames))) = {''};
catcolors = cell(1, length(groups));


for i = 1:length(groups)
    
    a = tbl_prb_daily{ig==i, end};
    a = reshape(a, NumHoles, length(a)/NumHoles);
    a = mat2cell(a,NumHoles,ones(size(a,2),1));
    a = cellfun(@(x) x(RotIdx), a, 'uniform',0); a = horzcat(a{:});
    MED(i) = {median(a,2)};    
    
    b = [median(a,2)+mad(a,1,2);...
        flipud(median(a,2))-flipud(mad(a,1,2))]; b = [b;b(1)];
    b = [b,[(1:NumHoles)';flipud((1:NumHoles)');1]];  
    Q(i) = {b};
    
    catcolors(i) = {cmap1(i,:)};
    
end


ymax = 60;

hold on
cellfun(@(x,y) patch(x(:,2),x(:,1),y,'FaceAlpha',.3,'EdgeColor','none'),...
    Q,catcolors)
h = cellfun(@(x,y) plot(1:NumHoles,x,'LineWidth',6), MED,catcolors);

% plot target pos
plot(ones(2,1).*find(RotIdx==goal),[0,100],'k--','LineWidth',2)

legend(h, groups, 'FontSize',18, 'Interpreter','none'); legend('boxoff')

xlabel('Angle from target', 'interpreter','none', 'FontSize',16, 'FontWeight','bold')
ylabel({'Time spent around each hole (s)'}, 'interpreter','none', 'FontSize',16, 'FontWeight','bold')

set(gca, 'XLim',[.5,NumHoles+.5], 'XTick',1:NumHoles, 'XTickLabel',xticknames,...
    'YLim',[0,ymax],'FontSize',18 ,'LineWidth',2)
axis square
box off

hold off

set(gcf,'color','w');

F = getframe(f); F = F.cdata;
imwrite(F,[parameters.PrjDir,'\cnv_prb_summ.png'])
close(f)