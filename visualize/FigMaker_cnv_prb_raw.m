function [] = FigMaker_cnv_prb_raw(tbl_prb_daily,parameters,options)

%%
NumHoles = length(unique(tbl_prb_daily.Holes));
[groups,~,ig] = unique(tbl_prb_daily.Group);


f = figure('Position', [10, 10, 900, 900]);


cmap1 = colormap(lines(length(groups)));

c = arrayfun(@(x)length(find(ig==x)), unique(ig), 'Uniform',false);
c = cell2mat(c);

rawdata = zeros(max(c), length(groups));
CatIdx = zeros(max(c), length(groups));

goal = 1;
RotIdx = circshift(1:NumHoles,5,2);


E = zeros(NumHoles,length(groups));

xticknames = cellfun(@num2str, num2cell(-150:30:180), 'uniform',0);
xticknames(strcmp(xticknames,'0')) = {'Target'};
xticknames(~ismember(1:length(xticknames), 2:2:length(xticknames))) = {''};
catcolors = cell(1, length(groups));


for i = 1:length(groups)
    
    a = tbl_prb_daily{ig==i, end};
    a = reshape(a, NumHoles, length(a)/NumHoles);
    a = mat2cell(a,NumHoles,ones(size(a,2),1));
    a = cellfun(@(x) x(RotIdx), a, 'uniform',0); a = horzcat(a{:});
    E(:,i) = median(a,2);
    
    rawdata(1:sum(ig==i),i) = a(:);
    CatIdx(1:sum(ig==i),i) = i;
    DistIdx(1:sum(ig==i),i) = tbl_prb_daily{ig==i, 'Holes'};
    catcolors(i) = {cmap1(i,:)};
    
end


CatIdx = CatIdx(:); DistIdx = DistIdx(:);

if any(CatIdx==0)
    catcolors = [{[0,0,0]}, catcolors];
    catmarkers = ['none', repmat({'.'},1,length(groups))];
else
    catmarkers = repmat({'.'},1,length(groups));
end


[~,pos] = plotSpread(rawdata, 'distributionIdx',DistIdx,...
    'categoryIdx',CatIdx, 'categoryColors',catcolors,...
    'categoryMarkers',catmarkers,...
    'showMM',0);

% plot raw value
cmap1 = colormap(lines(length(groups)));
cmap2 = cmap1+(1-cmap1).*0.9;
for j = 1:length(groups)
    posj = pos(CatIdx==j,:);
    xdot = reshape(posj(:,1), NumHoles, sum(CatIdx==j)/NumHoles);
    ydot = reshape(posj(:,2), NumHoles, sum(CatIdx==j)/NumHoles);
    
    plot(xdot,ydot, 'Color',cmap2(j,:));
    hold on
end

% plot representative value
h = line(1:NumHoles,E', 'LineStyle','-', 'LineWidth',2);
markeredgecolors = mat2cell(cmap1.*.7,ones(length(groups),1),3);
set(h, {'Color'},markeredgecolors);

% plot raw value
plotSpread(rawdata, 'distributionIdx',DistIdx,...
    'categoryIdx',CatIdx, 'categoryColors',catcolors,...
    'categoryMarkers',catmarkers,...
    'showMM',0);

% plot target pos
plot(ones(2,1).*find(RotIdx==goal),[0,100],'k--')


legend(h, groups, 'FontSize',18, 'Interpreter','none'); legend('boxoff')

xlabel('Angle from target', 'interpreter','latex', 'FontSize',16)
ylabel({'Time spent around each hole', '(s)'}, 'interpreter','latex', 'FontSize',16)

set(gca, 'XLim',[.5,NumHoles+.5], 'XTick',1:NumHoles, 'XTickLabel',xticknames,...
    'YLim',[0,50],'FontSize',18)
axis square
box off

hold off

F = getframe(f); F = F.cdata;
imwrite(F,[parameters.PrjDir,'\cnv_prb_raw.png'])
close(f)