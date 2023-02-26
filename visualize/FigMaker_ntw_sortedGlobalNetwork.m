function [] = FigMaker_ntw_sortedGlobalNetwork(tbl,parameters,options)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
labelflag = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


T = tbl;
step = options.step;
step = 1;

%% behavioral network
Days = unique(T.Day);
daynames = num2cell(Days);
daynames = cellfun(@num2str,daynames,'uniform',0);
groups = unique(T.Group);
phase = unique(T.Phase);

scrsz = get(groot,'ScreenSize');
figure('Position', scrsz.*[1,1,.9,.7]);

cmax = 30;
cmap = hot(cmax);
for i = 1:step:length(Days)
    for j = 1:length(groups)
        
        
        % read data
        Day = Days(i);
        group = groups(j);
        
        if ismember(phase, {'training','retraining','reversal_training'})
            reading = all([T{:,'Day'}==Day,...
                strcmp(T{:,'Group'},group),...
                ~cellfun(@isempty, T{:,'xy_2'})], 2);
        elseif ismember(phase, {'probe_test'})
            reading = all([strcmp(T{:,'Group'},group),...
                ~cellfun(@isempty, T{:,'xy_2'})], 2);
        end
        
        lnode_c = T{reading, 'xy_2'};
        NumNodes = cellfun(@size, vertcat(lnode_c(:)), 'UniformOutput',0);
        NumNodes = vertcat(NumNodes{:}); NumNodes = NumNodes(:,1); % local network のノード数
        
        lnode_xy = vertcat(lnode_c{:});
        PathBwTrials = cumsum(NumNodes); PathBwTrials(end) = [];
        
        a = cellfun(@iscell, lnode_c);
        b = [lnode_c{a}];
        lnode_c(a) = b;
        
        idnums = T{reading, 'SN'};
        
        idlabels = zeros(sum(NumNodes),1);
        for k = 1:length(NumNodes)
            idlabels(find(idlabels==0,1):find(idlabels==0,1)+NumNodes(k)-1) =...
                repmat(idnums(k),NumNodes(k),1);
        end
        
        
        % CCA
        [~,lnode2gnode,~] = CCA(lnode_xy, parameters.thresh4CCA);
        lnode2gnode = lnode2gnode';
        [gnode_idx,~,ind2gnode] = unique(lnode2gnode);
        gnode_idx = gnode_idx';
        
        N = histcounts(ind2gnode, .5:max(ind2gnode)+.5);
        N = N(ind2gnode); % number of nodes contained in each global node
        
        pwr = log(N).*10+4; % 対数変換
        nodepower = floor(N./sum(N).*100);
        nodepower(nodepower<=1) = 1;
        
        
        % 各使用経路の使用数
        % find crossing links between 2 successive trials
        Path = [lnode2gnode(1:end-1),lnode2gnode(2:end)];
        Path(PathBwTrials,:) = [];
        Path(Path(:,1)==Path(:,2),:) = [];
        
        [Path,~,ic] = unique(Path,'rows');
        pathusage = histcounts(ic,.5:max(ic)+.5);
        
        A = zeros(max(Path(:))); % adjacent matrix        
        linearInd = sub2ind(size(A), Path(:,1), Path(:,2));                       
        A(linearInd) = pathusage;
        A = A + A';
                       
        
        % rank of global degree
        B = sum(A>0,2);
        [C,I] = sort(B,'descend');
        C(C>cmax) = cmax;
        C = C(C>0);
        I = I(C>0);
        nnzA = A(I,I);
        [r,c] = find(nnzA);       
                    

        % make global nodes
        % http://goo.gl/n3DOvD
        NumPoints = length(C)+1; % Number of points making up the circle
        theta = linspace(0,2*pi,NumPoints); % 1000 evenly spaced points between 0 and 2pi
        radius_topo = 1; % Radius of the circle
        rho = ones(1,NumPoints)*radius_topo; % Radius should be 1 for all 100 points
        [x,y] = pol2cart(theta, rho);
        x = x(1:end-1); y = y(1:end-1);
        
        
        p = ceil(length(Days)/step)*(j-1)+(i-1)/step+1;
        subplot(length(groups), ceil(length(Days)/step), p);
        
        axis off
        axis square                
        
        plot([x(c);x(r)], [y(c);y(r)], 'k-')
        hold on
        scatter(x, y, 30, cmap(C,:), 'fill', 'MarkerEdgeColor','k')
        
        hold off
 
    end
end






%%

for i = 1:(length(groups)*ceil(length(Days)/step))
    subplot(length(groups), ceil(length(Days)/step), i);
    a = mod(i,ceil(length(Days)/step));
    b = step*(a+(a==0)*ceil(length(Days)/step))-step+1;
    title({groups{ceil(i/length(Days))},['Day ',daynames{b}]},...
        'FontSize',18, 'interpreter','none', 'FontWeight','bold')
end


xtick = linspace(-radius_topo*1.05,radius_topo*1.05,4);
xticklabel = '';

ytick = xtick;
yticklabel = xticklabel;


for i = 1:length(groups)*ceil(length(Days)/step)
    
    subplot(length(groups), ceil(length(Days)/step), i);
    axis on
    axis square
    
    
    set(gca, 'XLim',xtick([1,end]), 'XTick',xtick, 'XTickLabel',xticklabel,...
        'YLim',ytick([1,end]), 'YTick',ytick, 'YTickLabel',yticklabel,...
        'box','on', 'FontSize',18, 'LineWidth',2)
end

set(gcf,'color','w');
