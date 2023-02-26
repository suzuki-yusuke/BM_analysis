function [] = FigMaker_ntw_GlobalNetwork(tbl,parameters,options)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
labelflag = 0;
% step = options.step;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


T = tbl;
if strcmp(unique(T.Phase),'probe_test')
    step = 1;
else
    step = options.step;
    step = 1;
end


%% behavioral network
Days = unique(T.Day);
daynames = num2cell(Days);
daynames = cellfun(@num2str,daynames,'uniform',0);
groups = unique(T.Group);
phase = unique(T.Phase);

scrsz = get(groot,'ScreenSize');


polc = [255,151,163;...
    69,255,193;...
    255,255,160;...
    153,201,255;...
    198,228,255]./255;
% polc = [0.9110,1.0000,0.5238
%     1.0000,0.6713,0.5975
%     0.7834,0.5956,1.0000
%     0.5351,0.9999,1.0000
%     0.5486,1.0000,0.9456];

% polc = polcmap(5);
% imagesc(reshape(polc,5,1,3))

cmap_0 = hsv(360);

figure('Position', scrsz.*[1,1,.9,.7]);

for i = 1:step:length(Days)
    for j = 1:length(groups)
        
        % read data
        Day = Days(i);
        group = groups(j);
        
        if ismember(phase, {'training'})
            reading = all([T{:,'Day'}==Day,...
                strcmp(T{:,'Group'},group),...
                ~cellfun(@isempty, T{:,'xy_2'})], 2);
        elseif ismember(phase, {'probe_test'})
            reading = all([strcmp(T{:,'Group'},group),...
                ~cellfun(@isempty, T{:,'xy_2'})], 2);
        end
        
        local_node_xy = T{reading, 'xy_2'};
        NumNodes = cellfun(@size, vertcat(local_node_xy(:)), 'UniformOutput',0);
        NumNodes = vertcat(NumNodes{:}); NumNodes = NumNodes(:,1); % local network のノード数
        
        lnode_xy = vertcat(local_node_xy{:});
        PathBwTrials = cumsum(NumNodes); PathBwTrials(end) = [];
        
        a = cellfun(@iscell, local_node_xy);
        b = [local_node_xy{a}];
        local_node_xy(a) = b;
        
        idnums = T{reading, 'SN'};
        
        idlabels = zeros(sum(NumNodes),1);
        for k = 1:length(NumNodes)
            idlabels(find(idlabels==0,1):find(idlabels==0,1)+NumNodes(k)-1) =...
                repmat(idnums(k),NumNodes(k),1);
        end
        
        % CCA
        [global_node_xy,lnode2gnode] = CCA_visualize(lnode_xy, parameters.thresh4CCA);
%         lnode2gnode = lnode2gnode';
        [global_node_idx,~,ind2gnode] = unique(lnode2gnode);
        global_node_idx = global_node_idx';
        
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
        
        
        pathpower = log(A)+1; % 対数変換
        pathpower(isinf(pathpower)) = 1;

        
        % 各ノードのリンク数
        NumGlobalLinks = zeros(length(global_node_idx),1);
        for k = 1:length(global_node_idx)
            NumGlobalLinks(k) = sum((global_node_idx(k)==global_node_idx(:,1)));
        end
        linkpower = floor(NumGlobalLinks./sum(NumGlobalLinks).*100);
        linkpower(linkpower<=1) = 1;
        
        
        p = ceil(length(Days)/step)*(j-1)+(i-1)/step+1;
        subplot(length(groups), ceil(length(Days)/step), p);
        
        axis off
        axis square
        
        hold on
        
        
        % local links
        cellfun(@(x) plot(x(:,1),x(:,2),'-','Color',[.7,.7,.7]),local_node_xy)
        
        
        % global links
        for k = 1:size(Path,1)
            plot(global_node_xy([Path(k,1),Path(k,2)],1),...
                global_node_xy([Path(k,1),Path(k,2)],2),...
                '-', 'LineWidth',pathpower(Path(k,1),Path(k,2)), 'Color',[.5,.5,.5]) %
        end
        
        % local nodes        
        v1 = global_node_xy-repmat(parameters.O,size(global_node_xy,1),1);
        v1 = mat2cell(v1,ones(size(v1,1),1),2);
        theta = cellfun(@(x) direction(parameters.O+[1,0],x), v1);
        theta = floor(arrayfun(@(x) 180-sign(x)*(180-abs(x)),theta))+1;
        d = cellfun(@(x) norm(x),v1)./parameters.RAD;
        
        cmap_1 = cmap_0(theta,:).*d;
        
        a = cat(1,local_node_xy{:});
        scatter(a(:,1),a(:,2),20,cmap_1,'filled')
        
        % global nodes
        % http://goo.gl/n3DOvD
        numPoints = 100; % Number of points making up the circle
        theta = linspace(0, 2*pi, numPoints); % 100 evenly spaced points between 0 and 2pi
        for k = 1:length(global_node_idx)
            radius = pwr(global_node_idx(k)); % Radius of the circle
            rho = ones(1,numPoints) * radius; % Radius should be 1 for all 100 points
            [x2,y2] = pol2cart(theta, rho);
            H = patch(x2+global_node_xy(global_node_idx(k),1),...
                y2+global_node_xy(global_node_idx(k),2), 1); % Make a patch object
            set(H, 'FaceColor',cmap_1(global_node_idx(k),:),...
                'EdgeColor',cmap_1(global_node_idx(k),:).*.7,...
                'FaceAlpha',.3) % Set the color to grey
        end


        if labelflag
            % global node label
            cmap = colormap(hot(100));
            for k = 1:length(global_node_idx)
                patch([gnode_xy(k,1)-55,gnode_xy(k,1)+55,gnode_xy(k,1)+55,gnode_xy(k,1)-55],...
                    [gnode_xy(k,2)+20,gnode_xy(k,2)+20,gnode_xy(k,2)-20,gnode_xy(k,2)-20],...
                    [.8,.8,.8], 'FaceAlpha',.5, 'EdgeColor','none')
                text(gnode_xy(k,1),gnode_xy(k,2),...
                    ['(',num2str(global_node_idx(k)),',',...
                    [['\color[rgb]{', num2str(cmap(nodepower(k),:)),'}'], num2str(N(k))],',',...
                    [['\color[rgb]{', num2str(cmap(linkpower(k),:)),'}'], num2str(NumGlobalLinks(k))],')'],...
                    'FontWeight','bold', 'FontName','MigMix',...
                    'HorizontalAlignment','center',...
                    'VerticalAlignment','middle');
            end
        end
        
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

switch options.maze_dia
    case 1
        xtick = [0,250,500];
        xticklabel = {'0','0.5','1'};
    case 3
        xtick = linspace(0,parameters.pix,4);
        xticklabel = num2cell(linspace(0,parameters.WIDTH/1000,4));
        xticklabel = cellfun(@(x) num2str(x,2), xticklabel,'uniform',0);
end

ytick = xtick;
yticklabel = xticklabel;


for i = 1:length(groups)*ceil(length(Days)/step)
    

    subplot(length(groups), ceil(length(Days)/step), i);
    axis on
    axis square
        
    if mod(i-1,length(Days))==0
        ylabel('y (m)', 'Interpreter','none', 'FontWeight','bold')
        xlabel('x (m)', 'Interpreter','none', 'FontWeight','bold')
    end
    
    set(gca, 'XLim',xtick([1,end]), 'XTick',xtick, 'XTickLabel',xticklabel,...
        'YLim',ytick([1,end]), 'YTick',ytick, 'YTickLabel',yticklabel,...
        'box','on', 'FontSize',24, 'TickLength',[0.03,0.025], 'LineWidth',2)
    
end
set(gcf,'color','w');
