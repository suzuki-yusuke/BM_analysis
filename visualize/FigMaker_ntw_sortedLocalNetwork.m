function [] = FigMaker_ntw_sortedLocalNetwork(tbl,parameters,options)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
step = options.step; % step between days
champs = options.champs;
labelflag = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


T = tbl;


%% behavioral network
Days = unique(T.Day);
daynames = num2cell(Days);
daynames = cellfun(@num2str,daynames,'uniform',0);
groups = unique(T.Group);
phase = unique(T.Phase);

scrsz = get(groot,'ScreenSize');


if 1 % local network sample を描画する場合
    polc = [255,151,163;...
        69,255,193;...
        198,228,193;...
        153,201,255;...
        198,228,255]./255;
else
    polc = polcmap(4,.9);
end



a = 1:(length(groups)+1)*(length(Days)+1);
a = reshape(a,length(Days)+1,length(groups)+1)';
I = a; I(:,1) = []; I(1,:) = []; I = I'; ind = I(:);
G = a(2:end,1);
D = a(1,2:end);


for champ = 1:length(champs)
    
    figure('Position', scrsz.*[1,1,.9,.7]);
    
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
            reading = reading&(T.SN==champs(champ));
            
            
            localdata = T{reading, 'xy_2'};
            
            a = cellfun(@iscell, localdata);
            b = [localdata{a}];
            localdata(a) = b;
            map_local = vertcat(localdata{:});
            
            idnums = T{reading, 'SN'};
            NumNodes = cellfun(@size, vertcat(localdata(:)), 'UniformOutput',0);
            NumNodes = vertcat(NumNodes{:}); NumNodes = NumNodes(:,1); % local network のノード数
            
            PathBwTrials = cumsum(NumNodes); PathBwTrials(end) = [];
            
            idlabels = zeros(sum(NumNodes),1);
            for k = 1:length(NumNodes)
                idlabels(find(idlabels==0,1):find(idlabels==0,1)+NumNodes(k)-1) =...
                    repmat(idnums(k),NumNodes(k),1);
            end
            
            
            if sum(NumNodes)>1 
                % CCA
                [~,lnode2gnode,~] = CCA(map_local, parameters.thresh4CCA);
                lnode2gnode = lnode2gnode';
                [gnode_idx,~,ind2gnode] = unique(lnode2gnode);
                gnode_idx = gnode_idx';
            else
                lnode2gnode = 1;
                ind2gnode = 1;
                gnode_idx = 1;
            end
            
            N = histcounts(ind2gnode, .5:max(ind2gnode)+.5);
            N = N(ind2gnode); % number of nodes contained in each global node
            
            nodepower = floor(N./sum(N).*100);
            nodepower(nodepower<=1) = 1;
            
            
            % 各使用経路の使用数
            % find crossing links between 2 successive trials
            if sum(NumNodes)>1
                Path = [lnode2gnode(1:end-1),lnode2gnode(2:end)];
                Path(PathBwTrials,:) = [];
                Path = [Path;fliplr(Path)]; Path = unique(Path,'rows');
            else
                Path = [1,1];
            end

            pathusage = histcounts(Path(:,1),.5:max(Path(:,1))+.5);                                 
            pathpower = log(pathusage) + 1; % 対数変換
            
            % 各ノードのリンク数
            NumglobalLinks = zeros(length(gnode_idx),1);
            for k = 1:length(gnode_idx)
                NumglobalLinks(k) = sum((gnode_idx(k)==Path(:,1)));
            end
            linkpower = floor(NumglobalLinks./sum(NumglobalLinks).*100);
            linkpower(linkpower<=1) = 1;
            
            
            
            
            p = ceil(length(Days)/step)*(j-1)+(i-1)/step+1;
            subplot(length(groups), ceil(length(Days)/step), p);
            
            axis off
            axis square
            
            hold on
            
            
            [~,I] = sort(NumglobalLinks,'descend');
            [~,I] = sort(I); % link usage ranking
            
            % global nodes
            % http://goo.gl/n3DOvD
            numPoints = 1000; % Number of points making up the circle
            theta_topo = linspace(0,2*pi,numPoints); % 1000 evenly spaced points between 0 and 2pi
            radius_topo = 200; % Radius of the circle
            rho_topo = ones(1,numPoints) * radius_topo; % Radius should be 1 for all 100 points
            [topox,topoy] = pol2cart(theta_topo, rho_topo); Step = floor(numPoints/length(gnode_idx));
            topox = topox(1:Step:end); topoy = topoy(1:Step:end);
            toponodes = [topox',topoy'];
            
            
            % global links
            for k = 1:size(Path,1)-1
                plot(toponodes([I(ind2gnode(Path(k,1))),I(ind2gnode(Path(k,2)))], 1),... % from
                    toponodes([I(ind2gnode(Path(k,1))),I(ind2gnode(Path(k,2)))], 2),... % to
                    '-', 'LineWidth',1, 'Color',[.4,.4,.4])
            end
            
            cmap = colormap(hot(6));
            for k = 1:length(gnode_idx)
                H = plot(toponodes(I(k),1), toponodes(I(k),2), 'o'); % Make a patch object
                c = [NumglobalLinks(k)<3,...
                    and(3<=NumglobalLinks(k), NumglobalLinks(k)<6),...
                    and(6<=NumglobalLinks(k), NumglobalLinks(k)<9),...
                    and(9<=NumglobalLinks(k), NumglobalLinks(k)<12),...
                    and(12<=NumglobalLinks(k), NumglobalLinks(k)<15),...
                    15<=NumglobalLinks(k)];
                set(H, 'MarkerFaceColor',cmap(c,:),...
                    'MarkerEdgeColor',[.5,.5,.5],...
                    'MarkerSize',8)
            end
            
            
            if labelflag
                %global node label
                
                theta_label = linspace(0,2*pi,numPoints); % 1000 evenly spaced points between 0 and 2pi
                radius_label = 200; % Radius of the circle
                rho_label = ones(1,numPoints) * radius_label; % Radius should be 1 for all 100 points
                [labelx,labely] = pol2cart(theta_label, rho_label); Step = ceil(numPoints / length(gnode_idx));
                labelx = labelx(1:Step:end); labely = labely(1:Step:end);
                topolabel = [labelx',labely'];
                
                cmap = colormap(hot(100));
                for k = 1:length(gnode_idx)
                    patch([topolabel(k,1)-55,topolabel(k,1)+55, topolabel(k,1)+55,topolabel(k,1)-55],...
                        [topolabel(k,2)+20,topolabel(k,2)+20,topolabel(k,2)-20,topolabel(k,2)-20],...
                        [.8,.8,.8], 'FaceAlpha',.5, 'EdgeColor','none')
                    
                    text(topolabel(I(k),1),topolabel(I(k),2),...
                        ['(',num2str(gnode_idx(k)),',',...
                        [['\color[rgb]{', num2str(cmap(nodepower(k),:)),'}'], num2str(N(k))],',',...
                        [['\color[rgb]{', num2str(cmap(linkpower(k),:)),'}'], num2str(NumglobalLinks(k))],')'],...
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
        title(['Day ',daynames{b}], 'FontSize',18, 'interpreter','none')
    end
    
    
    xtick = linspace(-210,210,4);
    xticklabel = '';
    
    ytick = xtick;
    yticklabel = xticklabel;
    
    
    for i = 1:length(groups)*ceil(length(Days)/step)
        
        subplot(length(groups), ceil(length(Days)/step), i);
        axis on
        axis square
        
        if mod(i-1,length(Days))==0
            ylabel(groups{ceil(i/length(Days))}, 'Interpreter','none')
        end
        
        set(gca, 'XLim',xtick([1,end]), 'XTick',xtick, 'XTickLabel',xticklabel,...
            'YLim',ytick([1,end]), 'YTick',ytick, 'YTickLabel',yticklabel,...
            'box','on', 'FontSize',18)
    end
    
end

