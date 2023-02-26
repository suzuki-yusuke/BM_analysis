clear;
clc
close all



%% set parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
start_path = 'D:\documents\data';
options.maze_dia = 1; % BM type
options.offset = true; % offset; 32 in BM3, 41 in BM1
options.cutoff = 10; % cutoff duration (s) of micro-movements around goal supposed as goal entry
options.timeout = false;
options.view_network_graph = false;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




addpath(genpath(pwd))

switch options.maze_dia
    case 1
        options.step = 1; % x-axis view step (day)
        parameters = parameters_BM1; % make an object binding parameters
        parameters.cutoff_probe = 150; %60, 90, 150
    case 3
        options.step = 2; % x-axis view step (day)
        parameters = parameters_BM3; % make an object binding parameters
        parameters.cutoff_probe = 150; %240, 150, 90
end



[FileName,PathName] = uigetfile('*.mat', 'Select the LOG file', start_path);

k = strfind(PathName,'\');
PrjDir = PathName(1:k(end-1)-1);
TrnDir = dir(fullfile(PrjDir, 'training'));
PrbDir = dir(fullfile(PrjDir, 'probe_test'));
DatDir = [{TrnDir.name}, {PrbDir.name}];

load(fullfile(PathName,FileName));
parameters.PrjDir = PrjDir;
parameters.DatDir = DatDir;
parameters.groups = unique(tbl.Group);
parameters.SN = unique(str2double(tbl.SN));

N = size(tbl,1);

ConventionalFeatures = table;
ConventionalFeatures.xy = cell(N,1);
ConventionalFeatures.approach = cell(N,1);
ConventionalFeatures{:,{'no_of_errors','latency','travel_distance'}} = nan(N,3);


NetworkFeatures = table;
NetworkFeatures.xy_1 = cell(N,1);
NetworkFeatures{:,{'xy_2','A_2'}} = cell(N,2);
NetworkFeatures{:,{'o_2','n_2','l_2','m_2','CWS_2','rho_2','x_2','cl_2'}} = nan(N,8);


%% change duration of the probe test
% options.order = order_v0();
% 
% f = waitbar(0,'');
% for n = 1:N
%     XY = tbl.xy{n}(1:(parameters.cutoff_probe*parameters.fps),:);
%     phase = tbl.Phase{n};
%     tbl.xy{n} = XY;
%     ConventionalFeatures =...
%             feature_cnv_v2(parameters, options, ConventionalFeatures, n, phase, XY);    
%     tbl.approach{n} = ConventionalFeatures.approach{n};
%         
%     NetworkFeatures = feature_ntw_d_v0(parameters,NetworkFeatures,n,XY);
%     a = find(strcmp('xy_2',tbl.Properties.VariableNames));
%     b = find(strcmp('cl_2',tbl.Properties.VariableNames));
%     tbl(:,a:b) = NetworkFeatures(:,2:end);
%     
%     waitbar(n/N,f);
% end
% close(f)
% 
% %save([PathName,'\data_probe_test_cutoff.mat'],'tbl')



%% get number of errors in the probe test
threshold = parameters.bwHoles/2.5;        
f = waitbar(0,'');
for n = 1:N
    
    XY = tbl.xy{n};
    T = size(XY,1);
    
    if ~isempty(options.cutoff)
        
        target = parameters.holes(1,:); target = repmat(target, T,1);
        d = mat2cell(XY-target, ones(T,1),2);
        d = cellfun(@norm, d);
        
        fin = 0;
        initframe = find(d<threshold,1,'first'); endframe = 1;
                
        if and(initframe~=T, T>options.cutoff*parameters.fps)
            while ~fin
                while ~fin
                    fin = or(and(d(initframe+endframe,:)<threshold,...
                        endframe>options.cutoff*parameters.fps),...
                        initframe+endframe>=T-1);
                    endframe = endframe + 1;
                end
                initframe = initframe + endframe - 1;
                endframe = 1;
            end
        else
            initframe = T+1;
        end
        endframe = initframe-1;
        XY(endframe:end,:) = [];
    end
    
    
    switch options.maze_dia
        case 1
            if length(XY)>(parameters.cutoff_probe*parameters.fps)
              XY = XY(1:parameters.cutoff_probe*parameters.fps,:);
            end
        case 3
%             if length(XY)>(90*parameters.fps)
%                 XY = tbl.xy{n}(1:90*parameters.fps,:);
%             end
    end
    
    ConventionalFeatures =...
        feature_cnv_v2(parameters, options, ConventionalFeatures, n, 'training', XY);
    tbl.no_of_errors(n) = ConventionalFeatures.no_of_errors(n);
    tbl.latency(n) = ConventionalFeatures.latency(n);
    tbl.travel_distance(n) = ConventionalFeatures.travel_distance(n);
    
    ConventionalFeatures =...
        feature_cnv_v2(parameters, options, ConventionalFeatures, n, 'probe_test', XY);
    tbl.approach(n) = ConventionalFeatures.approach(n);
    
    NetworkFeatures = feature_ntw_d_v0(parameters,NetworkFeatures,n,XY);
    a = find(strcmp('xy_2',tbl.Properties.VariableNames));
    b = find(strcmp('cl_2',tbl.Properties.VariableNames));
    tbl(:,a:b) = NetworkFeatures(:,2:end);
    
    waitbar(n/N,f);
end
close(f)

writetable(tbl(:,{'Phase','Day','Group','SN','no_of_errors','latency','travel_distance'}),...
    [PathName,sprintf('\\BM%d_reanalyze_probetest.csv',options.maze_dia)])

save([PathName,sprintf('\\data_probe_test_BM%d_normalized.mat',options.maze_dia)],'tbl')