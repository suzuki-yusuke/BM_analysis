clear;
close all
clc


%% set parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
start_path = 'D:\documents\data';
options.maze_dia = 3; % BM type
options.offset = true; % offset; 32 in BM3, 41 in BM1
options.cutoff = 30; % cutoff = duration (s) till micro-movements around goals are supposed as goal entry
options.champs = []; % selet certain serial numbers to visualize whose local network
% options.order = order_v0();
options.timeout = false;
options.view_network_graph = false;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%% load data
addpath(genpath(pwd))

switch options.maze_dia
    case 1
        options.step = 1; % x-axis view step (day)
        parameters = parameters_BM1; % make an object binding parameters
    case 3
        options.step = 3; % x-axis view step (day)
        parameters = parameters_BM3; % make an object binding parameters
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

%%
N = size(tbl,1);
threshold = parameters.bwHoles/2.5;
for i = 1:N
    
    switch tbl.Group{i}
        case 'BM1'
            options.step = 1; % x-axis view step (day)
            parameters = parameters_BM1; % make an object binding parameters
        case 'BM3'
            options.step = 2; % x-axis view step (day)
            parameters = parameters_BM3; % make an object binding parameters
    end
    
    XY = tbl.xy{i};
    T = size(XY,1);
    
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
    
    T = size(XY,1);
    tbl.latency(i) = T/parameters.fps;
   
    
     % No. of errors
    threshold = parameters.bwHoles/2.5;
    XYC = mat2cell(XY,ones(T,1),2);
    
    d = cellfun(@(x) repmat(x,parameters.nHoles,1)-parameters.holes, XYC, 'uniform',0);
    d = cellfun(@(x) mat2cell(x,ones(parameters.nHoles,1),2),d,'uniform',0);
    d = cellfun(@(x) cellfun(@norm,x), d,'uniform',0);
    d = cellfun(@(x) find(x<threshold),d,'uniform',0);    
    d = cat(1,d{:});
    
    b = find(diff(d)~=0);
    
    if isempty(b)
        d = {[]};
    elseif length(b)==1
        d = mat2cell(d,[b,length(d)-b]);
    else        
        a = [1;b]; b = [b-1;length(d)];
        d = mat2cell(d,b-a+1,1);
    end    

    d(cellfun(@(x) length(x)<(parameters.fps*1),d)) = [];        
    d = cat(1,d{:});
    
    d(d==11) = [];
    
    l = diff(XY);
    l = mat2cell(l,ones(T-1,1),2);
    l = cellfun(@norm,l);
    l = sum((l.*parameters.spr)/10);
    if isempty(d)&&(l<=(parameters.RAD*parameters.spr)) % immobile around center
        d = 0;
    elseif isempty(d)&&(l>(parameters.RAD*parameters.spr))
        d = 0;
    else
        d = sum(diff(d)~=0,'omitnan');
    end
    
    if isnan(d)
        d = 0;
    end
    
    tbl.no_of_errors(i) = d;
    
end

X = tbl(:,{'Group','no_of_errors','latency'});
X = sortrows(X,{'Group'});
writetable(X,[PathName,'err&lat_prb.csv'],'delimiter',',')
save([PathName,'data_probe_test_errors_latency.mat'],'tbl')



%% compaarison of no. of erros and latentcy in the probe test
[~,~,ic] = unique(tbl.Group,'stable');
[Report,Level,M,N,F,p,MS,df,ES] = anova_1b([ic,tbl.no_of_errors])

[~,~,ic] = unique(tbl.Group,'stable');
[Report,Level,M,N,F,p,MS,df,ES] = anova_1b([ic,tbl.latency])

