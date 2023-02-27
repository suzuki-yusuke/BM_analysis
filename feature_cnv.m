function ConventionalFeatures = feature_cnv_v2(parameters,options,ConventionalFeatures,n,phase,XY)

T = size(XY,1);
XYC = mat2cell(XY,ones(T,1),2);

% xy coordinates
ConventionalFeatures.xy{n} = XY;


% distance travel
l = diff(XY);
l = mat2cell(l,ones(T-1,1),2);
l = cellfun(@norm,l);
l = sum((l.*parameters.spr)/10);
ConventionalFeatures.travel_distance(n) = l;


% No. of errors
threshold = parameters.bwHoles/2.5;

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

ConventionalFeatures.no_of_errors(n) = d;

if ~strcmp(phase, 'probe_test') % training, retraining, reversal_training �̏ꍇ        
    
    % latency
    ConventionalFeatures.latency(n) = T/parameters.fps;        
       
    
else % probe_test �̏ꍇ
        
    % latency
    cutoff = 10; % duration of continuous micro-movement around goal
    threshold = parameters.bwHoles/2.5;    
    target = parameters.holes(1,:); target = repmat(target, T,1);
    d = mat2cell(XY-target, ones(T,1),2);
    d = cellfun(@norm, d);
    
    fin = 0;
    initframe = find(d<threshold,1,'first'); endframe = 1;
    if and(initframe~=T, T>cutoff*parameters.fps)
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
    ConventionalFeatures.latency(n) = endframe/parameters.fps;
    
            
    % time spent around each hole
    threshold = parameters.bwHoles/2.5; % default = 40
    
    d = cellfun(@(x) repmat(x,parameters.nHoles,1)-parameters.holes, XYC, 'uniform',0);
    d = cellfun(@(x) mat2cell(x,ones(parameters.nHoles,1),2),d,'uniform',0);
    d = cellfun(@(x) ...
        cellfun(@norm,x),...
        d,'uniform',0);
    d = cellfun(@(x) find(x<threshold),d,'uniform',0);
    TF = cellfun(@isempty,d);
    d(TF) = [];
    d = vertcat(d{:});    
    
    visitcount = histcounts(d,.5:parameters.nHoles+.5);
    visitcount./parameters.fps;
    ConventionalFeatures.approach(n) = {visitcount./parameters.fps};
    
end