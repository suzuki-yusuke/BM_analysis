function ConventionalFeatures = feature_cnv_v2(parameters,options,ConventionalFeatures,n,phase,XY)

T = size(XY,1);
XYC = mat2cell(XY,ones(T,1),2);

% xy coordinates
ConventionalFeatures.xy{n} = XY;

if ~strcmp(phase, 'probe_test') % training, retraining, reversal_training �̏ꍇ        
    
    % latency
    ConventionalFeatures.latency(n) = T/parameters.fps;
    
    
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
    
else % probe_test �̏ꍇ
    
%     threshold = parameters.bwHoles/2.5; % default = 40
    threshold = parameters.bwHoles/2.5; % default = 40
%     threshold = 300/paramObj.spr; % default = 40
    
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