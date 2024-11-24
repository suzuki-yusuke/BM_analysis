function StrategyFeatures = feature_str(parameters,~,~,StrategyFeatures,n,~,XY)

%% Calculate Strategy Indices
%{
Nat Commun. 2014 Aug 19;5:4701. doi: 10.1038/ncomms5701.
The MK2/3 cascade regulates AMPAR trafficking and cognitive flexibility.
Eales KL, Palygin O, O'Loughlin T, Rasooli-Nejad S, Gaestel M, M?ller J, Collins DR, Pankratov Y, Corr?a SA.
%}

T = size(XY,1);
holes = parameters.holes;
startframe = 1;


%% discard initial trajectory
% switch options.maze_dia
%     case 1
%         threshold = 160/parameters.spr; % dia of the lift = 10 cm, body length = 80 mm
% %         threshold = parameters.bwHoles/2.5;
%     case 3
%         threshold = 160/parameters.spr; % dia of the lift is 11.6 cm
% %         threshold = parameters.bwHoles/2.5;
% end
% if OnlineTracking
%     d = 0;
%     while (d<=threshold)&&(startframe<T)
%         d = norm(XY(startframe,:) - parameters.O);
%         startframe = startframe + 1;
%     end
% end

XY = XY(startframe:end, :);
T = size(XY,1);



%% quadrant transition
xy = XY-ones(T,1)*parameters.O;
Q = zeros(T,1);
for i = 1:T
    Sign = [sign(xy(i,1)), sign(xy(i,2))];
    if (Sign(1)==1)&&(Sign(2)==1)
        Q(i) = 1;
    elseif (Sign(1)==1)&&(Sign(2)==-1)
        Q(i) = 4;
    elseif (Sign(1)==-1)&&(Sign(2)==1)
        Q(i) = 2;
    elseif (Sign(1)==-1)&&(Sign(2)==-1)
        Q(i) = 3;
    end
end
quadrant_visit = unique(Q(Q>0));



%% hole transition
threshold = 80/parameters.spr;
H = zeros(T,1);
for i = 1:parameters.nHoles
    Hi = zeros(T,1);
    for j = 1:T
        Hi(j) = norm(XY(j,:) - holes(i,:));
    end
    H(Hi<=threshold) = i;
end


% sequence of hole visit
if isempty(H)
    hole_visit = 0;
else
    k = [1;(diff(H)~=0)];
    hole_visit = H(logical(k));
    hole_visit(hole_visit==0) = [];
end




%% perimeter walk
%{
perimeter walk is a trajectory of outer area except 
the areas around the goal and the reversal goal.

distance to inner edge of outer area from the center
= 1980 mm = 3000 mm dia * 0.66 (BM3)
%}
threshold = parameters.WIDTH*0.66/2/parameters.spr;
P = false(T,1);
for i = 1:T
    P(i) = all([...
        norm(XY(i,:)-parameters.O)>threshold,...
        norm(XY(i,:)-holes(1,:))>(80/parameters.spr),...
        norm(XY(i,:)-holes(7,:))>(80/parameters.spr)]);
end



%% head angle
XY_smooth = [smooth(XY(:,1),0.25,'lowess'),smooth(XY(:,2),0.25,'lowess')];
angle = zeros(size(XY_smooth,1),1);
for i = 2:size(XY_smooth,1)
    dxy = XY_smooth(i,:)-XY_smooth(i-1,:);
    angle(i) = abs(atan2d(dxy(2), dxy(1)));
end




%%

spatial_flag = all([all(ismember(hole_visit,[12,1,2])),...
    sum(ismember(hole_visit,[12,2]))<3,...
    all(ismember(quadrant_visit,[1]))]);

serial_flag = all([all(ismember(abs(diff(hole_visit)),[1,2,10,11])),...
    all(ismember(quadrant_visit,[1,2,4]))]);


confirmatory_flag = all([any(ismember(hole_visit,[12,1,2])),...
    any(ismember(hole_visit,[6,7,8])),...
    (sum(ismember(H,[3,4,5,9,10,11]))/(sum(H>0)))<0.25,...
    (sum(ismember(Q,[3]))/(T-sum(H>1)))>0.25]);

perimeter_flag = all([sum(P)/(T-sum(H==1))>0.75,...
    length(unique(quadrant_visit))>=3]);




% Algorithm-based classification (cf. Garthe, et al., 2009)
if spatial_flag
    StrategyFeatures.Spatial(n) = 1;
elseif confirmatory_flag
    StrategyFeatures.Confirmatory(n) = 1;
elseif serial_flag
    StrategyFeatures.Serial(n) = 1;
elseif perimeter_flag
    StrategyFeatures.Perimeter(n) = 1;
else
    StrategyFeatures.Random(n) = 1;
end


