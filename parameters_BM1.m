function parameters = parameters_BM1()

% field parameters
WIDTH = 980; % mm
nHoles = 12; % 穴の数
nMidNodes = 12; % 中間ノードの数
nNodes = nHoles + nMidNodes + 1; % ノードの総数
PossibleLinks = nchoosek(1:nNodes, 2);
nPossibleLinks = 2*nchoosek(nNodes, 2);


% recording parameters
fps = 20;
pix = 500; % pixel 数
spr = WIDTH/pix; % 空間解像度，x mm / pix
DIA = WIDTH/spr; % フィールドの直径，pix表示
RAD = DIA/2; % フィールドの半径，pix表示
O = [RAD,RAD]; % フィールドの中心座標
dia = 40/spr; % 穴の直径，pix表示
diaInitPos = 130/spr;
rad = dia/2; % 穴の半径，pix表示
Edge2Holes = 90/spr; % フィールドの外周から穴の円周までの最短距離，pix表示
Center2Holes = RAD - (Edge2Holes + rad); % フィールドの中心から穴の中心までの距離，pix表示
Center2Holes_H = Center2Holes + rad; % フィールドの中心から穴の上までの距離，pix表示
Center2Holes_L = Center2Holes - rad; % フィールドの中心から穴の下までの距離，pix表示
diaStartBox = 130/spr; %

thresh4CCA = 80/spr;


%各穴の座標を計算
%{
GA = (250, 250)から(1, 1)へのベクトルを，
30°ずつ回転させた時のベクトルと座標と大きさを設定
x' = x*cos(theta) - y*sin(theta)
y' = x*sin(theta) + y*cos(theta)
%}
theta = (0:30:330)';
holes = zeros(length(theta), 2);
holes(:,1) = ((cosd(theta)-sind(theta))*...
    (Center2Holes/sqrt(2))+RAD)';
holes(:,2) = ((sind(theta)+cosd(theta))*...
    (Center2Holes/sqrt(2))+RAD)';
bwHoles = sqrt(Center2Holes^2-(Center2Holes*sind(75))^2) * 2; % distance b/w holes




%中間ノードの座標を計算
%{
GA = (250, 250)から(1, 1)へのベクトルを，
30°ずつ回転させた時のベクトルと座標と大きさを設定
x' = x*cos(theta) - y*sin(theta)
y' = x*sin(theta) + y*cos(theta)
%}
CenterEdge = (Center2Holes-(bwHoles/2))/2; % 中心ノードの円周
Center2MidNodes = CenterEdge + CenterEdge/2;
theta = (0:30:330)';
middle_nodes = zeros(length(theta), 2);
middle_nodes(:, 1) = ((cosd(theta)-sind(theta))*...
    (Center2MidNodes/sqrt(2))+RAD)';
middle_nodes(:, 2) = ((sind(theta)+cosd(theta))*...
    (Center2MidNodes/sqrt(2))+RAD)';
bwMidNodes = sqrt(Center2MidNodes^2-(Center2MidNodes*sind(75))^2) * 2;



% make structure of parameter
c = who;
params_cell = cell(length(c),1);
for i = 1:length(c)
    if strcmp(c{i}, 'c')
        continue
    end
    params_cell(i) = {eval(c{i})};
end
parameters = cell2struct(params_cell, c, 1);

