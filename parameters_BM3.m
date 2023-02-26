function parameters = parameters_BM3()

% field parameters
WIDTH = 3000; % mm
nHoles = 12; % 穴の数
nMidNodes = 12; % 中間ノードの数
nNodes = nHoles + nMidNodes + 1; % ノードの総数
PossibleLinks = nchoosek(1:nNodes, 2);
nPossibleLinks = 2*nchoosek(nNodes, 2);

fps = 20;
pix = 400; % pixel 数
spr = WIDTH / pix; % 空間解像度，x mm / pix
DIA = WIDTH / spr; % フィールドの直径，pix表示
RAD = DIA / 2; % フィールドの半径，pix表示
O = [RAD,RAD]; % フィールドの中心座標
dia = 160 / spr; % 穴の直径，pix表示
diaInitPos = 130 / spr;
rad = dia / 2; % 穴の半径，pix表示
Edge2Holes = 125 / spr; % フィールドの外周から穴の円周までの最短距離，pix表示
Center2Holes = RAD - (Edge2Holes + rad); % フィールドの中心から穴の中心までの距離，pix表示
Center2Holes_H = Center2Holes + rad; % フィールドの中心から穴の上までの距離，pix表示
Center2Holes_L = Center2Holes - rad; % フィールドの中心から穴の下までの距離，pix表示

thresh4CCA = 30; % default = 80/2/spr


%各穴の座標を計算
%{
LabVIEW の物体検出プログラムで抽出された各穴の座標を時計回りでソート
%}
listing = dir('*.txt');
listing = {listing.name}; listing = listing(end);
holes = dlmread(listing{:},'\t');
holes = circshift(holes,[-4,0]);
holes = holes(:,2:3);
w = repmat([0.86,0.86,0.86,0.87,...
    0.89,0.9,0.9,0.9,...
    0.89,0.88,0.87,0.86],2,1)';
p = (holes-O).*w+O;

bwHoles = sqrt(Center2Holes^2-(Center2Holes*sind(75))^2) * 2; % distance b/w holes

direct = [0,0, pix,0, pix,pix, 0,pix];

d = repmat(holes,1,4)-repmat(direct,size(holes,1),1);
d = mat2cell(d, ones(size(holes,1),1), ones(4,1).*2);
d = cellfun(@norm, d);
[~,GoalHoleNum] = min(d);


CenterEdge = (Center2Holes-(bwHoles/2))/2; % 中心ノードの円周

%中間ノードの座標を計算
%{
center から holes の a 倍の位置に生成されたノード
%}
a = 0.5;
middle_nodes = a.*holes;
middle_nodes = middle_nodes + ones(size(middle_nodes,1),1)*(O-mean(middle_nodes));


holes = p;



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
