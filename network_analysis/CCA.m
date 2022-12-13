function [unions,nb,na] = CCA(individuals, criterion)


% criterion = 20; % unifying criterion
unions = individuals;

na = 1:size(unions,1);
nb = na;
pairs = nchoosek(nb, 2); % ユニークノード間の全ての組み合わせ

X = 1;
ub = [];
uc = [];

while any(X)
    
    dist = unions(pairs(:,1), :) - unions(pairs(:,2), :);
    dist = mat2cell(dist, ones(size(dist,1),1), 2);
    dist = cellfun(@norm, dist);
    
    
    % ノード間の距離が閾値以下となるノードのペアを抽出
    X = and(dist<=criterion, dist>0);
    Y = pairs(X,:);
    
    % すでにいずれかのノードに含まれているノードを除去
    [Cy,~,~] = unique(Y(:,2), 'stable');
    Y(ismember(Y(:,1),Cy),:) = [];
    
    g = unique(Y(:,1)); % 統合ノードのインデックスを抽出
    NM = length(g); % 統合するノードの数
    
    % ノードのクラスターと重心座標を算出
    for k = 1:NM
        
        M = Y(Y(:,1)==g(k), 2); % 周辺ノード
        p = zeros(length(M)+1, 2);
        p(1,:) = unions(g(k), :);
        p(2:end,:) = unions(M, :);
        ua = [g(k); M];
        ub = [ub; ua]; % 中心ノードのインデックス
        uc = [uc; repmat(k,length(ua),1)]; % 中心ノードの修正されたインデックス
        
        nb(ua) = min(ua);
        na(ua) = k;
        unions(ua, :) = repmat(mean(p), [length(ua), 1]);
        
    end
    
end

[~, free, ~] = setxor(nb, ub); % 統合されていないノードのインデックス
[~, ia, ~] = unique(ub, 'last'); % 最終的な中心ノードのインデックス
uc = uc(ia); % 最終的な中心ノードの修正されたインデックス
free_indx = 1:size(individuals,1);
[~, free_indx, ~] = setxor(free_indx, uc); % 中心ノードのインデックスを除く
free_indx = free_indx(1:length(free));
na(free) = free_indx;