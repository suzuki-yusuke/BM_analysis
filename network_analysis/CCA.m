function [unions,nb,na] = CCA(individuals, criterion)


% criterion = 20; % unifying criterion
unions = individuals;

na = 1:size(unions,1);
nb = na;
pairs = nchoosek(nb, 2); % ���j�[�N�m�[�h�Ԃ̑S�Ă̑g�ݍ��킹

X = 1;
ub = [];
uc = [];

while any(X)
    
    dist = unions(pairs(:,1), :) - unions(pairs(:,2), :);
    dist = mat2cell(dist, ones(size(dist,1),1), 2);
    dist = cellfun(@norm, dist);
    
    
    % �m�[�h�Ԃ̋�����臒l�ȉ��ƂȂ�m�[�h�̃y�A�𒊏o
    X = and(dist<=criterion, dist>0);
    Y = pairs(X,:);
    
    % ���łɂ����ꂩ�̃m�[�h�Ɋ܂܂�Ă���m�[�h������
    [Cy,~,~] = unique(Y(:,2), 'stable');
    Y(ismember(Y(:,1),Cy),:) = [];
    
    g = unique(Y(:,1)); % �����m�[�h�̃C���f�b�N�X�𒊏o
    NM = length(g); % ��������m�[�h�̐�
    
    % �m�[�h�̃N���X�^�[�Əd�S���W���Z�o
    for k = 1:NM
        
        M = Y(Y(:,1)==g(k), 2); % ���Ӄm�[�h
        p = zeros(length(M)+1, 2);
        p(1,:) = unions(g(k), :);
        p(2:end,:) = unions(M, :);
        ua = [g(k); M];
        ub = [ub; ua]; % ���S�m�[�h�̃C���f�b�N�X
        uc = [uc; repmat(k,length(ua),1)]; % ���S�m�[�h�̏C�����ꂽ�C���f�b�N�X
        
        nb(ua) = min(ua);
        na(ua) = k;
        unions(ua, :) = repmat(mean(p), [length(ua), 1]);
        
    end
    
end

[~, free, ~] = setxor(nb, ub); % ��������Ă��Ȃ��m�[�h�̃C���f�b�N�X
[~, ia, ~] = unique(ub, 'last'); % �ŏI�I�Ȓ��S�m�[�h�̃C���f�b�N�X
uc = uc(ia); % �ŏI�I�Ȓ��S�m�[�h�̏C�����ꂽ�C���f�b�N�X
free_indx = 1:size(individuals,1);
[~, free_indx, ~] = setxor(free_indx, uc); % ���S�m�[�h�̃C���f�b�N�X������
free_indx = free_indx(1:length(free));
na(free) = free_indx;