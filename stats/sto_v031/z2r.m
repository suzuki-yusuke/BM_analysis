function [r] = z2r(z)

% z2r: フィッシャーのz変換の逆
% 
% [r] = z2r(z)
%    z   z値（-Inf 〜 +Inf）
%    r   変換された相関係数（-1 〜 +1）
% 
% 入力された z に対して、(exp(2*z)-1) / (exp(2*z)+1)
% を計算します。zは行列可です。
% 
% see also:   r2z
% 
% (2015/09/03, by R. NIIMI)

r = [];
if exp(2*z)==Inf
	r = 1;
else
	r = (exp(2*z)-1) ./ (exp(2*z)+1);
end

