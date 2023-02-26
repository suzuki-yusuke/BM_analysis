function P = direction(v0,v1)

% �Q�̃x�N�g���̂Ȃ��p
% v0 = [1, 0]; % �x�N�g�� 0
% v1 = [0, -1]; % �x�N�g�� 1

InnerProduct = v0(1)*v1(1) + v0(2)*v1(2);
CrossProduct = v0(1)*v1(2) - v0(2)*v1(1);

% Norm1 = norm(v1);
% Norm2 = norm(v0);
% cos = InnerProduct / (Norm1*Norm2);
% P = acos(cos)*180/pi; % 0:180
P = atan2(CrossProduct,InnerProduct)*180/pi; % -180:180
