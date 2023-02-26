function parameters = parameters_BM1()

% field parameters
WIDTH = 980; % mm
nHoles = 12; % ���̐�
nMidNodes = 12; % ���ԃm�[�h�̐�
nNodes = nHoles + nMidNodes + 1; % �m�[�h�̑���
PossibleLinks = nchoosek(1:nNodes, 2);
nPossibleLinks = 2*nchoosek(nNodes, 2);


% recording parameters
fps = 20;
pix = 500; % pixel ��
spr = WIDTH/pix; % ��ԉ𑜓x�Cx mm / pix
DIA = WIDTH/spr; % �t�B�[���h�̒��a�Cpix�\��
RAD = DIA/2; % �t�B�[���h�̔��a�Cpix�\��
O = [RAD,RAD]; % �t�B�[���h�̒��S���W
dia = 40/spr; % ���̒��a�Cpix�\��
diaInitPos = 130/spr;
rad = dia/2; % ���̔��a�Cpix�\��
Edge2Holes = 90/spr; % �t�B�[���h�̊O�����猊�̉~���܂ł̍ŒZ�����Cpix�\��
Center2Holes = RAD - (Edge2Holes + rad); % �t�B�[���h�̒��S���猊�̒��S�܂ł̋����Cpix�\��
Center2Holes_H = Center2Holes + rad; % �t�B�[���h�̒��S���猊�̏�܂ł̋����Cpix�\��
Center2Holes_L = Center2Holes - rad; % �t�B�[���h�̒��S���猊�̉��܂ł̋����Cpix�\��
diaStartBox = 130/spr; %

thresh4CCA = 80/spr;


%�e���̍��W���v�Z
%{
GA = (250, 250)����(1, 1)�ւ̃x�N�g�����C
30������]���������̃x�N�g���ƍ��W�Ƒ傫����ݒ�
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




%���ԃm�[�h�̍��W���v�Z
%{
GA = (250, 250)����(1, 1)�ւ̃x�N�g�����C
30������]���������̃x�N�g���ƍ��W�Ƒ傫����ݒ�
x' = x*cos(theta) - y*sin(theta)
y' = x*sin(theta) + y*cos(theta)
%}
CenterEdge = (Center2Holes-(bwHoles/2))/2; % ���S�m�[�h�̉~��
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

