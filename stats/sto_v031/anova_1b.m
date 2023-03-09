function [Report,Level,M,N,F,p,MS,df,ES] = anova_1b(Data)
% 
% anova_1b: 対応のない1要因分散分析
% 
% [Report,Level,M,N,F,p,MS,df,ES] = anova_1b(Data)
%    Data    対象となるデータ（2列の行列）
%    Report  分散分析表
%    Level   水準番号の一覧
%    M       各水準の平均値
%    N       各水準のサンプル数
%    F       F値
%    p       p値
%    MS      平均平方
%    df      自由度
%    ES      効果量（pEta^2:偏イータ2乗）
%
% 入力データ（Data）は、1列目に水準番号、2列目に測定値が書かれた行列。
% すなわち、1行が1サンプル。行の順序はランダムでかまいません。
% 3水準のデータの例：
%       Data = [
%          1    12.3
%          3    15.1
%          2     9.9
%          3    18.7
%          2    14.0
%          ...  ...
%       ];
% 水準番号は整数（負も可）。MやNは、水準番号の小さい順になっています。
% 
% see also:   anova_1w   anova_1wm
% 
% (2015/09/10, by R. NIIMI)

Report = ['']; Level = []; M = []; N = []; F = []; p = []; MS = []; df = []; ES = [];

if ~(size(Data,1)>1 & size(Data,2)==2)
    disp('  [anova_1b]:  Data must be a N by 2 matrix.');
    return;
end

Nall = size(Data,1);
Level = min(Data(:,1));
for k=min(Data(:,1))+1:max(Data(:,1))
   if sum(Data(:,1)==k) > 0
       Level = [Level k];
   end
end
Dat = [];
NL = length(Level);
for k=1:NL
    I = find(Data(:,1)==Level(k));
    M(k) = mean(Data(I,2));
    N(k) = length(I);
    Dat(1:N(k),k) = Data(I,2);
end

X = (sum(sum(Dat))^2)/Nall;
AS = sum(sum(Dat.^2));
A = sum((sum(Dat).^2)./N);

SS1 = A - X;
SSe = AS - A;
SSt = AS - X;
df1 = NL-1;
dfe = Nall-NL;
dft = Nall-1;
MS1 = SS1 / df1;
MSe = SSe / dfe;
F = MS1 / MSe;
p = p4F(F,df1,dfe);
ES = SS1/SSt;

temp = [
    length('SS') length('df') length('MS') length('F') length('p') length('pEta^2');
    length(num2str(SS1)) length(num2str(df1)) length(num2str(MS1)) length(num2str(F)) length(num2str(p)) length(num2str(ES));
    length(num2str(SSe)) length(num2str(dfe)) length(num2str(MSe)) 0 0 0 ;
    length(num2str(SSt)) length(num2str(dft)) 0 0 0 0 ;
];
temp = [max(temp); max(temp); max(temp); max(temp);] - temp;
sp = '   ';
Report(2,:) = ['   Source' sp 'SS' blanks(temp(1,1)) sp 'df' blanks(temp(1,2)) sp 'MS' blanks(temp(1,3)) sp 'F' blanks(temp(1,4)) sp 'p' blanks(temp(1,5)) sp 'pEta^2' blanks(temp(1,6))];
Report(1,:) = ['   ' char(ones(1,size(Report,2)-3)*45)];
Report(3,:) = Report(1,:);
Report(4,:) = ['   Factor' sp blanks(temp(2,1)) num2str(SS1) sp blanks(temp(2,2)) num2str(df1) sp blanks(temp(2,3)) num2str(MS1) sp blanks(temp(2,4)) num2str(F) sp blanks(temp(2,5)) num2str(p) sp blanks(temp(2,6)) num2str(ES)];
Report(5,:) = ['   Error ' sp blanks(temp(3,1)) num2str(SSe) sp blanks(temp(3,2)) num2str(dfe) sp blanks(temp(3,3)) num2str(MSe) sp blanks(temp(3,4)) sp blanks(temp(3,5)) sp blanks(temp(3,6))];
Report(6,:) = Report(1,:);
Report(7,:) = ['   Total ' sp blanks(temp(4,1)) num2str(SSt) sp blanks(temp(4,2)) num2str(dft) sp blanks(sum(temp(4,3:6))) sp sp sp];
Report(8,:) = Report(1,:);

SS = [SS1 SSe SSt];
df = [df1 dfe dft];
MS = [MS1 MSe];

