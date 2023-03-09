function [Report,HSD5,HSD1,Result5,Result1,Md] = tukey_1wm(Data)
% 
% tukey_1wm: 対応のある1要因の分散分析における、TukeyのHSD法による多重比較（データは集計表形式）
% 
% [Report,HSD5,HSD1,Result5,Result1,Md] = tukey_1wm(Data)
%    Data    対象となるデータ（2列の行列）
%    Report  検定結果の表（* p < .05, ** p < .01）
%    HSD5    有意水準5%でのHSD値
%    HSD1    有意水準1%でのHSD値
%    Result5 各比較における、有意水準5%での検定結果（1:有意, 0:n.s.）
%    Result1 各比較における、有意水準1%での検定結果（1:有意, 0:n.s.）
%    Md      平均値の差の表
% 
% 入力データ（Data）は、anova_1wm.mと同様、サンプル数×水準数の行列です。
% すなわち、同一行内の数値は対応したデータ。
% 4名の参加者による3水準のデータの例：
%       Data = [
%          12.3   18.9   15.5
%          15.1   21.2   16.1
%           9.9   23.5   12.0
%          11.3   20.9   14.7
%       ];
% 上記のようなデータに対して、左列から順に1,2,3...と水準番号が自動的に付さ
% れます。なお、Result5, Result1, Md は、水準数×水準数の行列になります。
% 
% この関数は、関数 anova_1wm, studragen_5, studrange_1 を利用します。
% 
% see also:   tukey_1b   tukey_1w   anova_1wm
% 
% (2015/09/03, by R. NIIMI)

Report = ['']; HSD5 = []; HSD1 = []; Result5 = []; Result1 = []; Md = [];

if ~(size(Data,1)>1 & size(Data,2)>1)
    disp('  [tukey_1wm]:  Data must be a matrix.');
    return;
end

[Report,Level,M,N,F,p,MS,df] = anova_1wm(Data);
if isempty(p) % failure on anova_1wm
    return;
end

NL = length(M); % 水準数
Md = (ones(NL,1) * M) - (M' * ones(1,NL));
q5 = studrange_5(NL,df(3));
q1 = studrange_1(NL,df(3));
HSD5 = sqrt(MS(3)/N) * q5;
HSD1 = sqrt(MS(3)/N) * q1;

Result5 = abs(Md)>HSD5;
Result1 = abs(Md)>HSD1;

Sp = char(ones(NL,3)*32);
Report = [Sp char(ones(NL,1)*91) num2str(Level') char(ones(NL,1)*93)];
Report = [blanks(size(Report,2)); Report];
Sp = [blanks(3); Sp];
for j=1:NL
    NowCol = num2str(Md(:,j));
    temp = size(NowCol,2) - (length(num2str(Level(j)))+2);
    if temp>=0
        NowCol = [blanks(temp) '[' num2str(Level(j)) ']'; NowCol];
    else
        NowCol = ['[' num2str(Level(j)) ']'; [char(ones(NL,-temp)*32) NowCol]]; 
    end
    temp = [0; Result5(:,j) + Result1(:,j)];
    for i=1:NL+1
        Ast(i,:) = [char(ones(1,temp(i))*42)  blanks(2-temp(i))];
    end
    Report = [Report Sp NowCol Ast];
end

