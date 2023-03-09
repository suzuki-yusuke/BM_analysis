function [Report,HSD5,HSD1,Result5,Result1,Md] = tukey_2wwm(Data,Fct)
% 
% tukey_2wwm:  2要因分散分析（2要因とも対応がある場合）での、TukeyのHSD法による多重比較（データは集計表形式）
% 
% [Report,HSD5,HSD1,Result5,Result1,Md] = tukey_2wwm(Data,Fct)
%    Data    対象となるデータ（集計表形式）
%    Fct     検定する要因はどちらかを、下のように指定します。
%             1 -> 要因Aの主効果についての多重比較
%             2 -> 要因Bの主効果についての多重比較
%    Report  検定結果の表（* p < .05, ** p < .01）
%    HSD5    各比較における、有意水準5%でのHSD値
%    HSD1    各比較における、有意水準1%でのHSD値
%    Result5 各比較における、有意水準5%での検定結果（1:有意, 0:n.s.）
%    Result1 各比較における、有意水準1%での検定結果（1:有意, 0:n.s.）
%    Md      平均値の差の表
% 
% 入力データ（Data）は、サンプル数×要因Bの水準数×要因Aの
% 水準数となる、3次元の行列です。同一行内の数値は対応したデ
% ータです（anova_2wwm.mに準じる）。
% 
% この関数は、関数 studragen_5, studrange_1 を利用します。
% 
% see also:   anova_2wwm   sme_2wwm  tukey_2ww
% 
% (2015/09/03, by R. NIIMI)

Report = ['']; HSD5 = []; HSD1 = []; Result5 = []; Result1 = []; Md = [];

if ~((size(Data,1)>1 & size(Data,2)>1) & size(Data,3)>1)
    disp('  [tukey_2wwm]:  Data must be a 3-dimensional matrix.');
    return;
end

%%%%%% ここからanova_2wwm.mと共通コード %%%%%%
N = size(Data,1);
NL = [size(Data,3) size(Data,2)]; % 水準数
Nall = N * NL(1) * NL(2);
for k=1:NL(1)
    for l=1:NL(2)
        M(k,l) = sum(Data(:,l,k));
    end
end

ABS = sum(sum(sum(Data.^2)));
X = sum(sum(sum(Data))).^2 / (N*NL(1)*NL(2));
A = sum(sum(M').^2) / (N*NL(2));
B = sum(sum(M).^2) / (N*NL(1));
AB = sum(sum(M.^2)) / N;
S = sum(sum(sum(Data,2),3).^2) / (NL(1)*NL(2));
AS = sum(sum(sum(Data,2).^2)) / NL(2);
BS = sum(sum(sum(Data,3).^2)) / NL(1);

%%%%%% ここまでanova_2wwm.mと共通コード %%%%%%

if Fct==1
    dfe = (NL(1)-1)*(N-1); % dfas of anova_2wwm
    MSe = (AS-A-S+X)/dfe; % MSas of anova_2wwm
    Level = [1:NL(1)];
    M = M';
elseif Fct==2
    dfe = (NL(2)-1)*(N-1); % dfbs of anova_2wwm
    MSe = (BS-B-S+X)/dfe; % MSbs of anova_2wwm
    Level = [1:NL(2)];
else
    disp('  [tukey_2wwm]:  Fct must be 1 or 2.');
    return;
end

M = mean(M/N);

Md = (ones(NL(Fct),1) * M) - (M' * ones(1,NL(Fct)));
q5 = studrange_5(NL(Fct),dfe);
q1 = studrange_1(NL(Fct),dfe);
HSD5 = sqrt(MSe/(N*NL(3-Fct))) * q5;
HSD1 = sqrt(MSe/(N*NL(3-Fct))) * q1;

Result5 = abs(Md)>HSD5;
Result1 = abs(Md)>HSD1;

NL = NL(Fct);
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

