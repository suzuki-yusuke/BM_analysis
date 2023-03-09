function [Report,HSD5,HSD1,Result5,Result1,Md] = tukey_2ww(Data,Fct)
% 
% tukey_2ww:  2要因分散分析（2要因とも対応がある場合）での、TukeyのHSD法による多重比較（データはリスト形式）
% 
% [Report,HSD5,HSD1,Result5,Result1,Md] = tukey_2ww(Data,Fct)
%    Data    対象となるデータ（3列の行列）
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
% 入力データ（Data）は、1列目に要因Aの水準番号、2列目に要因Bの
% 水準番号、そして3列目に測定値が書かれた行列（anova_2ww.mに準
% じる）。
% 
% この関数は、関数 studragen_5, studrange_1 を利用します。
% 
% see also:   anova_2ww   sme_2ww   tukey_2wwm
% 
% (2015/09/03, by R. NIIMI)

Report = ['']; HSD5 = []; HSD1 = []; Result5 = []; Result1 = []; Md = [];

if ~(size(Data,1)>1 & size(Data,2)==3)
    disp('  [tukey_2ww]:  Data must be a N by 3 matrix.');
    return;
end

%%%%%% ここからanova_2ww.mと共通コード %%%%%%
Nall = size(Data,1);
Level{1} = min(Data(:,1));
for k=min(Data(:,1))+1:max(Data(:,1))
   if sum(Data(:,1)==k) > 0
       Level{1} = [Level{1} k];
   end
end
Level{2} = min(Data(:,2));
for k=min(Data(:,2))+1:max(Data(:,2))
   if sum(Data(:,2)==k) > 0
       Level{2} = [Level{2} k];
   end
end
Dat = [];
NL = [length(Level{1}) length(Level{2})]; % 水準数
N = Nall / (NL(1) * NL(2));
for k=1:NL(1)
    for l=1:NL(2)
        I = find(Data(:,1)==Level{1}(k) & Data(:,2)==Level{2}(l));
        M(k,l) = sum(Data(I,3));
        Dat(:,l,k) = Data(I,3); % データを3次元（S*B*A）の集計表形式にする
    end
end

ABS = sum(sum(sum(Dat.^2)));
X = sum(sum(sum(Dat))).^2 / (N*NL(1)*NL(2));
A = sum(sum(M').^2) / (N*NL(2));
B = sum(sum(M).^2) / (N*NL(1));
AB = sum(sum(M.^2)) / N;
S = sum(sum(sum(Dat,2),3).^2) / (NL(1)*NL(2));
AS = sum(sum(sum(Dat,2).^2)) / NL(2);
BS = sum(sum(sum(Dat,3).^2)) / NL(1);

%%%%%% ここまでanova_2ww.mと共通コード %%%%%%

if Fct==1
    dfe = (NL(1)-1)*(N-1); % dfas of anova_2ww
    MSe = (AS-A-S+X)/dfe; % MSas of anova_2ww
    Level = Level{1};
    M = M';
elseif Fct==2
    dfe = (NL(2)-1)*(N-1); % dfbs of anova_2ww
    MSe = (BS-B-S+X)/dfe; % MSbs of anova_2ww
    Level = Level{2};
else
    disp('  [tukey_2ww]:  Fct must be 1 or 2.');
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

