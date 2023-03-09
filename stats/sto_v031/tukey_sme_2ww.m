function [Report,HSD5,HSD1,Result5,Result1,Md] = tukey_sme_2ww(Data,Fct,Lev)
% 
% tukey_sme_2ww:  2要因分散分析（2要因とも対応のある場合）での、単純主効果についてのTukeyのHSD法による多重比較（データはリスト形式）
% 
% [Report,HSD5,HSD1,Result5,Result1,Md] = tukey_sme_2ww(Data,Fct,Lev)
%    Data    対象となるデータ（3列の行列）
%    Fct     どちらの要因の単純主効果を検定するか、下のように指定します。
%             1 -> 要因Aの単純主効果についての多重比較
%             2 -> 要因Bの単純主効果についての多重比較
%    Lev     もう一方の要因のどの水準における検定かを、水準番号で指定します。
%    Report  検定結果の表（* p < .05, ** p < .01）
%    HSD5    各比較における、有意水準5%でのHSD値
%    HSD1    各比較における、有意水準1%でのHSD値
%    Result5 各比較における、有意水準5%での検定結果（1:有意, 0:n.s.）
%    Result1 各比較における、有意水準1%での検定結果（1:有意, 0:n.s.）
%    Md      平均値の差の表
% 
% 例：要因Aが2水準（水準番号は [1 0] ）、要因Bが3水準（水準番号は 
%  [1 2 3] ）の場合に、要因Aの水準0における要因Bの単純主効果につい
% ての多重比較をするには、tukey_sme_2ww(Data,2,0) とします。
% 
% 入力データ（Data）は、1列目に要因Aの水準番号、2列目に要因Bの水準
% 番号、そして3列目に測定値が書かれた行列（anova_2ww.mに準じる）。
% 
% この関数は、関数 studrange_5, studrange_1 を利用します。
% 
% see also:   anova_2ww   sme_2ww   tukey_sme_2bb   tukey_sme_2wwm
% 
% (2015/09/03, by R. NIIMI)

Report = ['']; HSD5 = []; HSD1 = []; Result5 = []; Result1 = []; Md = [];

if ~(size(Data,1)>1 & size(Data,2)==3)
    disp('  [tukey_sme_2ww]:  Data must be a N by 3 matrix.');
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

SS = [
    S - X; % SSs
    A - X; % SSa
    AS - A - S + X; % SSas
    B - X; % SSb
    BS - B - S + X; % SSbs
    AB - A - B + X; % SSab
    ABS - AB - AS - BS + A + B + S - X; % SSabs
    ABS - X; % SSt
];

dfs = N - 1;
dfa = NL(1) - 1;
dfas = dfa * dfs;
dfb = NL(2) - 1;
dfbs = dfb * dfs;
dfab = dfa * dfb;
dfabs = dfas * dfb;
dft = Nall - 1;
df = [dfs dfa dfas dfb dfbs dfab dfabs dft]';

%%%%%% ここまでanova_2ww.mと共通コード %%%%%%

M = M / N;
Report = [];
HSD5 = []; HSD1 = [];
Result5 = []; Result1 = [];
Md = [];

if Fct==1
    M = M';
    LevIndex = find(Level{2}==Lev);
    Level = Level{1};
elseif Fct==2
    LevIndex = find(Level{1}==Lev);
    Level = Level{2};
else 
    disp('  [tukey_sme_2ww]:  Fct must be 1 or 2.');
    return;
end
if isempty(LevIndex)
    if Fct==1
        disp('  [tukey_sme_2ww]:  No such level for Factor B.');
    else
        disp('  [tukey_sme_2ww]:  No such level for Factor A.');
    end
    return;
end

MS = SS./df;
MS = [MS(Fct*2) MS(Fct*2+1) MS(7)];
df = [df(Fct*2) df(Fct*2+1) df(7)];
dfpool = df(2)+df(3);
MSpool = (SS(Fct*2+1)+SS(7)) / dfpool;

M = M(LevIndex,:);
Md = (ones(NL(Fct),1) * M) - (M' * ones(1,NL(Fct)));

q5 = [studrange_5(NL(Fct),df(2)) studrange_5(NL(Fct),df(3))];
q5 = (q5(1)*MS(2) + q5(2)*MS(3)*(NL(3-Fct)-1)) / (MS(2) + MS(3)*(NL(3-Fct)-1));
q1 = [studrange_1(NL(Fct),df(2)) studrange_1(NL(Fct),df(3))];
q1 = (q1(1)*MS(2) + q1(2)*MS(3)*(NL(3-Fct)-1)) / (MS(2) + MS(3)*(NL(3-Fct)-1));
HSD5 = sqrt(MSpool/N) * q5;
HSD1 = sqrt(MSpool/N) * q1;

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
Header = ['ABba'];
Header = [blanks(3) Header(Fct) '(' Header(Fct+2) num2str(Lev) ')'];
Report = strvcat(Header,Report);

