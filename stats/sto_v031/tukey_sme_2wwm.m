function [Report,HSD5,HSD1,Result5,Result1,Md] = tukey_sme_2wwm(Data,Fct,Lev)
% 
% tukey_sme_2wwm:  2要因分散分析（2要因とも対応のある場合）での、単純主効果についてのTukeyのHSD法による多重比較（データは集計表形式）
% 
% [Report,HSD5,HSD1,Result5,Result1,Md] = tukey_sme_2wwm(Data,Fct,Lev)
%    Data    対象となるデータ（集計表形式）
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
% 例：要因Aが2水準（水準番号は [1 2] となります）、要因Bが3水準
% （水準番号は [1 2 3] となります）の場合に、要因Aの第1水準にお
% ける要因Bの単純主効果についての多重比較をするには、
%  tukey_sme_2wwm(Data,2,1) とします。
% 
% 入力データ（Data）は、サンプル数×要因Bの水準数×要因Aの水準数
% となる、3次元の行列です。同一行内の数値は対応したデータです（
% anova_2wwm.mに準じる）。
% 
% この関数は、関数 studrange_5, studrange_1 を利用します。
% 
% see also:   anova_2wwm   sme_2wwm   tukey_sme_2ww
% 
% (2015/09/03, by R. NIIMI)

Report = ['']; HSD5 = []; HSD1 = []; Result5 = []; Result1 = []; Md = [];

if ~((size(Data,1)>1 & size(Data,2)>1) & size(Data,3)>1)
    disp('  [tukey_sme_2wwm]:  Data must be a 3-dimensional matrix.');
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

%%%%%% ここまでanova_2wwm.mと共通コード %%%%%%

M = M / N;
Report = [];
HSD5 = []; HSD1 = [];
Result5 = []; Result1 = [];
Md = [];

if Fct==1
    M = M';
elseif Fct~=2
    disp('  [tukey_sme_2wwm]:  Fct must be 1 or 2.');
    return;
end
if sum(Lev==[1:NL(3-Fct)])==0
    if Fct==1
        disp('  [tukey_sme_2wwm]:  No such level for Factor B.');
    else
        disp('  [tukey_sme_2wwm]:  No such level for Factor A.');
    end
    return;
end

MS = SS./df;
MS = [MS(Fct*2) MS(Fct*2+1) MS(7)];
df = [df(Fct*2) df(Fct*2+1) df(7)];
dfpool = df(2)+df(3);
MSpool = (SS(Fct*2+1)+SS(7)) / dfpool;
Level = [1:NL(Fct)];

M = M(Lev,:);
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

