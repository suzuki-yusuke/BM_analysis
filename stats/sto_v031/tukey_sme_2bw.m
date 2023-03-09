function [Report,HSD5,HSD1,Result5,Result1,Md] = tukey_sme_2bw(Data,Fct,Lev)
% 
% tukey_sme_2bw:  2要因分散分析（混合計画の場合）での、単純主効果についてのTukeyのHSD法による多重比較
% 
% [Report,HSD5,HSD1,Result5,Result1,Md] = tukey_sme_2bw(Data,Fct,Lev)
%    Data    対象となるデータ（3列の行列）
%    Fct     どちらの要因の単純主効果を検定するか、下のように指定します。
%             1 -> 要因Aの単純主効果についての多重比較
%             2 -> 要因Bの単純主効果についての多重比較
%            （要因Aが対応のない要因、要因Bが対応のある要因です）
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
% ての多重比較をするには、tukey_sme_2bw(Data,2,0) とします。
% 
% 入力データ（Data）は、1列目に要因Aの水準番号、2列目に要因Bの水準
% 番号、そして3列目に測定値が書かれた行列（anova_2bw.mに準じる）。
% 
% この関数は、関数 studrange_5, studrange_1 を利用します。
% 
% see also:   anova_2bw   sme_2bw   tukey_sme_2bb   tukey_sme_2ww
% 
% (2015/09/03, by R. NIIMI)

Report = ['']; HSD5 = []; HSD1 = []; Result5 = []; Result1 = []; Md = [];

if ~(size(Data,1)>1 & size(Data,2)==3)
    disp('  [tukey_sme_2bw]:  Data must be a N by 3 matrix.');
    return;
end

%%%%%% ここからanova_2bw.mと共通コード %%%%%%

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
Alert = 0;
NL = [length(Level{1}) length(Level{2})]; % 水準数
for k=1:NL(2)
    I = find(Data(:,2)==Level{2}(k));
    if k>1
        if length(I)~=size(Dat,1) % 要因Bの各水準のデータ数が一致しているかチェック
            Alert = 1;
        else
            Dat(:,:,k) = Data(I,:);
            if size(Dat,1)~=sum(Dat(:,1,k-1)==Dat(:,1,k)) % 対応したデータ同士で要因Aの水準が一致しているかチェック
                Alert = 1;
            end
        end
    else
  %      temp = Data(I,1);
        Dat = Data(I,:);
    end        
    for l=1:NL(1)
        I = find(Data(:,1)==Level{1}(l) & Data(:,2)==Level{2}(k));
        N(l,k) = length(I);
        M(l,k) = sum(Data(I,3));
    end
end

if Alert
    disp('  [tukey_sme_2bw]:  Factor B must be a repated-measures factor.');
    return;
end

ABS = sum(Data(:,3).^2);
AS = sum(sum(Dat(:,3,:),3).^2)/NL(2);
AB = sum(sum((M.^2)./N));
A = sum((sum(M').^2) ./ sum(N'));
M = M ./ N; % 各セルの平均値
Gp = sum(sum(M));
Xp = Gp^2 / (NL(1)*NL(2));
Ap = sum(mean(M,2).^2) * NL(2);
Bp = sum(mean(M).^2) * NL(1);
ABp = sum(sum(M.^2));
Ntil = NL(1) / sum(1./mean(N'));

SS = [
    Ntil * (Ap - Xp); % SSa
    AS - A; % SSs(a)
    Ntil * (Bp - Xp); % SSb
    Ntil * (ABp - Ap - Bp + Xp); % SSab
    ABS - AB - AS + A; % SSbs(a)
];

df = [
    NL(1) - 1; % dfa
    sum(mean(N')) - NL(1); % dfs(a)
    NL(2) - 1; % dfb
    (NL(1)-1) * (NL(2)-1); % dfab
    (sum(mean(N')) - NL(1)) * (NL(2)-1); % dfbs(a)
];

MS = SS ./ df;

%%%%%% ここまでanova_2bw.mと共通コード %%%%%%

Report = [];
HSD5 = []; HSD1 = [];
Result5 = []; Result1 = [];
Md = [];

if Fct==1
    M = M';
    LevIndex = find(Level{2}==Lev);
    Level = Level{1};
    dfe = df(2)+df(5); % dfpool
    MSe = (SS(2)+SS(5))/dfe; % MSpool
    q5 = [studrange_5(NL(1),df(2)) studrange_5(NL(1),df(5))];
    q5 = (q5(1)*MS(2) + q5(2)*MS(5)*(NL(3-Fct)-1)) / (MS(2) + MS(5)*(NL(3-Fct)-1));
    q1 = [studrange_1(NL(1),df(2)) studrange_1(NL(1),df(5))];
    q1 = (q1(1)*MS(2) + q1(2)*MS(5)*(NL(3-Fct)-1)) / (MS(2) + MS(5)*(NL(3-Fct)-1));
elseif Fct==2
    LevIndex = find(Level{1}==Lev);
    Level = Level{2};
    dfe = df(5);
    MSe = MS(5);
    q5 = studrange_5(NL(2),dfe);
    q1 = studrange_1(NL(2),dfe);
else
    disp('  [tukey_sme_2bw]:  Fct must be 1 or 2.');
    return;
end
if isempty(LevIndex)
    if Fct==1
        disp('  [tukey_sme_2bw]:  No such level for Factor B.');
    else
        disp('  [tukey_sme_2bw]:  No such level for Factor A.');
    end
    return;
end

M = M(LevIndex,:);
Md = (ones(NL(Fct),1) * M) - (M' * ones(1,NL(Fct)));

HSD5 = sqrt(MSe/Ntil) * q5;
HSD1 = sqrt(MSe/Ntil) * q1;

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

