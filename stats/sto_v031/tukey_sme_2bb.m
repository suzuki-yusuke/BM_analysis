function [Report,HSD5,HSD1,Result5,Result1,Md] = tukey_sme_2bb(Data,Fct,Lev)
% 
% tukey_sme_2bb:  2要因分散分析（2要因とも対応のない場合）での、単純主効果についてのTukeyのHSD法による多重比較
% 
% [Report,HSD5,HSD1,Result5,Result1,Md] = tukey_sme_2bb(Data,Fct,Lev)
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
% ての多重比較をするには、tukey_sme_2bb(Data,2,0) とします。
% 
% 入力データ（Data）は、1列目に要因Aの水準番号、2列目に要因Bの水準
% 番号、そして3列目に測定値が書かれた行列（anova_2bb.mに準じる）。
% 
% この関数は、関数 sme_2bb, studrange_5, studrange_1 を利用します。
% 
% see also:   anova_2bb   sme_2bb   tukey_sme_2bw   tukey_sme_2ww
% 
% (2015/09/03, by R. NIIMI)

Report = ['']; HSD5 = []; HSD1 = []; Result5 = []; Result1 = []; Md = [];

if ~(size(Data,1)>1 & size(Data,2)==3)
    disp('  [tukey_sme_2bb]:  Data must be a N by 3 matrix.');
    return;
end

[Report,Level,M,N,F,p,MS,df] = sme_2bb(Data,Fct);
if isempty(p) % failure on sme_2bb
    return;
end

NL = [length(Level{1}) length(Level{2})];
Ntil = NL(1)*NL(2) / sum(sum(1./N));
MSwc = MS(length(MS));
dfwc = df(length(df));
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
    disp('  [tukey_sme_2bb]:  Fct must be 1 or 2.');
    return;
end
if isempty(LevIndex)
    if Fct==1
        disp('  [tukey_sme_2bb]:  No such level for Factor B.');
    else
        disp('  [tukey_sme_2bb]:  No such level for Factor A.');
    end
    return;
end

M = M(LevIndex,:);
Md = (ones(NL(Fct),1) * M) - (M' * ones(1,NL(Fct)));
q5 = studrange_5(NL(Fct),dfwc);
q1 = studrange_1(NL(Fct),dfwc);
HSD5 = sqrt(MSwc/Ntil) * q5;
HSD1 = sqrt(MSwc/Ntil) * q1;

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

