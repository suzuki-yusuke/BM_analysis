function [Report,HSD5,HSD1,Result5,Result1,Md] = tukey_2bb(Data,Fct)
% 
% tukey_2bb:  2要因分散分析（2要因とも対応のない場合）での、TukeyのHSD法による多重比較
% 
% [Report,HSD5,HSD1,Result5,Result1,Md] = tukey_2bb(Data,Fct)
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
% 水準番号、そして3列目に測定値が書かれた行列（anova_2bb.mに準
% じる）。
% 
% 多重比較の対象となる平均値は、各セルごとのサンプル数が異なって
% いてもそれに応じた重みづけをしていないもの（非加重平均）です。
% 例えば要因A（3水準）×要因B（2水準）の計画（6セル）で、要因Aの
% 効果についての多重比較をする場合、要因Aの第1水準の平均値は、セ
% ルごとに算出した6つの平均値のうち要因Aが第1水準である3つのセル
% の値を単純に平均したものです。
% 
% この関数は、関数 studragen_5, studrange_1 を利用します。
% 
% see also:   anova_2bb   sme_2bb   tukey_2ww   tukey_2bw
% 
% (2015/09/03, by R. NIIMI)

Report = ['']; HSD5 = []; HSD1 = []; Result5 = []; Result1 = []; Md = [];

if ~(size(Data,1)>1 & size(Data,2)==3)
    disp('  [tukey_2bb]:  Data must be a N by 3 matrix.');
    return;
end

%%%%%% ここからanova_2bb.mと共通コード %%%%%%
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
for k=1:NL(1)
    for l=1:NL(2)
        I = find(Data(:,1)==Level{1}(k) & Data(:,2)==Level{2}(l));
        N(k,l) = length(I);
        M(k,l) = sum(Data(I,3));
%        Dat(:,l,k) = Data(I,3); % データを3次元（S*B*A）の集計表形式にする
    end
end

ABS = sum(Data(:,3).^2);
AB = sum(sum((M.^2)./N));
M = M./N; % 各セルの平均値

%%%%%% ここまでanova_2bb.mと共通コード %%%%%%

Ntil = NL(1)*NL(2) / sum(sum(1./N));

dfa = NL(1) - 1;
dfb = NL(2) - 1;
dfwc = Nall - NL(1)*NL(2);
df = [dfa dfb dfwc]';
SSwc = ABS - AB;
MSwc = SSwc / dfwc;

if Fct==1
    M = M';
    Level = Level{1};
elseif Fct~=2
    disp('  [tukey_2bb]:  Fct must be 1 or 2.');
    return;
else
    Level = Level{2};
end

M = mean(M); % 非加重平均
Md = (ones(NL(Fct),1) * M) - (M' * ones(1,NL(Fct)));
q5 = studrange_5(NL(Fct),dfwc);
q1 = studrange_1(NL(Fct),dfwc);
HSD5 = sqrt(MSwc/(Ntil*NL(3-Fct))) * q5;
HSD1 = sqrt(MSwc/(Ntil*NL(3-Fct))) * q1;

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

