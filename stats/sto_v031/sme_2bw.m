function [Report,Level,M,N,F,p,MS,df] = sme_2bw(Data,Fct)
% 
% sme_2bw:  2要因分散分析（混合計画の場合）での単純主効果の検定
% 
% [Report,Level,M,N,F,p,MS,df] = sme_2bw(Data,Fct)
%    Data    対象となるデータ（3列の行列）
%    Fct     検定する要因はどちらかを、下のように指定します。
%             1 -> 要因Bの各水準における要因Aの単純主効果の検定
%             2 -> 要因Aの各水準における要因Bの単純主効果の検定
%            （要因Aが対応のない要因、要因Bが対応のある要因です）
%    Report  分析結果の一覧表
%    Level   水準番号の一覧（セルアレイ）
%    M       各条件の平均値
%    N       サンプル数
%    F       F値
%    p       p値
%    MS      平均平方
%    df      自由度
% 
% 入力データ（Data）は、1列目に要因Aの水準番号、2列目に要因Bの
% 水準番号、そして3列目に測定値が書かれた行列（anova_2bw.mに準
% じる）。
% 
% see also:   anova_2bw   sme_2bb   sme_2ww   tukey_2bw
% 
% (2008/12/25, by R. NIIMI)

Report = ['']; Level = []; M = []; N = []; F = []; p = []; MS = []; df = [];

if ~(size(Data,1)>1 & size(Data,2)==3)
    disp('  [sme_2bw]:  Data must be a N by 3 matrix.');
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
    disp('  [anova_2wb]:  Factor B must be a repated-measures factor.');
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

%%%%%% ここまでanova_2bw.mと共通コード %%%%%%

if Fct==1
    SSe = ABS - AB; % SSpool
    dfe = (sum(mean(N')) - NL(1)) + ((sum(mean(N')) - NL(1)) * (NL(2)-1)); % dfpool =  dfs(a) + dfbs(a)
else
    SSe = ABS - AB - AS + A;  % SSbs(a)
    dfe = (sum(mean(N')) - NL(1)) * (NL(2)-1); % dfbs(a)
end
MSe = SSe / dfe;
df = [NL(1)-1 NL(2)-1 dfe];
   
if Fct==2
    M = M';
elseif Fct~=1
    disp('  [sme_2bw]:  Fct must be 1 or 2.');
    return;
end
SS = Ntil * (sum(M.^2) - (sum(M).^2)/NL(Fct));
if Fct==2
    M = M';
end
MS = SS / df(Fct);
F = MS / MSe;
for k=1:NL(3-Fct)
    p(k) = p4F(F(k),df(Fct),dfe);
end
SS = [SS SSe];
MS = [MS MSe];
df = [ones(1,NL(3-Fct))*df(Fct) dfe];

sp = char(ones(NL(3-Fct)+2,3)*32);
if Fct==1
    F1 = 'A';
    F2 = 'b';
    temp = Level{2};
elseif Fct==2
    F1 = 'B';
    F2 = 'a';
    temp = Level{1};
end
Report = 'Source';
for k=1:NL(3-Fct)
    Report = strvcat(Report,[F1 '(' F2 num2str(temp(k)) ')']);
end
Report = strvcat(Report,'Error');
Report = [Report sp strvcat('SS',num2str(SS')) sp strvcat('df',num2str(df'))];
Report = [Report sp strvcat('MS',num2str(MS')) sp strvcat('F',num2str(F'), ' ') sp strvcat('p',num2str(p'), ' ')];
temp = char(ones(1,size(Report,2))*45);
Report = [temp; Report(1,:); temp; Report(2:NL(3-Fct)+1,:); temp; Report(NL(3-Fct)+2,:); temp];
Report = [char(ones(size(Report,1),3)*32) Report];

