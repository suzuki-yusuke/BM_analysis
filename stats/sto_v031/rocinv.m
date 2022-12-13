function [FA] = rocinv(HR,d,varargin)
% 
% rocinv: ROC曲線上でヒット率に対応する誤警報率の計算
% 
% [FA] = rocinv(HR,d,sdratio,wannaplot)
%    FA         誤警報率
%    HR         ヒット率
%    d          d'値（N分布の標準偏差を単位とする）
%    sdratio    N分布の標準偏差に対するSN分布の標準偏差の比
%    wannaplot  ROC曲線をプロットするには1を指定して下さい
% 
% 与えられたdとsdratioに基づいてROC曲線を決定し、その曲線において
% 与えられたヒット率HRに対応する誤警報率FAを計算します。逆に誤警報
% 率からヒット率を計算するにはrocを用いて下さい。
% 
% HRおよびdにはスカラーか等しい要素数のベクトルを指定します。ただ
% し、HRがベクトルでdがスカラーの場合、HRのすべての要素について
% 同じdの値に基づいて（つまり、単一のROC曲線に基づいて）対応する
% 誤警報率を計算します。HRとdがともにベクトルの場合、HRの第n要素に
% 対応するヒット率は、dの第n要素に基づいて計算されます。
% 
% sdratioはスカラーです。sdratioは省略された場合、1が指定されたも
% のとみなされます（つまり、両分布は等分散であると仮定します）。
% 
% wannaplotに1が指定されると、ROC曲線をグラフにプロットします。省
% 略されるとプロットしません。
% 
% see also:   roc   dprime   fitroc
% 
% (2009/02/23, by R. NIIMI)

FA = [];

if 1-isvector(d)
    disp('  [rocinv]:  d should be a scalar or vector.');
    return;
end
if 1-isvector(HR)
    disp('  [rocinv]:  HR should be a scalar or vector.');
    return;
end
Check = find((HR<0) | (HR>1));
if 1-isempty(Check)
    disp('  [rocinv]:  HR should be 0~1.');
    return;
end

sdratio = 1;
wannaplot = 0;
if length(varargin)>0
    sdratio = varargin{1};
    if length(sdratio)>1
        disp('  [rocinv]:  sdratio should be a scalar.');
        return;
    end
    if sdratio<=0
        disp('  [rocinv]:  sdratio should be larger than zero.');
        return;
    end
    if isempty(sdratio)
        sdratio = 1;
    end
    if length(varargin)>1
        wannaplot = varargin{2}(1);
    end
    if isempty(wannaplot)
        wannaplot = 0;
    end
end

L = length(HR);
if length(d)==1 & L>1
    d = ones(1,L)*d;
elseif length(d)~=L
    disp('  [rocinv]:  Numbers of elements of HR and d are not equal.');
    return;
end

Zsn = erfinv((1-HR)*2-1) * sqrt(2);
Zn = d + (Zsn * sdratio);
FA = 1 - (erf(Zn/sqrt(2))+1)/2;

if wannaplot==1
    plotFA  = [0.001 0.005 ([1:99]/100) 0.995 0.999];
    plotZn  = erfinv((1-plotFA)*2-1) * sqrt(2);
    plotZsn = d'*ones(1,103) - ones(length(d),1)*plotZn;
    plotZsn = plotZsn / sdratio;
    plotHR  = (erf(plotZsn/sqrt(2))+1)/2;
    plotFA = [0 plotFA 1];
    plotHR = [zeros(length(d),1) plotHR ones(length(d),1)];
    plot([0 1],[0 1],'k-'), hold on;
    plot(plotFA,plotHR,'r-'), hold on;
    plot(FA,HR,'b.'), hold on;
    axis square;
    xlim([0 1]); ylim([0 1]);
    xlabel('False Alarm'); ylabel('Hit');
end

