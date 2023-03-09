function [Report,U,p] = Utest(A,B,varargin)
% 
% Utest: マン・ホイットニーのU検定
% 
% [Report,U,p] = Utest(A,B,nop)
%    A        標本1（ベクトル）
%    B        標本2（ベクトル）
%    nop      p値を計算せずU値だけを計算するには1を指定して下さい
%    Report   検定結果
%    U        観測されたU値
%    p        観測されたp値（両側）
% 
% U統計量を用いて、両標本の中央値に差があるかを検定します。検定
% は両側です。引数nopに1が指定されると、p値の計算（長い時間を要
% することがあります）を行わず、U値だけを計算します。この場合、
% pには空配列が出力されます。nopは省略可能です。省略されるとp値
% の計算を行います。
% 
% U値は、標本1の要素が標本2の要素より大きい対の総数（すなわち、
% A(n)>B(n) の総数）と、その逆である対の総数（ B(n)>A(n) の総数）
% のうち、小さい方の値を報告します。両標本間で同値の対（結び）が
% ある場合には、0.5としてU値に計上しています。
% 
% A, Bともに要素数が15以下の場合には、帰無仮説の下で実際にどのよ
% うなU値が観測されうるかをすべて数え上げ、その結果に基づいて、
% 観測されたU値以下のU値が生起する確率を計算します。pはその確率
% の2倍の値です。なお、数え上げの際には与えられた標本に含まれて
% いる結びも考慮して数え上げています。
%
% A, Bどちらかの要素数が15を超える場合は、Uの分布が正規分布で近
% 似できることを利用して検定を行います。
% 
% この関数は、関数 p4ztest.m を利用します。
% 
% see also:   ranksum (Statistics Toolbox)
% 
% (2009/02/10, by R. NIIMI)

Report = ['']; U = []; p = [];

if ~(isvector(A) & isvector(B))
    disp('  [Utest]:  Please input vectors as data.');
    return;
end
nop = 0;
if ~isempty(varargin)
    if ~isempty(varargin{1})
        nop = varargin{1}(1);
    end
end

N = [length(A) length(B)];
Ns = sum(N);
Np = prod(N);

Am = A' * ones(1,N(2));
Bm = ones(N(1),1) * B;

Ueq = sum(sum(Am==Bm));
Ua = sum(sum(Am>Bm)) + Ueq/2;
Ub = Np-Ua;
% Ub = sum(sum(Bm>Am)) + Ueq/2;
U = min([Ua Ub]);

if nop~=1
    if sum(N>15)~=0
        if Ueq==0
            S = sqrt(Np*(Ns+1))/12;
        else
            temp = sort([A B]);
            temp = temp(1:Ns-1)==temp(2:Ns);
            temp = cumsum(temp).*temp;
            temp = [0 temp 0];
            temp = temp(find(temp(2:Ns+1)-temp(1:Ns)<0));
            temp = [0 temp];
            temp = temp(2:end) - temp(1:end-1) + 1;
            temp = sum(temp.^3 - temp);
            S = sqrt( Np/(Ns*(Ns-1)) * (Ns^3-Ns-temp)/12 );
        end
    %     Uaz = (Ua-(Np/2)) / S;
        Uaz = (abs(Ua-(Np/2))-0.5) / S;
        p = p4ztest(Uaz);
    else
        data = [A B];
        ListA = nchoosek([1:Ns],N(1));
        Ncase = size(ListA,1);
        ListB = [];
        UaAll = [];
        for c=1:Ncase
            NowA = ListA(c,:);
            I = sum((ones(N(1),1)*[1:Ns])==(NowA'*ones(1,Ns)));
            NowB = find(I==0);
            NowAm = data(NowA)' * ones(1,N(2));
            NowBm = ones(N(1),1) * data(NowB);
            NowUeq = sum(sum(NowAm==NowBm));
            NowUa = sum(sum(NowAm>NowBm)) + NowUeq/2;
    %        NowUb = sum(sum(NowAm<NowBm)) + NowUeq/2;
            UaAll(c) = NowUa;
    %        UbAll(c) = NowUb;
    %        ListB(c,:) = NowB;
        end
        NUa = [sum(UaAll<=Ua) sum(UaAll>=Ua)];
    %    UbAll = Np - UaAll;
    %    NUb = [sum(UbAll<=Ub) sum(UbAll>=Ub)];
        p = min(NUa)/Ncase;
        p = p*2;
    end
    Report = ['  [Utest]:  U = ' num2str(U) ',  p = ' num2str(p) ' (two-tailed)'];
else
    Report = ['  [Utest]:  U = ' num2str(U)];
end
Report = ['  [Utest]:  U = ' num2str(U)];
if nop~=1
    if sum(N>15)~=0
        Report = strvcat(Report,['            z = ' num2str(Uaz) ',  p = ' num2str(p) ' (two-tailed)']);
        Report = strvcat(Report,['            (Tested with approximation by normal distribution)']);
    else
        Report = [Report ',  p = ' num2str(p) ' (two-tailed)'];
    end
end

