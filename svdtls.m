clc,clear;
t=1:128;% 数据时间向量
N=length(t);% 数据个数
f=t/N;% 频率向量
w=2*pi*f;% 角频率向量
z=exp(j*w);
randn('state',sum(100*clock));% 每次计算给随机数产生设置不同的起点
wn=randn(size(t));% 功率为 1 的高斯白噪声
x=sqrt(20)*sin(2*pi*0.2*t)+sqrt(2)*sin(2*pi*0.213*t)+wn;% 观测数据
%估计自相关函数
R=xcorr(x);
%估计协方差函数
COR=xcov(x);
CORP=COR;
CORP(N)=COR(N)/2;
%用 SVD-TLS 法估计功率谱和谐波频率
peTLS=6;
qeTLS=10;%选取修正 Yule-Walker 方程的起点
%以下求增广矩阵
RTLS=R((N+qeTLS+1):-1:(N+qeTLS+1-peTLS));
for k=2:(peTLS+2)
RTLS=[RTLS;R((N+qeTLS+k):-1:(N+qeTLS+k-peTLS))];
end
%以下对 RTLS 做奇异值分解 SVD
[UTLS STLS VTLS]=svd(RTLS);
%奇异值归一化，并确定有效秩 p
STLS1=STLS/STLS(1,1);
pTLS=0;
k=peTLS;
while pTLS==0
if STLS1(k,k)>=0.005
pTLS=k;
else
k=k-1;
end
end
display('SVD-TLS 法确定的 AR 阶数为： ')
pTLS
%计算 Sp
SpTLS=zeros(pTLS+1,pTLS+1);
for i=1:pTLS
for k=1:(peTLS+1-pTLS)
SpTLS=SpTLS+STLS(i,i)*VTLS(k:k+pTLS,i)*(VTLS(k:k+pTLS,i))';
end
end
%AR 参数的估计值为：
SpnTLS=SpTLS^(-1);
aTLS=SpnTLS(:,1)/SpnTLS(1,1);
display('SVD-TLS 法确定的 AR 参数为： ')
aTLS
%采用 Cadzow 谱估计子进行功率谱估计
CCadzow=CORP(N:N+pTLS);
for i=1:pTLS
CCadzow=[CCadzow;CORP(N-i:N-i+pTLS)];
end
nCadzow=aTLS'*CCadzow;
Nz=zeros(1,N);
Az=zeros(1,N);
Nnz=zeros(1,N);
Anz=zeros(1,N);
for k=0:pTLS
Nz=Nz+nCadzow(k+1)*(z.^(-k));
Az=Az+aTLS(k+1)*(z.^(-k));
Nnz=Nnz+nCadzow(k+1)*(z.^k);
Anz=Anz+aTLS(k+1)*(z.^k);
end
PWxCadzow=Nz./Az+Nnz./Anz;
PWxCadzow=abs(PWxCadzow);
subplot(3,1,1);
stem(f,PWxCadzow,'filled');
title('Cadzow 谱估计子估计的功率谱 ');
%通过功率谱估计频率
display('Cadzow 谱估计子估计的频率： ')
fCadzow=find(PWxCadzow(1:(N/2-1))==max(PWxCadzow(1:(N/2-1))))/N
%MA 阶数和参数的辨识
Q=qeTLS;%Q=qe>q 
%以下求增广矩阵
RTLStemp=R((N+Q):-1:(N+Q-pTLS));
for i=1:(pTLS+5) % 使 RTLStemp 为超定矩阵
RTLStemp=[RTLStemp;R((N+Q+i):-1:(N+Q+i-pTLS))];
end
[UTLStemp STLStemp VTLStemp]=svd(R TLStemp);
OQ1=STLStemp(pTLS+1,pTLS+1);
flag=0;
while flag<=0.3
%以下求增广矩阵
Q=Q-1;
RTLStemp=R((N+Q):-1:(N+Q-pTLS));
for i=1:(pTLS+5) % 使 RTLStemp 为超定矩阵
RTLStemp=[RTLStemp;R((N+Q+i):-1:(N+Q+i-pTLS))];
end
[UTLStemp STLStemp VTLStemp]=svd(R TLStemp);
OQ2=STLStemp(pTLS+1,pTLS+1);
flag=abs((OQ1-OQ2)/OQ2);
qTLS=Q;
end
display(' zhang方法确定的 MA 阶数为： ')
qTLS
%Kaveh 估计子系数和 Newton-Raphson 算法估计 MA 系数
ck=zeros(qTLS+1,1);
for k=1:(qTLS+1)
for m=1:pTLS+1
for n=1:pTLS+1
ck(k)=ck(k)+aTLS(m)*conj(aTLS(n))*COR(N+k-1+n-m);
end
end
end
%Newton-Raphson 算法
bNewton=zeros(qTLS+1,1);
%以下是 MA 系数的初始值
biNewton=[ck(1)^(1/2);zeros(qTLS,1)];
while min(abs(bNewton))==0
%以下是第 i 次的拟合误差函数
I=find(abs(bNewton)~=0);
biNewton(I)=0;
biNewton=biNewton+bNewton;
fki=-ck;
for k=1:qTLS+1
for m=1:(qTLS+2-k)
fki(k)=fki(k)+biNewton(m)*biNewton(m+k-1);
end
end
Fi=hankel(biNewton)+rot90(hankel(biNewton(qTLS+1:-1:1)),3);
%以下是第 i+1 次跌代的 MA 系数
bi1=biNewton-Fi^(-1)*fki;
bNewton=bi1;
%判断第 i+1 次跌代的 MA 系数的收敛性
for k=1:qTLS+1
if abs((bi1(k)-biNewton(k))/bi1(k))<=0.05
bNewton(k)=bi1(k);
end
end
biNewton=bi1;
end
display('Newton-Raphson 算法确定的 MA 参数为： ')
bNewton
%用 ARMA 模型和 Kaveh 估计子进行功率谱估计
Bz=zeros(1,N);
Cz=zeros(1,N);
Bnz=zeros(1,N);
Cnz=zeros(1,N);
for k=0:qTLS
Bz=Bz+bNewton(k+1)*(z.^(-k));
Cz=Cz+ck(k+1)*(z.^(-k));
Bnz=Bnz+bNewton(k+1)*(z.^k);
Cnz=Cnz+ck(k+1)*(z.^k);
end
Cz(1)=Cz(1)/2;
Cnz(1)=Cnz(1)/2;
%ARMA 模型功率谱估计
PWxzARMA=(Bz.*Bnz)./(Az.*Anz);
PWx=abs(PWxzARMA);
subplot(3,1,2);
stem(f,PWx,'filled');
title('ARMA 谱估计的功率谱 ');
%通过功率谱估计频率
display('ARMA 谱估计的频率： ')
fARMAP=find(PWx(1:(N/2-1))==max(PWx(1:(N/2-1))))/N
%Kaveh 估计子功率谱估计
PWxzKaveh=(Cz+Cnz)./(Az.^2);
PWx=abs(PWxzKaveh);
subplot(3,1,3);
stem(f,PWx,'filled');
title('Kaveh 谱估计子估计的功率谱 ');
%通过功率谱估计频率
display('Kaveh 谱估计子估计的频率： ')
fKaveh=find(PWx(1:(N/2-1))==max(PWx(1:(N/2-1))))/N 