clc,clear;
t=1:128;% 数据时间向量
N=length(t);% 数据个数
randn('state',sum(100*clock));% 每次计算给随机数产生设置不同的起点
wn=randn(size(t));% 功率为 1 的高斯白噪声
x=sqrt(20)*sin(2*pi*0.2*t)+sqrt(2)*sin(2*pi*0.213*t)+wn;% 观测数据
%估计自相关函数
R=xcorr(x);
%Pisarenko 谐波分解及恢复
psin=2;% 若已知正弦波个数为 2
%以下求观测数据相关函数矩阵
RPisa=rot90(hankel(R((N+2*psin):-1:N),R(N:-1:(N-2*psin))));
[XPisa,DPisa]=eig(RPisa);
aPisa=XPisa(:,2*psin+1);
rPisa=roots(aPisa);
fPisa=atan(abs(imag(rPisa)./real(rPisa)))/(2*pi)
%ARMA 建模法谐波恢复
pe=6;
M=10;%M>>p
%以下求增广矩阵
RARMA=R((N+pe+1):-1:(N+1));
for i=2:M
RARMA=[RARMA;R((N+pe+i):-1:(N+i))];
end
%以下对 RARMA 做奇异值分解 SVD
[UA SA VA]=svd(RARMA); 
19
%奇异值归一化，并确定有效秩 p
SA1=SA/SA(1,1);
pARMA=0;
i=pe;
while pARMA==0
if SA1(i,i)>=0.01
pARMA=i;
else
i=i-2;
end
end
%计算 Sp
SpARMA=zeros(pARMA+1,pARMA+1);
for i=1:pARMA
for k=1:(pe-pARMA+1)
SpARMA=SpARMA+SA(i,i)*V A(k:k+pARMA,i)*(V A(k:k+pARMA,i))';
end
end
%AR 参数的估计值为：
SpnARMA=SpARMA^(-1);
aARMA=SpnARMA(:,1)/SpnARMA(1,1);
rARMA=roots(aARMA);
fARMA=atan(abs(imag(rARMA)./real(rARMA)))/(2*pi)
%MUSIC 方法实现谐波恢复
%以下求观测数据相关函数矩阵
RMUSIC=rot90(hankel(R((2*N-1):-1:N),R(N:-1:1)));
[XMUSIC,DMUSIC]=eig(RMUSIC);
%以下寻找谐波个数
pm=0;%谐波个数
k=1;
while pm==0
if (DMUSIC(k,k)-DMUSIC(k+1,k+1))/DMUSIC(k+1,k+1)>=1.2
pm=k;
else
k=k+1;
end
end
display(' 谐波个数为： ')
pm=4
S=XMUSIC(:,1:pm);
SHconj=S';
Ps=S*SHconj;
NMUSIC=256;% 分段数量
wiMUSIC=2*pi/NMUSIC;% 频率搜索步长
ww=1:NMUSIC;
f=ww/NMUSIC;% 频率向量
%以下构造 A 矩阵
NN=t'*ww;
awi=exp(j*wiMUSIC*NN);
awic=awi';
fMUSIC=0;
PWxMUSIC=zeros(1,NMUSIC);
for k=1:NMUSIC
PWxMUSIC(k)=1/(awic(k,:)*(eye(N,N)-Ps)*awi(:,k));
%搜索极大值
if k>=2&&k<NMUSIC/2
if PWxMUSIC(k)>fMUSIC
fMUSIC=k/NMUSIC;
end
end 
20
end
display('M USIC 方法估计的频率： ')
fMUSIC
figure(1);
stem(f,abs(PWxMUSIC),'filled');
title('MUSIC 方法得到的谱 ');
xlabel('f');
ylabel('PMUSIC');
%ESPRIT 方法实现谐波恢复
%第一步构造自相关矩阵 Rxx 和 Rxy
mESPRIT=4;
Rxx=rot90(hankel(R(N+mESPRIT-1:-1:N),conj(R(N:N+mESPRIT-1))));
Rxy=rot90(hankel(R(N-mESPRIT:N-1),R(N-1:N+mESPRIT-2)),3);
[XESPRIT,DESPRIT]=eig(Rxx);
Q2ESPRIT=min(diag(DESPRIT));
Cxx=Rxx-Q2ESPRIT*eye(size(Rxx));
Cxy=Rxy-Q2ESPRIT*diag(ones(1,mESPRIT-1),-1);
%LS-ESPRIT 算法
[XESPRIT,DESPRIT]=eig(Cxx,Cxy);
display('LS-ESPRIT 方法估计的频率： ')
fESPRIT_LS=atan(abs(imag(diag(DESPRIT))./real(diag(DESPRIT))))/(2*pi)
%TLS-ESPRIT 算法
[UESPRIT,SESPRIT,VESPRIT]=svd(Cxx);
%奇异值归一化，并确定有效秩 p
SESPRIT1=SESPRIT/SESPRIT(1,1);
pESPRIT=0;
k=mESPRIT;
while pESPRIT==0
if abs(SESPRIT1(k,k))>=0.002
pESPRIT=k;
else
k=k-1;
end
end
UESPRITp=UESPRIT(:,1:pESPRIT);
SESPRITp=SESPRIT(1:pESPRIT,1:pESPRIT);
VESPRITp=VESPRIT(:,1:pESPRIT);
Cxyp=UESPRITp'*Cxy*VESPRITp;
[XESPRITp,DESPRITp]=eig(SESPRITp,Cxyp);
display('T LS-ESPRIT 方法估计的频率： ')
fESPRIT_TLS=atan(abs(imag(diag(DESPRITp))./real(diag(DESPRITp))))/(2*pi) 