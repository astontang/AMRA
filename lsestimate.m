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
%用一般最小二乘法估计功率谱和谐波频率
pls=4;%AR 模型的阶数取 4
for qels=5:8;% 选取修正 Yule-Walker 方程的起点
%以下求相关函数估计值组成的 Hankel 矩阵
RLS4=R(qels:(qels+pls-1));
for k=1:(pls-1)
RLS4=[R((qels-k):(qels-k+pls-1));RLS4];
end
bls4=-R((qels+1):(qels+pls));%Y ule-Walker 方程右边的值
als4(:,(qels-pls))=(RLS4'*RLS4)^(-1)*RLS4'*(bls4)';% 最小二乘估计公式
end
display('AR 阶数为 4，qe 取 5,6,7,8 时，用一般最小二乘法估计的 AR 参数 :')
als4=als4(pls:-1:1,:)
%采用 Cadzow 谱估计子进行功率谱估计
als4=[ones(1,4);als4];
als4=als4';
CCadzowp4=CORP(N:N+pls);
for i=1:pls
CCadzowp4=[CCadzowp4;CORP(N-i:N-i+pls)];
end
nCadzowp4=als4*CCadzowp4;
Nzp4=zeros(4,N);
Azp4=zeros(4,N);
Nnzp4=zeros(4,N);
Anzp4=zeros(4,N);
for k=0:pls
Nzp4=Nzp4+nCadzowp4(:,k+1)*(z.^(-k));
Azp4=Azp4+als4(:,k+1)*(z.^(-k));
Nnzp4=Nnzp4+nCadzowp4(:,k+1)*(z.^k);
Anzp4=Anzp4+als4(:,k+1)*(z.^k);
end
PWxCadzowp4=Nzp4./Azp4+Nnzp4./Anzp4;
PWxCadzowp4=abs(PWxCadzowp4);
figure(1);
for i=1:4
subplot(4,1,i);
stem(f,PWxCadzowp4(i,:),'filled');
flsp4(i)=find(PWxCadzowp4(i,1:(N/2-1))==max(PWxCadzowp4(i,1:(N/2-1))))/N;
end
display('AR 阶数为 6，qe 取 7,6,7,8 时，用 Cadzow 谱估计子估计频率 :') 
15
flsp4
pls=6;%AR 模型的阶数取 6
for qels=7:10;% 选取修正 Yule-Walker 方程的起点
%以下求相关函数估计值组成的 Hankel 矩阵
RLS6=R(qels:(qels+pls-1));
for k=1:(pls-1)
RLS6=[R((qels-k):(qels-k+pls-1));RLS6];
end
bls6=-R((qels+1):(qels+pls));%Y ule-Walker 方程右边的值
als6(:,(qels-pls))=(RLS6'*RLS6)^(-1)*RLS6'*(bls6)';% 最小二乘估计公式
end
display('AR 阶数为 6，qe 取 7,8，9，10 时，用一般最小二乘法估计的 AR 参数 :')
als6=als6(pls:-1:1,:)
%采用 Cadzow 谱估计子进行功率谱估计
als6=[ones(1,4);als6];
als6=als6';
CCadzowp6=CORP(N:N+pls);
for i=1:pls
CCadzowp6=[CCadzowp6;CORP(N-i:N-i+pls)];
end
nCadzowp6=als6*CCadzowp6;
Nzp6=zeros(4,N);
Azp6=zeros(4,N);
Nnzp6=zeros(4,N);
Anzp6=zeros(4,N);
for k=0:pls
Nzp6=Nzp6+nCadzowp6(:,k+1)*(z.^(-k));
Azp6=Azp6+als6(:,k+1)*(z.^(-k));
Nnzp6=Nnzp6+nCadzowp6(:,k+1)*(z.^k);
Anzp6=Anzp6+als6(:,k+1)*(z.^k);
end
PWxCadzowp6=Nzp6./Azp6+Nnzp6./Anzp6;
PWxCadzowp6=abs(PWxCadzowp6);
figure(2);
for i=1:4
subplot(4,1,i);
stem(f,PWxCadzowp6(i,:),'filled');
flsp6(i)=find(PWxCadzowp6(i,1:(N/2-1))==max(PWxCadzowp6(i,1:(N/2-1))))/N;
end
display('AR 阶数为 6，qe 取 7,8，9，10 时，用 Cadzow 谱估计子估计频率 :')
flsp6 