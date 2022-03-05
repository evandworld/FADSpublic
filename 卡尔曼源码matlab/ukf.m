
clc;
clear;
T=1;%扫描周期
N=60/T;% 总的采样次数
X=zeros(4,N);% 目标真实位置、速度
X(:,1)=[-100,2,200,20];
Z=zeros(1,N);
delta_w=1e-4;
Q=delta_w*diag([0.5,1]);    %过程噪音均值
G=[T^2/2,0;T,0;0,T^2/2;0,T];
R=1;     %观测噪音方差
F=[1,T,0,0;0,1,0,0;0,0,1,T;0,0,0,1];
x0=200;
y0=300;
Xstation=[x0,y0];
v=sqrtm(R)*randn(1,N);
for t=2:N;
    X(:,t)=F*X(:,t-1)+G*sqrtm(Q)*randn(2,1);%%实际值
end
for t=1:N
    if length(Xstation)<=2
        d=sqrt((X(1,t)-Xstation(1))^2+(X(3,t)-Xstation(2))^2);
    else
        d=sqrt((X(1,t)-Xstation(1))^2+(X(3,t)-Xstation(3))^2);
    end
    Z(t)=d+v(t);      %%观测值
end
%%%%%%%%%%%%%%%%%
%UKF 滤波,UT变换
L=4;    %矩阵X的列数
alpha=1;
kalpha=0;
belta=2;
ramda=3-L;  %-1
for j=1:2*L+1    %j=1/2/3/4.../9  ，由比例修正对称采样策略，可得到2*L+1
                 %%wm和wc是随机的权重，加在一起后为3/2
    Wm(j)=1/(2*(L+ramda));    %均值为1/6
    Wc(j)=1/(2*(L+ramda));
end
Wm(1)=ramda/(L+ramda);  %-1/3
Wc(1)=ramda/(L+ramda)+1-alpha^2+belta; %权值计算
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Xukf=zeros(4,N);
Xukf(:,1)=X(:,1);%无迹Kalman滤波状态初始化
PO=eye(4);%协方差阵初始化
for t=2:N
    xestimate= Xukf(:,t-1);     %取上一轮滤波后的值
    P=PO;
    %第一步:获得一组 Sigma点集
    cho=(chol(P*(L+ramda)))';  %chol() 用于对矩阵进行Cholesky分解
                               %假定X是正定的R=chol(X)则产生一个上三角阵R
                               %使R'R=X.
    for k=1:L
        xgamaP1(:,k)=xestimate+cho(:,k);
        xgamaP2(:,k)=xestimate+cho(:,k);
    end
    Xsigma=[xestimate,xgamaP1,xgamaP2]; %Sigma点集（4行，3+3+3=9行的矩阵）
    %第二步:对 Sigma点集进行一步预测
    Xsigmapre=F*Xsigma; %矩阵4*9
    %第三步:利用第二步的结果计算均值和协方差
    Xpred=zeros(4,1); %均值
    for k=1:2*L+1
        Xpred=Xpred+Wm(k)*Xsigmapre(:,k);
    end
    Ppred=zeros(4,4);%协方差阵预测
    for k=1:2*L+1
        Ppred=Ppred+Wc(k)*(Xsigmapre(:,k)-Xpred)*(Xsigmapre(:,k)-Xpred)';
    end
    Ppred=Ppred+G*Q*G';
    %第四步:根据预测值，再一次使用UT变换，得到新的sigma点集
    chor=(chol((L+ramda)*Ppred))';
    for k=1:L
     XaugsigmaP1(:,k)=Xpred+chor(:,k);
     XaugsigmaP2(:,k)=Xpred-chor(:,k);
    end
    Xaugsigma=[Xpred XaugsigmaP1 XaugsigmaP2];
    %第五步:观测预测
    for k=1:2*L+1
        Zsigmapre(1,k)=sqrt((Xaugsigma(1,k)-Xaugsigma(1,k))^2+(Xaugsigma(3,k)-Xaugsigma(2,k))^2);
    end
    %第六步:计算观测预测均值和协方差
    %观测预测的均值
    Zpred=0;
    for k=1:2*L+1
        Zpred=Zpred+Wm(k)*Zsigmapre(1,k);
    end
    Pzz=0;
    for k=1:2*L+1
        Pzz=Pzz+Wc(k)*(Zsigmapre(1,k)-Zpred)*(Zsigmapre(1,k)-Zpred)';
    end
    Pzz=Pzz+R;%得到协方差Pzz
    Pxz=zeros(4,1);
    for k=1:2*L+1
        Pxz=Pxz+Wc(k)*(Xaugsigma(:,k)-Xpred)*(Zsigmapre(1,k)-Zpred)';
    end
    %第七步:计算Kalman增益
    %Kalman增益
    K=Pxz*inv(Pzz);
    %第八步:状态和方差更新
    xestimate=Xpred+K*(Z(t)-Zpred);%状态更新
    P=Ppred-K*Pzz*K';  %方差更新
    PO=P;
    Xukf(:,t)=xestimate;
end
%误差分析
for i=1:N
    if length(Xukf(2,i))<=2
        d=sqrt((X(1,i)-Xukf(1,i))^2+(X(3,i)-Xukf(2,i))^2);
    else
        d=sqrt((X(1,i)-Xukf(1,i))^2+(X(3,i)-Xukf(3,i))^2);
    end
    Err_KalmanFilter(i)=d;
%滤波后的误差
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%9
%画图
figure
hold on;box on;
plot(X(1,:),X(3,:),'-k.');
%真实轨迹
plot(Xukf(1,:),Xukf(3,:),'-r+');
% 无迹Kalman滤波轨迹
legend('真实轨迹','UKF轨迹')
figure
hold on;box on;
plot(Err_KalmanFilter,'-ks','MarkerFace','r');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % 子函数:求两点间的距离
% function d=Dist(X1,X2)
% if length(X2)<=2
%     d=sqrt((X1(1)-X2(1))^2+(X1(3)-X2(2))^2);
% else
%     d=sqrt((X1(1)-X2(1))^2+(X1(3)-X2(3))^2);
% end
% end
% %观测子函数:观测距离
% function [y]=hfun(x,xx)
% y=sqrt((x(1)-xx(1))^2+(x(3)-xx(2))^2);
% end

