function EKF_2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear;
T=1;%雷达扫描周期,
N=60/T;%总的采样次数
X=zeros(4,N);%目标真实位置、速度
X(:,1)=[-100,2,200,20];%目标初始位置、速度
Z=zeros(1,N);%传感器对位置的观测
w=1e-3;%如果增大这个参数，目标真实轨迹就是曲线了
Q=w*diag([0.5,1]);%过程噪声方差
G=[T^2/2,0;T,0;0,T^2/2;0,T];%过程噪声驱动矩阵
R=5;%观测噪声方差
F=[1,T,0,0;0,1,0,0;0,0,1,T;0,0,0,1];%状态转移矩阵
x0=200;
y0=300;
%观测站的位置，可以设为其他值
Xstation=[x0,y0];
for t=2:N
    X(:,t)=F*X(:,t-1)+G*sqrtm(Q)*randn(2,1); %目标真实轨迹
end
for t=1:N
    Z(t)=Dist(X(:,t),Xstation)+sqrtm(R)*randn;%对目标观测
end
% EKF滤波
Xekf=zeros(4,N);
Xekf(:,1)=X(:,1);% Kalman滤波状态初始化
PO=eye(4);%协方差阵初始化
D=zeros(N);
for i=2:N
    Xn=F*Xekf(:,i-1);%预测
    P1=F*PO*F'+G*Q*G';%预测误差协方差
    dd=Dist(Xn,Xstation);%观测预测
    D(i)=dd;
    %求雅可比矩阵H
    H=[(Xn(1,1)-x0)/dd,0,(Xn(3,1)-y0)/dd,0];%即为所求一阶近似
    K=P1*H'/(H*P1*H'+R);%增益
    Xekf(:,i)=Xn+K*(Z(:,i)-dd);%状态更新
    PO=(eye(4)-K*H)*P1;%滤波误差协方差更新
end
%误差分析
for i=1:N
    Err_KalmanFilter(i)=Dist(X(:,1),Xekf(:,i));%滤波后的误差
end
%%%%%%%%%%%%%%%%%%%%%%
%画图
figure
hold on;box on;
plot(X(1,:),X(3,:),'-k.');%真实轨迹
plot(Xekf(1,:),Xekf(3,:),'-r+');%扩展Kalman滤波轨迹

legend('真实轨迹','EKF轨迹');
xlabel('横坐标X/m');
ylabel('纵坐标Y/m');
figure
hold on;box on;
t=1:60;
plot(t,D,'-k');
plot(t,Z,'-b');
figure
hold on;box on;
plot(Err_KalmanFilter,'-ks','MarkerFace','r')
xlabel('时间/s');
ylabel('位置估计偏差/m');
%%子函数
function d=Dist(X1,X2)
if length(X2)<=2;
    d=sqrt((X1(1)-X2(1))^2+(X1(3)-X2(2))^2);
else
    d=sqrt((X1(1)-X2(1))^2+(X1(3)-X2(2))^2);
end
