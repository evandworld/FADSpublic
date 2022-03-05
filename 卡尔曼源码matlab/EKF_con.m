function EKF_con
%初始化
%总时间
T=50;
%Q的值改变，观察不同Q值时滤波结果
Q=10;
%测量噪声
R=1;
%产生过程噪声
w=sqrt(Q)*randn(1,T);
%产生观测噪声
v=sqrt(R)*randn(1,T);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%状态方程
x=zeros(1,T);
x(1)=0.1;
y=zeros(1,T);
y(1)=x(1)^2/20+v(1);
for k=2:T
    x(k)=0.5*x(k-1)+2.5*x(k-1)/(1+x(k-1)^2)+8*cos(1.2*k)+w(k-1);
    y(k)=x(k)^2/20+v(k);
end
% EKF滤波算法
Xekf=zeros(1,T);
Xekf(1)=y(1);
% Yekf=zeros(1,T);
% Yekf(1)=y(1);
PO=eye(1);
for k=2:T
    %状态预测
    Xn=0.5*Xekf(k-1)+2.5*Xekf(k-1)/(1+Xekf(k-1)^2)+8*cos(1.2*k);
    %观测预测
    Zn=Xn^2/20;
    %求状态矩阵F
    F=0.5+2.5*(1-Xn^2)/(1+Xn^2)^2;
    %求观测矩阵
    H=Xn/10;
    %协方差预测
    P=F*PO*F'+Q;
    %求Kalman增益
    K=P*H'/(H*P*H'+R);
    %状态更新
    Xekf(k)=Xn+K*(y(k)-Zn);
    %协方差阵更新
    PO=(eye(1)-K*H)*P;
end
%误差分析
Xstd=zeros(1,T);
for k=1:T
    Xstd(k)=Xekf(k)-x(k);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%画图
figure
hold on;box on;   %box on是将图表右边和上面加上边框
plot(x,'-ko','MarkerFace','g');
plot(Xekf,'-ks','MarkerFace','b');
plot(y,'-k.','MarkerFace','r');
plot(Xstd,'-r.','MarkerFace','r');
legend('真实值','Kalman滤波估计值','观测值','偏差')
xlabel('时间/s');
ylabel('状态值x');
%误差分析
figure
hold on;box on;
plot(Xstd,'-ko','MarkerFace','g');
xlabel('时间/s');
ylabel('状态估计偏差');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

