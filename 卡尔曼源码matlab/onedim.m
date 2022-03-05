clc;
clear all;
% function main
N=60;
%CON=25;
%Xexpect=CON*ones(1,N);
X=zeros(2,N);
P=zeros(2,2);
X(:,1)=(0,1);
P(1)=0;  %（协）方差初值，例程取0.01，这个值选多选少无所谓，不影响后面的精度
Z(1)=24.9;
Xkf(1)=Z(1);
Q=0.01;  %外界白噪声
R=0.25;  %观测误差
W=sqrt(Q)*randn(1,N);
V=sqrt(R)*randn(1,N);
F=1;
G=1;
H=1;
I=eye(1);
for k=2:N
    X(k)=F*X(k-1)+G*W(k-1);
    Z(k)=H*X(k)+V(k);
    X_pre=F*Xkf(k-1);
    P_pre=F*P(k-1)*F'+Q;
    Kg=H*P_pre*inv(H*P_pre*H'+R);   %H'是H的转置，inv是逆,可换成/
    e=Z(k)-H*X_pre;     %观测值与预测值*H后的差值
    Xkf(k)=X_pre+Kg*e;     %kalman滤波值
    P(k)=(I-Kg*H)*P_pre;   %计算下一个（协）方差
end
Err_Messure=zeros(1,N);
Err_Kalman=zeros(1,N);
for k=1:N
    Err_Messure(k)=abs(Z(k)-X(k));
    Err_Kalman(k)=abs(Xkf(k)-X(k));
end
t=1:N;
figure
%plot(t,Xexpect,'-b',t,X,'-r.',t,Z,'-ko',t,Xkf,'-g*');
plot(t,X,'-r.',t,Z,'-ko',t,Xkf,'-g*');
legend('真实值','观测值','滤波值');
xlabel('采样时间/s');
ylabel('温度值/℃');
figure
plot(t,Err_Messure,'-b.',t,Err_Kalman,'-k*');
legend('测量偏差','滤波偏差');
xlabel('采样时间/s');
ylabel('温度偏差值/℃');
%end