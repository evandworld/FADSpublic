clc;
clear;
% function main
N=100;
%CON=25;
%Xexpect=CON*ones(1,N);
X=zeros(2,N);
P=zeros(2,2);
X(:,1)=[0,1];
P0=[1,0;0,1];  %（协）方差初值，例程取0.01，这个值选多选少无所谓，不影响后面的精度
Z=zeros(1,N);

Q=[0.01,0;0,0.01];  %外界白噪声
R=1;  %观测误差
W=sqrt(Q)*randn(2,N);
V=sqrt(R)*randn(1,N);
F=[1,1;0,1];
G=1;
H=[1,0];
Z=H*X+V;
I=eye(2);
Xkf=zeros(2,N);
Xkf(:,1)=[Z(1),1];     %%%%%%此处有疑惑
for k=2:N
    X(:,k)=F*X(:,k-1)+G*W(:,k-1);   %实际值
    Z(k)=H*X(:,k)+V(k);         %观测值
    X_pre=F*Xkf(:,k-1);
    P_pre=F*P0*F'+Q;
    K=P_pre*H'/(H*P_pre*H'+R);   %H'是H的转置
    e=Z(k)-H*X_pre;     %观测值与预测值*H后的差值
    Xkf(:,k)=X_pre+K*e;     %kalman滤波值
    P0=(I-K*H)*P_pre;   %计算下一个（协）方差
end
Err_Messure=zeros(1,N);
Err_Kalman=zeros(1,N);
for k=1:N
    Err_Messure(k)=abs(Z(k)-X(1,k));
    Err_Kalman(k)=abs(Xkf(1,k)-X(1,k));
end
t=1:N;
figure
%plot(t,Xexpect,'-b',t,X,'-r.',t,Z,'-ko',t,Xkf,'-g*');
plot(t,X(1,:),'-r.',t,Z,'-ko',t,Xkf(1,:),'-g*',t,Err_Messure,'-b.',t,Err_Kalman,'-k*');
legend('位移真实值','位移观测值','位移滤波值','位置测量偏差','位置滤波偏差');
xlabel('采样时间/s');
ylabel('位移/m');

%end