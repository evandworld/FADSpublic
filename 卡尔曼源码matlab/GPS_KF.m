
function main

clc;
clear;

T=1;
N=80/T;
X=zeros(4,N);
X(:,1)=[-100,2,200,20];
Z=zeros(2,N);
Z(:,1)=[X(1,1),X(3,1)];
delta_w=1;
Q=delta_w*diag([0.5,1,0.5,1]);
R=100*eye(2);
F=[1,T,0,0;0,1,0,0;0,0,1,T;0,0,0,1];
H=[1,0,0,0;0,0,1,0];
for t=2:N;
    X(:,t)=F*X(:,t-1)+sqrtm(Q)*randn(4,1);
    Z(:,t)=H*X(:,t)+sqrtm(R)*randn(2,1);
end
Xkf=zeros(4,N);
Xkf(:,1)=X(:,1);
P0=eye(4);
for i=2:N
    Xn=F*Xkf(:,i-1);
    P1=F*P0*F'+Q;
    K=P1*H'/(H*P1*H'+R);
    Xkf(:,i)=Xn+K*(Z(:,i)-H*Xn);
    P0=(eye(4)-K*H)*P1;
end
for i=1:N
    Err_Observation(i)=RMS(X(:,i),Z(:,i));
    Err_KF(i)=RMS(X(:,i),Xkf(:,i));
end
figure
hold on;box on;
plot(X(1,:),X(3,:),'-k'); %��ʵ�켣
plot(Z(1,:),Z(2,:),'-b.');%�۲�켣
plot(Xkf(1,:),Xkf(3,:),'-r+');% Kalman�˲��켣
legend('��ʵ�켣','�۲�켣','�˲��켣');
xlabel('������X/m');
ylabel('������ Y/m');
figure
hold on;box on;
plot(Err_Observation,'-ko','MarkerFace','g')
plot(Err_KF,'-ks','MarkerFace','r')
legend('�˲�ǰ���','�˲������')
xlabel('�۲�ʱ��/s');
ylabel('���ֵ');
%%%%%%%%%%%%%%%%%%%%%%%%%
%����ŷ�Ͼ����Ӻ���
function dist=RMS(X1,X2);
if length(X2)<=2
    dist=sqrt((X1(1)-X2(1))^2+(X1(3)-X2(2))^2);
else
    dist=sqrt((X1(1)-X2(1))^2+(X1(3)-X2(3))^2);
end


