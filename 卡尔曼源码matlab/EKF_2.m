function EKF_2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear;
T=1;%�״�ɨ������,
N=60/T;%�ܵĲ�������
X=zeros(4,N);%Ŀ����ʵλ�á��ٶ�
X(:,1)=[-100,2,200,20];%Ŀ���ʼλ�á��ٶ�
Z=zeros(1,N);%��������λ�õĹ۲�
w=1e-3;%����������������Ŀ����ʵ�켣����������
Q=w*diag([0.5,1]);%������������
G=[T^2/2,0;T,0;0,T^2/2;0,T];%����������������
R=5;%�۲���������
F=[1,T,0,0;0,1,0,0;0,0,1,T;0,0,0,1];%״̬ת�ƾ���
x0=200;
y0=300;
%�۲�վ��λ�ã�������Ϊ����ֵ
Xstation=[x0,y0];
for t=2:N
    X(:,t)=F*X(:,t-1)+G*sqrtm(Q)*randn(2,1); %Ŀ����ʵ�켣
end
for t=1:N
    Z(t)=Dist(X(:,t),Xstation)+sqrtm(R)*randn;%��Ŀ��۲�
end
% EKF�˲�
Xekf=zeros(4,N);
Xekf(:,1)=X(:,1);% Kalman�˲�״̬��ʼ��
PO=eye(4);%Э�������ʼ��
D=zeros(N);
for i=2:N
    Xn=F*Xekf(:,i-1);%Ԥ��
    P1=F*PO*F'+G*Q*G';%Ԥ�����Э����
    dd=Dist(Xn,Xstation);%�۲�Ԥ��
    D(i)=dd;
    %���ſɱȾ���H
    H=[(Xn(1,1)-x0)/dd,0,(Xn(3,1)-y0)/dd,0];%��Ϊ����һ�׽���
    K=P1*H'/(H*P1*H'+R);%����
    Xekf(:,i)=Xn+K*(Z(:,i)-dd);%״̬����
    PO=(eye(4)-K*H)*P1;%�˲����Э�������
end
%������
for i=1:N
    Err_KalmanFilter(i)=Dist(X(:,1),Xekf(:,i));%�˲�������
end
%%%%%%%%%%%%%%%%%%%%%%
%��ͼ
figure
hold on;box on;
plot(X(1,:),X(3,:),'-k.');%��ʵ�켣
plot(Xekf(1,:),Xekf(3,:),'-r+');%��չKalman�˲��켣

legend('��ʵ�켣','EKF�켣');
xlabel('������X/m');
ylabel('������Y/m');
figure
hold on;box on;
t=1:60;
plot(t,D,'-k');
plot(t,Z,'-b');
figure
hold on;box on;
plot(Err_KalmanFilter,'-ks','MarkerFace','r')
xlabel('ʱ��/s');
ylabel('λ�ù���ƫ��/m');
%%�Ӻ���
function d=Dist(X1,X2)
if length(X2)<=2;
    d=sqrt((X1(1)-X2(1))^2+(X1(3)-X2(2))^2);
else
    d=sqrt((X1(1)-X2(1))^2+(X1(3)-X2(2))^2);
end
