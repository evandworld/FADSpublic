
clc;
clear;
T=1;%ɨ������
N=60/T;% �ܵĲ�������
X=zeros(4,N);% Ŀ����ʵλ�á��ٶ�
X(:,1)=[-100,2,200,20];
Z=zeros(1,N);
delta_w=1e-4;
Q=delta_w*diag([0.5,1]);    %����������ֵ
G=[T^2/2,0;T,0;0,T^2/2;0,T];
R=1;     %�۲���������
F=[1,T,0,0;0,1,0,0;0,0,1,T;0,0,0,1];
x0=200;
y0=300;
Xstation=[x0,y0];
v=sqrtm(R)*randn(1,N);
for t=2:N;
    X(:,t)=F*X(:,t-1)+G*sqrtm(Q)*randn(2,1);%%ʵ��ֵ
end
for t=1:N
    if length(Xstation)<=2
        d=sqrt((X(1,t)-Xstation(1))^2+(X(3,t)-Xstation(2))^2);
    else
        d=sqrt((X(1,t)-Xstation(1))^2+(X(3,t)-Xstation(3))^2);
    end
    Z(t)=d+v(t);      %%�۲�ֵ
end
%%%%%%%%%%%%%%%%%
%UKF �˲�,UT�任
L=4;    %����X������
alpha=1;
kalpha=0;
belta=2;
ramda=3-L;  %-1
for j=1:2*L+1    %j=1/2/3/4.../9  ���ɱ��������ԳƲ������ԣ��ɵõ�2*L+1
                 %%wm��wc�������Ȩ�أ�����һ���Ϊ3/2
    Wm(j)=1/(2*(L+ramda));    %��ֵΪ1/6
    Wc(j)=1/(2*(L+ramda));
end
Wm(1)=ramda/(L+ramda);  %-1/3
Wc(1)=ramda/(L+ramda)+1-alpha^2+belta; %Ȩֵ����
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Xukf=zeros(4,N);
Xukf(:,1)=X(:,1);%�޼�Kalman�˲�״̬��ʼ��
PO=eye(4);%Э�������ʼ��
for t=2:N
    xestimate= Xukf(:,t-1);     %ȡ��һ���˲����ֵ
    P=PO;
    %��һ��:���һ�� Sigma�㼯
    cho=(chol(P*(L+ramda)))';  %chol() ���ڶԾ������Cholesky�ֽ�
                               %�ٶ�X��������R=chol(X)�����һ����������R
                               %ʹR'R=X.
    for k=1:L
        xgamaP1(:,k)=xestimate+cho(:,k);
        xgamaP2(:,k)=xestimate+cho(:,k);
    end
    Xsigma=[xestimate,xgamaP1,xgamaP2]; %Sigma�㼯��4�У�3+3+3=9�еľ���
    %�ڶ���:�� Sigma�㼯����һ��Ԥ��
    Xsigmapre=F*Xsigma; %����4*9
    %������:���õڶ����Ľ�������ֵ��Э����
    Xpred=zeros(4,1); %��ֵ
    for k=1:2*L+1
        Xpred=Xpred+Wm(k)*Xsigmapre(:,k);
    end
    Ppred=zeros(4,4);%Э������Ԥ��
    for k=1:2*L+1
        Ppred=Ppred+Wc(k)*(Xsigmapre(:,k)-Xpred)*(Xsigmapre(:,k)-Xpred)';
    end
    Ppred=Ppred+G*Q*G';
    %���Ĳ�:����Ԥ��ֵ����һ��ʹ��UT�任���õ��µ�sigma�㼯
    chor=(chol((L+ramda)*Ppred))';
    for k=1:L
     XaugsigmaP1(:,k)=Xpred+chor(:,k);
     XaugsigmaP2(:,k)=Xpred-chor(:,k);
    end
    Xaugsigma=[Xpred XaugsigmaP1 XaugsigmaP2];
    %���岽:�۲�Ԥ��
    for k=1:2*L+1
        Zsigmapre(1,k)=sqrt((Xaugsigma(1,k)-Xaugsigma(1,k))^2+(Xaugsigma(3,k)-Xaugsigma(2,k))^2);
    end
    %������:����۲�Ԥ���ֵ��Э����
    %�۲�Ԥ��ľ�ֵ
    Zpred=0;
    for k=1:2*L+1
        Zpred=Zpred+Wm(k)*Zsigmapre(1,k);
    end
    Pzz=0;
    for k=1:2*L+1
        Pzz=Pzz+Wc(k)*(Zsigmapre(1,k)-Zpred)*(Zsigmapre(1,k)-Zpred)';
    end
    Pzz=Pzz+R;%�õ�Э����Pzz
    Pxz=zeros(4,1);
    for k=1:2*L+1
        Pxz=Pxz+Wc(k)*(Xaugsigma(:,k)-Xpred)*(Zsigmapre(1,k)-Zpred)';
    end
    %���߲�:����Kalman����
    %Kalman����
    K=Pxz*inv(Pzz);
    %�ڰ˲�:״̬�ͷ������
    xestimate=Xpred+K*(Z(t)-Zpred);%״̬����
    P=Ppred-K*Pzz*K';  %�������
    PO=P;
    Xukf(:,t)=xestimate;
end
%������
for i=1:N
    if length(Xukf(2,i))<=2
        d=sqrt((X(1,i)-Xukf(1,i))^2+(X(3,i)-Xukf(2,i))^2);
    else
        d=sqrt((X(1,i)-Xukf(1,i))^2+(X(3,i)-Xukf(3,i))^2);
    end
    Err_KalmanFilter(i)=d;
%�˲�������
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%9
%��ͼ
figure
hold on;box on;
plot(X(1,:),X(3,:),'-k.');
%��ʵ�켣
plot(Xukf(1,:),Xukf(3,:),'-r+');
% �޼�Kalman�˲��켣
legend('��ʵ�켣','UKF�켣')
figure
hold on;box on;
plot(Err_KalmanFilter,'-ks','MarkerFace','r');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % �Ӻ���:�������ľ���
% function d=Dist(X1,X2)
% if length(X2)<=2
%     d=sqrt((X1(1)-X2(1))^2+(X1(3)-X2(2))^2);
% else
%     d=sqrt((X1(1)-X2(1))^2+(X1(3)-X2(3))^2);
% end
% end
% %�۲��Ӻ���:�۲����
% function [y]=hfun(x,xx)
% y=sqrt((x(1)-xx(1))^2+(x(3)-xx(2))^2);
% end

