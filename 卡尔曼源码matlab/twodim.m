clc;
clear all;
N=100;
Q=[0,0;0,0];
R=1;
W=sqrt(Q)*randn(2,N);  %��ʵȫ��0
V=sqrt(R)*randn(1,N);
A=[1,1;0,1];%״̬ת�ƾ���
%������
B=[0.5;1];
U=-9.8;   %������-1�������������ٶ�Ӧ����9.8
%�۲����
H=[1,0];
%��ʼ��
X=zeros(2,N);%������ʵ״̬
X(:,1)=[95;1];%��ʼλ�ƺ��ٶ�
P0=[10,0;0,1];%��ʼ���
Z=zeros(1,N);
Z(1)=H*X(:,1)+V(1);%��ʼ�۲�ֵ
Xkf=zeros(2,N);% Kalman����״̬��ʼ��
Xkf(:,1)=Z(1);
err_P=zeros(N,2);
err_P(1,1)=P0(1,1);
err_P(1,2)=P0(2,2);
I=eye(2);
%��άϵͳ
for k=2:N
    %�������䣬��״̬���̵�����
    X(:,k)=A*X(:,k-1)+B*U+W(k);
    %λ�ƴ�������Ŀ����й۲�
    Z(k)=H*X(:,k)+V(k);
    % Kalman�˲�
    X_pre=A*Xkf(:,k-1)+B*U;%״̬Ԥ��
    P_pre=A*P0*A'+Q;%Э����Ԥ��
    Kg=P_pre*H'/(H*P_pre*H'+R);%����Kalman����
    Xkf(:,k)=X_pre+Kg*(Z(k)-H*X_pre);%״̬����
    P0=(I-Kg*H)*P_pre;%�������
    %������ֵ
    err_P(k,1)=P0(1,1);
    errP(k,2)=P0(2,2);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ������
messure_err_x=zeros(1,N);%λ�ƵĲ������
kalman_err_x=zeros(1,N);% Kalman���Ƶ�λ��ƫ��
kalman_err_v=zeros(1,N);%Kalman���Ƶ��ٶ�ƫ��
for k=1:N
    messure_err_x(k)=Z(k)-X(1,k);
    kalman_err_x(k)=Xkf(1,k)-X(1,k);
    kalman_err_v(k)=Xkf(2,k)-X(2,k);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%��ͼ���
%����ͼ
figure
plot(W);
xlabel('����ʱ��/s');
ylabel('��������');
figure
plot(V);
xlabel('����ʱ��/s');
ylabel('��������');
%λ��ƫ��
figure
hold on,box on;
plot(messure_err_x,'-r.');%������λ�����
plot(kalman_err_x,'-g.');%kalman����λ�����
legend('����','Kalman����')
xlabel('����ʱ��/s');
ylabel('λ��ƫ��/m');
%Kalman�ٶ�ƫ��
figure
plot(kalman_err_v);
xlabel('����ʱ��/s');
ylabel('�ٶ�ƫ��/(m/s)');
%����ֵ
figure
plot(err_P(:,1));
xlabel('����ʱ��/s');
ylabel('λ��������ֵ');
figure
plot(err_P(:,1));
xlabel('����ʱ��/s');
ylabel('�ٶ�������ֵ');
for k=1:20
    X_x(k)=X(1,k);
    Xkf_x(k)=Xkf(1,k);
end
t=1:20;
figure
%plot(t,X(1,:),'-g',t,Xkf(1,:),'*b');
plot(t,X_x,'-g.',t,Xkf_x,'-ko');
legend('ʵ��λ��','Kalman����λ��')
xlabel('����ʱ��/s');
ylabel('λ��/m');