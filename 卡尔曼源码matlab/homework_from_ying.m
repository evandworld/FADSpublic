clc;
clear;
% function main
N=100;
%CON=25;
%Xexpect=CON*ones(1,N);
X=zeros(2,N);
P=zeros(2,2);
X(:,1)=[0,1];
P0=[1,0;0,1];  %��Э�������ֵ������ȡ0.01�����ֵѡ��ѡ������ν����Ӱ�����ľ���
Z=zeros(1,N);

Q=[0.01,0;0,0.01];  %��������
R=1;  %�۲����
W=sqrt(Q)*randn(2,N);
V=sqrt(R)*randn(1,N);
F=[1,1;0,1];
G=1;
H=[1,0];
Z=H*X+V;
I=eye(2);
Xkf=zeros(2,N);
Xkf(:,1)=[Z(1),1];     %%%%%%�˴����ɻ�
for k=2:N
    X(:,k)=F*X(:,k-1)+G*W(:,k-1);   %ʵ��ֵ
    Z(k)=H*X(:,k)+V(k);         %�۲�ֵ
    X_pre=F*Xkf(:,k-1);
    P_pre=F*P0*F'+Q;
    K=P_pre*H'/(H*P_pre*H'+R);   %H'��H��ת��
    e=Z(k)-H*X_pre;     %�۲�ֵ��Ԥ��ֵ*H��Ĳ�ֵ
    Xkf(:,k)=X_pre+K*e;     %kalman�˲�ֵ
    P0=(I-K*H)*P_pre;   %������һ����Э������
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
legend('λ����ʵֵ','λ�ƹ۲�ֵ','λ���˲�ֵ','λ�ò���ƫ��','λ���˲�ƫ��');
xlabel('����ʱ��/s');
ylabel('λ��/m');

%end