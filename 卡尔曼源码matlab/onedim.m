clc;
clear all;
% function main
N=60;
%CON=25;
%Xexpect=CON*ones(1,N);
X=zeros(2,N);
P=zeros(2,2);
X(:,1)=(0,1);
P(1)=0;  %��Э�������ֵ������ȡ0.01�����ֵѡ��ѡ������ν����Ӱ�����ľ���
Z(1)=24.9;
Xkf(1)=Z(1);
Q=0.01;  %��������
R=0.25;  %�۲����
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
    Kg=H*P_pre*inv(H*P_pre*H'+R);   %H'��H��ת�ã�inv����,�ɻ���/
    e=Z(k)-H*X_pre;     %�۲�ֵ��Ԥ��ֵ*H��Ĳ�ֵ
    Xkf(k)=X_pre+Kg*e;     %kalman�˲�ֵ
    P(k)=(I-Kg*H)*P_pre;   %������һ����Э������
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
legend('��ʵֵ','�۲�ֵ','�˲�ֵ');
xlabel('����ʱ��/s');
ylabel('�¶�ֵ/��');
figure
plot(t,Err_Messure,'-b.',t,Err_Kalman,'-k*');
legend('����ƫ��','�˲�ƫ��');
xlabel('����ʱ��/s');
ylabel('�¶�ƫ��ֵ/��');
%end