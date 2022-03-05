hold on
% plot(simout1,'.k');
plot(simout2,'-r');    %改进后模糊
plot(simout3,'-.');   %最后模糊
legend('原模糊控制','改进后的模糊控制');
title('Step Response');
xlabel('t/s');
ylabel('amplitude');