%% Plot for Amdahl
clear all; close all; clc
%% Define data
data_all = importdata('v1_11.txt');
TP   = data_all(5:10:end,3);
P    = unique(data_all(:,1));

% Ideal
x = 1:1:24;
y = x;

T1   = TP(1);
SP   = T1./TP;
f    = (1-TP/T1)./(1-(1./P));
f(1) = 1;
f4   = f(3);
f8   = f(4);

% Estimated performance w. 4 processors
E_TP4 = (f4./P)*T1+(1-f4)*T1;
E_SP4 = T1./E_TP4;
E_SP4(1) = 1;

% Estimated performance w. 8 processors
E_TP8 = (f8./P)*T1+(1-f8)*T1;
E_SP8 = T1./E_TP8;
E_SP8(1) = 1;

figure(1)
plot(P,SP,'-o'); hold on
plot(x,y,'-o')
plot(P,E_SP4,'-.')
plot(P,E_SP8,'-.')
legend('S(P)','Ideal','Est._4 f=0.64','Est._8 f=0.61','Location','Northwest')
xlabel('Processors'); ylabel('Speed-up'); grid on
title('Amdahls Law');
set(gca,'FontSize',12)
set(gcf, 'Position',  [100, 100, 800, 500])


