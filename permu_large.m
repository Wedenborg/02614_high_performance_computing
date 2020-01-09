%% PLOT
clear all; clc; close all;
T = (readtable('temp.txt'));
M = table2array(T(:,1));
memory = ceil(unique(M))';
speed  = table2array(T(:,2));

i = 1:6:size(speed,1);

for j =1:size(memory,2)
mnk(j) = speed(i(j));
mkn(j) = speed(i(j)+1);
kmn(j) = speed(i(j)+2);
knm(j) = speed(i(j)+3);
nkm(j) = speed(i(j)+4);
nmk(j) = speed(i(j)+5);
end 

figure
semilogx(memory,mnk,'-o'); hold on
plot(memory,mkn,'-o');
plot(memory,nmk,'-o');
plot(memory,nkm,'-o');
plot(memory,kmn,'-o');
plot(memory,knm,'-o'); hold off
xlabel('Memory usage [kBytes]');
ylabel('Performance [Mflop/s]')
legend('mnk','mkn','nmk','nkm','kmn','knm','Location','east');
set(gca,'FontSize',10)





