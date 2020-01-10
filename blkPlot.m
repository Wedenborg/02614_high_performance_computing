%% PLOT
clear all; clc; close all;

T1 = (readtable('blk_NonOptim.txt'));
Size1 = table2array(T1(:,5));
speed1  = table2array(T1(:,2));

T2 = (readtable('blk_Optim.txt'));
Size2 = table2array(T2(:,5));
speed2  = table2array(T2(:,2));

figure
semilogx(Size1,speed1,'-o'); 
xlabel('Size of block');
ylabel('Performance [Mflop/s]')
legend('Non Optim');
set(gca,'FontSize',10)

figure
semilogx(Size2, speed2, '-or')
xlabel('Size of block');
ylabel('Performance [Mflop/s]')
legend( 'Optim');
set(gca,'FontSize',10)





