%% PLOT FOR MATRIX MULTIPLICATION
clear all; clc; close all;
%% For mm_batch.nat.txt
T = readtable('mm_batch_nat.txt');
M = table2array(T(:,1));
memory = round(unique(M))';    % kB
speed  = table2array(T(:,2));  % Mflops/s

% Cache model
L1 = 32; L2 = 256; L3 = 30720; % kB

% Index for 6 permutations
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
plot(memory,knm,'-o'); 
l1 = line([L1 L1], [0 600],'Color','k','LineWidth',0.1,'LineStyle','-');
l2 = line([L1+L2 L1+L2], [0 600],'Color','k','LineWidth',0.1,'LineStyle','-');
l3 = line([L1+L2+L3 L1+L2+L3], [0 600],'Color','k','LineWidth',0.1,'LineStyle','-');
hold off
xlabel('Memory usage [kB]');
ylabel('Performance [Mflop/s]')
%axis([-inf inf 0 6000]);
legend('mnk','mkn','nmk','nkm','kmn','knm','Location','Southwest');
title('Performance (Native)')
set(gca,'FontSize',12)
set(gcf, 'Position',  [100, 100, 800, 400])

%% For mm_batch_O3.txt
T = readtable('mm_batch_O3.txt');
M = table2array(T(:,1));
memory = round(unique(M))';    % kB
speed  = table2array(T(:,2));  % Mflops/s

% Index for 6 permutations
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
xlabel('Memory usage [kB]');
ylabel('Performance [Mflop/s]')
legend('mnk','mkn','nmk','nkm','kmn','knm','Location','Southwest');
title('Performance (-O3)')
set(gca,'FontSize',12)
set(gcf, 'Position',  [100, 100, 800, 400])

%% For mm_batch_O3_fast.txt
T = readtable('mm_batch_O3_fast.txt');
M = table2array(T(:,1));
memory = round(unique(M))';          % kB
memory = memory(1:length(memory)-1); % NB: last 1 must be excluded
speed  = table2array(T(:,2));        % Mflops/s
speed  = speed(1:length(speed)-6);   % NB: last 6 must be excluded

% Index for 6 permutations
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
xlabel('Memory usage [kB]');
ylabel('Performance [Mflop/s]')
legend('mnk','mkn','nmk','nkm','kmn','knm','Location','Southwest');
title('Performance (-O3 -fast)')
set(gca,'FontSize',12)
set(gcf, 'Position',  [100, 100, 800, 400])

%% For mm_batch_O3_fast_unroll.txt
T = readtable('mm_batch_O3_fast_unroll.txt');
M = table2array(T(:,1));
memory = round(unique(M))';    % kB
speed  = table2array(T(:,2));  % Mflops/s

% Index for 6 permutations
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
xlabel('Memory usage [kB]');
ylabel('Performance [Mflop/s]')
legend('mnk','mkn','nmk','nkm','kmn','knm','Location','Southwest');
title('Performance (-O3 -fast -unroll)')
set(gca,'FontSize',12)
set(gcf, 'Position',  [100, 100, 800, 400])
