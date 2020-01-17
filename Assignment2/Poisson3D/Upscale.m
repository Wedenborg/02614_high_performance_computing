%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot for upscale
clear all; close all; clc
%% Define data
v1       = importdata('v1_11.txt');
N        = unique(v1(:,2));   % N

mem = N.^3*3*64;              % bits
mem = mem*1/8*0.001;          % kB

thread  = unique(v1(:,1));    % Number of threads
memory  = v1(:,2).^3*3*64;    % Memory in bits
memory  = memory*(1/8)*0.001;   % Memory in kB
time    = v1(:,3);            % Time [sec]
iter    = 2000./time;         % Iterations = 2000, Performance [iter/sec]
g       = N.^3;
g       = repmat(g,length(thread));
g       = g(:,1);

per     = (g*12*2000*10^(-6))./time; %Mflops/s

%% Plot
z = 1:length(N):length(memory);
s = length(N)-1;
%z = z(1:7);    % adjust number of visualized threads


figure(1)
for i=1:length(z)
    semilogx(memory(z(i):z(i)+s),per(z(i):z(i)+s),'-o','Linewidth',0.8); hold on
end
grid on
xlabel('Memory [kB]'); ylabel('Performance [Mflops/s]');
axis([5000 500000 1500 6500]);
%xticks([logspace(3,6,10)])
legend(' 1  Thread',' 2  Threads',' 4  Threads',' 8  Threads','12 Threads',...
    '16 Threads','20 Threads','24 Threads','Location','Southwest')
title('Upscale: Naïve parallel Jacobi')
set(gca,'FontSize',12)
set(gcf, 'Position',  [100, 100, 800, 500])
