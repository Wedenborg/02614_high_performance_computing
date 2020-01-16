size = [ 75 100 125 150];

mem = size.^3*3*64;% bits

mem = mem*1/8*0.001; % kB

v1 = importdata('v1UpScale.txt')

v1(:,3)=2000./v1(:,3)