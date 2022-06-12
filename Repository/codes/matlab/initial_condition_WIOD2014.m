    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Generate dataset initial conditions  -- 50-sector model
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


cd('Repository/simulation/input')


close all
clear all
clc

N = 44;
J = 50;

% expenditure shares
alphas=importdata('alphas.txt');

% value added cost shares
B=importdata('B.txt');

% input-output matrix
G=importdata('G.txt');

% initial trade share matrix
Din=importdata('Din2014.txt');

% value added
VAn=importdata('VA_n.txt');

% trade imbalance (net exports) + net inventory stock increase
SIn=importdata('SI_n.txt');

% trade imbalance (net exports) only
Sn=importdata('S_n.txt');

% expenditure
X0=importdata('X0.txt');

% trade tariff weighted shares
F=importdata('F.txt');

% inventory data
Inv=importdata('inventory2014.txt');


% actual input-output matrix for Leontief inverse
A=csvread('ioc2014_2.txt',0,0);

% sectoral demand for inventory changes
VP=importdata('VP.txt');

tau=importdata('tariffs2014.txt');                       % tariffs 2014

save '../MATLAB/initial_condition_2014'



% initial VA flows (VA originating in sector s in country o that is consumed
% in country d) -- based on original IO coefficients

clearvars -except workdir

beta=importdata('beta4va.txt');
fexp=importdata('fexp4va.txt');
A=importdata('io4va.txt');

va = diag(beta)*inv(eye(size(A,1))-A)*fexp;
dlmwrite('../Results/va_baseline.txt',va)





