%%%%%%%%%%%%%%%%%%%%% BASELINE 3-SECTOR MODEL  
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Generate dataset initial conditions 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd('Repository/simulation/input_3s')


close all
clear all
clc

N = 44;
J = 3;

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


% sectoral demand for inventory changes
VP=importdata('VP.txt');

tau=importdata('tariffs2014.txt');                       % tariffs 2014


% actual input-output matrix for Leontief inverse
A=csvread('ioc2014_2.txt',0,0);

save ../MATLAB/initial_condition_2014_3sec










