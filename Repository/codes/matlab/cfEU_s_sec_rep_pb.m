
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Brexit & Post-Brexit EU collapse: BOOTSTRAP 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clear all
clc

workdir='Repository/simulation';
cd(workdir)



load MATLAB/initial_condition_2014

% import estimated thetas 
Tmat=importdata('input/Bdraws/sectoral/bltau-sec-nnew.txt');
Tmat=Tmat.data;
he = size(Tmat,1);
Tmat = -ones(he,J)./Tmat;
Tmat=Tmat';
% (each columns corresponds to one draw)

% import estimated PTA coefficients
hel=importdata('input/Bdraws/sectoral/bschengen-sec-nnew.txt');
delta_schengen_mat = hel.data;
delta_schengen_mat=-delta_schengen_mat'.*Tmat;
%delta_schengen_mat(~isfinite(delta_schengen_mat))=0; % setting cost changes to zero when trade elasticity is zero


hel=importdata('input/Bdraws/sectoral/bbotheea-sec-nnew.txt');
delta_Smarket_mat = hel.data;
delta_Smarket_mat=-delta_Smarket_mat'.*Tmat;
%delta_Smarket_mat(~isfinite(delta_Smarket_mat))=0; % setting cost changes to zero when trade elasticity is zero

hel=importdata('input/Bdraws/sectoral/bbotheuro-sec-nnew.txt');
delta_euro_mat = hel.data;
delta_euro_mat=-delta_euro_mat'.*Tmat;
%delta_euro_mat(~isfinite(delta_euro_mat))=0; % setting cost changes to zero when trade elasticity is zero


hel=importdata('input/Bdraws/sectoral/beukor-sec-nnew.txt');
delta_eukor_mat = hel.data;
delta_eukor_mat=-delta_eukor_mat'.*Tmat;
%delta_eukor_mat(~isfinite(delta_eukor_mat))=0; % setting cost changes to zero when trade elasticity is zero


hel=importdata('input/Bdraws/sectoral/bbothrta-sec-nnew.txt');
delta_oRTA_mat = hel.data;
delta_oRTA_mat=-delta_oRTA_mat'.*Tmat;
%delta_oRTA_mat(~isfinite(delta_oRTA_mat))=0; % setting cost changes to zero when trade elasticity is zero

clear hel



%convert T from cif to fob

Tmat=1./Tmat;
Tmat = -Tmat + ones(size(Tmat));
Tmat=-1./Tmat;

%set trade elasticity to a value near zero if zero, i.e. 1/theta is set to
%a big number
Tmat(~isfinite(Tmat))=1000; 

% BREXIT

fta_hat_oRTA_brexit=importdata('input/brexit_dummy_diff_oRTA.txt'); 
fta_hat_eukor_brexit=importdata('input/brexit_dummy_diff_eukor.txt'); 
fta_hat_Smarket_brexit=importdata('input/brexit_dummy_diff_Smarket.txt'); 
fta_hat_schengen_brexit=importdata('input/brexit_dummy_diff_Schengen.txt'); % #borders_old - # borders_new
fta_hat_euro_brexit=importdata('input/brexit_dummy_diff_euro.txt'); 

%counterfactual tariffs
taupBR=importdata('input/taup_brexit2014.txt'); 


% POST-BREXIT

fta_hat_oRTA=importdata('input/pb_dummy_diff_oRTA.txt');  %post brexit dummy changes
fta_hat_eukor=importdata('input/pb_dummy_diff_eukor.txt');  %post brexit dummy changes
fta_hat_Smarket=importdata('input/pb_dummy_diff_Smarket.txt'); 
fta_hat_schengen=importdata('input/pb_dummy_diff_Schengen.txt'); % #borders_old - # borders_new
fta_hat_euro=importdata('input/pb_dummy_diff_euro.txt'); 

%counterfactual tariffs
taupEU=importdata('input/taup_eu2014.txt'); 
%counterfactual tariffs
taupoRTA=importdata('input/taup_oRTA2014.txt'); 
taupallEU=importdata('input/taup_allEU2014.txt'); 



save 'MATLAB/initial_condition_brexit_reps'




%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

clearvars -except workdir
cd(strcat(workdir,'/MATLAB'))
Reps=1000;


for reps=1:Reps
    
    load 'initial_condition_brexit_reps'
    
    time=datestr(now)
reps    

   T(:,1)=Tmat(:,reps);
   delta_schengen(:,1)=delta_schengen_mat(:,reps);
   delta_Smarket(:,1)=delta_Smarket_mat(:,reps);
   delta_euro(:,1)=delta_euro_mat(:,reps);
   delta_oRTA(:,1)=delta_oRTA_mat(:,reps);
   delta_eukor(:,1)=delta_eukor_mat(:,reps);

D_schengen=[];
D_euro=[];
D_Smarket=[];
D_oRTA=[];
D_eukor=[];


for j=1:J
    for i=1:N
        D_schengen = [D_schengen ; delta_schengen(j)];
        D_euro = [D_euro ; delta_euro(j)];
        D_Smarket = [D_Smarket ; delta_Smarket(j)];
        D_oRTA = [D_oRTA ; delta_oRTA(j)];
        D_eukor = [D_eukor ; delta_eukor(j)];
    end
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  generate inventory free equilibrium with new algorithm %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% calculate surplus as share of income.
% ignore inventories, as those will be set to zero

R=(X0'*(1-F)).*eye(N,N)*ones(N,1);
sd = Sn./(VAn+R);
VP=zeros(J,N);
SIn = Sn;
Iv = zeros(N,1);
tau_hat=tau./tau;
taup=tau;

vfactor  = -.1;
tol      = 1E-6;
maxit    = 200;
wstart   = ones(N,1);
pstart   = ones(J,N);


   
[wf0_all pf0_all PQ_all Fp_all Dinp_all emax] = equi_s_emax(tau_hat,taup,alphas,T,B,G,sd,Din,J,N,maxit,tol,VAn,Iv,VP,vfactor,wstart,pstart); %add wstart for use with equilibrium_LC
% this version delivers the max iteration




% if no convergence, choose smaller factor
emax2=[];  maxit2    = 200;
if (emax == maxit)
    vf2=-.01;
    tol2      = 1E-7;
  
[wf0_all pf0_all PQ_all Fp_all Dinp_all emax2] = equi_s_emax(tau_hat,taup,alphas,T,B,G,sd,Din,J,N,maxit2,tol2,VAn,Iv,VP,vf2,wstart,pstart); %add wstart for use with equilibrium_LC
end


% if no convergence, choose smaller factor
if (emax2 == maxit2)
    vf2=-.001;
    tol2      = 1E-8;
    maxit2    = 500;
[wf0_all pf0_all PQ_all Fp_all Dinp_all emax] = equi_s_emax(tau_hat,taup,alphas,T,B,G,sd,Din,J,N,maxit2,tol2,VAn,Iv,VP,vf2,wstart,pstart); %add wstart for use with equilibrium_LC
end




% new io-matrix

%      s1    s2
%s1 x1 i1 i2 i1 i2
%   x2
%s2 x1
%   x2

% new baseline
Din = Dinp_all;
F = Fp_all;
X0 = PQ_all;
VAn=wf0_all.*VAn;
Inv = zeros(J,N);

clear Dinp_all Fp_all pf0_all PQ_all SIn wf0_all 



%%% sim Brexit

DDD_oRTA=D_oRTA.*fta_hat_oRTA_brexit;
DDD_eukor=D_eukor.*fta_hat_eukor_brexit;
DDD_schengen=D_schengen.*fta_hat_schengen_brexit;
DDD_euro=D_euro.*fta_hat_euro_brexit;
DDD_Smarket=D_Smarket.*fta_hat_Smarket_brexit;
taup=taupBR;
tau_hat=taupBR./tau .* exp(DDD_oRTA).*exp(DDD_eukor).* exp(DDD_schengen).* exp(DDD_euro).* exp(DDD_Smarket);  


vfactor  = -.1;
tol      = 1E-6;
maxit    = 200;
wstart   = ones(N,1);
pstart   = ones(J,N);

[wf0_all pf0_all PQ_all Fp_all Dinp_all] = equi_s(tau_hat,taup,alphas,T,B,G,sd,Din,J,N,maxit,tol,VAn,Iv,VP,vfactor,wstart,pstart); %add wstart for use with equilibrium_LC




% if no convergence, choose smaller factor
emax2=[];  maxit2    = 200;
if (emax == maxit)
    vf2=-.01;
    tol2      = 1E-7;
  
[wf0_all pf0_all PQ_all Fp_all Dinp_all emax2] = equi_s_emax(tau_hat,taup,alphas,T,B,G,sd,Din,J,N,maxit2,tol2,VAn,Iv,VP,vf2,wstart,pstart); %add wstart for use with equilibrium_LC
end


% if no convergence, choose smaller factor
if (emax2 == maxit2)
    vf2=-.001;
    tol2      = 1E-8;
    maxit2    = 500;
[wf0_all pf0_all PQ_all Fp_all Dinp_all emax] = equi_s_emax(tau_hat,taup,alphas,T,B,G,sd,Din,J,N,maxit2,tol2,VAn,Iv,VP,vf2,wstart,pstart); %add wstart for use with equilibrium_LC
end



I=(VAn+(X0'*(1-F)).*eye(N,N)*ones(N,1)).*(1-sd);
Ip=(diag(wf0_all)*VAn+(PQ_all'*(1-Fp_all)).*eye(N,N)*ones(N,1)).*(1-sd);

% income change
Ihat=Ip./I;
% change in price index
P=pf0_all.^alphas;
P0_all=prod(P,1);
P0_all=P0_all';
what=Ihat./P0_all;



if (reps==1)
    dlmwrite(strcat(workdir,'/Results/Brexit1000_secB/what.txt'),what)

end

if (reps>1)
    dlmwrite(strcat(workdir,'/Results/Brexit1000_secB/what.txt'),what,'-append')

end



% new baseline
Din = Dinp_all;
F = Fp_all;
X0 = PQ_all;
VAn=wf0_all.*VAn;
Inv = zeros(J,N);
tau=taup;



clear Dinp_all Fp_all pf0_all PQ_all SIn wf0_all


%% ALL EU


DDD_oRTA=D_oRTA.*fta_hat_oRTA; 
DDD_Smarket=D_Smarket.*fta_hat_Smarket;                                           
DDD_euro=D_euro.*fta_hat_euro;                   
DDD_schengen=D_schengen.*fta_hat_schengen;                   
DDD_eukor=D_eukor.*fta_hat_eukor;                   

taup=taupallEU;  
tau_hat=taup./tau .* exp(DDD_oRTA).* exp(DDD_schengen).* exp(DDD_euro).* exp(DDD_Smarket).* exp(DDD_eukor);     

%load starting values from scenario with parameter means
load(strcat(workdir,'/Results/allEU_smdB/allEU_smd_pb'), 'pf0_all', 'wf0_all')                 
wstart = wf0_all;
pstart = pf0_all;
clear wf0_all pf0_all

[wf0_all pf0_all PQ_all Fp_all Dinp_all emax] = equi_s_emax(tau_hat,taup,alphas,T,B,G,sd,Din,J,N,maxit,tol,VAn,Iv,VP,vfactor,wstart,pstart); %add wstart for use with equilibrium_LC


% if no convergence, choose smaller factor
emax2=[];  maxit2    = 200;
if (emax == maxit)
    vf2=-.01;
    tol2      = 1E-7;
  
[wf0_all pf0_all PQ_all Fp_all Dinp_all emax2] = equi_s_emax(tau_hat,taup,alphas,T,B,G,sd,Din,J,N,maxit2,tol2,VAn,Iv,VP,vf2,wstart,pstart); %add wstart for use with equilibrium_LC
end


% if no convergence, choose smaller factor
if (emax2 == maxit2)
    vf2=-.001;
    tol2      = 1E-8;
    maxit2    = 500;
[wf0_all pf0_all PQ_all Fp_all Dinp_all emax] = equi_s_emax(tau_hat,taup,alphas,T,B,G,sd,Din,J,N,maxit2,tol2,VAn,Iv,VP,vf2,wstart,pstart); %add wstart for use with equilibrium_LC
end


% results output

I=(VAn+(X0'*(1-F)).*eye(N,N)*ones(N,1)).*(1-sd);
Ip=(diag(wf0_all)*VAn+(PQ_all'*(1-Fp_all)).*eye(N,N)*ones(N,1)).*(1-sd);

% income change
Ihat=Ip./I;

% change in price index
P=pf0_all.^alphas;
P0_all=prod(P,1);
P0_all=P0_all';
what=Ihat./P0_all;


clear DDD Dinp_all Fp_all PQ_all pf0_all wf0_all

if (reps==1)
dlmwrite(strcat(workdir,'/Results/allEU1000_secB/what_pb.txt'),what)                  
end

if (reps>1)
dlmwrite(strcat(workdir,'/Results/allEU1000_secB/what_pb.txt'),what,'-append')                  
end



end
