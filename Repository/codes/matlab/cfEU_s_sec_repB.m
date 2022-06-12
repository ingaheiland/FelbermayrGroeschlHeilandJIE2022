% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Main Counterfactuals BOOTSTRAP   %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



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


hel=importdata('input/Bdraws/sectoral/bbotheea-sec-nnew.txt');
delta_Smarket_mat = hel.data;
delta_Smarket_mat=-delta_Smarket_mat'.*Tmat;

hel=importdata('input/Bdraws/sectoral/bbotheuro-sec-nnew.txt');
delta_euro_mat = hel.data;
delta_euro_mat=-delta_euro_mat'.*Tmat;


hel=importdata('input/Bdraws/sectoral/beukor-sec-nnew.txt');
delta_eukor_mat = hel.data;
delta_eukor_mat=-delta_eukor_mat'.*Tmat;


hel=importdata('input/Bdraws/sectoral/bbothrta-sec-nnew.txt');
delta_oRTA_mat = hel.data;
delta_oRTA_mat=-delta_oRTA_mat'.*Tmat;

clear hel



%convert T from cif to fob

Tmat=1./Tmat;
Tmat = -Tmat + ones(size(Tmat));
Tmat=-1./Tmat;

%set trade elasticity to a value near zero if zero, i.e. 1/theta is set to
%a big number
Tmat(~isfinite(Tmat))=1000; 

fta_hat_oRTA=importdata('input/dummy_diff_oRTA.txt'); 
fta_hat_Smarket=importdata('input/dummy_diff_Smarket.txt'); 
fta_hat_schengen=importdata('input/dummy_diff_Schengen.txt'); % #borders_old - # borders_new
fta_hat_euro=importdata('input/dummy_diff_euro.txt'); 
fta_hat_eukor=importdata('input/dummy_diff_eukor.txt'); 


%counterfactual tariffs
taupEU=importdata('input/taup_eu2014.txt'); 
%counterfactual tariffs
taupoRTA=importdata('input/taup_oRTA2014.txt'); 
taupallEU=importdata('input/taup_allEU2014.txt'); 


save 'MATLAB/initial_condition_eu_repsB'




%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

%% SECTORAL EST. 1000 REPETITIONS


clearvars -except workdir
cd(strcat(workdir,'/MATLAB'))
Reps=1000;


for reps=1:Reps
    
    load 'initial_condition_eu_repsB'
    
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

DDD_oRTA=D_oRTA.*fta_hat_oRTA; 
DDD_Smarket=D_Smarket.*fta_hat_Smarket;                                           
DDD_euro=D_euro.*fta_hat_euro;                   
DDD_schengen=D_schengen.*fta_hat_schengen;                   
DDD_eukor=D_eukor.*fta_hat_eukor;                   




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

vfactor  = -.5;
tol      = 1E-7;
maxit    = 200;
wstart   = ones(N,1);
pstart   = ones(J,N);


   
[wf0_all pf0_all PQ_all Fp_all Dinp_all emax] = equi_s_emax(tau_hat,taup,alphas,T,B,G,sd,Din,J,N,maxit,tol,VAn,Iv,VP,vfactor,wstart,pstart); %add wstart for use with equilibrium_LC
% this version delivers the max iteration



% if no convergence, choose smaller factor
emax2=[];  maxit2    = 200;
if (emax == maxit)
    vf2=-.1;
    tol2      = 1E-8;
  
[wf0_all pf0_all PQ_all Fp_all Dinp_all emax2] = equi_s_emax(tau_hat,taup,alphas,T,B,G,sd,Din,J,N,maxit2,tol2,VAn,Iv,VP,vf2,wstart,pstart); %add wstart for use with equilibrium_LC
end


% if no convergence, choose smaller factor
if (emax2 == maxit2)
    vf2=-.01;
    tol2      = 1E-9;
    maxit2    = 500;
[wf0_all pf0_all PQ_all Fp_all Dinp_all emax] = equi_s_emax(tau_hat,taup,alphas,T,B,G,sd,Din,J,N,maxit2,tol2,VAn,Iv,VP,vf2,wstart,pstart); %add wstart for use with equilibrium_LC
end


% new io-matrix

%      s1    s2
%s1 x1 i1 i2 i1 i2
%   x2
%s2 x1
%   x2

% new baseline io tab
pihat=Din./Dinp_all;
for j=1:J
    Ares((j-1)*N+1:j*N,:)=pihat((j-1)*N+1:j*N,1:N)';
end
Ares(isnan(Ares))=1;
AA=kron(ones(1,J),Ares).*A;


% new baseline
Din = Dinp_all;
F = Fp_all;
X0 = PQ_all;
VAn=wf0_all.*VAn;
Inv = zeros(J,N);

clear Dinp_all Fp_all pf0_all PQ_all SIn wf0_all Ares A pihat



%%%% simulate scenarios starting from inventory-free baseline


%% ALL EU

taup=taupallEU;  
tau_hat=taup./tau .* exp(DDD_oRTA).* exp(DDD_schengen).* exp(DDD_euro).* exp(DDD_Smarket).* exp(DDD_eukor);     

%load starting values from scenario with parameter means
load(strcat(workdir,'/Results/allEU_smdB/allEU_smd'), 'pf0_all', 'wf0_all')                 
wstart = wf0_all;
pstart = pf0_all;
clear wf0_all pf0_all

[wf0_all pf0_all PQ_all Fp_all Dinp_all emax] = equi_s_emax(tau_hat,taup,alphas,T,B,G,sd,Din,J,N,maxit,tol,VAn,Iv,VP,vfactor,wstart,pstart); %add wstart for use with equilibrium_LC
% if no convergence, choose smaller factor
emax2=[];  maxit2    = 200;
if (emax == maxit)
    vf2=-.1;
    tol2      = 1E-8;
  
[wf0_all pf0_all PQ_all Fp_all Dinp_all emax2] = equi_s_emax(tau_hat,taup,alphas,T,B,G,sd,Din,J,N,maxit2,tol2,VAn,Iv,VP,vf2,wstart,pstart); %add wstart for use with equilibrium_LC
end


% if no convergence, choose smaller factor
if (emax2 == maxit2)
    vf2=-.01;
    tol2      = 1E-9;
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
wrhat=wf0_all./P0_all;

% exports - same structure as trade shares
%           exp1 exp2
%sec1 imp1
%     imp2
%sec2 imp1

EXf = (Din./tau).*reshape(X0',[N*J,1]);
EXfp = (Dinp_all./taup).*reshape(PQ_all',[N*J,1]);

% value added exports
% final demand
C = alphas.*I'; C=reshape(C',[N*J,1]);
Cp = alphas.*Ip'; Cp=reshape(Cp',[N*J,1]);
r = reshape(B',[N*J,1]);
ahat=Dinp_all./Din.*tau./taup;
for j=1:J
    Ares((j-1)*N+1:j*N,:)=ahat((j-1)*N+1:j*N,1:N)';
end
Ares(isnan(Ares))=1;
AAp=kron(ones(1,J),Ares).*AA;

FEXf=(Din./tau).*C;
FEXfp=(Dinp_all./taup).*Cp;

for j=1:J
    RFEXf((j-1)*N+1:j*N,:)=FEXf((j-1)*N+1:j*N,1:N)';
    RFEXfp((j-1)*N+1:j*N,:)=FEXfp((j-1)*N+1:j*N,1:N)';
end
RVA = diag(r)*inv(eye(size(AA,1))-AA)*RFEXf;
RVAp = diag(r)*inv(eye(size(AAp,1))-AAp)*RFEXfp;


clear DDD Dinp_all Fp_all PQ_all pf0_all wf0_all AAp Ares ahat

if (reps==1)
dlmwrite(strcat(workdir,'/Results/allEU1000_secB/rva.txt'),RVA)                
dlmwrite(strcat(workdir,'/Results/allEU1000_secB/rvap.txt'),RVAp)                
dlmwrite(strcat(workdir,'/Results/allEU1000_secB/exf.txt'),EXf)                
dlmwrite(strcat(workdir,'/Results/allEU1000_secB/exfp.txt'),EXfp)                
dlmwrite(strcat(workdir,'/Results/allEU1000_secB/wr.txt'),wrhat)                  
dlmwrite(strcat(workdir,'/Results/allEU1000_secB/what.txt'),what)                  
end

if (reps>1)
dlmwrite(strcat(workdir,'/Results/allEU1000_secB/rva.txt'),RVA,'-append')                
dlmwrite(strcat(workdir,'/Results/allEU1000_secB/rvap.txt'),RVAp,'-append')                
dlmwrite(strcat(workdir,'/Results/allEU1000_secB/exf.txt'),EXf,'-append')                
dlmwrite(strcat(workdir,'/Results/allEU1000_secB/exfp.txt'),EXfp,'-append')                
dlmwrite(strcat(workdir,'/Results/allEU1000_secB/wr.txt'),wrhat,'-append')                  
dlmwrite(strcat(workdir,'/Results/allEU1000_secB/what.txt'),what,'-append')                  
end




%% EURO

taup=tau;
tau_hat=taup./tau .*  exp(DDD_euro);     

%load starting values from scenario with parameter means
load(strcat(workdir,'/Results/Euro_smdB/Euro_smd'), 'pf0_all', 'wf0_all')                 
wstart = wf0_all;
pstart = pf0_all;
clear wf0_all pf0_all

[wf0_all pf0_all PQ_all Fp_all Dinp_all emax] = equi_s_emax(tau_hat,taup,alphas,T,B,G,sd,Din,J,N,maxit,tol,VAn,Iv,VP,vfactor,wstart,pstart); %add wstart for use with equilibrium_LC
% if no convergence, choose smaller factor
emax2=[];  maxit2    = 200;
if (emax == maxit)
    vf2=-.1;
    tol2      = 1E-8;
  
[wf0_all pf0_all PQ_all Fp_all Dinp_all emax2] = equi_s_emax(tau_hat,taup,alphas,T,B,G,sd,Din,J,N,maxit2,tol2,VAn,Iv,VP,vf2,wstart,pstart); %add wstart for use with equilibrium_LC
end


% if no convergence, choose smaller factor
if (emax2 == maxit2)
    vf2=-.01;
    tol2      = 1E-9;
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
wrhat=wf0_all./P0_all;


% exports - same structure as trade shares
%           exp1 exp2
%sec1 imp1
%     imp2
%sec2 imp1

EXf = (Din./tau).*reshape(X0',[N*J,1]);
EXfp = (Dinp_all./taup).*reshape(PQ_all',[N*J,1]);

% value added exports
% final demand
C = alphas.*I'; C=reshape(C',[N*J,1]);
Cp = alphas.*Ip'; Cp=reshape(Cp',[N*J,1]);
r = reshape(B',[N*J,1]);
ahat=Dinp_all./Din.*tau./taup;
for j=1:J
    Ares((j-1)*N+1:j*N,:)=ahat((j-1)*N+1:j*N,1:N)';
end
Ares(isnan(Ares))=1;
AAp=kron(ones(1,J),Ares).*AA;

FEXf=(Din./tau).*C;
FEXfp=(Dinp_all./taup).*Cp;

for j=1:J
    RFEXf((j-1)*N+1:j*N,:)=FEXf((j-1)*N+1:j*N,1:N)';
    RFEXfp((j-1)*N+1:j*N,:)=FEXfp((j-1)*N+1:j*N,1:N)';
end
RVA = diag(r)*inv(eye(size(AA,1))-AA)*RFEXf;
RVAp = diag(r)*inv(eye(size(AAp,1))-AAp)*RFEXfp;

clear DDD Dinp_all Fp_all PQ_all pf0_all wf0_all AAp Ares ahat

if (reps==1)
dlmwrite(strcat(workdir,'/Results/Euro1000_secB/rva.txt'),RVA)                
dlmwrite(strcat(workdir,'/Results/Euro1000_secB/rvap.txt'),RVAp)                
dlmwrite(strcat(workdir,'/Results/Euro1000_secB/exf.txt'),EXf)                
dlmwrite(strcat(workdir,'/Results/Euro1000_secB/exfp.txt'),EXfp)                
dlmwrite(strcat(workdir,'/Results/Euro1000_secB/wr.txt'),wrhat)                  
dlmwrite(strcat(workdir,'/Results/Euro1000_secB/what.txt'),what)                  
end

if (reps>1)
dlmwrite(strcat(workdir,'/Results/Euro1000_secB/rva.txt'),RVA,'-append')                
dlmwrite(strcat(workdir,'/Results/Euro1000_secB/rvap.txt'),RVAp,'-append')                
dlmwrite(strcat(workdir,'/Results/Euro1000_secB/exf.txt'),EXf,'-append')                
dlmwrite(strcat(workdir,'/Results/Euro1000_secB/exfp.txt'),EXfp,'-append')                
dlmwrite(strcat(workdir,'/Results/Euro1000_secB/wr.txt'),wrhat,'-append')                  
dlmwrite(strcat(workdir,'/Results/Euro1000_secB/what.txt'),what,'-append')                  
end












%% SCHENGEN


taup=tau;
tau_hat=taup./tau .*  exp(DDD_schengen);     

%load starting values from scenario with parameter means
load(strcat(workdir,'/Results/Schengen_smdB/Schengen_smd'), 'pf0_all', 'wf0_all')                 
wstart = wf0_all;
pstart = pf0_all;
clear wf0_all pf0_all

[wf0_all pf0_all PQ_all Fp_all Dinp_all emax] = equi_s_emax(tau_hat,taup,alphas,T,B,G,sd,Din,J,N,maxit,tol,VAn,Iv,VP,vfactor,wstart,pstart); %add wstart for use with equilibrium_LC
% if no convergence, choose smaller factor
emax2=[];  maxit2    = 200;
if (emax == maxit)
    vf2=-.1;
    tol2      = 1E-8;
  
[wf0_all pf0_all PQ_all Fp_all Dinp_all emax2] = equi_s_emax(tau_hat,taup,alphas,T,B,G,sd,Din,J,N,maxit2,tol2,VAn,Iv,VP,vf2,wstart,pstart); %add wstart for use with equilibrium_LC
end


% if no convergence, choose smaller factor
if (emax2 == maxit2)
    vf2=-.01;
    tol2      = 1E-9;
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
wrhat=wf0_all./P0_all;

% exports - same structure as trade shares
%           exp1 exp2
%sec1 imp1
%     imp2
%sec2 imp1

EXf = (Din./tau).*reshape(X0',[N*J,1]);
EXfp = (Dinp_all./taup).*reshape(PQ_all',[N*J,1]);

% value added exports
% final demand
C = alphas.*I'; C=reshape(C',[N*J,1]);
Cp = alphas.*Ip'; Cp=reshape(Cp',[N*J,1]);
r = reshape(B',[N*J,1]);
ahat=Dinp_all./Din.*tau./taup;
for j=1:J
    Ares((j-1)*N+1:j*N,:)=ahat((j-1)*N+1:j*N,1:N)';
end
Ares(isnan(Ares))=1;
AAp=kron(ones(1,J),Ares).*AA;
FEXf=(Din./tau).*C;
FEXfp=(Dinp_all./taup).*Cp;

for j=1:J
    RFEXf((j-1)*N+1:j*N,:)=FEXf((j-1)*N+1:j*N,1:N)';
    RFEXfp((j-1)*N+1:j*N,:)=FEXfp((j-1)*N+1:j*N,1:N)';
end
RVA = diag(r)*inv(eye(size(AA,1))-AA)*RFEXf;
RVAp = diag(r)*inv(eye(size(AAp,1))-AAp)*RFEXfp;

clear DDD Dinp_all Fp_all PQ_all pf0_all wf0_all AAp Ares ahat

if (reps==1)
dlmwrite(strcat(workdir,'/Results/Schengen1000_secB/rva.txt'),RVA)                
dlmwrite(strcat(workdir,'/Results/Schengen1000_secB/rvap.txt'),RVAp)                
dlmwrite(strcat(workdir,'/Results/Schengen1000_secB/exf.txt'),EXf)                
dlmwrite(strcat(workdir,'/Results/Schengen1000_secB/exfp.txt'),EXfp)                
dlmwrite(strcat(workdir,'/Results/Schengen1000_secB/wr.txt'),wrhat)                  
dlmwrite(strcat(workdir,'/Results/Schengen1000_secB/what.txt'),what)                  
end

if (reps>1)
dlmwrite(strcat(workdir,'/Results/Schengen1000_secB/rva.txt'),RVA,'-append')                
dlmwrite(strcat(workdir,'/Results/Schengen1000_secB/rvap.txt'),RVAp,'-append')                
dlmwrite(strcat(workdir,'/Results/Schengen1000_secB/exf.txt'),EXf,'-append')                
dlmwrite(strcat(workdir,'/Results/Schengen1000_secB/exfp.txt'),EXfp,'-append')                
dlmwrite(strcat(workdir,'/Results/Schengen1000_secB/wr.txt'),wrhat,'-append')                  
dlmwrite(strcat(workdir,'/Results/Schengen1000_secB/what.txt'),what,'-append')                  
end






%% SMARKET


taup=tau;
tau_hat=taup./tau .*  exp(DDD_Smarket);     

%load starting values from scenario with parameter means
load(strcat(workdir,'/Results/Smarket_smdB/Smarket_smd'), 'pf0_all', 'wf0_all')                 
wstart = wf0_all;
pstart = pf0_all;
clear wf0_all pf0_all

[wf0_all pf0_all PQ_all Fp_all Dinp_all emax] = equi_s_emax(tau_hat,taup,alphas,T,B,G,sd,Din,J,N,maxit,tol,VAn,Iv,VP,vfactor,wstart,pstart); %add wstart for use with equilibrium_LC
% if no convergence, choose smaller factor
emax2=[];  maxit2    = 200;
if (emax == maxit)
    vf2=-.1;
    tol2      = 1E-8;
  
[wf0_all pf0_all PQ_all Fp_all Dinp_all emax2] = equi_s_emax(tau_hat,taup,alphas,T,B,G,sd,Din,J,N,maxit2,tol2,VAn,Iv,VP,vf2,wstart,pstart); %add wstart for use with equilibrium_LC
end


% if no convergence, choose smaller factor
if (emax2 == maxit2)
    vf2=-.01;
    tol2      = 1E-9;
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
wrhat=wf0_all./P0_all;


% exports - same structure as trade shares
%           exp1 exp2
%sec1 imp1
%     imp2
%sec2 imp1

EXf = (Din./tau).*reshape(X0',[N*J,1]);
EXfp = (Dinp_all./taup).*reshape(PQ_all',[N*J,1]);

% value added exports
% final demand
C = alphas.*I'; C=reshape(C',[N*J,1]);
Cp = alphas.*Ip'; Cp=reshape(Cp',[N*J,1]);
r = reshape(B',[N*J,1]);
ahat=Dinp_all./Din.*tau./taup;
for j=1:J
    Ares((j-1)*N+1:j*N,:)=ahat((j-1)*N+1:j*N,1:N)';
end
Ares(isnan(Ares))=1;
AAp=kron(ones(1,J),Ares).*AA;
FEXf=(Din./tau).*C;
FEXfp=(Dinp_all./taup).*Cp;

for j=1:J
    RFEXf((j-1)*N+1:j*N,:)=FEXf((j-1)*N+1:j*N,1:N)';
    RFEXfp((j-1)*N+1:j*N,:)=FEXfp((j-1)*N+1:j*N,1:N)';
end
RVA = diag(r)*inv(eye(size(AA,1))-AA)*RFEXf;
RVAp = diag(r)*inv(eye(size(AAp,1))-AAp)*RFEXfp;

clear DDD Dinp_all Fp_all PQ_all pf0_all wf0_all AAp Ares ahat

if (reps==1)
dlmwrite(strcat(workdir,'/Results/Smarket1000_secB/rva.txt'),RVA)                
dlmwrite(strcat(workdir,'/Results/Smarket1000_secB/rvap.txt'),RVAp)                
dlmwrite(strcat(workdir,'/Results/Smarket1000_secB/exf.txt'),EXf)                
dlmwrite(strcat(workdir,'/Results/Smarket1000_secB/exfp.txt'),EXfp)                
dlmwrite(strcat(workdir,'/Results/Smarket1000_secB/wr.txt'),wrhat)                  
dlmwrite(strcat(workdir,'/Results/Smarket1000_secB/what.txt'),what)                  
end

if (reps>1)
dlmwrite(strcat(workdir,'/Results/Smarket1000_secB/rva.txt'),RVA,'-append')                
dlmwrite(strcat(workdir,'/Results/Smarket1000_secB/rvap.txt'),RVAp,'-append')                
dlmwrite(strcat(workdir,'/Results/Smarket1000_secB/exf.txt'),EXf,'-append')                
dlmwrite(strcat(workdir,'/Results/Smarket1000_secB/exfp.txt'),EXfp,'-append')                
dlmwrite(strcat(workdir,'/Results/Smarket1000_secB/wr.txt'),wrhat,'-append')                  
dlmwrite(strcat(workdir,'/Results/Smarket1000_secB/what.txt'),what,'-append')                  
end





%% Other RTA

taup=taupoRTA;
DDD_oRTA=D_oRTA.*fta_hat_oRTA;
DDD_eukor=D_eukor.*fta_hat_eukor;

tau_hat=taup./tau .* exp(DDD_oRTA).*exp(DDD_eukor); 

%load starting values from scenario with parameter means
load(strcat(workdir,'/Results/oRTA_smdB/oRTA_smd'), 'pf0_all', 'wf0_all')                 
wstart = wf0_all;
pstart = pf0_all;
clear wf0_all pf0_all

[wf0_all pf0_all PQ_all Fp_all Dinp_all emax] = equi_s_emax(tau_hat,taup,alphas,T,B,G,sd,Din,J,N,maxit,tol,VAn,Iv,VP,vfactor,wstart,pstart); %add wstart for use with equilibrium_LC
% if no convergence, choose smaller factor
emax2=[];  maxit2    = 200;
if (emax == maxit)
    vf2=-.1;
    tol2      = 1E-8;
  
[wf0_all pf0_all PQ_all Fp_all Dinp_all emax2] = equi_s_emax(tau_hat,taup,alphas,T,B,G,sd,Din,J,N,maxit2,tol2,VAn,Iv,VP,vf2,wstart,pstart); %add wstart for use with equilibrium_LC
end


% if no convergence, choose smaller factor
if (emax2 == maxit2)
    vf2=-.01;
    tol2      = 1E-9;
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
wrhat=wf0_all./P0_all;

% exports - same structure as trade shares
%           exp1 exp2
%sec1 imp1
%     imp2
%sec2 imp1

EXf = (Din./tau).*reshape(X0',[N*J,1]);
EXfp = (Dinp_all./taup).*reshape(PQ_all',[N*J,1]);

% value added exports
% final demand
C = alphas.*I'; C=reshape(C',[N*J,1]);
Cp = alphas.*Ip'; Cp=reshape(Cp',[N*J,1]);
r = reshape(B',[N*J,1]);
ahat=Dinp_all./Din.*tau./taup;
for j=1:J
    Ares((j-1)*N+1:j*N,:)=ahat((j-1)*N+1:j*N,1:N)';
end
Ares(isnan(Ares))=1;
AAp=kron(ones(1,J),Ares).*AA;
FEXf=(Din./tau).*C;
FEXfp=(Dinp_all./taup).*Cp;

for j=1:J
    RFEXf((j-1)*N+1:j*N,:)=FEXf((j-1)*N+1:j*N,1:N)';
    RFEXfp((j-1)*N+1:j*N,:)=FEXfp((j-1)*N+1:j*N,1:N)';
end
RVA = diag(r)*inv(eye(size(AA,1))-AA)*RFEXf;
RVAp = diag(r)*inv(eye(size(AAp,1))-AAp)*RFEXfp;

clear DDD Dinp_all Fp_all PQ_all pf0_all wf0_all AAp Ares ahat

if (reps==1)
dlmwrite(strcat(workdir,'/Results/oRTA1000_secB/rva.txt'),RVA)                
dlmwrite(strcat(workdir,'/Results/oRTA1000_secB/rvap.txt'),RVAp)                
dlmwrite(strcat(workdir,'/Results/oRTA1000_secB/exf.txt'),EXf)                
dlmwrite(strcat(workdir,'/Results/oRTA1000_secB/exfp.txt'),EXfp)                
dlmwrite(strcat(workdir,'/Results/oRTA1000_secB/wr.txt'),wrhat)                  
dlmwrite(strcat(workdir,'/Results/oRTA1000_secB/what.txt'),what)                  
end

if (reps>1)
dlmwrite(strcat(workdir,'/Results/oRTA1000_secB/rva.txt'),RVA,'-append')                
dlmwrite(strcat(workdir,'/Results/oRTA1000_secB/rvap.txt'),RVAp,'-append')                
dlmwrite(strcat(workdir,'/Results/oRTA1000_secB/exf.txt'),EXf,'-append')                
dlmwrite(strcat(workdir,'/Results/oRTA1000_secB/exfp.txt'),EXfp,'-append')                
dlmwrite(strcat(workdir,'/Results/oRTA1000_secB/wr.txt'),wrhat,'-append')                  
dlmwrite(strcat(workdir,'/Results/oRTA1000_secB/what.txt'),what,'-append')                  
end









%% Customs Union


tau_hat=taupEU./tau;  
taup=taupEU;

%load starting values from scenario with parameter means
load(strcat(workdir,'/Results/Cunion_smdB/Cunion_smd'), 'pf0_all', 'wf0_all')                 
wstart = wf0_all;
pstart = pf0_all;
clear wf0_all pf0_all

[wf0_all pf0_all PQ_all Fp_all Dinp_all emax] = equi_s_emax(tau_hat,taup,alphas,T,B,G,sd,Din,J,N,maxit,tol,VAn,Iv,VP,vfactor,wstart,pstart); %add wstart for use with equilibrium_LC
% if no convergence, choose smaller factor
emax2=[];  maxit2    = 200;
if (emax == maxit)
    vf2=-.1;
    tol2      = 1E-8;
  
[wf0_all pf0_all PQ_all Fp_all Dinp_all emax2] = equi_s_emax(tau_hat,taup,alphas,T,B,G,sd,Din,J,N,maxit2,tol2,VAn,Iv,VP,vf2,wstart,pstart); %add wstart for use with equilibrium_LC
end


% if no convergence, choose smaller factor
if (emax2 == maxit2)
    vf2=-.01;
    tol2      = 1E-9;
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
wrhat=wf0_all./P0_all;

% exports - same structure as trade shares
%           exp1 exp2
%sec1 imp1
%     imp2
%sec2 imp1

EXf = (Din./tau).*reshape(X0',[N*J,1]);
EXfp = (Dinp_all./taup).*reshape(PQ_all',[N*J,1]);

% value added exports
% final demand
C = alphas.*I'; C=reshape(C',[N*J,1]);
Cp = alphas.*Ip'; Cp=reshape(Cp',[N*J,1]);
r = reshape(B',[N*J,1]);
ahat=Dinp_all./Din.*tau./taup;
for j=1:J
    Ares((j-1)*N+1:j*N,:)=ahat((j-1)*N+1:j*N,1:N)';
end
Ares(isnan(Ares))=1;
AAp=kron(ones(1,J),Ares).*AA;
FEXf=(Din./tau).*C;
FEXfp=(Dinp_all./taup).*Cp;

for j=1:J
    RFEXf((j-1)*N+1:j*N,:)=FEXf((j-1)*N+1:j*N,1:N)';
    RFEXfp((j-1)*N+1:j*N,:)=FEXfp((j-1)*N+1:j*N,1:N)';
end
RVA = diag(r)*inv(eye(size(AA,1))-AA)*RFEXf;
RVAp = diag(r)*inv(eye(size(AAp,1))-AAp)*RFEXfp;

clear DDD Dinp_all Fp_all PQ_all pf0_all wf0_all AAp Ares ahat

if (reps==1)
dlmwrite(strcat(workdir,'/Results/Cunion1000_secB/rva.txt'),RVA)                
dlmwrite(strcat(workdir,'/Results/Cunion1000_secB/rvap.txt'),RVAp)                
dlmwrite(strcat(workdir,'/Results/Cunion1000_secB/exf.txt'),EXf)                
dlmwrite(strcat(workdir,'/Results/Cunion1000_secB/exfp.txt'),EXfp)                
dlmwrite(strcat(workdir,'/Results/Cunion1000_secB/wr.txt'),wrhat)                  
dlmwrite(strcat(workdir,'/Results/Cunion1000_secB/what.txt'),what)                  
end

if (reps>1)
dlmwrite(strcat(workdir,'/Results/Cunion1000_secB/rva.txt'),RVA,'-append')                
dlmwrite(strcat(workdir,'/Results/Cunion1000_secB/rvap.txt'),RVAp,'-append')                
dlmwrite(strcat(workdir,'/Results/Cunion1000_secB/exf.txt'),EXf,'-append')                
dlmwrite(strcat(workdir,'/Results/Cunion1000_secB/exfp.txt'),EXfp,'-append')                
dlmwrite(strcat(workdir,'/Results/Cunion1000_secB/wr.txt'),wrhat,'-append')                  
dlmwrite(strcat(workdir,'/Results/Cunion1000_secB/what.txt'),what,'-append')                  
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%      allEU w transfers;   SECTORAL ESTIMATIONS 

S_prime = importdata(strcat(workdir,'/input/S_n_prime.txt'));
sd_prime=sd.*S_prime./Sn;


DDD_oRTA=D_oRTA.*fta_hat_oRTA;
DDD_eukor=D_eukor.*fta_hat_eukor;

DDD_schengen=D_schengen.*fta_hat_schengen;
DDD_euro=D_euro.*fta_hat_euro;
DDD_Smarket=D_Smarket.*fta_hat_Smarket;
taup=taupallEU;
tau_hat=taup./tau .* exp(DDD_oRTA).* exp(DDD_eukor).* exp(DDD_schengen).* exp(DDD_euro).* exp(DDD_Smarket);  


%load starting values from scenario with parameter means
load(strcat(workdir,'/Results/allEUtransfer_smdB/allEUtransfer_smd'), 'pf0_all', 'wf0_all')                 
wstart = wf0_all;
pstart = pf0_all;
clear wf0_all pf0_all

[wf0_all pf0_all PQ_all Fp_all Dinp_all emax] = equi_s_emax(tau_hat,taup,alphas,T,B,G,sd_prime,Din,J,N,maxit,tol,VAn,Iv,VP,vfactor,wstart,pstart); %add wstart for use with equilibrium_LC
% if no convergence, choose smaller factor
emax2=[];  maxit2    = 200;
if (emax == maxit)
    vf2=-.1;
    tol2      = 1E-8;
  
[wf0_all pf0_all PQ_all Fp_all Dinp_all emax2] = equi_s_emax(tau_hat,taup,alphas,T,B,G,sd_prime,Din,J,N,maxit2,tol2,VAn,Iv,VP,vf2,wstart,pstart); %add wstart for use with equilibrium_LC
end


% if no convergence, choose smaller factor
if (emax2 == maxit2)
    vf2=-.01;
    tol2      = 1E-9;
    maxit2    = 500;
[wf0_all pf0_all PQ_all Fp_all Dinp_all emax] = equi_s_emax(tau_hat,taup,alphas,T,B,G,sd_prime,Din,J,N,maxit2,tol2,VAn,Iv,VP,vf2,wstart,pstart); %add wstart for use with equilibrium_LC
end

% results output

I=(VAn+(X0'*(1-F)).*eye(N,N)*ones(N,1)).*(1-sd);
Ip=(diag(wf0_all)*VAn+(PQ_all'*(1-Fp_all)).*eye(N,N)*ones(N,1)).*(1-sd_prime);

% income change
Ihat=Ip./I;

% change in price index
P=pf0_all.^alphas;
P0_all=prod(P,1);
P0_all=P0_all';
what=Ihat./P0_all;
wrhat=wf0_all./P0_all;

% exports - same structure as trade shares
%           exp1 exp2
%sec1 imp1
%     imp2
%sec2 imp1

EXf = (Din./tau).*reshape(X0',[N*J,1]);
EXfp = (Dinp_all./taup).*reshape(PQ_all',[N*J,1]);

% value added exports
% final demand
C = alphas.*I'; C=reshape(C',[N*J,1]);
Cp = alphas.*Ip'; Cp=reshape(Cp',[N*J,1]);
r = reshape(B',[N*J,1]);
ahat=Dinp_all./Din.*tau./taup;
for j=1:J
    Ares((j-1)*N+1:j*N,:)=ahat((j-1)*N+1:j*N,1:N)';
end
Ares(isnan(Ares))=1;
AAp=kron(ones(1,J),Ares).*AA;
FEXf=(Din./tau).*C;
FEXfp=(Dinp_all./taup).*Cp;

for j=1:J
    RFEXf((j-1)*N+1:j*N,:)=FEXf((j-1)*N+1:j*N,1:N)';
    RFEXfp((j-1)*N+1:j*N,:)=FEXfp((j-1)*N+1:j*N,1:N)';
end
RVA = diag(r)*inv(eye(size(AA,1))-AA)*RFEXf;
RVAp = diag(r)*inv(eye(size(AAp,1))-AAp)*RFEXfp;

clear DDD Dinp_all Fp_all PQ_all pf0_all wf0_all AAp Ares ahat

if (reps==1)
dlmwrite(strcat(workdir,'/Results/allEUtransfer1000_secB/rva.txt'),RVA)                
dlmwrite(strcat(workdir,'/Results/allEUtransfer1000_secB/rvap.txt'),RVAp)                
dlmwrite(strcat(workdir,'/Results/allEUtransfer1000_secB/exf.txt'),EXf)                
dlmwrite(strcat(workdir,'/Results/allEUtransfer1000_secB/exfp.txt'),EXfp)                
dlmwrite(strcat(workdir,'/Results/allEUtransfer1000_secB/wr.txt'),wrhat)                  
dlmwrite(strcat(workdir,'/Results/allEUtransfer1000_secB/what.txt'),what)                  
end

if (reps>1)
dlmwrite(strcat(workdir,'/Results/allEUtransfer1000_secB/RVA.txt'),RVA,'-append')                
dlmwrite(strcat(workdir,'/Results/allEUtransfer1000_secB/RVAp.txt'),RVAp,'-append')                
dlmwrite(strcat(workdir,'/Results/allEUtransfer1000_secB/exf.txt'),EXf,'-append')                
dlmwrite(strcat(workdir,'/Results/allEUtransfer1000_secB/exfp.txt'),EXfp,'-append')                
dlmwrite(strcat(workdir,'/Results/allEUtransfer1000_secB/wr.txt'),wrhat,'-append')                  
dlmwrite(strcat(workdir,'/Results/allEUtransfer1000_secB/what.txt'),what,'-append')                 
end



end



