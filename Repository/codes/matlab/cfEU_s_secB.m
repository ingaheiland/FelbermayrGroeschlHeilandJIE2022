% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Main Counterfactuals: AT MEANS OF PARAMETER ESTIMATES   %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clear all
clc


workdir='Repository/simulation';

cd(strcat(workdir,'/MATLAB'))

vfactor  = -.1;
tol      = 1E-13;
maxit    = 1000;

load initial_condition_2014



% calculate surplus as share of income.
% ignore inventories, as those will be set to zero

R=(X0'*(1-F)).*eye(N,N)*ones(N,1);
sd = Sn./(VAn+R);


% inventories
Iv = SIn-Sn;


wstart      =  ones(N,1); 
pstart      =  ones(J,N);

T=importdata(strcat(workdir,'/input','/bltau_smdB.txt'));
T=-T-ones(J,1); % manufacturing thetas estimated with fob data, Egger et al's elasticity for sectors is also FOB
T = 1./T; % code is for -1/(cif trade elasticity) 





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Placebo test %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

taup=tau;
tau_hat=taup./tau;

[wf0_all pf0_all PQ_all Fp_all Dinp_all] = equi_s(tau_hat,taup,alphas,T,B,G,sd,Din,J,N,maxit,tol,VAn,Iv,VP,vfactor,wstart,pstart); %add wstart for use with equilibrium_LC


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  generate inventory free equilibrium with new algorithm %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

VP=zeros(J,N);
SIn = Sn;
Iv = zeros(N,1);
taup=tau;
tau_hat=taup./tau;


[wf0_all pf0_all PQ_all Fp_all Dinp_all] = equi_s(tau_hat,taup,alphas,T,B,G,sd,Din,J,N,maxit,tol,VAn,Iv,VP,vfactor,wstart,pstart); %add wstart for use with equilibrium_LC

%save(strcat(workdir,'/Results/NewBaseline_smd/IvZero_secB'))

%load(strcat(workdir,'/Results/NewBaseline_smd/IvZero_secB'))


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

I=VAn+(X0'*(1-F)).*eye(N,N)*ones(N,1).*(1-sd);
%dlmwrite(strcat(workdir,'/Results/NewBaseline_incomeB.txt'),I)


clear Dinp_all Fp_all pf0_all PQ_all SIn wf0_all T Ares pihat A
save(strcat(workdir,'/Results/NewBaseline_smd/NewInitialConditionsSecB'))




load(strcat(workdir,'/Results/NewBaseline_smd/NewInitialConditionsSecB'))

fta_hat_oRTA=importdata(strcat(workdir,'/input/dummy_diff_oRTA.txt')); 
fta_hat_eukor=importdata(strcat(workdir,'/input/dummy_diff_eukor.txt')); 
fta_hat_Smarket=importdata(strcat(workdir,'/input/dummy_diff_Smarket.txt')); 
fta_hat_schengen=importdata(strcat(workdir,'/input/dummy_diff_Schengen.txt')); % #borders_old - # borders_new
fta_hat_euro=importdata(strcat(workdir,'/input/dummy_diff_euro.txt')); 

%counterfactual tariffs
taupEU=importdata(strcat(workdir,'/input/taup_eu2014.txt')); 
taupoRTA=importdata(strcat(workdir,'/input/taup_oRTA2014.txt')); 
taupallEU=importdata(strcat(workdir,'/input/taup_allEU2014.txt')); 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%             SECTORAL ESTIMATIONS 
%%%%%%%%%                   at means of 1000 draws
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

T=importdata(strcat(workdir,'/input','/bltau_smdB.txt'));
T = -ones(J,1)./T; 

%schengen
delta_schengen = importdata(strcat(workdir,'/input/bschengen_smdB.txt')); % read in dummy coefficients (theoretically -delta/theta); delta will be
delta_schengen=-delta_schengen.*T;
D_schengen=[];
for j=1:J
    for i=1:N
        D_schengen = [D_schengen ; delta_schengen(j)];
    end
end       

%euro
delta_euro = importdata(strcat(workdir,'/input/bbotheuro_smdB.txt')); % read in dummy coefficients (theoretically -delta/theta); delta will be
delta_euro=-delta_euro.*T;
D_euro=[];
for j=1:J
    for i=1:N
        D_euro = [D_euro ; delta_euro(j)];
    end
end       


%single market
delta_Smarket = importdata(strcat(workdir,'/input/bbotheea_smdB.txt')); % read in dummy coefficients (theoretically -delta/theta); delta will be
delta_Smarket=-delta_Smarket.*T;
D_Smarket=[];
for j=1:J
    for i=1:N
        D_Smarket = [D_Smarket ; delta_Smarket(j)];
    end
end       



%other RTAs
delta_oRTA = importdata(strcat(workdir,'/input/bbothrta_smdB.txt')); % read in dummy coefficients (theoretically -delta/theta); delta will be
delta_oRTA=-delta_oRTA.*T;
D_oRTA=[];
for j=1:J
    for i=1:N
        D_oRTA = [D_oRTA ; delta_oRTA(j)];
    end
end       


%other RTAs
delta_eukor = importdata(strcat(workdir,'/input/beukor_smdB.txt')); % read in dummy coefficients (theoretically -delta/theta); delta will be
delta_eukor=-delta_eukor.*T;
D_eukor=[];
for j=1:J
    for i=1:N
        D_eukor = [D_eukor ; delta_eukor(j)];
    end
end       

% convert to cif elasticity
T=1./T;
T=T-1;
T = ones(J,1)./T; %code is written for 1/elasticity


     






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%            allEU;   SECTORAL ESTIMATIONS 


DDD_oRTA=D_oRTA.*fta_hat_oRTA;
DDD_eukor=D_eukor.*fta_hat_eukor;
DDD_schengen=D_schengen.*fta_hat_schengen;
DDD_euro=D_euro.*fta_hat_euro;
DDD_Smarket=D_Smarket.*fta_hat_Smarket;
taup=taupallEU;
tau_hat=taup./tau .* exp(DDD_oRTA).* exp(DDD_eukor).* exp(DDD_schengen).* exp(DDD_euro).* exp(DDD_Smarket);  


[wf0_all pf0_all PQ_all Fp_all Dinp_all] = equi_s(tau_hat,taup,alphas,T,B,G,sd,Din,J,N,maxit,tol,VAn,Iv,VP,vfactor,wstart,pstart); %add wstart for use with equilibrium_LC


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

dlmwrite(strcat(workdir,'/Results/allEU_smdB/cf_rva.txt'),RVA)
dlmwrite(strcat(workdir,'/Results/allEU_smdB/cf_rvap.txt'),RVAp)
dlmwrite(strcat(workdir,'/Results/allEU_smdB/cf_what.txt'),what)
dlmwrite(strcat(workdir,'/Results/allEU_smdB/cf_wrhat.txt'),wrhat)
dlmwrite(strcat(workdir,'/Results/allEU_smdB/cf_EXf.txt'),EXf)
dlmwrite(strcat(workdir,'/Results/allEU_smdB/cf_EXfp.txt'),EXfp)


save(strcat(workdir,'/Results/allEU_smdB/allEU_smd'))


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% SCHENGEN; SECTORAL ESTIMATIONS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


taup=tau;
DDD=D_schengen.*fta_hat_schengen;
tau_hat=taup./tau .* exp(DDD);  


[wf0_all pf0_all PQ_all Fp_all Dinp_all] = equi_s(tau_hat,taup,alphas,T,B,G,sd,Din,J,N,maxit,tol,VAn,Iv,VP,vfactor,wstart,pstart); %add wstart for use with equilibrium_LC


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


dlmwrite(strcat(workdir,'/Results/Schengen_smdB/cf_what.txt'),what)
dlmwrite(strcat(workdir,'/Results/Schengen_smdB/cf_wrhat.txt'),wrhat)
dlmwrite(strcat(workdir,'/Results/Schengen_smdB/cf_EXf.txt'),EXf)
dlmwrite(strcat(workdir,'/Results/Schengen_smdB/cf_EXfp.txt'),EXfp)
dlmwrite(strcat(workdir,'/Results/Schengen_smdB/cf_rva.txt'),RVA)
dlmwrite(strcat(workdir,'/Results/Schengen_smdB/cf_rvap.txt'),RVAp)
save(strcat(workdir,'/Results/Schengen_smdB/Schengen_smd'))


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%            EURO; SECTORAL ESTIMATIONS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


taup=tau;

DDD=D_euro.*fta_hat_euro;
tau_hat=taup./tau .* exp(DDD);  


[wf0_all pf0_all PQ_all Fp_all Dinp_all] = equi_s(tau_hat,taup,alphas,T,B,G,sd,Din,J,N,maxit,tol,VAn,Iv,VP,vfactor,wstart,pstart); %add wstart for use with equilibrium_LC



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

dlmwrite(strcat(workdir,'/Results/Euro_smdB/cf_what.txt'),what)
dlmwrite(strcat(workdir,'/Results/Euro_smdB/cf_wrhat.txt'),wrhat)
dlmwrite(strcat(workdir,'/Results/Euro_smdB/cf_EXf.txt'),EXf)
dlmwrite(strcat(workdir,'/Results/Euro_smdB/cf_EXfp.txt'),EXfp)
dlmwrite(strcat(workdir,'/Results/Euro_smdB/cf_rva.txt'),RVA)
dlmwrite(strcat(workdir,'/Results/Euro_smdB/cf_rvap.txt'),RVAp)
save(strcat(workdir,'/Results/Euro_smdB/Euro_smd'))




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%            SINGLE MARKET;   SECTORAL ESTIMATIONS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

taup=tau;
DDD=D_Smarket.*fta_hat_Smarket;
tau_hat=taup./tau .* exp(DDD);  


[wf0_all pf0_all PQ_all Fp_all Dinp_all] = equi_s(tau_hat,taup,alphas,T,B,G,sd,Din,J,N,maxit,tol,VAn,Iv,VP,vfactor,wstart,pstart); %add wstart for use with equilibrium_LC



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


dlmwrite(strcat(workdir,'/Results/Smarket_smdB/cf_what.txt'),what)
dlmwrite(strcat(workdir,'/Results/Smarket_smdB/cf_wrhat.txt'),wrhat)
dlmwrite(strcat(workdir,'/Results/Smarket_smdB/cf_EXf.txt'),EXf)
dlmwrite(strcat(workdir,'/Results/Smarket_smdB/cf_EXfp.txt'),EXfp)
dlmwrite(strcat(workdir,'/Results/Smarket_smdB/cf_rva.txt'),RVA)
dlmwrite(strcat(workdir,'/Results/Smarket_smdB/cf_rvap.txt'),RVAp)
save(strcat(workdir,'/Results/Smarket_smdB/Smarket_smd'))



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%            other RTAs;   SECTORAL ESTIMATIONS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

taup=taupoRTA;
DDD_oRTA=D_oRTA.*fta_hat_oRTA;
DDD_eukor=D_eukor.*fta_hat_eukor;

tau_hat=taup./tau .* exp(DDD_oRTA).*exp(DDD_eukor);  


[wf0_all pf0_all PQ_all Fp_all Dinp_all] = equi_s(tau_hat,taup,alphas,T,B,G,sd,Din,J,N,maxit,tol,VAn,Iv,VP,vfactor,wstart,pstart); %add wstart for use with equilibrium_LC



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

dlmwrite(strcat(workdir,'/Results/oRTA_smdB/cf_what.txt'),what)
dlmwrite(strcat(workdir,'/Results/oRTA_smdB/cf_wrhat.txt'),wrhat)
dlmwrite(strcat(workdir,'/Results/oRTA_smdB/cf_EXf.txt'),EXf)
dlmwrite(strcat(workdir,'/Results/oRTA_smdB/cf_EXfp.txt'),EXfp)
dlmwrite(strcat(workdir,'/Results/oRTA_smdB/cf_rva.txt'),RVA)
dlmwrite(strcat(workdir,'/Results/oRTA_smdB/cf_rvap.txt'),RVAp)
save(strcat(workdir,'/Results/oRTA_smdB/oRTA_smd'))




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%            TARIFFS;   SECTORAL ESTIMATIONS 


tau_hat=taupEU./tau;  
taup=taupEU;

[wf0_all pf0_all PQ_all Fp_all Dinp_all] = equi_s(tau_hat,taup,alphas,T,B,G,sd,Din,J,N,maxit,tol,VAn,Iv,VP,vfactor,wstart,pstart); %add wstart for use with equilibrium_LC


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



dlmwrite(strcat(workdir,'/Results/Cunion_smdB/cf_what.txt'),what)
dlmwrite(strcat(workdir,'/Results/Cunion_smdB/cf_wrhat.txt'),wrhat)
dlmwrite(strcat(workdir,'/Results/Cunion_smdB/cf_EXf.txt'),EXf)
dlmwrite(strcat(workdir,'/Results/Cunion_smdB/cf_EXfp.txt'),EXfp)
dlmwrite(strcat(workdir,'/Results/Cunion_smdB/cf_rva.txt'),RVA)
dlmwrite(strcat(workdir,'/Results/Cunion_smdB/cf_rvap.txt'),RVAp)
save(strcat(workdir,'/Results/Cunion_smdB/CUnion_smd'))








%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%      allEU w transfers;   SECTORAL ESTIMATIONS 

S_prime = importdata(strcat(workdir,'/input/S_n_prime.txt'));
sd_prime=sd.*S_prime./Sn;
SD=horzcat(sd,sd_prime);
dlmwrite(strcat(workdir,'/Results/allEUtransfer_smdB/SD.txt'),SD)


DDD_oRTA=D_oRTA.*fta_hat_oRTA;
DDD_eukor=D_eukor.*fta_hat_eukor;

DDD_schengen=D_schengen.*fta_hat_schengen;
DDD_euro=D_euro.*fta_hat_euro;
DDD_Smarket=D_Smarket.*fta_hat_Smarket;
taup=taupallEU;
tau_hat=taup./tau .* exp(DDD_oRTA).* exp(DDD_eukor).* exp(DDD_schengen).* exp(DDD_euro).* exp(DDD_Smarket);  


[wf0_all pf0_all PQ_all Fp_all Dinp_all] = equi_s(tau_hat,taup,alphas,T,B,G,sd_prime,Din,J,N,maxit,tol,VAn,Iv,VP,vfactor,wstart,pstart); %add wstart for use with equilibrium_LC




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


dlmwrite(strcat(workdir,'/Results/allEUtransfer_smdB/cf_what.txt'),what)
dlmwrite(strcat(workdir,'/Results/allEUtransfer_smdB/cf_wrhat.txt'),wrhat)
dlmwrite(strcat(workdir,'/Results/allEUtransfer_smdB/cf_EXf.txt'),EXf)
dlmwrite(strcat(workdir,'/Results/allEUtransfer_smdB/cf_EXfp.txt'),EXfp)
dlmwrite(strcat(workdir,'/Results/allEUtransfer_smdB/cf_rva.txt'),RVA)
dlmwrite(strcat(workdir,'/Results/allEUtransfer_smdB/cf_rvap.txt'),RVAp)
save(strcat(workdir,'/Results/allEUtransfer_smdB/allEUtransfer_smd'))




