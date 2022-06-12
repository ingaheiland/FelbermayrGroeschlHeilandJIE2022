% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Brexit & Post-Brexit EU collapse: AT MEANS OF SECTORAL PARAMETERS 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



close all
clear all
clc



workdir='Repository/simulation';

cd(strcat(workdir,'/MATLAB'))

vfactor  = -.1;
tol      = 1E-13;
maxit    = 1000;

load(strcat(workdir,'/Results/NewBaseline_smd/NewInitialConditionsSecB'))


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


     



%%%%%%%%%%%%%%%%%%% GENERATE NEW BENCHMARK %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%         hard Brexit   SECTORAL ESTIMATIONS 


% BREXIT

fta_hat_oRTA_brexit=importdata('../input/brexit_dummy_diff_oRTA.txt'); 
fta_hat_eukor_brexit=importdata('../input/brexit_dummy_diff_eukor.txt'); 
fta_hat_Smarket_brexit=importdata('../input/brexit_dummy_diff_Smarket.txt'); 
fta_hat_schengen_brexit=importdata('../input/brexit_dummy_diff_Schengen.txt'); % #borders_old - # borders_new
fta_hat_euro_brexit=importdata('../input/brexit_dummy_diff_euro.txt'); 

%counterfactual tariffs
taupBR=importdata('../input/taup_brexit2014.txt'); 




DDD_oRTA=D_oRTA.*fta_hat_oRTA_brexit;
DDD_eukor=D_eukor.*fta_hat_eukor_brexit;
DDD_schengen=D_schengen.*fta_hat_schengen_brexit;
DDD_euro=D_euro.*fta_hat_euro_brexit;
DDD_Smarket=D_Smarket.*fta_hat_Smarket_brexit;
taup=taupBR;
tau_hat=taupBR./tau .* exp(DDD_oRTA).*exp(DDD_eukor).* exp(DDD_schengen).* exp(DDD_euro).* exp(DDD_Smarket);  


[wf0_all pf0_all PQ_all Fp_all Dinp_all] = equi_s(tau_hat,taup,alphas,T,B,G,sd,Din,J,N,maxit,tol,VAn,Iv,VP,vfactor,wstart,pstart); %add wstart for use with equilibrium_LC


I=VAn+(X0'*(1-F)).*eye(N,N)*ones(N,1).*(1-sd);
Ip=diag(wf0_all)*VAn+(PQ_all'*(1-Fp_all)).*eye(N,N)*ones(N,1).*(1-sd);
% income change
Ihat=Ip./I;
% change in price index
P=pf0_all.^alphas;
P0_all=prod(P,1);
P0_all=P0_all';
what=Ihat./P0_all;

dlmwrite(strcat(workdir,'/Results/Brexit_smdB/what.txt'),what)


% new baseline
Din = Dinp_all;
F = Fp_all;
X0 = PQ_all;
VAn=wf0_all.*VAn;
Inv = zeros(J,N);
tau = taup;


clear Dinp_all Fp_all pf0_all PQ_all SIn wf0_all
save(strcat(workdir,'/Results/Brexit_smdB/pbBaselineB'))








%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% post brexit scenarios
%%%%%%%%%%%%%%%%%%%%%%


%tau=importdata('tariffs2014.txt');                       % tariffs 2014
fta_hat_oRTA=importdata(strcat(workdir,'/input/pb_dummy_diff_oRTA.txt')); 
fta_hat_eukor=importdata(strcat(workdir,'/input/pb_dummy_diff_eukor.txt')); 
fta_hat_Smarket=importdata(strcat(workdir,'/input/pb_dummy_diff_Smarket.txt')); 
fta_hat_schengen=importdata(strcat(workdir,'/input/pb_dummy_diff_Schengen.txt')); % #borders_old - # borders_new
fta_hat_euro=importdata(strcat(workdir,'/input/pb_dummy_diff_euro.txt')); 

%counterfactual tariffs
taupEU=importdata(strcat(workdir,'/input/taup_eu2014.txt')); 
taupoRTA=importdata(strcat(workdir,'/input/taup_oRTA2014.txt')); 
taupallEU=importdata(strcat(workdir,'/input/taup_allEU2014.txt')); 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%             SECTORAL ESTIMATIONS 
%%%%%%%%%                   at means of 1000 draws
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



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


dlmwrite(strcat(workdir,'/Results/allEU_smdB/cf_what_pb.txt'),what)
save(strcat(workdir,'/Results/allEU_smdB/allEU_smd_pb'))










