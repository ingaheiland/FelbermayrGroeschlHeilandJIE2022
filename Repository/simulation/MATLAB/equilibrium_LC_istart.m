% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%  Main Program Counterfactuals
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [wf0 pf0 PQ Fp Dinp krit] = equilibrium_LC_istart(tau_hat,taup,alphas,T,B,G,Din,J,N,maxit,tol,VAn,SIn,Sn,VP,vfactor,wstart,pstart) %note that originally ***SIn*** was not in this place and VP is also new

% user-supplied starting values
 wf0      =  wstart;  pf0      =  pstart;

%wf0      =  ones(N,1); 
%pf0      =  ones(J,N);

krit=1;
 
 
 wfmax = 1; e = 1;
while (e <= maxit) && (wfmax > tol);

[pf0 c] = PH(wf0,tau_hat,T,B,G,Din,J,N,maxit,tol);
       
% Calculating trade shares
Dinp = Dinprime(Din,tau_hat,c,T,J,N);
Dinp_om = Dinp./taup;

for j   = 1:1:J
irow    = 1+N*(j-1):1:N*j;
Fp(j,:) = sum((Dinp(irow,:)./taup(irow,:))');
end
 
% % Expenditure MATRIX
PQ = expenditure_i(alphas,B,G,Dinp,taup,Fp,VAn,wf0,SIn,VP,J,N); %note that originally there was ***Sn*** in this place
% expenditures Xji in long vector: PQ_vec=(X11 X12 X13...)' 
PQ_vec   = reshape(PQ',1,J*N)'; 

for n = 1:1:N
    DP(:,n)  = Dinp_om(:,n).*PQ_vec; 
end
LHS = sum(DP)'; %exports

% calculating RHS (Imports) trade balance
PF = PQ.*Fp;
RHS = sum(PF)'; %imports

% excess function (trade balance)
ZW2 = -(RHS - LHS + Sn)./(abs(VAn)); % originally we had VA to scale down the error and Sn in place of SIn
%Iteration factor prices
wf1 =wf0.*(1-vfactor*ZW2./wf0);

wfmax=sum(abs(wf1-wf0));

disp(['wage tolerance = ' num2str(sum(abs(wf1-wf0)))])
%disp(wf0(7,1))
%disp(wf0(26,1))
wfmax0 = wfmax;
wf0=(wf1);

%vector that collects the corresponding tolerance (to find the minimum
%later)
%krit=horzcat(krit,wfmax);
krit=vertcat(e,wfmax);
disp(e);

e=e+1;
end

