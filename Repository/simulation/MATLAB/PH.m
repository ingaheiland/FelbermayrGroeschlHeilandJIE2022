 function [pf0 c]= PH(wages_N,tau_hat,T,B,G,Din,J,N,maxit,tol)

% reformating theta vector 
for   j    = 1:1:J;
      idx  = 1+(j-1)*N:1:N*j;
 LT(idx,1) = ones(N,1)*T(j);
end


% initialize vectors of ex-post wage and price factors 
wf0      =  wages_N;
pf0      =  ones(J,N);

pfmax = 1;
 it       = 1;     

while (it <= maxit) && (pfmax > tol)
lw       =  log(wf0);
lp       =  log(pf0);

% calculating log cost
for          i = 1:1:N;
       lc(:,i) = B(:,i)*lw(i) + (1-B(:,i)).*(G(1+(i-1)*J:J*i,:)'*lp(:,i));
end
             c = exp(lc);
        Din_om = Din.*( tau_hat.^(-1./(LT*ones(1,N))));

% calculating phat
for          j = 1:1:J
         for n = 1:1:N
    phat(j , n) = Din_om( n + (j-1)*N , :)*( c( j , :).^( - 1 / T( j )) )';

% this happens because of sectors with zero VA
% Note that we do not want logs of zero
    if    phat(j , n) == 0
          phat(j , n) = 1;
    else
          phat(j , n) = phat(j , n)^( - T( j ) );
    end
          end
end
      pfdev    = abs(phat - pf0); % Checking tolerance
      pf0      = phat;
      pfmax    = max(max(pfdev));
      it       = it + 1;
end
