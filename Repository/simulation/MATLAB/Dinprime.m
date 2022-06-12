function Dinp = Dinprime(Din,tau_hat,c,T,J,N)

% reformating theta vector
for   j    = 1:1:J;
      idx  = 1+(j-1)*N:1:N*j;
 LT(idx,1) = ones(N,1)*T(j);
end

for n=1:1:N
    cp(:,n) = c(:,n).^( - 1./T );
end

 Din_om = Din.*( tau_hat.^(-1./(LT*ones(1,N))));
  
 for n = 1:1:N
 idx = n:N:length(Din)-(N-n);
 DD(idx,:) = Din_om(idx,:).*cp;
 end
 
 phat = sum(DD')'.^-LT;
 
 for n = 1:1:N
 Dinp(:,n) = DD(:,n).*(phat.^(1./LT));
 end
 
 
 