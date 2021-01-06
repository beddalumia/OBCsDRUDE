function [w,res,nMax,mMax] = obcKubo(E,psi,nF,a,L,dx,step,basisDIM,cutoff)
    
    % global a                       % Physical parameters
    % global dx step basisDIM cutoff % Computational parameters
    % -> global variables are deprecated within parfor loops
    
  %% First of all, we apply the velocity operator to the eigenstates
   % $\hat{v}\psi_n(x) = i \frac{\mathrm{d}}{\mathrm{d}x} \psi_n(x)$
     x = linspace(0,L*a,step*L*a);
     vpsi = zeros(length(x),basisDIM);
     dpsi_dx = zeros(1,length(x));
     for n = 1:basisDIM
         dpsi_dx(1:(length(x)-1)) = diff(psi(:,n))/dx;
         dpsi_dx(length(x)) = dpsi_dx(length(x)-1); % Border Trick
         vpsi(:,n) = 1i * dpsi_dx;
     end

  %% Journey through all possible* excited (many-body) states
   % *compatible with Pauli principle
     n = 1;
     w = zeros(basisDIM,basisDIM);
     res = zeros(basisDIM,basisDIM);
     while 1
         if n == nF+1
             break
         end
         m = nF+1;
         while 1
             dw = E(m)-E(n);
             if dw > cutoff || m > basisDIM
                 break
             end
             w(n,m) = dw; % Pole frequency
             braket_nm = conj(psi(:,n)).*vpsi(:,m);
             braket_mn = conj(vpsi(:,m)).*psi(:,n); % Redundant
             M1 = sum(braket_nm)*dx;
             M2 = sum(braket_mn)*dx; % Redundant -> But necessary for d > 1
             if abs(M1) < 0.001 || abs(M2) < 0.001
                 res(n,m) = 0; % Helps avoiding numerical noise
             else
                  res(n,m) = 4*M1*M2/(dw*L*a);              
             end
              m = m+1;
         end
         n = n+1;
     end
     nMax = n; % An ugly workaround, to be
     mMax = m; % polished at some point...
end

