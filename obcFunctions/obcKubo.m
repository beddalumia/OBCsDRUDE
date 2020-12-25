function [w,res] = obcKubo(E,psi,nF,L)

    global dx basisDIM cutoff % Computational parameters

  %% First of all, we apply the velocity operator to the eigenstates
   % $\hat{v}\psi_n(x) = i \frac{\mathrm{d}}{\mathrm{d}x} \psi_n(x)$

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
             braket_mn = conj(vpsi(:,m)).*psi(:,n);
             M1 = sum(braket_nm)*dx;
             M2 = sum(braket_mn)*dx;
             if abs(M1) < 0.001 || abs(M2) < 0.001
                 res(n,m) = 0;
             else
                  res(n,m) = 4*M1*M2/(dw*L);              
             end
              m = m+1;
         end
         n = n+1;
     end
     nMax = n;
     mMax = m;
     
  %% Reshaffling: (NO IDIOT, poles&residues are in the main scope!)
%      deltaL = L-Lmin;
%      nPoles = zeros(1,mMax-1);
%      nResidues = zeros(1,mMax-1);
%      for n = 1:nF
%          for m_n = 1:(mMax-n)
%              m = n + m_n; % i.e. m-n = m_n ;)
%              nPoles(m_n) = w(n,m);
%              nResidues(m_n) = res(n,m);
%          end
%          poles(n, deltaL/Lstep+1, 1:length(nPoles)) = nPoles;
%          residues(n, deltaL/Lstep+1, 1:length(nResidues)) = nResidues;
%      end
end

