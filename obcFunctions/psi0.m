function psi0_n = psi0(x,n,L)
%% Imperturbed OBC eigenfunctions
   A = sqrt(2/L); % Normalization
   psi0_n = A*sin(pi/L * n * x);
end

