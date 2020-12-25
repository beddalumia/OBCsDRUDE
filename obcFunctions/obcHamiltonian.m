function H = obcHamiltonian(L)
%%  Hamiltonian matrix for the OBC system
%   We build the crystal hamiltonian computing matrix elements
%   $\{\langle \psi^0_n | K + V(x) | \psi^0_m\rangle\}$    

    global V0 W % Physical parameters
    global step dx basisDIM % Computational parameters

    H = zeros(basisDIM,basisDIM);
    x = linspace(0,L,step*L);
    Potential = V(x,L,V0,W);
    for n = 1:basisDIM
        for m = 1:basisDIM

            %% Kinetic (diagonal) part %%
            if n == m
                H(n,n) = eigenE0(n,L);
            end

            %% Periodic potential part 
            braUket = conj(psi0(x,n,L)).*psi0(x,m,L).*Potential;
            H(n,m) = H(n,m) + sum(braUket)*dx;
        end
    end
end

