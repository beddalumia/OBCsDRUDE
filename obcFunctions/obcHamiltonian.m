function H = obcHamiltonian(L,shape,csym)
%%  Hamiltonian matrix for the OBC system
%   We build the crystal hamiltonian computing matrix elements
%   $\{\langle \psi^0_n | K + V(x) | \psi^0_m\rangle\}$    

    global a V0 W % Physical parameters
    global step dx basisDIM % Computational parameters

    H = zeros(basisDIM,basisDIM);
    x = linspace(0,L*a,step*L*a);
    Potential = V(x,a,L,V0,W,shape,csym);
    for n = 1:basisDIM
        for m = 1:basisDIM

            %% Kinetic (diagonal) part %%
            if n == m
                H(n,n) = eigenE0(n,L*a);
            end

            %% Periodic potential part 
            braUket = conj(psi0(x,n,L*a)).*psi0(x,m,L*a).*Potential;
            H(n,m) = H(n,m) + sum(braUket)*dx;
        end
    end
end

