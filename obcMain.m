global particleDensity V0 W     % Physical parameters
global step dx basisDIM cutoff  % Computational parameters

%% First of all we define the physical parameters of our system:

particleDensity = 1 % n = N/L
cutoff = 55         % In ``Hartree''
V0 = 20             % In ``Hartree''
W = 0.2             % Width of V(x)

%% And some computational parameter:

step = 100;     % Accuracy of numerical integrals
dx   = 1/step;
basisDIM = 700; % Basis dimension

%% Thermodynamic Limit:
%  $L\to\infty$, $N\to\infty$ such that $N/L = n$ is constant. 
%  So we need a cycle over increasing L values:

Lmin  = 002;         %
Lmax  = 162;         % Only even values <=> Spinless electrons
Lstep = 020;         %

poles = zeros(step,step,step);
residues = zeros(step,step,step);

for L = Lmin:Lstep:Lmax 
    
    N = particleDensity*L;
    
 %% Building Crystallite's Hamiltonian
  % TO DO: check if an Hamiltonian is already saved...and else save it!
    H = obcHamiltonian(L);
    [c,E] = eig(H,'vector');
    
 %% Computing Fermi Energy
    if mod(N,2) == 0
        nF = N/2;
    else
        nF = (N+1)/2;
    end
    EF = E(nF); EF
    
 %% Geometrical Drude Weight
  % $D = 2v_\mathrm{F}$, with some care on what $v_\mathrm{F}$ is:
    if nF > 2
     % OBC Fermi velocity
        vF = (E(nF+1)-E(nF-1))/(2*(pi/L));
     % OBC Geometrical Drude Weight
        deltaL = L-Lmin;
        Dw(deltaL/Lstep + 1) = 2*vF/pi;
    end
    
  % [Plot Energy Scheme]
  
 %% Building up wavefunctions:
     x = linspace(0,L,step*L);
     psi = zeros(length(x),basisDIM);
     for n = 1:basisDIM
          psi_x = zeros(1,length(x));
          for m = 1:basisDIM
              psi_x = psi_x + psi0(x,m,L) * c(m,n);
          end
          psi(:,n) = psi_x;
     end
     
   % [Plot Wavefunctions]
   
   
 %% Kubo Poles and Residues
  % Now we proced with linear response evaluation (Kubo Formulae)
  [w,res] = obcKubo(E,psi,nF,L);
end