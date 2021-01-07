addpath ../abcFunctions
addpath ../obcFunctions

% First of all we define the physical parameters of our system:

a = 0.2;             % Lattice parameter; In bohr.
particleDensity = 1; % IN UNITS OF LATTICE PARAMETER!
V0 = 20;             % Strength of V(x); In hartree.
W = 0.2;             % Width of V(x); In bohr.
shape = 'HV';        % Functional form of V(x): 'HV' | 'MT' | 'KP'
csym = true;         % Do we want centrosymmetry in V(x)?  true | false

% And some computational parameter:

step = 100;         % Accuracy of numerical integrals
dx   = 1/step;      % " " " " " " " " " " " " " " " " 
basisDIM = 700;     % Basis dimension (for Hamiltonian diagonalization)
cutoff = 55;        % In ``Hartree''

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Thermodynamic Limit Calculation:
%  $L\to\infty$, $N\to\infty$ such that $N/L = n$ is constant. 
%  So we need a cycle over increasing L values:

Lmin  = 002;        %
Lmax  = 122;        % Only even values <=> Spinless electrons
Lstep = 002;        %

fprintf('###########################################\n');
fprintf('Lattice parameter: %f bohr\n',a);
fprintf('Strength of V(x): %f hartree\n',V0);
fprintf('Width of V(x): %f bohr\n',W);
if csym == true
fprintf('V(x) is a centrosymmetric %s-type potential\n',shape);
else
fprintf('V(x) is a noncentrosymmetric %s-type potential\n',shape);   
end
fprintf('###########################################\n');
fprintf('Calculation with %d electrons per unit cell\n',particleDensity);
fprintf('Using %d FE-states as a diagonalization basis\n',basisDIM);
fprintf('And a cutoff of %f hartree for the Kubo series\n',cutoff);

physID = sprintf('N%dA%.1fV%.1fW%.1f%s%d',particleDensity,a,V0,W,shape,csym);
diagID = [physID,sprintf('ACC%dDIM%d',step,basisDIM)];
kuboID = [diagID,sprintf('cut@%.1f',cutoff)];

Lvalues = Lmin:Lstep:Lmax;      % IN UNITS OF LATTICE PARAMETER HERE!

for i = 1:length(Lvalues)      %-> For sequential execution
%parfor i = 1:length(Lvalues)    %-> For parallel execution
    
    L = Lvalues(i);
    
    fprintf('~~~~~~~~~~~~~\n',L);
    fprintf('# Cells = %d\n',L);
    fprintf('~~~~~~~~~~~~~\n',L);
    
    N = particleDensity*L;
    
  % Building (loading) crystallite's Hamiltonian
    hamiltonianID = ['../matFiles/H_',sprintf('L%d',L),diagID,'.mat'];
    if isfile(hamiltonianID)
        fprintf('Loading the Hamiltonian..\n');
        WS = load(hamiltonianID); 
        H = WS.H; fprintf('.DONE!\n');
    else
        fprintf('Building up the Hamiltonian..');
        H = obcHamiltonian(a,V0,W,L,shape,csym,step,dx,basisDIM); 
        parsave(hamiltonianID,H); fprintf('.DONE!\n');
    end
    fprintf('Diagonalizing Hamiltonian..');
    [c,E] = eig(H,'vector'); fprintf('.DONE!\n');
    
  % Computing Fermi Energy
    if mod(N,2) == 0
        nF = N/2;
    else
        nF = (N+1)/2;
    end
    EF = E(nF);
    
%   % Geometrical Drude Weight
%   % $D = 2v_\mathrm{F}$, with some care on what $v_\mathrm{F}$ is:
%     if nF > 2
%      % OBC Fermi velocity
%         vF = (E(nF+1)-E(nF-1))/(2*(pi/(L*a)));
%      % OBC Geometrical Drude Weight
%         deltaL = L-Lmin;
%         gDw(deltaL/Lstep + 1) = 2*vF/pi;
%     end
    
  % [Plot Energy Scheme]
  
  % Building up wavefunctions:
  fprintf('Building up the wavefunctions..');
     x = linspace(0,L*a,step*L*a);
     psi = zeros(length(x),basisDIM);
     for n = 1:basisDIM
          psi_x = zeros(1,length(x));
          for m = 1:basisDIM
              psi_x = psi_x + psi0(x,m,L*a) * c(m,n);
          end
          psi(:,n) = psi_x;
     end; fprintf('.DONE!\n');
     
   % [Plot Wavefunctions]
   
   
 %% Kubo Poles and Residues
    fprintf('Kubo linear response..');
    [w,res,nMax,mMax] = obcKubo(E,psi,nF,a,L,dx,step,basisDIM,cutoff); 
    polesID = ['../matFiles/w&Res_',sprintf('L%d',L),kuboID,'.mat'];
    parsave(polesID,w,res); fprintf('.DONE!\n');
  
  % Reshaffling (ugly and unuseful in the end, I think)
%     deltaL = L-Lmin;
%     nPoles = zeros(1,mMax-1);
%     nResidues = zeros(1,mMax-1);
%     for n = 1:nF
%         for m_n = 1:(mMax-n)
%             m = n + m_n; % i.e. m-n = m_n ;)
%             nPoles(m_n) = w(n,m);
%             nResidues(m_n) = res(n,m);
%         end
%      poles(n, deltaL/Lstep+1, 1:length(nPoles)) = nPoles;
%      residues(n, deltaL/Lstep+1, 1:length(nResidues)) = nResidues;
%     end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Post Calculation Analysis
%  Sum-Rules and Conductivity Spectra

%% Drude Weight and f-sum:

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% All the plotting

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Basic routines for parpools
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function parsave(fname,varargin)
  numvars=numel(varargin);
  for i=1:numvars
     eval([inputname(i+1),'=varargin{i};']); 
     % This should allow giving it a variable number of arguments and the 
     % function recognizes their names.
  end
  save('-mat',fname,inputname(2));
  for i = 2:numvars    
    save('-mat',fname,inputname(i+1),'-append');
  end
end