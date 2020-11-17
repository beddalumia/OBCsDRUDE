function U = V(x,a,L,V0,W,shape,csym)
%% Periodic potential function: V(x+a) = V(x)
%  ------------------------------------------------------------------------
%% Input variables
%  > x is the whole domain of positions (a discretization of [0,L])
%  > a is the lattice parameters (bohr)
%  > L is the total length of the crystal (in units of a)
%  > V0 is the strenght of the potential (hartree) [applies to all shapes]
%  > W is the FWHM of the potential (units of a) [HV & KP shapes only]
%  > shape is a string identyfing the functional form of the potential:
%  ->'HV' for He-Vanderbilt Gaussian barriers (Phys.Rev.Lett.86,5341)
%  ->'MT' for Matthieu single-sine potential (Phys.Rev.87,807)
%  ->'KP' for Kronig-Penney square barriers (Proc.R.Soc.Lond.A.130,499–513)
%  > csym is a flag for centrosymmetry (is the crystal centrosymmetric?)
%  ------------------------------------------------------------------------

%% Centrosymmetry
%  If csym is set to TRUE we'll have centrosymmetric periodic potential, 
%  which is the "default result for dx = 0; If FALSE we set dx = a/4 and
%  do the trick x -> x-dx. See the implementations below for clarification.
   if csym == true
                    dx = 0;
   else
                    dx = a/4; % arbitrary "cell-polarization" (*)
   end
%
% *It may be interesting to make dx tunable, but for now let's stay here,
%  with the value that shifts of pi/4 the Mathieu sine potential.
 
%% He-Vanderbilt Gaussian atomic barriers (cf. Phys.Rev.Lett.86,5341)
   if shape == 'HV'
       U = 0;
       for x_at = 0:a:L*a
           U_at = V0.*exp(-(x-x_at-dx).^2 ./ (W*a)^2);
           U = U + U_at;
       end
%% Mathieu single-cosine periodic potential (cf. Phys.Rev.87,807)     
   elseif shape == 'MT'
       U = V0/2 * (1 + sin( 2*pi/a*(x+a/4-dx)) );
%                                     -csym-       
%% Kronig-Penney square atomic barriers (cf. Proc.R.Soc.Lond.A.130,499–513)  
   elseif shape == 'KP'
       U = V0/2 * ( 1 + square(2*pi/a*(x+W/2-dx),W*100) );
%                                        -csym-  -duty-
%% Error handling if shape is not recognized       
   else
       ErrorText = sprintf('Periodic potential shape is invalid!\n');
       Details = sprintf('> please enter either HV, MT or KP...');
       error([ErrorText,Details]);
   end
end
%% Notes 
%
%  1. He-Vanderbilt potential is implemented through the straightforward 
%     atomic recipe:
%                       U = 0;
%                       for x_at = 0:L*a
%                           U = U + U_at(x-x_at-dx)
%                       end
%     where U_at is a generic "inter-atomic" barrier function.
%     This way the obtained periodic potential satisfies always PBCs. 
%     Furthermore the crystal has inversion symmetry if and only if dx = 0.
%     Unfortunately the atomic decomposition makes for a quite inefficient
%     implementation within matlab language. Hamiltonian matrix build will
%     be severely slowed down with respect to 'MT' and 'KP' options.
%
%  2. Mathieu single-sine potential is of course implemented directly,
%     without any atomic decomposition: sin(x) is a primitive function.
%     We notice furthermore that the harmonic decomposition of a single sin
%     potential is indeed short-tailed, which is crucial regarding sparsity
%     of reciprocal-space (PBCs) and trigonometric (OBCs) hamiltonians. 
%     These facts imply that the 'MT' option makes indeed for a way faster
%     hamiltonian-matrix evaluation. Given that MT hamiltonians gives also
%     very flat lowest bands the 'MT' option could give in general the best
%     performance for a fixed target value of the Drude weight.
%
%  3. Kronig-Penney square barriers potential has been also implemented
%     through a native matlab function: square(t,duty), which gives a
%     square-wave signal with period 2pi and duty-cycle given in percentage 
%     by "duty". That's why we put duty = W*100, being W the total width of
%     the square barrier, given in units of the lattice period a.
%     Thanks to this native implementation the computation of the square
%     potential is quite fast (very similar to the Mathieu sine). 
%     Nevertheless the Fourier decomposition (either plane-wave or
%     trigonometric) will be for sure quite long tailed, so we expect the
%     hamiltonian matrix build to be not so efficient.
%
%  So, as far as speed is concerned, the preferred option should be 'MT'.
%  But being a sine a very irrealistic model for a crystalline potential we
%  argue that the other two options may be quite useful. Indeed in the main
%  publication produced with this code the 'HV' option has been used (with
%  optimized parameters...), seeking a balance between 'realism' and speed.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% MIT License 
% Copyright (c) 2020 Gabriele Bellomia  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
