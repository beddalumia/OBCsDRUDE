%% Basic-Test of the function

addpath ../abcFunctions

a  = 0.2       ;% Lattice parameter
L  = 3         ;% Number of cells [units of a]
V0 = 1         ;% Strenght of V_at(x) [hartree]
W  = 0.1       ;% FWHM of V_at(x) [units of a]
shape = 'HV'   ;% Functional form: 'HV' | 'MT' | 'KP'
csym  = 1      ;% Do we want centrosymmetry? 1 | 0

x = linspace(0,L*a,100*L);
y = V(x,a,L,V0,W,shape,csym);

figName = sprintf('Test generation of the %s periodic potential',shape);
figure("Name",figName);
area(x,y,'FaceColor', [0.9 0.8 0.9]); hold on
xline(L*a/2); 
xlim([0,L*a]); hold off

%% Speed-Benchmark of the three shapes (expected MT > KP >> HV)

a = 1; V0 = 1; W = 0.2; csym = 1;

shape = {'MT','HV','KP'}; L = round(logspace(0,3,25)); 
t = zeros(length(shape),length(L));
for i = 1:length(shape)
   for j = 1:length(L)
       x = linspace(0,L(j)*a,10*L(j));
       tic;
       y = V(x,a,L(j),V0,W,shape{i},csym);
       t(i,j) = toc;
   end
end

figure("Name",'Speed benchmark for the three potential shapes')
for i = 1:3
   plot(L,t(i,:),'LineWidth',1); hold on 
end
set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');
legend(shape,'Location','best');
title('Speed benchmark for the three potential shapes');
xlabel('Number of crystal cells');
ylabel('Total time to compute V(x)');
hold off
