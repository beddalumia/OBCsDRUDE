particleDensity = 0.2; % n = N\L

x = 10:10:500; i = 1;
yr = cell(1,length(x));
yf = cell(1,length(x));
for L = x
    fprintf('L=%d\n',L);
    freq = [];
    res = [];
    
    N = particleDensity*L;
    fprintf('N=%d\n',N);
    nF = N/2;
    EF = eigenE0(nF,L);
    cut = 0.01;
    
    n = 1;
    while eigenE0(n,L) <= EF
        m = nF+1;
        dw = eigenE0(m,L) -  eigenE0(n,L);
        while dw < cut
            Jnm = sqrt(1/(2*pi))*(cos((m-n)*pi)/(m-n) - cos((m+n)*pi)/(m+n)); 
            Res = 4*(Jnm)^2/(dw*L); 
            freq = [freq,dw];
            res = [res,Res];
            m = m+1;
           dw = eigenE0(m,L) -  eigenE0(n,L);
        end
        n = n+1;
    end
   
   %To scatter cool poles with radius ~ res -> Almost visual art!
   %scatter(L*ones(1,length(freq)),freq,res*50,'k','filled'); hold on
   %To scatter f-sum computed from pole-residue sums -> Dirty but works
   %scatter(L,sum(res)); hold on
   [yf{i}, idx] = sort(freq);
   yr{i} = res(idx);
   i = i+1;
end

l = x;
for i = 1:length(x)
    l(i) = length(yr{i});
end
lmax = max(l);

yPol = zeros(length(x),lmax);
yRes = zeros(length(x),lmax);
for i = 1:length(x)
   j = 1;
   l = length(yr{i});
   while j <= l
       yRes(i,j) = yr{i}(j);
       yPol(i,j) = yf{i}(j);
       j = j+1; 
   end
end

figure("Name",'Poles')
for j = 1:lmax
   y = yPol(:,j);
   plot(x(y~=0),y(y~=0),'k'); hold on 
end
figure("Name",'Residues')
for j = 1:lmax
   y = yRes(:,j);
   plot(x(y~=0),y(y~=0),'k'); hold on 
end
set(gca, 'YScale', 'log')

function En = eigenE0(n,L)
         En = (pi^2 * n^2)/(2 * L^2);
end