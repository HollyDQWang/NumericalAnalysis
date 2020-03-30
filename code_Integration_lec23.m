%% Integration, trapezium & Simpson
%f = @(x) exp(x); 
f = @(x) 1./(25*x.^2+1); 

% exact value
ff = chebfun(@(x)f(x)); I = sum(ff);

n = 2; %1: trap, 2: simpson
x = linspace(-1,1,n+1).'; 
F = f(x);
fNC = chebfun.interp1(x,F); % Interpolation in Newton-Cotes
chebfunsetting , clf
plot(ff,LW,lw), hold on
plot(x,F,'k.','markersize',14,LW,lw), grid on
plot(fNC,'.--',LW,lw),shg
title(['exact = ',num2str(I),', approx =',num2str(sum(fNC)),', error =',num2str(sum(fNC)-I)])

%% Composite rules: trap
N = 20; % #sample points
x = linspace(-1,1,N).'; 
F = f(x);
clf, plot(ff,LW,lw), hold on
Iquad = 0; 
for ii = 1:N-1
xnow = x(ii:ii+1); Fnow = F(ii:ii+1);
fNC = chebfun.interp1(xnow,Fnow); % Interpolation in Newton-Cotes
plot(x,F,'k.','markersize',14,LW,lw), grid on
plot(fNC,':',LW,lw),shg
%Iquad = Iquad + sum(fNC); 

% exact integral of interpolant
h = diff(xnow); Iquad = Iquad + h/2*(Fnow(1)+Fnow(2)); 
end
title(['exact = ',num2str(I),', approx =',num2str(Iquad),', error =',num2str(I-Iquad)])

%% Composite rules: Simpson's 
N = 21; % #sample points, 2k+1
x = linspace(-1,1,N).'; 
F = f(x);
clf, plot(ff,LW,lw), hold on
Iquad = 0; 
h = x(2)-x(1);
for ii = 1:(N-1)/2;
xnow = x(2*ii-1:2*ii+1); Fnow = F(2*ii-1:2*ii+1);
fNC = chebfun.interp1(xnow,Fnow); % Interpolation in Newton-Cotes
plot(x,F,'k.','markersize',14,LW,lw), grid on
plot(fNC,':',LW,lw),shg
Iquad = Iquad + h/3*(Fnow(1)+4*Fnow(2)+Fnow(3)); 
%Iquad = Iquad + sum(fNC); 
end
title(['exact = ',num2str(I),', approx =',num2str(Iquad),', error =',num2str(I-Iquad)])


%% convergence 

it = 0;
ns = 2.^(1:10);
err = zeros(numel(ns),1);
for n = ns
it = it+1;
x = linspace(-1,1,n).'; 
F = f(x);
Iquad = 0; 
for ii = 1:n-1
xnow = x(ii:ii+1); Fnow = F(ii:ii+1);
% exact integral of interpolant
h = diff(xnow); Iquad = Iquad + h/2*(Fnow(1)+Fnow(2)); 
end
err(it) = I-Iquad; 
end
clf, loglog(ns,err,'.--',MS,ms),shg
text(ns(end-2),abs(err(end-2)),'Trap',FS,fs)
grid on
xlabel('#samples',FS,fs)
ylabel('Integration error',FS,fs)

%% simpson's 
ns = ns+1; % #sample points, 2k+1
errS = zeros(numel(ns),1);
it = 0;
for n = ns
it = it+1;    
x = linspace(-1,1,n).'; 
F = f(x);
Iquad = 0; 
h = x(2)-x(1);
for ii = 1:(n-1)/2;
xnow = x(2*ii-1:2*ii+1); Fnow = F(2*ii-1:2*ii+1);
Iquad = Iquad + h/3*(Fnow(1)+4*Fnow(2)+Fnow(3)); 
end
errS(it) = I-Iquad; 
end
hold on, loglog(ns,abs(errS),'.--',MS,ms),shg
text(ns(end-2),errS(end-2),'Simpson',FS,fs)
grid on
xlabel('#samples',FS,fs)
ylabel('Integration error',FS,fs)

%% Gauss quadrature
errG = zeros(numel(ns),1);
it = 0;
for n = ns
it = it+1;    
[x,w] = legpts(n);
Iquad = w*ff(x); % n-point Gauss quadrature
errG(it) = I-Iquad; 
end
loglog(ns,abs(errG)+eps,'.-',MS,ms),shg
text(ns(end-2),errG(end-2)+eps,'Gauss',FS,fs)

%% Gauss converges exponentially
clf
semilogy(ns,abs(err),'.-',MS,ms),hold on
semilogy(ns,abs(errS),'.-',MS,ms)
semilogy(ns,abs(errG)+eps,'.-',MS,ms)
grid on
% xlim([1 200])

%% adaptive_simpson
a = -1; b = 1; tol = 1e-5; nmax = 50; 
adaptive_simpson(ff, a, b, tol, nmax)
