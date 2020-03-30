%% some of the codes require that 'Chebfun' is installed (google should find it without problem)
%% 1st example
clf
 %lagrange([1,1.2,1.3,1.4],[4,3.5,3,7]);
 lagrange([1,1.2,1.3,1.4],[4,3.5,3,0]);
shg, grid on, yy = ylim;
%% Lagrange polynomials
 x = [1,1.2,1.3,1.4]; 
 p1 = lagrange(x,[4,0,0,0]);
 p2 = lagrange(x,[0,3.5,0,0]);
 p3 = lagrange(x,[0,0,3,0]);
 p4 = lagrange(x,[0,0,0,0]);

 lagrange([1,1.2,1.3,1.4],[4,3.5,3,0]);shg
 grid on,hold on
 plot(x,0*x,'k.','markersize',15)
 xx = xlim; 
 xx = linspace(xx(1),xx(2),1000);
 plot(xx,polyval(p1,xx))
 plot(xx,polyval(p2,xx))
 plot(xx,polyval(p3,xx))
 plot(xx,polyval(p4,xx))
 ylim([yy])

 %%
 lagrange([0,2.3,3.5,3.6,4.7,5.9],[0,0,0,1,1,1]);
 
%% convergence (hopefully)
% f = @(x)exp(x); % exponential
f = @(x) 1./(25*x.^2+1); % Runge
x = linspace(-1,1,20) ; % try 5,10,20
y = f(x);
p = lagrange(x,y),shg, grid on, hold on
figure
xx = linspace(-1,1,1000);
plot(xx,polyval(p,xx)-f(xx)), grid on, hold on
plot(x,0*x,'mx','markersize',18,'linewidth',3)
title('error'),shg


%% poly for Hermite interpolation
clear
n = 4; 
x = linspace(-1,1,n+1);
k = 2; 
xk = x(k);
xknow = x.'; xknow(k) = [];
Lk = @(x) prod(x-xknow); 
Lk = @(x) Lk(x)/Lk(xk); % Lagrange
dLk = diff(chebfun(@(x)Lk(x))); 
Hk = @(x) ((Lk(x)).^2).*(1-2*(x-xk).*dLk(xk));
clf,fplot(Hk,[-1,1],'linewidth',2),shg, grid on, hold on
plot(x,Hk(x),'k.','markersize',15)
plot(xk,Hk(xk),'r.','markersize',20)
xkk = (2*xk+x(min(k+1,end)))/3;
text(xkk,Hk(xkk),'Hk','fontsize',20,'color','b')

Kk = @(x) ((Lk(x)).^2).*(x-xk);
hold on,fplot(Kk,[-1,1],'r-','linewidth',2),shg, grid on, hold on
plot(x,Kk(x),'k.','markersize',15)
plot(xk,Kk(xk),'r.','markersize',20)
text(xkk,1.3*Kk(xkk),'Kk','fontsize',20,'color','r')

%% Hermite interpolation

% f = chebfun(@(x)exp(x)); 
f = chebfun(@(x)1./(25*x.^2+1)); 
df = diff(f); 

n = 6; 
x = linspace(-1,1,n+1);
F = f(x);
dF = df(x); 
p = @(x)0; 
for k = 1:n+1 
xk = x(k);
xknow = x.'; xknow(k) = [];
Lk = @(x) prod(x-xknow); 
Lk = @(x) Lk(x)/Lk(xk); % Lagrange
dLk = diff(chebfun(@(x)Lk(x))); 
Hk = @(x) ((Lk(x)).^2).*(1-2*(x-xk).*dLk(xk));
Kk = @(x) ((Lk(x)).^2).*(x-xk);
p = @(x) p(x) + Hk(x)*F(k) + Kk(x)*dF(k);
end

close all
plot(f,'linewidth',2), hold on, 
fplot(p,[-1 1],'--','linewidth',2),shg, grid on
figure
fplot(@(x)f(x)-p(x),[-1 1],'linewidth',2)
title('error'), grid on 
alignfigs

%% Chebfun demo
f = @(x) 1./(25*x.^2+1); 
%x = linspace(-1,1,15); % sample points 
x = cos(linspace(0,pi,50)); % better sample points
clf, subplot(2,1,1)
xx = linspace(-1,1,1000); % test points
p = chebfun.interp1(x,f(x)); 
plot(xx,f(xx),'k--'),hold on
title('functions')
plot(p,'linewidth',2),grid on,shg
hold on, plot(x,f(x),'r.','markersize',13)
subplot(2,1,2)
plot(p(xx)-f(xx),'linewidth',2)
title('error'), grid on

