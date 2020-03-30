% Enjoy Gerschgorin disks
n = 3; 
% n = 10; 
rng(2)
%rng(5)
A = randn(n)+1i*randn(n);
A(1,:) = A(1,:)*0.3;
A(n,n) = 3+1i;

D = diag(diag(A)); 
F = A-D; 

ts = linspace(0,1,30);
co = {'b','r','k','m'};
for t = ts
Anow = D+t*F;    
for ii = 1:n
c = Anow(ii,ii); % center
R = sum(abs(Anow(ii,:)))-abs(Anow(ii,ii)); % radius
plot(chebfun(@(t)c+R*exp(2i*pi*t),[0 1]),co{mod(ii,4)+1}) % plot circle
%circle([real(c),imag(c)],R,'color',co{ii}) % plot circle, other way
hold on
end
axis equal, grid on
plot(eig(Anow),'r.','markersize',18)
shg, snapnow
pause(0.2)
hold off
end