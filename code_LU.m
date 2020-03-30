% want to show pivots are important for stability; though often no pivot is fine too
n = 100;
A = randn(n);
%[L,U] = lu(A); % this does pivots!
[L,U] = lu_nopivot(A); % this does no pivots!
subplot(1,2,1), spy(L), subplot(1,2,2), spy(abs(U)>1e-10), shg
norm(L*U-A)/norm(A)

b = randn(n,1);
y = L\(b);
x = U\y;
norm(A*x-b)/norm(b) % we want this to be O(1e-15)
%% with pivots
A = randn(n);

[L,U] = lu(A); % PA=LU, note spy(L)
subplot(1,2,1), spy(L), title('L'), subplot(1,2,2), spy(abs(U)>1e-10), shg
pause
[L,U,P] = lu(A); 
subplot(1,2,1), spy(L), title('PL'), subplot(1,2,2), spy(abs(U)>1e-10), shg

y = L\(P*b);
x = (U\y);
norm(A*x-b)/norm(b)
