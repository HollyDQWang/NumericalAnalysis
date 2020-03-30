%% Householder QR 
m = 30; n = 20; 
A = randn(m,n);
A = gallery('randsvd',[m n],1e10);

Anow = A; Q = eye(m); R = []; AA = A; 
spy(abs(AA)>1e-14), snapnow
for ii = 1:n    
 pause    
% pause(0.3)
q = Anow(:,1)-norm(Anow(:,1))*eye(m-ii+1,1);
q = q/norm(q);
Anow = Anow-2*q*(q'*Anow); % Householder reflector Q*a
Q = Q*blkdiag(eye(m-length(q)),eye(length(q))-2*q*q');
R = [R;zeros(1,ii-1) Anow(1,:)];
AA(ii:end,ii:end) = Anow; 
spy(abs(AA)>1e-14), snapnow, 
Anow = Anow(2:end,2:end);
end
[norm(A-Q(:,1:n)*R) norm(Q'*Q -eye(m))/sqrt(m)] % residual and orthogonality


%% QR factorization, gram-schmidt
m = 30; n = 20; 
%A = randn(m,n);
A = gallery('randsvd',[m n],1e6);

R = norm(A(:,1));
Q = A(:,1)/R;

% classical GS 
for ii = 2:n
R12 = Q'*A(:,ii); % ith column of R (above diag)
Anow = A(:,ii)-Q*R12; % ith column of A orthogonal to previous Q 
R22 = norm(Anow); % (i,i) element of R 
R = [[R;zeros(1,size(R,2))] [R12;R22]];
Q = [Q Anow/norm(Anow)];
end
spy(R)
Q'*Q 
norm(A-Q*R)

%% modified GS 
%m = 4; n = 3; A = gallery('randsvd',[m n],1e6);
%A = randn(m,n);

R = norm(A(:,1));
Q = A(:,1)/R;

for ii = 2:n
%R12 = Q'*A(:,ii); Anow = A(:,ii)-Q*R12; % CGS
R12 = zeros(ii-1,1); Anow = A(:,ii);
    for jj = 1:size(R12,1) % only difference from CGS
        R12(jj) = Q(:,jj)'*Anow;  
        Anow = Anow - Q(:,jj)*R12(jj);
    end
R22 = norm(Anow);
R = [[R;zeros(1,size(R,2))] [R12;R22]];
Q = [Q Anow/norm(Anow)];
end
[norm(Q'*Q-eye(n)) norm(Q*R-A)]

