function [V, H] = myArnoldi_mgs(A, b, n)
% ARNOLDI METHOD
%   [V,H] = myArnoldi(A,b,n) produces an orthonormal basis V of the Krylov
%   subspace K_n(A,b) and the upper Hessenberg matrix H such that
%       A*V(:,1:n) = V(:,1:n)*H(1:n,1:n) + H(n+1,n)*V(:,n+1) * [0 ... 0 1]
%       or equivalently
%       A*V(:,1:n) = V*H
%
% Inputs:
%   A - square matrix (m x m)
%   b - initial vector (length m)
%   n - number of iterations
%
% Outputs:
%   V - matrix with orthonormal columns (m x (n+1))
%   H - upper Hessenberg matrix ((n+1) x n)

tol = 1e-14;

m = size(A,1);
V = zeros(m, n+1);
H = zeros(n+1, n);

% Normalize the starting vector
V(:,1) = b / norm(b);

for j = 1:n
    w = A * V(:,j);
    % Modified Gram-Schmidt orthogonalization
    for i = 1:j
        H(i,j) = V(:,i)' * w;
        w = w - H(i,j) * V(:,i);
    end
    H(j+1,j) = norm(w);
    if H(j+1,j) < tol
        fprintf('Breakdown. Arnoldi terminated early at step %d\n', j);
        return
    end
    V(:,j+1) = w / H(j+1,j);
end
end
