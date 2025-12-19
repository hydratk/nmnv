function [V, T] = myLanczos(A, b, n)
% LANCZOS   Lanczos iteration for symmetric matrix A.
%   [V,T] = myLanczos(A,b,n) computes an orthonormal basis V and
%   tridiagonal rectangular matrix T such that
%       A*V(:,1:n) = V(:,1:n)*T(1:n,1:n) + T(n+1,n)*V(:,n+1) * [0 ... 0 1]
%       or equivalently
%       A*V(:,1:n) = V*T
%
% Inputs:
%   A - symmetric matrix (m x m)
%   b - initial vector (length m)
%   n - number of iterations
%
% Outputs:
%   V - orthonormal basis (m x (n+1))
%   T - tridiagonal matrix ((n+1) x n)

tol = 1e-14;

m = size(A,1);
V = zeros(m, n+1);
T = zeros(n+1, n);

% Normalize the starting vector
V(:,1) = b / norm(b);
w = A * V(:,1);
alpha = V(:,1)' * w;
w = w - alpha * V(:,1);
T(1,1) = alpha;

for j = 1:n-1
    beta = norm(w);
    if beta < tol
        fprintf('Lanczos terminated early at step %d\n', j);
        return
    end
    V(:,j+1) = w / beta;
    T(j+1,j) = beta;
    T(j,j+1) = beta;
    T(j,j) = alpha;

    w = A * V(:,j+1) - beta * V(:,j);
    alpha = V(:,j+1)' * w;
    T(j+1,j+1) = alpha;
    w = w - alpha * V(:,j+1);

end
beta = norm(w);
V(:,n+1) = w / beta;    % SP: this were missing in the previous code
T(n+1,n) = beta;
end
