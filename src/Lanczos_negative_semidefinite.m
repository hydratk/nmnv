% Illustration of the error bound given by Theorem 2 using Lanczos method,
% since Arnoldi and Lanczos approximations here coincide.

clear all;
close all;

%% Construction of matrix A and vector v from article
m = 1001;
A = diag(linspace(-40, 0, m));

v = rand(1001,1);
v = v/norm(v);

%% convergence cruve for exp(A)*v and (I-A)^-1 *v

exact_exp_times_v = expmv(A,v);
exact_inv_times_v = (eye(m) - A) \ v;

maxit = 80;
[V,T] = myLanczos(A,v,maxit); % function myLanczos() from Lab6

error_exp_times_v = zeros(1,maxit);
error_inv_times_v = zeros(1,maxit);
for j = 1:maxit
    V_j = V(:, 1:j);
    T_j = T(1:j, 1:j);
    e_1 = eye(j,1);
    Krylov_exp_times_v = V_j * expmv(T_j,e_1);
    Krylov_inv_times_v = V_j * (eye(j)-T_j)^-1 * e_1;
    error_exp_times_v(j) = norm(Krylov_exp_times_v - exact_exp_times_v);
    error_inv_times_v(j) = norm(Krylov_inv_times_v - exact_inv_times_v);
end

%% tabulate estimates

rho = 10;
tau = 1;
e = exp(1);

Krylov_estimate = zeros(1,maxit);
prefactor = 10 * (rho*tau)^(-1) * exp(-rho*tau);
for m = ceil(sqrt(4*rho*tau)) : ceil(2*rho*tau)-1
    Krylov_estimate(m) = 10 * exp((-m^2)/(5*rho*tau));
end
for m = ceil(2*rho*tau):maxit
    Krylov_estimate(m) = prefactor * (e*rho*tau/m)^m;
end

Taylor_estimate = zeros(1,maxit);
normA = norm(0.5*A);
for m = 1 : maxit
    Taylor_estimate(m) = 2 * normA^m / factorial(m);
end

%% plot figure

fig = figure;
semilogy(error_exp_times_v) 
hold on
semilogy(error_inv_times_v)
semilogy(Krylov_estimate)
semilogy(Taylor_estimate)
ylim([1e-10, 1e1])
title('Lanczos method - negative semidefinite matrix')
legend('e^A v', '(I-A)^{-1} v', 'estimate (3.1), (3.2)', 'Taylor estimate', 'Location', 'northeast')
xlabel("number of iterations")
ylabel("relative error")
hold off
exportgraphics(fig ,'Lanczos_negative_semidefinite.png')