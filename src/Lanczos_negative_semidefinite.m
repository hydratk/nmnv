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
run('../styles/style.m')

fig = figure;
ax1 = subplot(1, 2, 1); 
semilogy(1:length(error_exp_times_v), error_exp_times_v, '-o', 'DisplayName', '$e^A v$')
hold on
semilogy(1:length(error_inv_times_v), error_inv_times_v, '-o', 'DisplayName', '$(I-A)^{-1} v$')
semilogy(1:length(Krylov_estimate), Krylov_estimate, '-o', 'DisplayName', 'Estimate from article')
semilogy(1:length(Taylor_estimate), Taylor_estimate, '-o', 'DisplayName', 'Taylor estimate')

% Axes settings
ax1.YScale = 'log';
ax1.YMinorTick = 'off';
ax1.YMinorGrid = 'off';
ax1.YTick = 10.^(-10:2:0);
ax1.YLim = [1e-10 1e1];

title(ax1, 'Lanczos Negative Semidefinite Matrix')
xlabel(ax1, 'Number of iterations')
ylabel(ax1, 'Relative error')
legend(ax1, 'show', 'Location', 'northeast', 'Interpreter', 'latex');
hold off

ax2 = subplot(1, 2, 2);

% Capture handles (h1, h2) to explicitly pass to legend later
semilogy(1:length(error_exp_times_v), error_exp_times_v, '-o', 'DisplayName', '$e^A v$')
hold on
semilogy(1:length(error_inv_times_v), error_inv_times_v, '-o', 'DisplayName', '$(I-A)^{-1} v$')
semilogy(1:length(Krylov_estimate), Krylov_estimate, '-o', 'DisplayName', 'Estimate from article')
semilogy(1:length(Taylor_estimate), Taylor_estimate, '-o', 'DisplayName', 'Taylor estimate')

% Axes settings
ax2.YScale = 'log';
ax2.YMinorTick = 'off';
ax2.YMinorGrid = 'off';
% Adjusted limits to ensure data is visible. 
% If error is > 1e-4, the previous limits [1e-16 1e-4] would show a blank plot.
ax2.YLim = [min(error_exp_times_v) max(Taylor_estimate)]; 
ax2.YTick = 10.^(-16:2:ceil(log10(max(Taylor_estimate))));

title(ax2, 'Lanczos Negative Semidefinite Matrix (Zoomed)')
xlabel(ax2, 'Number of iterations')
ylabel(ax2, 'Relative error')

% Pass the handles [h1, h2] to force legend to recognize these lines
legend(ax2, 'Location', 'northeast', 'Interpreter', 'latex');
hold off

% Export
exportgraphics(fig, '../graphics/lanczos_negative_semidefinite.pdf', 'ContentType', 'vector');
