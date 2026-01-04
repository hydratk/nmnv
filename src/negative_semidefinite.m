clear all;
close all;
run('../styles/style.m')

%% construction of the matrix A and vector v from the article
A = zeros(1001);
for j = 1:1001
    A(j,j) = -40 + (j-1) * 40/1000;
end
v = rand(1001,1);
v = v/norm(v);

%% convergence curve for exp(A)*v and (I-A)^-1 *v
exact_exp_times_v = expmv(A,v);
exact_inv_times_v = (eye(1001) - A) \ v;

maxit = 80;
[V,H] = myArnoldi_mgs(A,v,maxit); 

error_exp_times_v = zeros(1,maxit);
error_inv_times_v = zeros(1,maxit);
for j = 1:maxit
    V_j = V(:, 1:j);
    H_j = H(1:j, 1:j);
    e_1 = zeros(j,1);
    e_1(1,1) = 1;
    Krylov_exp_times_v = V_j * expmv(H_j,e_1);
    Krylov_inv_times_v = V_j * (eye(j)-H_j)^-1 * e_1;
    error_exp_times_v(j) = norm(Krylov_exp_times_v - exact_exp_times_v);
    error_inv_times_v (j) = norm(Krylov_inv_times_v - exact_inv_times_v);
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
% Make figure wider so subplots have room for legends/labels
% fig.Position = [100, 100, 1200, 500]; 

% --- Subplot 1: Full Range with Estimates ---
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

title(ax1, 'Negative Semidefinite Matrix')
xlabel(ax1, 'Number of iterations')
ylabel(ax1, 'Relative error')
legend(ax1, 'show', 'Location', 'northeast', 'Interpreter', 'latex');
hold off


% --- Subplot 2: Zoomed View (Errors Only) ---
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

title(ax2, 'Negative Semidefinite Matrix (Zoomed)')
xlabel(ax2, 'Number of iterations')
ylabel(ax2, 'Relative error')

% Pass the handles [h1, h2] to force legend to recognize these lines
legend(ax2, 'Location', 'northeast', 'Interpreter', 'latex');
hold off

% Export
exportgraphics(fig, '../graphics/negative_semidefinite.pdf', 'ContentType', 'vector');
