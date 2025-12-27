% Illustration of the error bound given by Theorem 6 (Sectorial Matrix)
% using Arnoldi method.

clear all;
close all;

%% 1. Construction of matrix A and vector v
% We construct a matrix with eigenvalues in a sector of angle theta*pi.
% Theorem 6 assumes numerical range in rho * S_theta.
% S_theta has opening angle theta*pi at the origin.

m = 1001;
theta = 0.7; % The parameter theta from the paper (0 < theta < 1)
rho = 10;
tau = 1;

% Generate eigenvalues distributed inside the sector
% We distribute them typically up to a magnitude, e.g., 4*rho like previous examples
max_magnitude = 4 * rho; 
r_vals = linspace(0, max_magnitude, m)';

% Distribute angles linearly between -theta*pi/2 and +theta*pi/2
% This fills the "wedge"
max_angle = (theta * pi) / 2; 
angles = linspace(-max_angle, max_angle, m)';

% Shuffle angles to avoid structure correlation with magnitude
rng(1); % Fixed seed for reproducibility
angles = angles(randperm(m));

% Construct diagonal A (Normal matrix => Numerical Range = Convex Hull of Spectrum)
% Note: Eigenvalues are negative/left-half plane, so we use -r * exp(i*angle)
A = diag(-r_vals .* exp(1i * angles));

v = rand(m,1);
v = v/norm(v);

%% 2. Convergence curve for exp(A)*v and (I-A)^-1 *v

% Note: If you do not have 'expmv', replace with 'expm(A)*v'
% exact_exp_times_v = expmv(A,v); 
exact_exp_times_v = expm(A)*v; 

exact_inv_times_v = (eye(m) - A) \ v;

maxit = 80;
[V,H] = myArnoldi_mgs(A, v, maxit); 

error_exp_times_v = zeros(1,maxit);
error_inv_times_v = zeros(1,maxit);

for j = 1:maxit
    V_j = V(:, 1:j);
    H_j = H(1:j, 1:j);
    e_1 = zeros(j,1); e_1(1) = 1;
    
    % Krylov approximations
    % Note: Using expm for small H_j is fast
    Krylov_exp_times_v = V_j * expm(H_j) * e_1; 
    Krylov_inv_times_v = V_j * ((eye(j)-H_j) \ e_1);
    
    error_exp_times_v(j) = norm(Krylov_exp_times_v - exact_exp_times_v);
    error_inv_times_v(j) = norm(Krylov_inv_times_v - exact_inv_times_v);
end

%% 3. Tabulate Estimates (Theorem 6)

Krylov_estimate = zeros(1,maxit);
e = exp(1);

% We implement the bound from Eq (3.9) for large m:
% eps <= C * m * exp(-(r - psi(r))*rho*tau) * (e*rho*tau/m)^m
% where r = m / (rho*tau)
% and psi(w) = (1 - 1/w)^(2-theta) * w

% Constant C is generic in the paper, we set it to align roughly with the curve
% for visualization purposes (scaling factor).
C_const = 10; 

start_idx = ceil(2 * rho * tau);

for m_iter = start_idx : maxit
    r = m_iter / (rho * tau);
    
    % Conformal map psi(r) for real r > 1
    % Eq just before Theorem 6
    psi_r = (1 - 1/r)^(2-theta) * r;
    
    % Term 1: Exponential decay related to geometry
    decay_geom = exp(-(r - psi_r) * rho * tau);
    
    % Term 2: Superlinear decay term
    decay_super = (e * rho * tau / m_iter)^m_iter;
    
    Krylov_estimate(m_iter) = C_const * m_iter * decay_geom * decay_super;
end

Taylor_estimate = zeros(1,maxit);
normA = norm(A); % Or norm(0.5*A) if checking specific bounds
for k = 1 : maxit
    Taylor_estimate(k) = 2 * normA^k / factorial(k);
end

%% 4. Plot Figure

fig = figure;
semilogy(error_exp_times_v, 'LineWidth', 1.5) 
hold on
semilogy(error_inv_times_v, 'LineWidth', 1.5)
% Only plot estimate where valid (m >= 2*rho*tau)
semilogy(start_idx:maxit, Krylov_estimate(start_idx:end), '--', 'LineWidth', 1.5)
semilogy(Taylor_estimate, ':', 'LineWidth', 1.5)

ylim([1e-12, 1e1])
title(['Sectorial Matrix (\theta = ' num2str(theta) ')'])
legend('e^{\tau A} v', '(I-\tau A)^{-1} v', 'Thm 6 Estimate (Eq 3.9)', 'Taylor estimate', 'Location', 'southwest')
xlabel("number of iterations")
ylabel("relative error")
grid on;
hold off

exportgraphics(fig, 'Sectorial_Matrix.png');