%% physical params defined by the project
m = 1.587; %kg
m_p = 0.25; %kg
g = 9.8065; %m/s^2
A = 0.0625; %m^2
rho = 1.225; %kg/m^3
R_rotor = 0.127; %m
A_r = pi * R_rotor^2; %m^2
L = 0.243; %m
l = 0.5; %m
k = 3.17 * 10^-5; %Ns^2
b = 7.69 * 10^-7; %Nms^2
I_xx = 0.0213; %kgm^2
I_yy = 0.0213; %kgm^2
I_zz = 0.0282; %kgm^2
k_D = 0.25; %kg/s
M = m+m_p;

%% actuator 
u2omega = (1/(4*b*L*k)) * [b*L, 2*b, 0, L*k; b*L, 0, 2*b, -L*k; b*L, -2*b, 0, L*k; b*L, 0, -2*b, -L*k];
omega2u = [k, k, k, k; L*k, 0, -L*k, 0; 0, L*k, 0, -L*k; b, -b, b, -b];

%% A linearized state matrix (phase 1)
Avv = diag([-k_D/m, -k_D/m, -k_D/M]);
Av_eta = [ 0, (M*g) / m,  0;
    -(M*g)/m,   0,   0;
    0,      0,   0];
Av_alpha = [-(m_p * g)/m;
    0;
    0];
a_alpha_v   = [k_D / (M * l), 0, 0];
a_alpha_eta = [0, -g/l, 0];

A_sys = zeros(14, 14);
A_sys(1:3, 4:6) = eye(3);
A_sys(4:6, 4:6)   = Avv;
A_sys(4:6, 7:9)   = Av_eta;
A_sys(4:6, 13)    = Av_alpha;
A_sys(7:9, 10:12) = eye(3);
A_sys(13, 14) = 1;
A_sys(14, 4:6) = a_alpha_v;
A_sys(14, 7:9) = a_alpha_eta;
A_sys(14, 13)  = g / l;

%% B matrix
B_v = [0,   0, 0, 0;
       0,   0, 0, 0;
       1/M, 0, 0, 0];

B_nu = [0, 1/I_xx, 0,      0;
        0, 0,      1/I_yy, 0;
        0, 0,      0,      1/I_zz];

B_sys = zeros(14, 4);
B_sys(4:6, :)   = B_v;
B_sys(10:12, :) = B_nu;

%% C & D matrix
measured_indices = [1, 2, 3, 7, 8, 9]; % 6 measured states

C_sys = zeros(6, 14);
for idx = 1:6
    C_sys(idx, measured_indices(idx)) = 1;
end
D_sys = zeros(6, 4);

%% controllability and rank
rank_ctrb = rank(ctrb(A_sys, B_sys));
rank_obsv = rank(obsv(A_sys, C_sys));

disp('Controllability Rank:');
disp(rank_ctrb)
disp('Observability Rank:');
disp(rank_obsv);

%% LQR for linearized controller design, kalman for observer (phase 2)
Q = diag([25.8,  25.8,  41.2,   ...  % x, y, z
    3.2,   3.2,  44.6,    ...  % x_dot, y_dot, z_dot
    4.86,  4.86,  26.3,   ...  % phi, theta, psi
    0.06,  0.12,  0.37,    ...  % p, q, r
    18.1, 2880         ...  % alpha, alpha_dot
    ]); % tunable!

R = diag([1.74, 5460, 3.09, 1310]); % delta_u_T, u_phi, u_theta, u_psi

% "optimal"
Q = diag([44.44   44.44   47.79   3.607   3.607   50.28   7.335   7.335    26.3 0.08601   0.172  0.5304   60.61    6678]);
R = diag([1.99   563 5.567  1310]);

K = lqr(A_sys, B_sys, Q, R);

G = eye(14);
B_aug = [B_sys, G];
D_aug = [D_sys, zeros(6, 14)];

Q_obs = diag([ ...
    10,  10,  10,   ...  % x, y, z
    100,  10,  10,   ...  % x_dot, y_dot, z_dot
    10,  10,  10,   ...  % phi, theta, psi
    10,  100,  10,   ...  % p, q, r
    300, 400        ...  % alpha, alpha_dot
    ]);
Q_obs = 1.42 * Q_obs;

R_obs = 0.01 * eye(6);

%"optimal"
Q_obs = diag([16.63 16.63 16.63 166.3 16.63 16.63 16.63 16.63 16.63 16.63 166.3 16.63 498.9 665.2]);

sys_augmented = ss(A_sys, B_aug, C_sys, D_aug);

sensors = 1:6;
known_inputs = 1:4;
[kalmf_linear, L_obs, P_linear] = kalman(sys_augmented, Q_obs, R_obs, [], sensors, known_inputs);
%% check
cl_ctrl_poles_c = eig(A_sys - B_sys * K);
cl_obs_poles_c  = eig(A_sys - L_obs * C_sys);

fprintf('Max real part of controller poles: %.4f \n', max(real(cl_ctrl_poles_c)));
fprintf('Max real part of observer poles:   %.4f \n', max(real(cl_obs_poles_c)));

obs_eigs = eig(A_sys - L_obs * C_sys);
fprintf('Slowest observer pole (real): %.4f (Target: -2 to -5)\n', max(real(obs_eigs)));
fprintf('Fastest observer pole (real): %.4f (Target: > -30)\n', min(real(obs_eigs)));

A_cl = [A_sys - B_sys*K,       B_sys*K;
    zeros(14, 14),         A_sys - L_obs*C_sys];
x0_aug = [zeros(12,1); 0.05; 0;   zeros(12,1); 0.05; 0]; % [x(0); e(0)]
[t_lin, x_lin] = ode45(@(t,z) A_cl*z, [0 5], x0_aug);
u_lin = (-K * (x_lin(:, 1:14) - x_lin(:, 15:28))')';
fprintf('Max linear pitch torque |tau_theta|: %.4f N*m\n', max(abs(u_lin(:,3))));
[t_lin, z_lin] = ode45(@(t,z) A_cl*z, [0 5], x0_aug);
x_true_lin = z_lin(:, 1:14);
x_hat_lin  = x_true_lin - z_lin(:, 15:28);
u_tot_lin  = (-K * x_hat_lin')' + [M*g, 0, 0, 0];
omega2_lin = (u2omega * u_tot_lin')';

fprintf('Max squared rotor speed in linear closed-loop: %.0f (Limit: 250000)\n', max(omega2_lin, [], 'all'));
fprintf('Min squared rotor speed in linear closed-loop: %.0f (Limit: 0)\n', min(omega2_lin, [], 'all'));
%% discretized controller and observer (phase 3)
Ts = 0.01;

sys_c = ss(A_sys, B_sys, C_sys, D_sys);
sys_d = c2d(sys_c, Ts, 'zoh');

A_d = sys_d.A;
B_d = sys_d.B;
C_d = sys_d.C;
D_d = sys_d.D;

sys_cl_c = ss(A_sys - B_sys * K, B_sys, C_sys, D_sys);
sys_cl_d = c2d(sys_cl_c, Ts, 'zoh');
K_discrete = dlqr(A_d,B_d,Q*Ts,R*Ts);

G_d = eye(14);
B_aug_d = [B_d, G_d];
D_aug_d = [D_d, zeros(6, 14)];
sys_augmented_discrete = ss(A_d, B_aug_d, C_d, D_aug_d, Ts);

Q_obs_d = Q_obs * Ts;
R_obs_d = R_obs / Ts;

[kalmf_discrete, L_discrete, P_discrete, M_discrete] = kalman(sys_augmented_discrete, Q_obs_d, R_obs_d, [], sensors, known_inputs);

%% --- Stability Sanity Checks ---
cl_ctrl_poles = eig(A_d - B_d * K_discrete);
cl_obs_poles  = eig(A_d - L_discrete * C_d);

fprintf('\n--- Phase 3 Discrete Stability Verification ---\n');
fprintf('Max discrete controller pole magnitude: %.4f \n', max(abs(cl_ctrl_poles)));
fprintf('Max discrete observer pole magnitude:   %.4f \n', max(abs(cl_obs_poles)));

assert(max(abs(cl_ctrl_poles)) < 1, 'Discrete controller is unstable!');
assert(max(abs(cl_obs_poles)) < 1, 'Discrete observer is unstable!');