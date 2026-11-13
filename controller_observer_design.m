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
Q = diag([10,  10,  20,   ...  % x, y, z
          2,   2,  2,    ...  % x_dot, y_dot, z_dot
          50,  50,  10,   ...  % phi, theta, psi
          5,  5,  10,    ...  % p, q, r
          150, 20         ...  % alpha, alpha_dot
          ]); % tunable!

% R penalizes control inputs, tunable!
R = diag([0.1, 1, 1, 1]); % delta_u_T, u_phi, u_theta, u_psi

K = lqr(A_sys, B_sys, Q, R);

G = eye(14);
B_aug = [B_sys, G];
D_aug = [D_sys, zeros(6, 14)];
Q_obs = eye(14);
Q_obs(13, 13) = 100;
Q_obs(14, 14) = 50;
R_obs = 0.01 * eye(6);

sys_augmented = ss(A_sys, B_aug, C_sys, D_aug);

sensors = 1:6;
known_inputs = 1:4;
[kalmf_linear, L_obs, P_linear] = kalman(sys_augmented, Q_obs, R_obs, [], sensors, known_inputs);

%% check
cl_ctrl_poles_c = eig(A_sys - B_sys * K);
cl_obs_poles_c  = eig(A_sys - L_obs * C_sys);

fprintf('Max real part of controller poles: %.4f (Must be < 0)\n', max(real(cl_ctrl_poles_c)));
fprintf('Max real part of observer poles:   %.4f (Must be < 0)\n', max(real(cl_obs_poles_c)));

%% discretized controller and observer (phase 3)
Ts = 0.02;
% 1. Discretize the physical plant via ZOH
sys_c = ss(A_sys, B_sys, C_sys, D_sys);
sys_d = c2d(sys_c, Ts, 'zoh');

A_d = sys_d.A;
B_d = sys_d.B;
C_d = sys_d.C;
D_d = sys_d.D;

% 2. Discretize the Phase 2 Continuous Controller K via ZOH step-invariance
%    (Matches the ZOH discrete closed-loop transition to c2d of the continuous closed-loop)
sys_cl_c = ss(A_sys - B_sys * K, B_sys, C_sys, D_sys);
sys_cl_d = c2d(sys_cl_c, Ts, 'zoh');
K_discrete = pinv(B_d) * (A_d - sys_cl_d.A);

% 3. Discrete Observer via kalman()
G_d = eye(14);
B_aug_d = [B_d, G_d];
D_aug_d = [D_d, zeros(6, 14)];
sys_augmented_discrete = ss(A_d, B_aug_d, C_d, D_aug_d, Ts);

Q_obs_d = eye(14);
Q_obs_d(13, 13) = 1e4;
Q_obs_d(14, 14) = 1e4;
R_obs_d = 1e-4 * eye(6);

[kalmf_discrete, L_discrete, P_discrete, M_discrete] = kalman(sys_augmented_discrete, Q_obs_d, R_obs_d, [], sensors, known_inputs);

%% --- Stability Sanity Checks ---
cl_ctrl_poles = eig(A_d - B_d * K_discrete);
cl_obs_poles  = eig(A_d - L_discrete * C_d);

fprintf('\n--- Phase 3 Discrete Stability Verification ---\n');
fprintf('Max discrete controller pole magnitude: %.4f (Must be < 1)\n', max(abs(cl_ctrl_poles)));
fprintf('Max discrete observer pole magnitude:   %.4f (Must be < 1)\n', max(abs(cl_obs_poles)));

assert(max(abs(cl_ctrl_poles)) < 1, 'Discrete controller is unstable!');
assert(max(abs(cl_obs_poles)) < 1, 'Discrete observer is unstable!');
disp('SUCCESS: Both ZOH-discretized controller and Kalman observer are stable!');