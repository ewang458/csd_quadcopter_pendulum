%% physical params defined by the project
m = 1.587; %kg
m_p = 0.25; %kg
g = 9.8065; %m/s^2
A = 0.0625; %m^2
rho = 1.225; %kg/m^3
R = 0.127; %m
A_r = pi * R^2; %m^2
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
measured_indices = [1, 2, 3, 7, 8, 9, 10, 11, 12, 13];
C_sys = eye(14);
D_sys = 0;

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
D_aug = [D_sys, zeros(10, 14)];

Q_obs = eye(14);
Q_obs(13, 13) = 100;
Q_obs(14, 14) = 50;
R_obs = 0.01 * eye(10);
sys_augmented = ss(A_sys, B_aug, C_sys, D_aug);

[kalmf_linear, L_obs, P_linear] = kalman(sys_augmented, Q_obs, R_obs, [], 1:10, 1:4);

%% discretized controller and observer (phase 3)
Ts = 0.02; % 50hz as per pdf
G = eye(14);
sys_augmented_discrete = c2d(sys_augmented, Ts, 'zoh');

A_d = sys_augmented_discrete.A;
B_d = sys_augmented_discrete.B(:, 1:4);
C_d = sys_augmented_discrete.C;

K_discrete = dlqr(A_d, B_d, Q, R);
[kalmf_discrete, L_discrete, P_discrete, M_discrete] = kalman(sys_augmented_discrete, Q_obs, R_obs);
