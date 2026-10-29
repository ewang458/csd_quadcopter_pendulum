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

states = [x; 
          y; 
          z; 
          x_dot;
          y_dot;
          z_dot;
          phi;
          theta;
          psi;
          p;
          q;
          r;
          alpha;
          alpha_dot];

inputs = [u_T;
          u_phi;
          u_theta;
          u_psi];

%% 15.2.1 coord frames and generalized coordinates
xi = [x; 
      y; 
      z]; % positions

eta = [phi;
       theta;
       psi]; % roll, pitch, yaw

nu = [p;
      q;
      r]; % body angular velocity

%% 15.2.2 rotation and rate relations
W = [1,  0,          -sin(thetalim);
     0,  cos(theta), sin(phi)*cos(theta);
     0, -sin(theta), cos(phi)*cos(theta)];

W_inv = [1,  sin(phi)*tan(theta), cos(phi)*tan(theta);
         0,  cos(phi),            -sin(phi);
         0,  sin(phi)/cos(theta), cos(phi)/cos(theta)];

eta_dot = W_inv*nu;

nu_b = [u;
       v;
       w]; % linear velocity of quadcopter

xi_dot = R_body_to_inertial * nu_b;

%% 15.2.3 Force Model
F_g = [0;
      0;
      -mg];

F_thrust = [0;
            0;
            u_T];

F_thrust_inertial_coords = R_body_to_inertial * F_thrust;

F_D = -k_D * [x_dot;y_dot;z_dot]; %linear drag

%% 15.2.4
tau = [u_phi;
       u_theta;
       u_psi];

omega_1;
omega_2;
omega_3;
omega_4;

u_thrust = k*(omega_1^2 + omega_2^2 + omega_3^2 + omega_4^2);

%% 15.2.5

%% 15.2.6 Pendulum Dynamics
alpha = 0
[x_p; z_p] = [x; z] + l * [sin(alpha); c_alpha];
[x_p_dotdot; z_p_dotdot] = [x_dotdot; z_dotdot] + l * [-sin(alpha)*alpha_dot^2 + cos(alpha)*alpha_dotdot; -cos(alpha)*alpha_dot^2 - sin(alpha)*alpha_dotdot]
m_p * [x_p_dotdot; z_p_dotdot] = [0; -m_p * g] - [T*sin(alpha); T*cos(alpha)];
l * alpha_dotdot = -1 * x_p_dotdot*cos(alpha) - z_p_dotdot*sin(alpha) + g*sin(alpha);
%% 15.2.7
F_p = m_p * (l*alpha_dot^2 - x_p_dotdot*sin(alpha) - z_p_dotdot*cos(alpha) - g*cos(alpha)) * [sin(alpha); 0; cos(alpha)];

%% 15.2.8
