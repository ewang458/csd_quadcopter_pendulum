function state_dot = drone(states, inputs, parameters)
    m     = parameters.m;
    m_p    = parameters.m_p;
    M     = m + m_p;
    g     = parameters.g;
    l     = parameters.l;
    k_D    = parameters.k_D;
    I_xx  = parameters.I_xx;
    I_yy  = parameters.I_yy;
    I_zz  = parameters.I_zz;

    xi = states(1:3);
    xi_dot = states(4:6);
    eta = states(7:9);
    nu = states(10:12);
    alpha = states(13);
    alpha_dot = states(14);

    u_T = inputs(1);
    tau = inputs(2:4);

    x_dot = xi_dot(1);
    y_dot = xi_dot(2);
    z_dot = xi_dot(3);

    % kinematics
    phi = eta(1);
    theta = eta(2);
    psi = eta(3);
    W_inv = [1,  sin(phi)*tan(theta), cos(phi)*tan(theta);
        0,  cos(phi),            -sin(phi);
        0,  sin(phi)/cos(theta), cos(phi)/cos(theta)];

    eta_dot = W_inv*nu;

    %rotational dynamics
    p_dot = (tau(1)/I_xx) + ((I_yy-I_zz)/I_xx)*nu(2)*nu(3);
    q_dot = (tau(2)/I_yy) + ((I_zz-I_xx)/I_yy)*nu(1)*nu(3);
    r_dot = (tau(3)/I_zz) + ((I_xx-I_yy)/I_zz)*nu(1)*nu(2);
    nu_dot = [p_dot; q_dot; r_dot];

    %translational dynamics
    c_alpha = cos(alpha);
    s_alpha = sin(alpha);
    c_phi = cos(phi); s_phi = sin(phi);
    c_theta = cos(theta); s_theta = sin(theta);
    c_psi = cos(psi); s_psi = sin(psi);

    Fx = (c_psi*s_theta*c_phi + s_psi*s_phi) * u_T;
    Fy = (s_psi*s_theta*c_phi - c_psi*s_phi) * u_T;
    Fz = c_theta*c_phi * u_T;

    x_dotdot = (1/(m*M)) * ( (m + m_p*c_alpha^2)*(Fx - k_D*x_dot) - m_p*s_alpha*c_alpha*(Fz - k_D*z_dot) + m*m_p*l*alpha_dot^2*s_alpha );
    y_dotdot = (Fy - k_D*y_dot) / m;
    z_dotdot = (1/(m*M)) * ( (m + m_p*s_alpha^2)*(Fz - k_D*z_dot) - m_p*s_alpha*c_alpha*(Fx - k_D*x_dot) + m*m_p*l*alpha_dot^2*c_alpha ) - g;
    xi_dotdot = [x_dotdot; y_dotdot; z_dotdot];

    % pendulum dynamics
    alpha_dotdot = (x_dotdot*c_alpha - z_dotdot*s_alpha + g*s_alpha) / l;

    % assemble state_dot
    state_dot = [xi_dot; xi_dotdot; eta_dot; nu_dot; alpha_dot; alpha_dotdot];
end