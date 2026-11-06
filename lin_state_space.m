%% non-linear state space
function state_dot = lin_state_space(states, inputs, parameters)
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

    phi = eta(1);
    theta = eta(2);
    psi = eta(3);

    u_phi = tau(1);
    u_theta = tau(2);
    u_psi = tau(3);
    p = nu(1);
    q = nu(2);
    r = nu(3);
    
    x_dot = xi_dot(1);
    y_dot = xi_dot(2);
    z_dot = xi_dot(3);

    x_dotdot = (theta*M*g - k_D*x_dot - m_p*g*alpha)/m;
    y_dotdot = (-phi*M*g - k_D*y_dot)/m;
    z_dotdot = ((u_T - M*g) - k_D*z_dot)/M;

    phi_dot = p;
    theta_dot = q;
    psi_dot = r;

    p_dot = u_phi/I_xx;
    q_dot = u_theta/I_yy;
    r_dot = u_psi/I_zz;

    eta_dot = [phi_dot; theta_dot; psi_dot];
    nu_dot = [p_dot; q_dot; r_dot];
    xi_dotdot = [x_dotdot; y_dotdot; z_dotdot];
    alpha_dotdot = (1/(M*l)) * (M*g*(alpha - theta) + k_D*x_dot);

    % assemble state_dot
    state_dot = [xi_dot; xi_dotdot; eta_dot; nu_dot; alpha_dot; alpha_dotdot];

end