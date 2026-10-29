function state_dot = drone(states, inputs, parameters)
    xi = states(1:3);
    xi_dot = state(4:6);
    eta = state(7:9);
    nu = state(10:12);
    alpha = state(13);
    alpha_dot = state(14);

    u_T = inputs(1);
    tau = inputs(2,4);

    % kinematics
    W_inv = [1,  sin(phi)*tan(theta), cos(phi)*tan(theta);
        0,  cos(phi),            -sin(phi);
        0,  sin(phi)/cos(theta), cos(phi)/cos(theta)];

    eta_dot = W_inv*nu;

    p_dot = (tau(1)/I_xx) + ((I_yy-I_zz)/I_xx)*nu(2)*nu(3);
    q_dot = (tau(2)/I_yy) + ((I_zz-I_xx)/I_yy)*nu(1)*nu(3);
    r_dot = (tau(3)/I_zz) + ((I_xx-I_yy)/I_zz)*nu(1)*nu(2);
    nu_dot = [p_dot; q_dot; r_dot];



    state_dot = [xi_dot, xi_dotdot, eta_dot, nu_dot, alpha_dot, alpha_dotdot];
end