function decision_vars = setup_final(mode, opti, N, parameters_nlp_values, fCasJfun, fCasJJfun)
    import casadi.*

    nZ = 2;
    n_zeta = 4;

    % initial conditions
    theta_0 = 0;
    % terminal constraints
    if strcmp(mode,'rampup') || strcmp(mode,'rampdown')
        theta_final = 4*pi;
    else 
        theta_final = 2*pi;
    end
    
    % ---- decision variables ---------
    Z = opti.variable(nZ,N+1); % states of virtual system
    theta = Z(1,:);
    theta_dot = Z(2,:);
    v = opti.variable(N)';   % virtual control input 
    T = opti.variable();      % final time
    z_N = opti.variable(n_zeta,N); 

    q_d = fCasJfun(theta);
    q_dd = fCasJJfun(theta);
    
    decision_vars = struct( ...
        'Z',Z,      ... % virtual state [m, m/s]
        'v',v,      ... % virtual control [m/s^2]
        'T',T,      ... % time [s]
        'z_N',z_N   ... % virtual variables for sin/cos terms.
    );
    
    % ---- objective ------------------
    opti.minimize(T+parameters_nlp_values.c_par/N*sum(v.^2)); %minimal time + regularization term

    % ---- dynamic constraints --------
    f1 = @(x,v) [x(2); v]; % ODE virtual system 
    
    dt = T/N;                                                              
    for k=1:N % loop over control intervals
        % Runge-Kutta 4 integration
        k1 = f1(Z(:,k),         v(k));
        k2 = f1(Z(:,k)+dt/2*k1, v(k));
        k3 = f1(Z(:,k)+dt/2*k2, v(k));
        k4 = f1(Z(:,k)+dt*k3,   v(k));                                      %FRAGE: Warum nicht v(:,k+1)?
        
        z_next = Z(:,k) + dt/6*(k1+2*k2+2*k3+k4); 
        con=Z(:,k+1)==z_next;
        opti.subject_to(con); % close the gaps

        % --- !!! constrains on the z_N terms !!! --- 
        % 1) Constrain relevant z_N's to make inversion feasible:
        opti.subject_to(-parameters_nlp_values.arc_lim <= z_N(1:4,k) <= parameters_nlp_values.arc_lim); 
        
        % 2) trigonometric Pythagoras to give some knowledge about the
        % problem structure to the optimizer:
        opti.subject_to(z_N(2,k)^2 + z_N(3,k)^2 == 1)

        % --- !!! Particle dynamics as constraints !!! ---      
        dyn_lhs = parameters_nlp_values.m*(q_dd(:,k).*theta_dot(k)^2 + q_d(:,k).*v(k)); 
        con6 = dyn_lhs(2:3) - [parameters_nlp_values.A_r*z_N(1,k)*z_N(2,k); parameters_nlp_values.A_z*z_N(3,k)*z_N(4,k)];
        opti.subject_to(con6(:) == 0);
    end
    
    % ---- boundary conditions --------
    % theta position
    opti.subject_to(theta(1)==theta_0);
    opti.subject_to(theta(N+1)==theta_final); 
    % theta velocity                                                                       
    if strcmp(mode,'periodic')
        opti.subject_to(theta_dot(1)==theta_dot(N+1)); 
    elseif strcmp(mode,'rampup')
        opti.subject_to(theta(ceil(2/3*N))==2*pi);
        opti.subject_to(theta_dot(1)==0); 
        opti.subject_to(theta_dot(ceil(2/3*N))==theta_dot(N+1)); 
    elseif strcmp(mode,'rampdown')
        opti.subject_to(theta(ceil(1/3*N))==2*pi);
        opti.subject_to(theta_dot(N+1)==0); 
        opti.subject_to(theta_dot(1)==theta_dot(ceil(1/3*N))); 
    else
        opti.subject_to(theta_dot(1)==0); % from stand-still 
        opti.subject_to(theta_dot(N+1)==0); % to stand-still 
    end

    % ---- misc. constraints  ----------
    opti.subject_to(T>=0); % time must be positive
end