function [bead_sim_forces_pos, bead_sim_forces_vel] = simulation_force(opt_T, N, forces, params, initial_state)
    deltaT = opt_T / N;         % time discretization
    bead_sim_forces_pos = zeros(3, N+1); % bead simulated position with forces from the NLP solver
    bead_sim_forces_vel = zeros(3, N+1); % bead simulated velocity with forces from the NLP solver
    bead_sim_forces_pos(:,1) = initial_state(1:3); % initial condition added
    bead_sim_forces_vel(:,1) = initial_state(4:6); % initial condition added
    sim_state = [bead_sim_forces_pos(:,1); bead_sim_forces_vel(:,1)]; % state = bead pos + vel
    for k=1:N
        tmp = sim_state;
        k1 = [tmp(4:6); forces(:,k) / params.m];
        tmp = sim_state + deltaT/2*k1;
        k2 = [tmp(4:6); forces(:,k) / params.m];
        tmp = sim_state + deltaT/2*k2;
        k3 = [tmp(4:6); forces(:,k) / params.m];
        tmp = sim_state + deltaT*k3;
        k4 = [tmp(4:6); forces(:,k) / params.m];
        sim_state = sim_state + deltaT/6*(k1+2*k2+2*k3+k4);
%         sim_state = sim_state + deltaT*k1;
        bead_sim_forces_pos(:,k+1) = sim_state(1:3);
        bead_sim_forces_vel(:,k+1) = sim_state(4:6);
    end
end