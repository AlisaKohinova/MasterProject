%% Minimum Time Path for a Levitating Particle in an a Acoustic Field

% ---- initialize ----------------
clear vars; close all; clc;
cd 'C:\Users\chris\Documents\GitHub\MasterProject\optitrap_alg_final' %CHANGE

addpath('C:\Users\chris\Documents\GitHub\MasterProject\casadi-3.6.3-windows64-matlab2018b') %CHANGE This is the solver
import casadi.*
opti = casadi.Opti();

addpath(genpath('shapes'))  % add path to matlab files containing the shapes
addpath(genpath('helperfunctions'))  % add path to matlab files containing helper functions

set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

save_yes = true;
load_yes = false;

%1. CONFIGURE INPUT
% in this case, we already have a csv file with periodic data in it.
% we can skip the pre-processing step by directly setting bead_nlp_pos to
% the positions in our csv files


% ---- post-processing ------------
% 1a) Preliminaries: extract particle position, velocity, and acceleration
bead_nlp_pos = full(fCasfun(theta_sol));
bead_nlp_vel = full(fCasJfun(theta_sol)) .* theta_dot_sol;
bead_nlp_acc = full(fCasJJfun(theta_sol)) .* theta_dot_sol.^2 + full(fCasJfun(theta_sol)) .* [v_sol v_sol(1)]; % AF: dim(v_sol) = dim(theta_sol) - 1

% 1b) Preliminaries: extract forces
force_nlp_unscaled = [zeros(size(z_N(1,:))); ...
                          z_N_sol(1,:).*z_N_sol(2,:); ...
                          z_N_sol(3,:).*z_N_sol(4,:)];

force_nlp = [parameters_nlp_values.A_r; parameters_nlp_values.A_r; parameters_nlp_values.A_z] .* force_nlp_unscaled; 

% 3) Perform system inversion
disp('(4/6) System inversion')
[dd_num, force_num_unscaled] = sys_inv_num(N, parameters_simulation_values, force_nlp_unscaled);
force_num = [parameters_nlp_values.A_r; parameters_nlp_values.A_r; parameters_nlp_values.A_z] .* force_num_unscaled; 

% plot deltas
figure
plot(dd_num','-o')
legend('dy-num','dz-num');

% plot forces
figure
plot(force_num','-')
title('Acoustic Force (F_x, F_y, F_z) after numerical system inversion')
legend('Fx-num', 'Fy-num', 'Fz-num')

% plot force differences
figure
plot(abs(force_num' - force_nlp'))
title('|F_{num} - F_{nlp}| (NOT scaled with A_r, A_z)')
legend('Fx-diff', 'Fy-diff', 'Fz-diff')

% 4) Calculate, save, and plot trap positions -- Assumption based on Diego's Excel: delta = u - p
% 4a) Calculate trap positions
disp('(5/6) Calculate and save trap positions')
traps_num = zeros(3,N);
traps_num(2:3, :) = bead_nlp_pos(2:3, 1:end-1) + dd_num;

% 4b) Save trap positions
csvwrite(sprintf('data/u_trap_solutions_%s_fsolve.csv', shape), traps_num)
traps = traps_num./1000;
traps(4,:) = 1;
traps_line = reshape(traps,[],1)'/1;
filename = sprintf('data/%s_%s_%d_ready_to_play.csv', shape, mode, N);
csvwrite(filename,round(traps_line,5))

if save_yes
    filename = sprintf('data/ocp_valid_%s_%d_%s.mat', shape, N, mode);
    save(filename,'v_sol','Z_sol','z_N_sol','opt_T', 'traps_num');
end

% 4c) Plot trap positions
scaleFactor = 1e0;
figure
plot3(scaleFactor*traps_num(1,:), scaleFactor*traps_num(2,:), scaleFactor*(traps_num(3,:)), '-o')
hold on
plot3(scaleFactor*bead_nlp_pos(1,:), scaleFactor*bead_nlp_pos(2,:), scaleFactor*(bead_nlp_pos(3,:)),'-o')
grid on
xlabel('x')
ylabel('y')
zlabel('z')
for i = 1:10:length(traps_num(1,:))
    text(scaleFactor*traps_num(1,i), scaleFactor*(traps_num(2,i)+0.0001), scaleFactor*(traps_num(3,i)+0.0001), num2str(i));
end
text(scaleFactor*traps_num(1,N), scaleFactor*(traps_num(2,N)+0.0001), scaleFactor*(traps_num(3,N)+0.0001), num2str(N));
hold off
title('trap and bead [in mm]')
legend('trap','bead')
view(90,0)
saveas(gcf,'Bead-Trap-Num.png')

% 5) Simulation with obtained trap position and corresponding plots
% 5a) Simulation
disp('(6/6) Simulation with obtained trap positions')
scaleFactor = 1e0;
deltaT = opt_T / N;         % time discretization
force_trap = zeros(3, N);       % Forces based on traps + simulated bead positions
force_trap_nlp = zeros(3, N);   % Forces based on traps + nlp bead positions
bead_sim_trap_pos = zeros(3, N+1);
bead_sim_trap_vel = zeros(3, N+1);
bead_sim_trap_pos(:,1) = scaleFactor * bead_nlp_pos(:,1); % initial condition added
bead_sim_trap_vel(:,1) = bead_nlp_vel(:,1); % initial condition added (CAREFUL: scaleFactor not implemented!)
trap_u = traps_num;% - [0;0;1e3*0.1194];
sim_state = [bead_sim_trap_pos(:,1); bead_sim_trap_vel(:,1)];
% - dynamics: see function F_acoustic
for k=1:N
    force_trap(:,k) = F_acoustic(sim_state(1:3), trap_u(:,k), parameters_simulation_values);
    force_trap_nlp(:,k) = F_acoustic(bead_nlp_pos(:,k), trap_u(:,k), parameters_simulation_values);

    k1 = ode_rhs_new((k-1)*deltaT, deltaT, sim_state,               trap_u, parameters_simulation_values);
    k1(1:3) = k1(1:3) / scaleFactor;
    k2 = ode_rhs_new((k-1)*deltaT, deltaT, sim_state + deltaT/2*k1, trap_u, parameters_simulation_values);
    k2(1:3) = k2(1:3) / scaleFactor;
    k3 = ode_rhs_new((k-1)*deltaT, deltaT, sim_state + deltaT/2*k2, trap_u, parameters_simulation_values);
    k3(1:3) = k3(1:3) / scaleFactor;
    k4 = ode_rhs_new((k-1)*deltaT, deltaT, sim_state + deltaT*k3,   trap_u, parameters_simulation_values);

    sim_state = sim_state + deltaT/6*(k1+2*k2+2*k3+k4);
    bead_sim_trap_pos(:,k+1) = sim_state(1:3);
    bead_sim_trap_vel(:,k+1) = sim_state(4:6);
end

% 5b) Plots of the simulation
figure
quiver(bead_sim_trap_pos(2,:), bead_sim_trap_pos(3,:), bead_sim_trap_vel(2,:), bead_sim_trap_vel(3,:), 0)
title('Particle position and velocity (trap simulation)')

figure()
scaleFactor = 1;
plot3(scaleFactor*trap_u(1,:), scaleFactor*trap_u(2,:), scaleFactor*(trap_u(3,:)), '-o')
hold on
plot3(scaleFactor*bead_nlp_pos(1,:), scaleFactor*bead_nlp_pos(2,:), scaleFactor*(bead_nlp_pos(3,:)),'-o')
hold on
plot3(scaleFactor*bead_sim_trap_pos(1,1:end), scaleFactor*bead_sim_trap_pos(2,1:end), scaleFactor*(bead_sim_trap_pos(3,1:end)),'-o')
grid on
xlabel('x')
ylabel('y')
zlabel('z')
hold off
title(sprintf('[SIMULATION] trap and bead [in mm] - %s, N=%d traps, optT=%dms',shape, N, round(opt_T,3,'significant')))
legend('trap','bead','bead (sim)')
view(90,0)
saveas(gcf,'Bead-Trap-Sim.png')