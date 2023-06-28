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

%1. CONFIGURE INPUT (Shape, speed, size)
% ---- pick a shape ----------------
POV_freq = 10;
num_shape = 2;                                                     

switch num_shape
    case 1
        shape = "circle"; 
        size_coeff = 0.035;
        POV_freq = 15;
    case 2
        shape = "cardioid"; 
        size_coeff = 0.035;
    case 3
        shape = "squircle";
        size_coeff = 0.0265;
    case 4
        shape = "fish";
        size_coeff = 0.0425;
    case 5
        shape = "flower";
        size_coeff = 0.036;
    case 6
        shape = "cat";
        size_coeff = 3e-5;
    case 7
        shape = "dolphin";
        size_coeff = 2e-5;
    case 8
        shape = "heart";
        size_coeff = 0.0017;
    otherwise
        warning('Unexpected number.')
end

% ---- pick a mode ----------------
num_mode = 1;

switch num_mode
    case 0
        mode = 'rest_to_rest';
    case 1
        mode = 'periodic';    
    case 2
        mode = 'rampup';  
    case 3
        mode = 'rampdown';  
end

% ---- parameters----------------
if strcmp(mode,'periodic')
     N = ceil(10000/POV_freq);
else 
     N = ceil(10000/POV_freq)*3;
end

%2. CONFIGURE SYSTEM (Forces, agressiveness)
%2.a Parameters to colve the NLP: We need this to determine timing 
%Parameters defining the non-linear problem: force model, regularization
%term in cost function (remove kinks), arc_lim inversion parameters (closer
%to 1 more agressive, but less stable).

parameters_nlp_values = struct( ...
    'm',        0.07  ,... % mass of the particle [mg]
    'A_r',      0.0243, ... % maximum horizontal force [mN]
    'A_z',      0.042, ... % maximum vertical force [mN]
    'c_par',    1e8,... % weight of the control cost term in the cost function [no unit]
    'arc_lim',  0.85    ... % threshold for the virtual sin/cos terms for easier inversion later [no unit]; the closer to 1, the more numerical errors in fsolve/inversion
);

%2.b. Parameters for inversion (i.e compute trap placement once we get timing right).   
parameters_simulation_values = struct( ...
    'm',        parameters_nlp_values.m, ...
    'A_r',      parameters_nlp_values.A_r, ...
    'A_z',      parameters_nlp_values.A_z, ...
    'c_z1',     13.23667588,   ...  % 2*pi/c_z1 = 474.68
    'c_z2',     004.901845302,  ...  % 2*pi/c_z2 = 1281.8
    'c_z5',     004.901845302,  ...  % 2*pi/c_z5 = 1281.8 [same as c_z2]
    'c_z6',     21.56428358 ... % 2*pi/c_z6 = 291.37
);
% ---- parameters END -------------

%3. CREATE FUNCTIONS DEFINING THE SHAPE
% ---- compute theta derivatives of the reference path
shape_name = str2func(sprintf('OT%s',shape));
thSym = SX.sym('th'); %Create a symbolic variable (thSym) 

fCas = 1e3*shape_name(thSym, size_coeff); %Function defining shape
fCasJ = jacobian(fCas,thSym);             %First and second derivatives wrt "thSym" (not time)  
fCasJJ = jacobian(fCasJ,thSym);              
fCasfun = Function('fCasfun',{thSym},{fCas});   %Callable functions: e.g. fCasFun(0.1,0.5)
fCasJfun = Function('fCasJfun',{thSym},{fCasJ});
fCasJJfun = Function('fCasJJfun',{thSym},{fCasJJ});

% ---- setup NLP- ALL DECLARATIVE; NO COMPUTATIONG YET--------------------
disp('(1/6) Setting up NLP')
tic
decision_vars = setup_final(mode, opti, N, parameters_nlp_values, fCasJfun, fCasJJfun);
toc

% decision variables
Z = decision_vars.Z;        % states of virtual system
v = decision_vars.v;        % virtual control input 
T = decision_vars.T;        % final time
z_N = decision_vars.z_N;    % auxiliary variables for the force computation

% ---- initial guess for solver ---- %AF: If 'Infeasible Problem', try
% changing the initial value for the *last* z value. Do not change the
% formula for the *second to last* z value!

if load_yes %start optimization from a previous valid solution
    filename = sprintf('data/ocp_valid_%s_%d_%s.mat', shape, N, mode);
    load(filename);
    opti.set_initial(Z, Z_sol);
    opti.set_initial(v, v_sol);
    opti.set_initial(T, opt_T);
    opti.set_initial(z_N, z_N_sol);
else
    opti.set_initial(Z,zeros(size(Z)));
    opti.set_initial(v,zeros(size(v)));
    opti.set_initial(T, 500);

    tmpinit = min(0.5,parameters_nlp_values.arc_lim);
    opti.set_initial(z_N(1,:), zeros(1,N));
    opti.set_initial(z_N(2,:), tmpinit*ones(1,N)); 
    opti.set_initial(z_N(3,:), sqrt(1-(tmpinit*ones(1,N)).^2));
    opti.set_initial(z_N(4,:), zeros(1,N));
end

% ---- solve NLP -------------------
disp('(2/6) Solving NLP')
tic
p_opts = struct('expand',true);
s_opts = struct('max_iter',50000); % change max iteration 
opti.solver('ipopt',p_opts, s_opts); % set numerical backend
sol = opti.solve();
toc

% ---- Get and save solution ---------
Z_sol = sol.value(Z);
theta_sol = Z_sol(1,:);
theta_dot_sol = Z_sol(2,:);
opt_T = sol.value(T);
v_sol = sol.value(v);
z_N_sol = sol.value(z_N);
                                                                           
% if save_yes
%     filename = sprintf('data/ocp_valid_%s_%d_%s.mat', shape, N, mode);
%     save(filename,'v_sol','Z_sol','z_N_sol','opt_T');
% end

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