function [deltas_num, force_num] = sys_inv_num(N, params, force_nlp_unscaled) % numerically solving for delta_x, delta_y, delta_z by solving 0=M(...). 
    flaags = 10*ones(1,N);
    options2 = optimoptions('fsolve','Display','off','Algorithm','levenberg-marquardt','OptimalityTolerance',1e-12,'FunctionTolerance',1e-8);
    options = optimoptions('fsolve','Display','final-detailed','Algorithm','trust-region-dogleg','OptimalityTolerance',1e-12,'FunctionTolerance',1e-8);

    force_num = zeros(3,N);
        deltas_num = zeros(2,N);
        [deltas_num(:,1),fval,flaags(1),ouput] = fsolve(@(dd)eqs_2d_ver(dd, params, force_nlp_unscaled(:,1)), 0.5*ones(2,1), options);
        force_num(2:3,1) = f_2d_ver(deltas_num(:,1), params);
        for k=2:N
            [deltas_num(:,k),fval,flaags(k),ouput] = fsolve(@(dd)eqs_2d_ver(dd, params, force_nlp_unscaled(:,k)), 0.5*ones(2,1), options);
            if flaags(k) < 1
                [deltas_num(:,k),fval,flaags(k),ouput] = fsolve(@(dd)eqs_2d_ver(dd, params, force_nlp_unscaled(:,k)), deltas_num(:,k-1), options2);
            end
            force_num(2:3,k) = f_2d_ver(deltas_num(:,k), params);
        end
end

function eqs = eqs_2d_ver(dd, params, force_nlp_current)
    force = f_2d_ver(dd, params);
    eqs = force - force_nlp_current(2:3);
end

function force_unscaled = f_2d_ver(dd, params)
    force_unscaled = [sin((2*pi)/params.c_z1*sqrt(dd(1)^2))*cos((2*pi)/params.c_z2*dd(2))*dd(1)/sqrt(dd(1)^2); 
                    sin((2*pi)/params.c_z2*dd(2))*cos((2*pi)/params.c_z6*sqrt(dd(1)^2))];
end