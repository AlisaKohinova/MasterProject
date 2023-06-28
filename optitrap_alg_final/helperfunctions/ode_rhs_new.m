function result = ode_rhs_new(t, trap_dt, state, trap_u, params)
    scaleTemp = 1e0;
    trap_number = round(t/trap_dt)+1;
    if trap_number == length(trap_u(1,:))+1
        trap_number = length(trap_u(1,:));
    end
    result = [state(4:6); scaleTemp*F_acoustic(state(1:3)/scaleTemp, trap_u(:,trap_number)/scaleTemp, params) / params.m];
end