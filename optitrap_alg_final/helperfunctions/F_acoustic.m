function force = F_acoustic(bead_p, trap_u, params)
    dr = sqrt((trap_u(1) - bead_p(1))^2 + (trap_u(2) - bead_p(2))^2);
    dz = trap_u(3) - bead_p(3);
    F_z = params.A_z * sin(dz * (2*pi)/params.c_z2) * cos(dr * (2*pi)/params.c_z6);
    if abs(dr) <= 1e-16
        force = [0; 0; F_z];
    else
        F_r = params.A_r * sin(dr * (2*pi)/params.c_z1) * cos(dz * (2*pi)/params.c_z2);
        force = [F_r * (trap_u(1) - bead_p(1))/dr; F_r * (trap_u(2) - bead_p(2))/dr; F_z];
    end
end
