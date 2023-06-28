function q = OTheart(theta,r)
h_center = 0.1194;

q_y = 16*sin(theta).^3;
q_z = 13*cos(theta)-5*cos(2*theta)-2*cos(3*theta)-cos(4*theta);
q=[0*theta; r*q_y; r*q_z+h_center];
end