function q = OTflower(theta,r)
h_center = 0.1194;
sc = 0.5;

q_y = sin(3*sc*theta).*cos(sc*theta);
q_z = sin(sc*theta).*sin(3*sc*theta);
q=[0*theta; r*q_y; r*q_z+h_center];
end