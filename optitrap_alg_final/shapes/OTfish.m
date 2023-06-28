function q = OTfish(theta,r)
h_center = 0.1194;

q_y=cos(theta)-(sin(theta).^2/sqrt(2));
q_z=sin(theta).*cos(theta);
q=[0*theta; r*q_y; r*q_z+h_center];
end