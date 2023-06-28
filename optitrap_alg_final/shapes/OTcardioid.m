function q=OTcardioid(theta,c)
h_center = 0.1194;
q=[0*theta;c*cos(theta).*(1+cos(theta))-c;c*sin(theta).*(1+cos(theta))+h_center];
end