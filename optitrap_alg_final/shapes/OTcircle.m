function q = OTcircle(theta,r)
h_center = 0.1194;
theta=theta+pi/2;
q=[0*theta; r*sin(theta); r*(1-cos(theta))+h_center-r];
end