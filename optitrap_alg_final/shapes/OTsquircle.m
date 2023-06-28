function q=OTsquircle(theta,r)
s = 0.9;
h_center = 0.1194;
q=[0*theta; r*cos(theta); r*sin(theta)./(sqrt(1-s*cos(theta).^2))+h_center];
end