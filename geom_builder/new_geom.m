function [r,d,d2] = new_geom(t,d,A)
% parameterization of funky boundary with period d and amplitude A
omega = 2*pi/d;
r = [t(:) + 1/omega*sin(2*omega*t(:)), A*cos(omega*t(:))].';
d = [ones(length(t),1) + 2*cos(2*omega*t(:)), -omega*A*sin(omega*t(:))].';
d2 = [-4*omega*sin(2*omega*t(:)), -omega^2*A*cos(omega*t(:))].';
end