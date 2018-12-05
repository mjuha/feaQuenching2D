function [N,dN,jac,inv_jac] = shape(gp,xe)

% local coordinate
r = gp(1);
s = gp(2);

% shape functions
N = [ 1-r-s, r, s ];
N_r = [ -1, 1, 0 ];
N_s = [-1, 0, 1];

dN = zeros(3,2);

x_r = N_r * xe(:,1);
x_s = N_s * xe(:,1);
y_r = N_r * xe(:,2);
y_s = N_s * xe(:,2);

jacobian = [x_r, x_s;  y_r, y_s];
jac = det(jacobian);

% check jacobian
if jac < 1.0e-14
    error('Negative jacobian, element too distorted!');
end

inv_jac = inv(jacobian);
% Note: inv_jac = [r_x, r_y; s_x, s_y]

for i=1:3
    dN(i,:) = [N_r(i), N_s(i)] * inv_jac; %#ok<MINV>
end

end
