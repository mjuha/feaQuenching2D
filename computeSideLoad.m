function [ke,fe] = computeSideLoad(el,xe,de,isFluxLoad)

global convectionLoad fluxLoad HTableData isNBCTempDependent STATE

% Gauss - Legendre rule. 1 point in [0,1]
gp = 0.5;
w = 1.0;

% compute residual: side loads
if ~isFluxLoad
    edge = convectionLoad(el,2);
    if ~isNBCTempDependent
        h = convectionLoad(el,5); % coefficient
    end
    Ta = convectionLoad(el,6); % ambient temperature
else
    edge = fluxLoad(el,2);
    q = fluxLoad(el,5); % heat flux
end
fe = zeros(3,1);
ke = zeros(3,3);
if edge == 1 % local nodes 1-2
    r = gp;
    % s = 0
    Nshape = [ 1-r, r, 0 ];
    N_r = [ -1, 1, 0 ];
elseif edge == 2 % local nodes 2-3
    r = gp;
    % s = 1-r
    Nshape = [ 0, r, 1-r ];
    N_r = [ 0 ,1, -1 ];
elseif edge == 3 % local nodes 3-1
    % r = 0
    s = gp;
    Nshape = [ 1-s, 0, s ];
    N_r = [ -1, 0, 1 ];
else
    error('Wrong edge, check input!');
end

% compute jacobian
x_r = N_r * xe(:,1);
y_r = N_r * xe(:,2);
jac = sqrt(x_r^2 + y_r^2);
if jac < 1.0e-12
    error('Jacobian less than zero, check input or element too distorted!');
end
    
% surface temperature
Ts = Nshape * de;

if isNBCTempDependent
    h = interp1(HTableData(:,1), HTableData(:,2), Ts, 'linear','extrap');
end
if h < 0.0
    error('Film coefficient is negative!')
end

if ~isFluxLoad
    if STATE(1) == 1 % axisymmetric
        % interpolate radial coordinate
        r = Nshape * xe(:,1);
        fe = fe - Nshape' * ( h * ( Ta - Ts ) ) * r * w * jac;
        ke = ke + Nshape' * h * Nshape * r * w * jac;
    else
        fe = fe - Nshape' * ( h * ( Ta - Ts ) ) * w * jac;
        ke = ke + Nshape' * h * Nshape * w * jac;
    end
else
    if STATE(1) == 1 % axisymmetric
        % interpolate radial coordinate
        r = N * xe(:,1);
        fe = fe - Nshape' * q * r * w * jac;
    else
        fe = fe - Nshape' * q * w * jac;
    end
    %ke = zeros(4,4);
end

end