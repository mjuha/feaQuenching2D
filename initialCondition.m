function initialCondition(IC)

global U nn Phase PhaseOld scheilRule

% solution vector
U = zeros(2,nn);

% Phase content
Phase = zeros(3,nn);
Phase(1,:) = 0.99; % austenite
%
PhaseOld = zeros(3,nn);
PhaseOld(1,:) = 0.99; % austenite
%
% Scheil's rule
scheilRule = zeros(1,nn);

name = IC{1};

if strcmp(name,'Constant')
    value = IC{2};
end

for i = 1:nn
    U(1,i) = value;
end


end