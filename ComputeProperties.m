function [ rho, kappa, cp ] = ComputeProperties( T, prop, phi )

% get properties
%tmp = prop{1};
% tempData = tmp{8};

% Ta = tempData(1); % starting temperature for transformation of austenite
% Tb = tempData(2); % starting temperature for transformation of perlite

% Austenite
rA = prop{1};
rhoA = interp1(rA(:,1), rA(:,2),T,'linear','extrap');
cA = prop{3};
cpA = interp1(cA(:,1), cA(:,2),T,'linear','extrap');
kA = prop{5};
kappaA = interp1(kA(:,1), kA(:,2),T,'linear','extrap');

% perlite
rA = prop{2};
rhoP = interp1(rA(:,1), rA(:,2),T,'linear','extrap');
cA = prop{4};
cpP = interp1(cA(:,1), cA(:,2),T,'linear','extrap');
kA = prop{6};
kappaP = interp1(kA(:,1), kA(:,2),T,'linear','extrap');

% martensite
%rA = tmp{1};
%rhoM = interp1(rA(:,1), rA(:,2),T,'linear','extrap');
rhoM = rhoA;
% cA = tmp{3};
%cpM = interp1(cA(:,1), cA(:,2),T,'linear','extrap');
cpM = cpA;
kA = prop{7};
kappaM = interp1(kA(:,1), kA(:,2),T,'linear','extrap');

% compute properties based on rule of mixture
% density
rho = rhoA * phi(1) + rhoP * phi(2) + rhoM * phi(3);
% heat capacity
cp = cpA * phi(1) + cpP * phi(2) + cpM * phi(3);
% conductivity
kappa = kappaA * phi(1) + kappaP * phi(2) + kappaM * phi(3);

% if T > Ta % austenite
%     rhoA = tmp{1};
%     rho = interp1(rhoA(:,1), rhoA(:,2),T,'linear','extrap');
% %     if isnan(rho)
% %         error('Temperature is aout of range, please check input data!')
% %     end
%     cpA = tmp{3};
%     cp = interp1(cpA(:,1), cpA(:,2),T,'linear','extrap');
% %     if isnan(cp)
% %         error('Temperature is aout of range, please check input data!')
% %     end
%     kA = tmp{5};
%     kappa = interp1(kA(:,1), kA(:,2),T,'linear','extrap');
% %     if isnan(kappa)
% %         error('Temperature is aout of range, please check input data!')
% %     end
% elseif T > Tb
%     rhoA = tmp{2};
%     rho = interp1(rhoA(:,1), rhoA(:,2),T,'linear','extrap');
% %     if isnan(rho)
% %         error('Temperature is aout of range, please check input data!')
% %     end
%     cpA = tmp{4};
%     cp = interp1(cpA(:,1), cpA(:,2),T,'linear','extrap');
% %     if isnan(cp)
% %         error('Temperature is aout of range, please check input data!')
% %     end
%     kA = tmp{6};
%     kappa = interp1(kA(:,1), kA(:,2),T,'linear','extrap');
%     if isnan(kappa)
%         error('Temperature is aout of range, please check input data!')
%     end
% else % martensite
%     rhoA = tmp{1};
%     rho = interp1(rhoA(:,1), rhoA(:,2),T,'linear','extrap');
% %     if isnan(rho)
% %         error('Temperature is aout of range, please check input data!')
% %     end
%     cpA = tmp{3};
%     cp = interp1(cpA(:,1), cpA(:,2),T,'linear','extrap');
% %     if isnan(cp)
% %         error('Temperature is aout of range, please check input data!')
% %     end
%     kA = tmp{7};
%     kappa = interp1(kA(:,1), kA(:,2),T,'linear','extrap');
% %     if isnan(kappa)
% %         error('Temperature is aout of range, please check input data!')
% %     end
% end


end

