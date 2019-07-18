function [fe,me,ke] = weakform(el,xe,de,deOld,ae,pe,peOld,~)

global convectionLoad MAT fluxLoad STATE TS radiationLoad

dt = TS{1};

% 1 point formula - degree of precision 1
gp =  [ 1/3, 1/3];
w = 0.5;

% get material properties
%prop = cell2mat(MAT(matNum));
% only one material
prop = MAT(1,:);

ke = zeros(3,3);
me = zeros(3,3);
fe = zeros(3,1);
for i = 1:length(w)
    % stress-strain displacement matrix
    B = zeros(2,3);
    % loop over gauss points
    [N,dN,jac] = shape(gp(i,:),xe);
    % compute temperature
    T = N * de;
    % compute tempertature t
    T1 = N * deOld(1,:)';
    % compute tempertature t-dt
    T2 = N * deOld(2,:)';
    % interpolate phases
    phi = [ dot(N,pe(1,:)); dot(N,pe(2,:)); dot(N,pe(3,:)) ];
    phiOld = [ dot(N,peOld(1,:)); dot(N,peOld(2,:)); dot(N,peOld(3,:)) ];
    % compute properties
    [ rho, k, cp ] = ComputeProperties( T, prop, phi );
    for j=1:3 % loop over local nodes
        B(:,j) = dN(j,:);
    end
    if STATE(1) == 1 % axisymmetric
        % interpolate radial coordinate
        r = N * xe(:,1);
        ke = ke + B' * k * B * r * w(i) * jac;
        %
        me = me + N' * (rho*cp) * N * r * w(i) * jac;
        % latent heat
        [H1] = LatentHeat(T1);
        [H2] = LatentHeat(T2);
        Hp1 = H1(1);
        Hp2 = H2(1);
        Hm = H1(2);
        %
        Q = ( (1.5*Hp1 - 0.5*Hp2) * ( phi(2) - phiOld(2) ) + ...
            Hm * ( phi(3) - phiOld(3) ) ) / dt;
        
        fe = fe - N' * Q * r * w(i) * jac; 
    else
        ke = ke + B' * k * B * w(i) * jac;
        %
        me = me + N' * (rho*cp) * N * w(i) * jac;
        % latent heat
        [H1] = LatentHeat(T1);
        [H2] = LatentHeat(T2);
        Hp1 = H1(1);
        Hp2 = H2(1);
        Hm = H1(2);
        %
        Q = ( (1.5*Hp1 - 0.5*Hp2) * ( phi(2) - phiOld(2) ) + ...
            Hm * ( phi(3) - phiOld(3) ) ) / dt;
        fe = fe - N' * Q * w(i) * jac; 
    end
end

% add contribution from temperature gradient
fe = fe + ke * de;

if size(convectionLoad,1) > 0
    index = find(convectionLoad(:,1)==el,1); % 1 face
    flag = size(index,1);
    % compute side load
    if flag > 0
        [fe1] = computeSideLoad(index,xe,deOld(1,:)','convection');
        fe = fe + fe1;
        %ke = ke + ke1;
    end
end

% now flux load
if size(fluxLoad,1) > 0
    index = find(fluxLoad(:,1)==el,1); % 1 face
    flag = size(index,1);
    % compute side load
    if flag > 0
        [fe1] = computeSideLoad(index,xe,de,'flux');
        fe = fe + fe1;
    end
end

% now radiation load
if size(radiationLoad,1) > 0
    index = find(radiationLoad(:,1)==el,1); % 1 face
    flag = size(index,1);
    % compute side load
    if flag > 0
        [fe2] = computeSideLoad(index,xe,de,'radiation');
        fe = fe + fe2;
    end
end

%fe = fe - ke * de;
fe = fe + me*ae;
%
%me = me + alpha*dt*ke;

end
