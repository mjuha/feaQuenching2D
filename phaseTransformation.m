function phaseTransformation(dt)
global MAT Phase scheilRule U nn

% get properties, only one material
% tauS
tauS = MAT{9};
%tauF
tauF = MAT{10};
%
Fs = 0.01;
Ff = 0.99;
% Austenite, Perlite, Bainite and Martensite
%vF = [1.0, 0, 0, 0];
% temperature data
tempData = MAT{8};
% Ta = tempData(1); % starting temperature for transformation of austenite
% Tb = tempData(2); % starting temperature for transformation of perlite
Ms = tempData(3); % starting temperature for transformation of martensite

for i = 1:nn % all nodes in domain
    temp = U(1,i); % temperature node i
    vF = Phase(:,i);
    sRe = scheilRule(i);
    tau_s = interp1(tauS(:,2),tauS(:,1),temp,'linear',Inf);
    tau_f = interp1(tauF(:,2),tauF(:,1),temp,'linear',Inf);
    if tau_s ~= Inf
        sRe = sRe + dt/tau_s;
    end
    % phase transformation started
    if ( sRe >= 1 ) || (temp <= Ms)
        % perlite
        if ( tau_f ~= Inf ) && ( tau_s ~= Inf )
            n = ( log(log(Fs)) - log(log(Ff)) ) / ( log(tau_f) - log(tau_s) );
            %
            a = -log(Ff) * tau_s^(-n);
            % compute ficticious time
            tfict = ( -log( 1 - vF(2) ) / a )^(1/n);
            % compute actual time
            tj = dt + tfict;
            %
            % compute volume phase (perlite)
            vF(2) = ( 1 - exp( -a * tj^n ) );
            vF(1) = 0.99 - vF(2); % (austenite)
        else % martensite
            if temp <= Ms
                old = (0.99 - vF(2));
                Fm = ( 1 - exp( -0.011 * (Ms - temp) ) ) * old;
                if Fm < vF(3) % it is an irreversible process
                    Fm = vF(3);
                end
                % autenite
                tmp = old - Fm;
                if tmp >= 0
                    vF(1) = tmp;
                    vF(3) = Fm;
                else
                    vF(1) = 0;
                    vF(3) = old;
                end
            end
        end
    end
    % update global variables
    scheilRule(i) = sRe;
    Phase(:,i) = vF;
end
end