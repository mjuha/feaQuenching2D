function [H] = LatentHeat(T)

Hp = 1.56e9 - 1.5e6 * T;

% Perlite , Martensite, Bainite
H = [Hp , 640e6, Hp];

end

