function [ normF ] = computeF
global neq nel coordinates elements LM U Phase PhaseOld Uold

F = zeros(neq,1);
for i=1:nel
    xe = coordinates(elements(i,2:4),:);
    de = U(1,elements(i,2:4));
    ae = U(2,elements(i,2:4));
    deOld = Uold(:,elements(i,2:4));
    % get phases on elements nodes
    pe = Phase(:,elements(i,2:4));
    peOld = PhaseOld(:,elements(i,2:4));
    matNum = elements(i,1); % element material number
    [fe,~,~] = weakform(i,xe,de',deOld,ae',pe,peOld,matNum);
    for k=1:3
        i_index = LM(k,i);
        if (i_index > 0)
            F(i_index) = F(i_index) + fe(k);
        end
    end
end

normF = norm(F);

end

