function [ normF ] = computeF
global neq nel coordinates elements LM U Phase PhaseOld Uold
global convectionLoad fluxLoad radiationLoad

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
% add convection load
conv_elem = size(convectionLoad,1);
if conv_elem > 0
    for i=1:conv_elem
        el = convectionLoad(i,1);
        xe = coordinates(elements(el,2:4),:);
        deOld = Uold(:,elements(el,2:4));
        [fe] = computeSideLoad(i,xe,deOld(1,:)','convection');
        for k=1:3
            i_index = LM(k,el);
            if (i_index > 0)
                F(i_index) = F(i_index) + fe(k);
            end
        end
    end
end
% now flux load
flux_elem = size(fluxLoad,1);
if flux_elem > 0
    for i=1:flux_elem
        el = fluxLoad(i,1);
        xe = coordinates(elements(el,2:4),:);
        deOld = Uold(:,elements(el,2:4));
        [fe] = computeSideLoad(i,xe,deOld(1,:)','flux');
        for k=1:3
            i_index = LM(k,el);
            if (i_index > 0)
                F(i_index) = F(i_index) + fe(k);
            end
        end
    end
end
% add radiation load
rad_elem = size(radiationLoad,1);
if rad_elem > 0
    for i=1:rad_elem
        el = radiationLoad(i,1);
        xe = coordinates(elements(el,2:4),:);
        deOld = Uold(:,elements(el,2:4));
        [fe] = computeSideLoad(i,xe,deOld(1,:)','radiation');
        for k=1:3
            i_index = LM(k,el);
            if (i_index > 0)
                F(i_index) = F(i_index) + fe(k);
            end
        end
    end
        end

normF = norm(F);

end

