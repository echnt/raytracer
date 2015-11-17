function beam = raytrace(OO,beam,maxInteractionCount)
    beam = propagateBeam(OO,beam,maxInteractionCount);
end

function beamProp = propagateBeam(OO,beam,maxInteractionCount)
    %propagtes (vectorial) Beamclass trough optical setup of OO
    %opticalObject until maxInteractionCount is reached. Returns a Beamclass.
    
    beamN = beam;
    for k=1:maxInteractionCount %hard limit number of calculations to prevent physical endless loops
        I = [];
        S = [];
        for j=1:length(OO) %compute distances to all objects         
            I = [I,ones(1,length(OO(j).optElements))*j]; %index lookup
            S = [S,1:length(OO(j).optElements)];
            if(j == 1)
                intersection = intersectOpticalObject(OO(j),beam);
            else
                intersection = [intersection,intersectOpticalObject(OO(j),beam)];
            end
        end
        intersection_OO = getMinTfromArray(intersection);
        
        sel = intersection_OO.isIntersect;
        beamN.t(sel) = intersection_OO.t(sel);
        beamN.ID(sel) = intersection_OO.ID(sel);
        if(k == 1)
            beamProp = beamN;
        else
            beamProp = mergeBeams(beamProp, beamN);
        end
        
        L = unique(intersection_OO.index); %group handling to single opt. surfaces
        
        beamN = Beamclass();
        for j=1:length(L)            
            sel = (L(j) == intersection_OO.index & intersection_OO.isIntersect);
            if(sum(sel) > 0)
                X = unique(I(intersection_OO.index(sel)));
                Y = unique(S(intersection_OO.index(sel)));
                beamN(j) = OpticalSurf.interactOpticalSurf(...
                    OO(X).optElements(Y),...
                    selectBeamByIndex(beam, sel),...
                    selectISbyIndex(intersection_OO,sel));
            else
                beamN(j) = Beamclass();
            end
        end

        beamN = mergeBeams(beamN);      
        beam = beamN;
        if(isempty(beamN.r))
            fprintf('Nothing to compute, stop raytrace after %d iterations\n',k);
            break; %stop if there is no more beam            
        end
    end
end

