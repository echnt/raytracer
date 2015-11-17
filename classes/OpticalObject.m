classdef OpticalObject
    %Class to store array of optElements to group Elements like mirrors,
    %lenses, etc.
    properties
        bbox; %box with limited surfaces to determine if intersect has to run
        optElements; %array of optical surfs
        center; %center which is added to all optical surfs
        phi; %rotation vector for all optical surfs to (0,0,0) center
        name; %descriptive name
    end
    methods
        function obj = OpticalObject(bbox, optElements, name,phi,center)
            if(nargin > 0)
                obj.bbox = bbox;
                obj.optElements = optElements;
                obj.name = name;
                obj.phi = [0;0;0];
                obj.center = [0;0;0];
                switch nargin
                    case 3
                        phi = [0;0;0];
                        center = [0;0;0];
                end
                obj = OpticalObject.updateRotateShiftOpticalObject(obj,phi,center);
            else
                obj.bbox = [];
                obj.optElements = [];
                obj.name = [];
                obj.phi = [];
                obj.center = [];
            end
            
        end
        
        function intersection = intersectOpticalObject(obj,beam)
            %returns list of beam intersects with the first element hit in
            %Optical Object
            
            %at this point it should check if bbox is hit first...                        
            intersection(length(obj.optElements)) = ISclass();
            for k=1:length(obj.optElements)
                intersection(k) = ISclass.intersectOpticalSurf(obj.optElements(k),beam);
                intersection(k).ID = obj.optElements(k).ID*ones(1,size(beam.v,2));
            end
        end
        
    end
    
    methods(Static)   
        
        function beam = interactOpticalObject(obj,beam,intersection)
            %calculate new beam at intersection point
            beam(length(obj.optElements)) = Beamclass();
            L = unique(intersection.index);
            for k=1:length(L)
                S = L(k) == intersection.index;
                beam(k) = OpticalSurf.interactOpticalSurf(obj.optElements(k),...
                    selectBeamByIndex(beam,S),...
                    selectISbyIndex(intersection,S));
            end
            beam = mergeBeams(beam);
        end
        
        function obj = updateRotateShiftOpticalObject(obj,phi,center)
            %transform array of optical objects by phi, center replacing the previous used transformation
            Mold = roti(obj.phi(1),'x')*roti(obj.phi(2),'y')*roti(obj.phi(3),'z');
            Mnew = roti(phi(1),'x')*roti(phi(2),'y')*roti(phi(3),'z');
            for k=1:length(obj.optElements)
                obj.optElements(k).geometry.n = Mnew * (Mold\obj.optElements(k).geometry.n);
                obj.optElements(k).geometry.center = Mnew * ( Mold\(obj.optElements(k).geometry.center - obj.center) ) + center;
                if(strcmp(obj.optElements(k).interactionType,OpticalSurf.interactionTypes{4})) %diffraction
                    obj.optElements(k).interaction.g = Mnew * ( Mold\obj.optElements(k).interaction.g );
                end
                if(strcmp(obj.optElements(k).geometry.type,OpticalSurf.geometryTypes{2})) %rect
                    obj.optElements(k).geometry.height = Mnew * ( Mold\obj.optElements(k).geometry.height);
                end
                if(strcmp(obj.optElements(k).geometry.type,OpticalSurf.geometryTypes{3})) %quad surf
                    %transform xi'*A*xi + 2*b*x' + J(=1) = 0; with M'*A*M
                    %and M*b
                    obj.optElements(k).geometry.M = Mnew;
                    A0 = obj.optElements(k).geometry.A;
                    if(size(obj.optElements(k).geometry.A,1) == 3 && size(obj.optElements(k).geometry.A,2) == 3)                       
                        A0(4,4) = 1;
                    end
                    %extract A and b
                    b = Mnew * (Mold \ A0(1:3,4)); %bi
                    A = A0(1:3,1:3);
                    A = Mnew'*(Mold'\A/Mold)*Mnew; %ai
                    A(4,4) = A0(4,4);
                    A(4,1:3) = b;
                    A(1:3,4) = b;
                    obj.optElements(k).geometry.A = A;
                end                
            end
            obj.phi = phi;
            obj.center = center;
        end
    end
end
