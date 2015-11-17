function OO = mirror_template(w,d,h,n,opticalInteractionType,phi,center)
%returns 6 plane object with n2 inside and n2 outside
%width, height, thickness, n2,n1,refl/refr, rotation, center
%plane rect: center, n, width,height,thickness

switch length(n)
    case 1
        interaction_front(1).n = 1;
        interaction_front(2).n = n;
        interaction_front(1).refractiveIndexType = 'constant';
        interaction_front(2).refractiveIndexType = 'constant';
        interaction_back = interaction_front;
    case 2
        interaction_front(1).n = n(2);
        interaction_front(2).n = n(1);
        interaction_front(1).refractiveIndexType = 'constant';
        interaction_front(2).refractiveIndexType = 'constant';        
        interaction_back = interaction_front;
    case 3
        interaction_front(1).n = n(2);
        interaction_front(2).n = n(1);
        interaction_front(1).refractiveIndexType = 'constant';
        interaction_front(2).refractiveIndexType = 'constant';    
        
        interaction_back(1).n = n(3);
        interaction_back(2).n = n(1);
        interaction_back(1).refractiveIndexType = 'constant';
        interaction_back(2).refractiveIndexType = 'constant';    
end

efficiency_mode = 0;
efficiency_R = 1;
efficiency_T = 1;
isOpticalActive = 1;

switch length(opticalInteractionType)
    case 1
        oit_front = opticalInteractionType{1};
        oit_back = opticalInteractionType{1};
        oit_side = opticalInteractionType{1};
    case 2
        oit_front = opticalInteractionType{1};
        oit_back = opticalInteractionType{2};
        oit_side = opticalInteractionType{2};
    case 3
        oit_front = opticalInteractionType{1};
        oit_back = opticalInteractionType{2};
        oit_side = opticalInteractionType{3};
end

%front
opticalInteractionType = oit_front;
geometry.width = w;
geometry.height = [0;0;h];
geometry.center = -[d/2;0;0];
geometry.type = 'planeRect';
geometry.n = -[1;0;0]; 
optElements(1) = OpticalSurf(geometry,opticalInteractionType,interaction_front,efficiency_mode,efficiency_R,efficiency_T,isOpticalActive);

%back
opticalInteractionType = oit_back;
geometry.center = [+d/2;0;0];
geometry.n = [1;0;0]; 
optElements(2) = OpticalSurf(geometry,opticalInteractionType,interaction_back,efficiency_mode,efficiency_R,efficiency_T,isOpticalActive);

%sides

%left
opticalInteractionType = oit_side;
geometry.center = [0;w/2;0];
geometry.n = [0;1;0]; 
geometry.width = d;
geometry.height = [0;0;h];
optElements(3) = OpticalSurf(geometry,opticalInteractionType,interaction_front,efficiency_mode,efficiency_R,efficiency_T,isOpticalActive);

%right
geometry.center = -[0;w/2;0];
geometry.n = -[0;1;0]; 
geometry.height = [0;0;h];
optElements(4) = OpticalSurf(geometry,opticalInteractionType,interaction_front,efficiency_mode,efficiency_R,efficiency_T,isOpticalActive);

%up
geometry.center = [0;0;h/2];
geometry.n = [0;0;1];
geometry.width = w;
geometry.height = [d;0;0];
optElements(5) = OpticalSurf(geometry,opticalInteractionType,interaction_front,efficiency_mode,efficiency_R,efficiency_T,isOpticalActive);

%down
geometry.center = -[0;0;h/2];
geometry.n = -[0;0;-1]; 
optElements(6) = OpticalSurf(geometry,opticalInteractionType,interaction_front,efficiency_mode,efficiency_R,efficiency_T,isOpticalActive);

bbox = [];
name = 'Mirror';
OO = OpticalObject(bbox, optElements, name,phi,center);
