function OO = prism_template(L,w,n,opticalInteractionType,phi,center)
%L leg length [L1,L2] or L12, 
%w width perpendicular
%n = [n2(inside),n1(outside)],
%rotation, center

switch length(L) 
    case 1 %
        L1 = L;
        L2 = L;
        L3 = L;
        h = sqrt(3)/2*L;
    case 2
        L1 = L(1);
        L2 = L(1);
        L3 = L(2);
        h = sqrt(L1^2 - (L2/2)^2);
    case 3
        L1 = L(1);
        L2 = L(2);
        L3 = L(3);
        s = sum(L)/2;
        F = sqrt(s*(s-L1)*(s-L2)*(s-L3));
        h = 2*F/L3;
end

p = sqrt(L1^2 - h^2);
q = sqrt(L2^2 - h^2);
alpha1 = asind(h/L1);
alpha2 = asind(h/L2);
if(iscell(n))
    interaction(1).materialType = n{1};
    interaction(2).materialType = n{2};
    interaction(1).refractiveIndexType = 'sellmeier';
    interaction(2).refractiveIndexType = 'sellmeier';
else
    switch length(n)
        case 1
            interaction(1).n = 1;
            interaction(2).n = n;
            interaction(1).refractiveIndexType = 'constant';
            interaction(2).refractiveIndexType = 'constant';
        case 2
            interaction(2).n = n(1);
            interaction(1).n = n(2);
            interaction(1).refractiveIndexType = 'constant';
            interaction(2).refractiveIndexType = 'constant';
    end
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
M = roti(0,'x')*roti(90-alpha1,'y')*roti(0,'z');
geometry.width = w;
geometry.height = M*[0;0;L1];
geometry.center = [-q/2;0;h/2];
geometry.type = 'planeRect';
geometry.n = M*[-1;0;0];
optElements(1) = OpticalSurf(geometry,opticalInteractionType,interaction,efficiency_mode,efficiency_R,efficiency_T,isOpticalActive);


% %back
opticalInteractionType = oit_back;
M = roti(0,'x')*roti(-(90-alpha2),'y')*roti(0,'z');
geometry.width = w;
geometry.height = M*[0;0;L2];
geometry.center = [p/2;0;h/2];
geometry.type = 'planeRect';
geometry.n = M*[1;0;0];
optElements(2) = OpticalSurf(geometry,opticalInteractionType,interaction,efficiency_mode,efficiency_R,efficiency_T,isOpticalActive);

% %bottom
opticalInteractionType = oit_side;
geometry.width = w;
geometry.height = [L3;0;0];
geometry.center = [0;0;0];
geometry.type = 'planeRect';
geometry.n = [0;0;-1];
optElements(3) = OpticalSurf(geometry,opticalInteractionType,interaction,efficiency_mode,efficiency_R,efficiency_T,isOpticalActive);

%sides
%not implemented

bbox = [];
name = 'Mirror';
OO = OpticalObject(bbox, optElements, name,phi,center);
