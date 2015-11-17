function obj = lens_template(nm,R,D,w,opticalInteractionType,phi,center)
%nm refractive index of medium n1|n2|n3
%R radius of curvature R1,R2, can be inf for plano convex
%D diameter of lense apterture
%w width between apex(concave)/plane(convex) of left side and right side
%optical interaction type: cell array with opticalInteractionTypes
%center position of center of spacing of width w
%phi orientation of objects

%% first/"left" lense
if(length(w) == 1)
    w(1:2) = w/2;
end

switch length(nm)
    case 1
        nm = [1;nm;1];
    case 2
        nm = [nm(1);nm(2);nm(1)];
end

if(size(R,1) == 3)
    geometry.radius = R(:,1);
else
    geometry.radius = R(1)*[1;1;1];
end

if(length(opticalInteractionType) == 1)
    opticalInteractionType{1:2} = opticalInteractionType{1};
end

if(max(isinf(geometry.radius)))
    geometry.radius = D/2;
    geometry.center = [-w(1)/2;0;0];
    geometry.type = 'planeDisc';
    geometry.n = [-1;0;0];
else
    if(max(geometry.radius < 0))
        geometry.center = [-sqrt(-(D/2)^2 + (geometry.radius(1))^2 ) - w(1)/2;0;0];
    else
        geometry.center = [+sqrt(-(D/2)^2 + (geometry.radius(1))^2 ) - w(1)/2;0;0];
    end
    geometry.D = D;
    geometry.phi = [0;0;0]; %degree
    geometry.M = roti(geometry.phi(1),'x4').*roti(geometry.phi(2),'y4').*roti(geometry.phi(3),'z4');
    geometry.A = zeros(4,4);
    geometry.A(sub2ind([4,4],1:3,1:3)) = 1./geometry.radius.^2;
    geometry.A(4,4) = -1;
    geometry.A = geometry.M*geometry.A;
    geometry.type = 'quadsurf';
    geometry.n = [-1;0;0];
end

geometry.n = geometry.n/norm(geometry.n);

efficiency_mode = 0;
efficiency_R = 1;
efficiency_T = 1;
isOpticalActive = 1;
    
interaction.n_normal = nm(1);
interaction.n_back = nm(2);

interaction(1).n = nm(1);
interaction(2).n = nm(2);
interaction(1).refractiveIndexType = 'constant';
interaction(2).refractiveIndexType = 'constant';

optElements(1) = OpticalSurf(geometry,opticalInteractionType{1},interaction,efficiency_mode,efficiency_R,efficiency_T,isOpticalActive);

%% second/"right" lense
if(size(R,1) == 3)
    geometry.radius = R(:,2);
else
    geometry.radius = R(2)*[1;1;1];
end

if(max(isinf(geometry.radius)))
    geometry.radius = D/2;
    geometry.center = [w(2)/2;0;0];
    geometry.type = 'planeDisc';
    geometry.n = [1;0;0];
else
    if(max(geometry.radius < 0))
        geometry.center = [- sqrt(-(D/2)^2 + (geometry.radius(1))^2 ) + w(2)/2;0;0];
    else 
        geometry.center = [- sqrt(-(D/2)^2 + (geometry.radius(1))^2 ) + w(2)/2;0;0];
    end
    geometry.D = D;
    geometry.phi = [0,0,0]; %degree
    geometry.M = roti(geometry.phi(1),'x4').*roti(geometry.phi(2),'y4').*roti(geometry.phi(3),'z4');
    geometry.A = zeros(4,4);
    geometry.A(sub2ind([4,4],1:3,1:3)) = 1./geometry.radius.^2;
    geometry.A(4,4) = -1;
    geometry.A = geometry.M*geometry.A;
    geometry.type = 'quadsurf';
    geometry.n = [1;0;0];
end

geometry.n = geometry.n/norm(geometry.n);

interaction(1).n = nm(3);
interaction(2).n = nm(2);
interaction(1).refractiveIndexType = 'constant';
interaction(2).refractiveIndexType = 'constant';

efficiency_mode = 0;
efficiency_R = 1;
efficiency_T = 1;
isOpticalActive = 1;

optElements(2) = OpticalSurf(geometry,opticalInteractionType{2},interaction,efficiency_mode,efficiency_R,efficiency_T,isOpticalActive);

bbox = [];
name = sprintf('Lense f=%3.4f m',1 / ((nm(2)-nm(1)) * (1/R(1) - -1/R(2) + (nm(2)-nm(1))*sum(w)/nm(2)/R(1)/-R(2) )) );
obj = OpticalObject(bbox, optElements, name, phi,center);
