function obj = parabol_template(R,D,opticalInteractionType,phi,center)
%R parameter vector Ax^2 + By^2 - 2cz = 0
%D diameter of apterture
%optical interaction type: cell array with opticalInteractionTypes
%center position of center of spacing of width w
%phi orientation of objects

geometry.center = [0;0;0];
geometry.D = D;
geometry.phi = [0;0;0]; %degree
geometry.M = roti(geometry.phi(1),'x4').*roti(geometry.phi(2),'y4').*roti(geometry.phi(3),'z4');
geometry.A = zeros(4,4);
geometry.A(1,1) = R(1); %x^2
geometry.A(2,2) = R(2); %y^2
geometry.A(3,4) = R(3)/2; %z
geometry.A(4,3) = R(3)/2; %z
geometry.A = geometry.M*geometry.A;
geometry.type = 'quadsurf';
geometry.n = [0;0;1];

geometry.n = geometry.n/norm(geometry.n);

efficiency_mode = 0;
efficiency_R = 1;
efficiency_T = 1;
isOpticalActive = 1;

interaction(1).n = 1;
interaction(2).n = 1.5;
interaction(1).refractiveIndexType = 'constant';
interaction(2).refractiveIndexType = 'constant';

optElements(1) = OpticalSurf(geometry,opticalInteractionType{1},interaction,efficiency_mode,efficiency_R,efficiency_T,isOpticalActive);

bbox = [];
name = sprintf('Parabol f=%1.3f m',-R(3)/2);

obj = OpticalObject(bbox, optElements, name, phi,center);
