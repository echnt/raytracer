function prism_ex()
%This example shows the possibility to spectrally seperate a beam with the
%help of material dispersion in glass. The rotation of the prism is scanned
%to find the minimal deviation angle. Fullfilling the minimum deviation
%condition allows to easily estimate the refrative index if the wavelength
%is known. The measurement is performed for the central wavelength and
%deviation from the ideal measurement is displayed
%for an experimental description see e.g.
%http://www.physik.uni-jena.de/Versuch_406.html

%% beam
v = [1;0;0];
v = v/norm(v);
r = [-40;0;0]*1e-3;
beamDiameter = 3e-3;
N = [1,1];
E = 1;
pol = 2;
l0 = 580e-9;
dl = 180e-9;
Nl = 50;
wavelength = linspace(l0-dl,l0+dl,Nl);
type = 'circular';
intensity_function = 'gauss';
beam = beam_template(v,r,beamDiameter,N,E,pol,wavelength,type,intensity_function);

%% detection
geometry.radius = 40e-3;
geometry.center = [20;0;-40]*1e-3;
geometry.type = 'planeDisc';
geometry.n = [1;0;0]; 
opticalInteractionType = 'refraction';

interaction(1).n = 1;
interaction(2).n = 1;
interaction(1).refractiveIndexType = 'constant';
interaction(2).refractiveIndexType = 'constant';

efficiency_mode = 0;
efficiency_R = 0;
efficiency_T = 0;
isOpticalActive = 1;

optElements(1) = OpticalSurf(geometry,opticalInteractionType,interaction,efficiency_mode,efficiency_R,efficiency_T,isOpticalActive);

bbox = [];
name = 'Imaging Plane';
OO(1) = OpticalObject(bbox, optElements, name);

%% optics
L = [20,20,20]*1e-3;
delta = 60;
w = 10e-3;
nm{1} = 'Vacuum';
nm{2} = 'SF10';
center = [0;0;-L(1)/3-4e-3];
opticalInteractionType = {'reflrefr','reflrefr'};

phi = linspace(24,50,80);
theta = zeros(length(phi),length(wavelength));
% I_t0_col = 0.4;
% t_t0_len = 10e-3;
% t_tmax_len = 1;

for k=1:length(phi)
    OO(2) = prism_template(L,w,nm,opticalInteractionType,[0;phi(k);0],center);
    maxInteractionCount = 3;
    beamN = raytrace(OO,beam,maxInteractionCount);
    
%     figure(1);
%     clf;
%     plotRaytracing(OO,beamN,I_t0_col,t_t0_len,t_tmax_len);
%     axis equal;
%     view(0,0);
%     drawnow;
    
    S = OO(1).optElements(1).ID == beamN.ID & beamN.t > 0;
    beamN = selectBeamByIndex(beamN,S);
    devAngle = acosd(dot(beamN.v,repmat(v,1,size(beamN.v,2))));
    theta(k,1:length(devAngle)) = devAngle;
end

%% evaluate aquired data
wavelength =  beamN.wavelength;

figure(2);
[x,y] = meshgrid(phi,wavelength*1e9);
surf(x,y,theta','edgecolor','none');
xlabel('prisma rotation angle /deg');
ylabel('wavelength /nm');
zlabel('deviation angle \delta /deg');
view(0,0);

figure(3);
mt = min(theta);
plot(wavelength*1e9,mt);
xlabel('wavelength /nm');
ylabel('minimum deviation angle \delta /deg');
grid on;

figure(4);
lctr = 532e-9;
[~,S] = min(abs(wavelength - lctr));
mt = min(theta(:,S),[],2);
plot(phi,mt);
xlabel('prisma rotation angle /deg');
ylabel(sprintf('deviation angle at %3.0f nm/deg',lctr*1e9));
grid on;

%inverse calculation
figure(5);
interaction(1).materialType = nm{1};
interaction(2).materialType = nm{2};
interaction(1).refractiveIndexType = 'sellmeier';
interaction(2).refractiveIndexType = 'sellmeier';
ni = OpticalSurf.calcRefractiveIndex(interaction(2),wavelength);

[~,S] = min(theta(:,S)); %select lctr minimum
devAngle = theta(S,:);
n = sind( (delta + devAngle)/2 )./sind(delta/2);
plot(wavelength*1e9,n,'x',wavelength*1e9,ni,'--');
xlabel('wavelength /nm');
ylabel(sprintf('refractive index (prisma %3.0f nm position)',lctr*1e9));
grid on;
legend('measurement n=sin( [\delta + \epsilon]/2) / sin(\epsilon/2)',sprintf('sellmeier input %s to %s',nm{1},nm{2}));
