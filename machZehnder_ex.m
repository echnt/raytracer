function machZehnder_ex()
% Mach Zehnder Interferometer

addpath('./template/','./classes/');

%% beam
v = [1;0;0];
v = v/norm(v);
r = [-0.2;0;0];
beamDiameter = 18e-3;
N = [100,2];
E = 1;
pol = 0;
l0 = 800e-9;
dl = 130e-9;
Nl = 1;
wavelength = linspace(l0-dl,l0+dl,Nl);
type = 'circular';
intensity_function = 'gauss';
beam = beam_template(v,r,beamDiameter,N,E,pol,wavelength,type,intensity_function);

%% mirrors
w = 150e-3;
h = 50e-3;
d = 20e-3;
n_in = 1.5;
n_out = 1;
% opticalInteractionType = {'reflrefr'};
opticalInteractionType = {'reflrefr','refraction'};
phi = [0;0;45];
center = [0;0;0]*1e-3;
OO(1) = mirror_template(w,d,h,[n_in,n_out],opticalInteractionType,phi,center);
opticalInteractionType = {'refraction','reflrefr'};
center = [0;-0.5;0];
phi = [0;0;45];
OO(2) = mirror_template(w,d,h,[n_in,n_out],opticalInteractionType,phi,center);
opticalInteractionType = {'reflrefr','refraction'};
center = [0.5;0;0];
phi = [0;0;45];
OO(3) = mirror_template(w,d,h,[n_in,n_out],opticalInteractionType,phi,center);
opticalInteractionType = {'refraction','reflrefr'};
center = [0.5;-0.5;0];
phi = [0;0;45];
OO(4) = mirror_template(w,d,h,[n_in,n_out],opticalInteractionType,phi,center);

%% delay glass

d = 1e-3;
opticalInteractionType = {'refraction','refraction'};
center = [0.2;0;0];
phi = [0;0;0];
OO(5) = mirror_template(w,d,h,[n_in,n_out],opticalInteractionType,phi,center);

%% detection
geometry.radius = 50e-3;
geometry.center = [0.6;-0.48;0];
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
OO(length(OO)+1) = OpticalObject(bbox, optElements, name);

%% calc and display
maxInteractionCount = 100;
beam = raytrace(OO,beam,maxInteractionCount);

figure(1);
clf;
plotRaytracing(OO,beam);
axis equal;
view(0,90);

figure(2);
S = OO(length(OO)).optElements(1).ID == beam.ID & beam.t > 0;
plot(beam.dist(S)-min(beam.dist(S)),'.');
xlabel('beam index');
ylabel('optical delay by glass /m');
title(sprintf('Path Difference with a %3.1fmm glass',d*1e3));
grid on;
drawnow;
