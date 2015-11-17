function machZehnder_ex()
% Mach Zehnder Interferometer

addpath('./template/','./classes/');

%% beam
v = [1;0;0];
v = v/norm(v);
r = [-0.2;0;0];
beamDiameter = 18e-3;
Nr = 100;
Np = 2;
N = [Nr,Np];
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
phi = [0;1;45];
OO(4) = mirror_template(w,d,h,[n_in,n_out],opticalInteractionType,phi,center);

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
subplot(2,1,1);
S = OO(length(OO)).optElements(1).ID == beam.ID & beam.t > 0;
D = beam.dist(S) + MV3norm(beam.r(:,S) + beam.v(:,S).*repmat(beam.t(S),3,1));
r1 = D(1:Nr*Np);
r2 = D(Nr*Np+1:2*Nr*Np);
y = beam.r(2,S);
y = y(1:Nr*Np);
y = y-mean(y);
wl = beam.wavelength(S);
wl = wl(1:Nr*Np);
E0 = beam.E0(S);
E0 = E0(1:Nr*Np);

plot(y*1e3,(r1-r2)*1e6,'.');
xlabel('r(beam line out) /mm');
ylabel('optical path difference /mu m');
title(sprintf('Path Difference with misaligned beam splitter'));
grid on;
drawnow;

subplot(2,1,2);
E = E0.*exp(1i * 2*pi./wl .* (r1-r2));
[y,I] = sort(y);
E = E(I);
E0 = E0(I);

plot(y*1e3,(E+E0).*conj(E+E0),'.-',y*1e3,E0.*conj(E0));
xlabel('r(beam line out) /mm');
ylabel('Intensity');
grid on;
legend('Interference','Initial beam');
