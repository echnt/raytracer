function fresnel_ex
%Display fresnel equation T and R coefficients as function of angle

addpath('./template/','./classes/');

%% beam

v = [-1;0;0];
v = v/norm(v);
r = [-0.1;0;0];
beamDiameter = 248e-3;
N = [100,2];
E = 1;
pol = 2;
l0 = 800e-9;
dl = 130e-9;
Nl = 1;
wavelength = linspace(l0-dl,l0+dl,Nl);
type = 'circular';
intensity_function = 'flat';
beam = beam_template(v,r,beamDiameter,N,E,pol,wavelength,type,intensity_function);

%% optics

%fresnel type reflection
w = 150e-3;
h = 50e-3;
d = 20e-3;
n_in = 1.5;
n_out = 1;
opticalInteractionType = {'reflrefr','refraction'};
phi = [0;0;45];
center = [1;1;0]*d/2/sqrt(2);
OO(1) = mirror_template(w,d,h,[n_in,n_out],opticalInteractionType,phi,center);
OO(1).optElements(1).efficiency_mode = 1;
OO(1).optElements(2).efficiency_mode = 1;

%init optics for suitable beam
D = 1.5*beamDiameter;
opticalInteractionType = {'reflection'};
center = r+[-0.05;0;0];
f = -center(1);
R = [1,1,-4*f]; %x^2+y^2 = 2*(2f)*z

phi = [0,90,0];
OO(2) = parabol_template(R,D,opticalInteractionType,phi,center);

%% detection
geometry.radius = 50e-3;
geometry.center = [0.05;0;0];
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
OO(3) = OpticalObject(bbox, optElements, name);

geometry.n = [0;1;0]; 
geometry.center = [0;-0.05;0];
optElements(1) = OpticalSurf(geometry,opticalInteractionType,interaction,efficiency_mode,efficiency_R,efficiency_T,isOpticalActive);
OO(4) = OpticalObject(bbox, optElements, name);

%% calculate beam path
maxInteractionCount = 5;

%P,2
beamP = raytrace(OO,beam,maxInteractionCount);
%S,1
beam.polarization = beam.polarization - 1;
beamS = raytrace(OO,beam,maxInteractionCount);


%% Plot beam path and objects
I_t0_col = 0;
t_t0_len = 0;
t_tmax_len = 1;
plotRaytracing(OO,beamS,I_t0_col,t_t0_len,t_tmax_len);

axis equal;
axis([-0.05,0.05,-0.05,0.05,-0.05,0.05]);
view(0,90);

%% extract data and plot fresnel R,T

Spt = OO(3).optElements(1).ID == beamP.ID & beamP.t > 0 & beamP.polarization == 2;
Spr = OO(4).optElements(1).ID == beamP.ID & beamP.t > 0 & beamP.polarization == 2;
Sst = OO(3).optElements(1).ID == beamS.ID & beamS.t > 0 & beamS.polarization == 1;
Ssr = OO(4).optElements(1).ID == beamS.ID & beamS.t > 0 & beamS.polarization == 1;
Spp = OO(1).optElements(1).ID == beamP.ID;
Sps = OO(1).optElements(1).ID == beamS.ID;

%get angle of incidence beam and output beam
n = repmat(beamP.v(:,N(1)*N(2)+1),1,size(beamP.v(:,Spp),2));
v = beamP.v(:,Spp);
phi = acosd(-dot(v,n));
phiP = phi - mean(phi) + 45;

n = repmat(beamS.v(:,N(1)*N(2)+1),1,size(beamS.v(:,Sps),2));
v = beamS.v(:,Sps);
phi = acosd(-dot(v,n));
phiS = phi - mean(phi) + 45;

figure(2);
clf;
plot(phiS,abs(beamS.E0(Sst)).^2,'.',phiS,abs(beamS.E0(Ssr)).^2,'.',phiP,abs(beamP.E0(Spt)).^2,'.',phiP,abs(beamP.E0(Spr)).^2,'.');
xlabel('angle of first incidence /deg');
ylabel('intensity /normalized to input');
grid on;
legend('T_s','R_s','T_p','R_p','Location','West');
