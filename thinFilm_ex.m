function thinFilm_ex
%prints an example picture for "thin"film interference

%% beam

v = [1;1;0];
v = v/norm(v);
r = [-0.1;-0.1;0]/2;
beamDiameter = 3e-3;
N = [100,2];
E = 1;
pol = 2;
l0 = 800e-9;
dl = 130e-9;
Nl = 1;
wavelength = linspace(l0-dl,l0+dl,Nl);
type = 'circular';
intensity_function = 'gauss';
beam = beam_template(v,r,beamDiameter,N,E,pol,wavelength,type,intensity_function);

%% optics

%fresnel type reflection
w = 150e-3;
h = 50e-3;
d = 20e-3;
n_back = 1.8;
n_in = 1.5;
n_out = 1;
opticalInteractionType = {'reflrefr','reflrefr'};
phi = [0;0;0];
center = [1;3;0]*d/2;
OO(1) = mirror_template(w,d,h,[n_in,n_out,n_back],opticalInteractionType,phi,center);
OO(1).optElements(1).efficiency_mode = 1;
OO(1).optElements(2).efficiency_mode = 1;

%% detection
geometry.radius = 60e-3;
geometry.center = [40;70;0]*1e-3;
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
geometry.center = [-40;70;0]*1e-3;
optElements(2) = OpticalSurf(geometry,opticalInteractionType,interaction,efficiency_mode,efficiency_R,efficiency_T,isOpticalActive);

bbox = [];
name = 'Imaging Plane';
OO(2) = OpticalObject(bbox, optElements, name);


maxInteractionCount = 11;
tic;
beam = raytrace(OO,beam,maxInteractionCount);
toc;

tic;
fHandle = figure(1);
I_t0_col = 0.4;
t_t0_len = 0.05;
t_tmax_len = 1;
plotRaytracing(OO,beam,I_t0_col,t_t0_len,t_tmax_len);
toc;
axis equal;
view(90,90);
axis off;

text(-0.02,0.11,'n_1');
text(+0.01,0.11,'n_2');
text(+0.04,0.11,'n_3');
text(0.01,-0.04,'d');
text(-0.02,-0.01,'\alpha');
hold on;
plot3([-0.05,0.00],[0,0],[0,0]);
axis([-40e-3,40e-3,-30e-3,200e-3,-100e-3,100e-3]);

set(gca,'LooseInset',get(gca,'TightInset'))

set(fHandle,'units','centimeters');
set(fHandle,'position',[5,5,12,12]);

set(findall(fHandle,'-property','FontSize'),'FontSize',14);
set(findall(fHandle,'-property','LineWidth'),'LineWidth',1);
set(fHandle,'Renderer','opengl');

set(fHandle, 'PaperUnits', 'centimeters', 'PaperPosition', [5,5,12,12]);
print(fHandle,'-dpng','-r300','./thinFilm.png');
