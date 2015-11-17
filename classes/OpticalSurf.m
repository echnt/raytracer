classdef OpticalSurf
    %collection of parameters describing simple optical surfaces. Helper
    %function calculate amplitude efficiencies and beam direction after
    %interaction with optical surf
    properties
        geometry; %shape information, additional geometry like grating 
                  %orientation
                  %
                  %property of all versions:
                  %type: implemented geometryTypes
                  %
                  %type specific information:
                  %
                  %plane disc: center, radius, n(normal vector of shape)
                  %
                  %plane rect: center, n, width,height(vector in height
                  %direction)
                  %
                  %ellispoid section: center, radius(3x1),
                  %phi(3x1),M(rotation),A(full solution matrix),
                  %n(direction of used aperture), D(diameter following r -> r-sqrt(
                  %-(D/2)^2 + R^2)
        
        interactionType; %type of optical interation handling: 
                                %relfection/refraction/diffraction 
                                %and relf/refr     
        interaction; %interaction parameters: 
                     %
                     %Array:
                     %(1) normal(refractive index on normal vector side),
                     %(2) back(other side)
                     %
                     %for refl/refr:
                     %refractiveIndexType: constant|sellmeier
                     %
                     %constant:
                     %n
                     %
                     %sellmeier:
                     %glassType: 'custom'|'BK7', 'UVFS' etc. to implemented glass type
                     %custom: array of B/C coefficients (no dim, mu m^2)
                     %
                     %diffraction: 
                     %g vector (grating orientation,normalized), 
                     %d grating line separation, 
                     %m diffraction orders to compute
                                
        efficiency_mode; %use fresnel equation or hard code input of R/T
        efficiency_R; %I = I*R_efficiency, intensity here because of measurement, internally handled as E because of theory!
        efficiency_T; %I = I*T_efficiency        
        isOpticalActive; %do optical interaction on surface
        ID; %identifier to track interaction
    end
    
    properties(Constant)
        geometryTypes = {'planeDisc','planeRect','quadsurf','cone','planePolygon'};
        efficiencyMode = {'iso'}; %this can be used to include non isotropic scattering mechanisms of surfaces
        interactionTypes = {'reflection','refraction','reflrefr','diffraction'}; %most common (linear) light-matter interactions
    end
    
    methods
        function obj = OpticalSurf(geometry,interactionType,interaction,efficiency_mode,efficiency_R,efficiency_T,isOpticalActive)
            obj.geometry = geometry;
            obj.interactionType = interactionType;
            obj.interaction = interaction;
            obj.efficiency_mode = efficiency_mode;
            obj.efficiency_R = efficiency_R;
            obj.efficiency_T = efficiency_T;
            obj.isOpticalActive = isOpticalActive;
            obj.ID = rand(1); %char(java.util.UUID.randomUUID);
        end
    end
    methods(Static)
        
        function [n_normal,n_back] = refractiveIndex(obj,beam)
            %returns front and backside refractive index
            n_normal = OpticalSurf.calcRefractiveIndex(obj.interaction(1),beam.wavelength);
            n_back = OpticalSurf.calcRefractiveIndex(obj.interaction(2),beam.wavelength);
        end
        
        function ni = calcRefractiveIndex(interaction,wl)
            %interaction: substruct of OpticalSurf
            %wl: wavelength in nm
            %returns refractive index. This is done by either returning a
            %const or calculating the t and r coefficients of the fresnel
            %equations. t and r can be negative indicating a phase jump at
            %the surface!
            switch interaction.refractiveIndexType
                case 'constant'
                    ni = repmat(interaction.n,1,length(wl));
                case 'sellmeier'
                    switch interaction.materialType
                        case 'custom'
                            B = interaction.B;
                            C = interaction.C;
                            p = interaction.p;
                        case 'BK7' %Schott BK7 glass
                            %http://refractiveindex.info/?shelf=glass&book=BK7&page=SCHOTT
                            %M. N. Polyanskiy. Refractive index database. Available at http://refractiveindex.info (accessed 11/2015)
                            B = [1.03961212     0.231792344     1.01046945];
                            C = [0.00600069867  0.0200179144    103.560653];
                            p = 0.5;
                        case 'UVFS' %Fused Silica
                            %http://refractiveindex.info/?shelf=glass&book=fused_silica&page=Malitson
                            %M. N. Polyanskiy. Refractive index database. Available at http://refractiveindex.info (accessed 11/2015)
                            B = [0.6961663  0.4079426   0.8974794];
                            C = [0.0684043  0.1162414   9.896161];
                            p = 0.5;
                        case 'SF10' %Dense Flint
                            %http://refractiveindex.info/?shelf=glass&book=SF10&page=SCHOTT
                            %M. N. Polyanskiy. Refractive index database. Available at http://refractiveindex.info (accessed 11/2015)
                            B = [1.62153902     0.256287842     1.64447552];
                            C = [0.0122241457   0.0595736775    147.468793];
                            p = 0.5;
                        case 'CaF2' %Calcium Fluoride
                            %http://refractiveindex.info/?shelf=main&book=CaF2&page=Li
                            %M. N. Polyanskiy. Refractive index database. Available at http://refractiveindex.info (accessed 11/2015)
                            B = [0.33973    0.69913     0.11994     4.35181];
                            C = [0          0.09374     21.18       38.46];
                            p = 0.5;
                        case 'MgF2' %Mangan Fluoride, ordinary
                            %http://refractiveindex.info/?shelf=main&book=MgF2&page=Li-o
                            %M. N. Polyanskiy. Refractive index database. Available at http://refractiveindex.info (accessed 11/2015)
                            B = [0.27620    0.60967     0.0080  2.14973];
                            C = [0          0.08636     18.0    25.0];
                            p = 0.5;
                        case 'Air'
                            %http://refractiveindex.info/?shelf=other&book=air&page=Ciddor
                            %M. N. Polyanskiy. Refractive index database. Available at http://refractiveindex.info (accessed 11/2015)
                            %warning: highly dependent on ambient
                            %parameters! For more complex cases use
                            %Ciddor/Edlen Formula e.g. on NIST
                            B = [0.05792105/238.0185    0.00167917/57.362];
                            C = [1/238.0185             1/57.362];                            
                            p = 1;
                        case 'Vacuum'
                            B = 0;
                            C = 0;
                            p = 1;
                    end
                    %see https://en.wikipedia.org/wiki/Sellmeier_equation
                    %use p=0.5, C=0 and B=n^2-1 for equivalent 'constant' mode
                    wl = wl*1e6; %to mu m
                    ni = (1 +  sum( (repmat(B',1,length(wl)) .* repmat(wl.^2,length(B),1))./(repmat(wl.^2,length(C),1) - repmat(C',1,length(wl))) )).^p;
            end
        end
                
        function [E_refr,E_refl] = efficiency(obj,beam,eta,c12)
            %calculates field amplitudes of beam hitting obj either by
            %0: constant value 
            %1: fresnel equation
            %this implementation should, but does not use jones calculus
            %e.g. see https://en.wikipedia.org/wiki/Jones_calculus

            E_refr = zeros(size(beam.polarization));
            E_refl = zeros(size(beam.polarization));
            
            %constant values
            C = beam.polarization == 0 | obj.efficiency_mode == 0;
            E_refr(C) = sqrt(obj.efficiency_R);
            E_refl(C) = sqrt(obj.efficiency_T);
            %fresnel equations https://en.wikipedia.org/wiki/Fresnel_equations
            %eta = n_in/n_out
            %c12 = cos(theta_in)/cos(theta_out)
            S = beam.polarization == 1 & obj.efficiency_mode == 1;
            E_refr(S) = 2 .* eta(S) .* c12(S) ./ (eta(S) .* c12(S) + 1);%vertical/S pol
            E_refl(S) = (eta(S).*c12(S) - 1) ./(eta(S) .* c12(S) + 1);%vertical/S pol
                        
            P = beam.polarization == 2 & obj.efficiency_mode == 1;
            E_refr(P) = 2 .* eta(P) .* c12(P) ./ (eta(P) + c12(P));%parallel polarization
            E_refl(P) = (-eta(P) + c12(P))./(eta(P) + c12(P));%parallel polarization
            
%             E_refl = E_refl;
            E_refr = E_refr./sqrt(eta.*c12); %correct for medium change
            
            I = (obj.efficiency_mode == 0) & true(1,length(beam.polarization));
            E_refr(I) = obj.efficiency_T;
            E_refl(I) = obj.efficiency_R;

            %critical angle handling
            I = imag(c12) ~= 0;
            E_refr(I) = 0;
            E_refl(I) = 1;
        end
        
        function vm = diffract(interaction,n,beam)
            %returns linear theory diffraction orders
            %ISclass interaction 
            %3x1 or 3xs n (s=number of beams) vector
            %Beamclass beam
            
            %e.g. see
            %https://en.wikipedia.org/wiki/Diffraction_grating for grating
            %equation
            vp = beam.v;
            
            if(size(n,2) == 1)
                nd = repmat(n,1,size(beam.v,2));
            end
            
            cos_t1 = -dot(nd,beam.v);             
            invS = ((cos_t1 > 0) - 0.5)*2 ; %angle of incidence > pi/2 -> beam arrives from backward direction
            nd = repmat(invS,3,1).*nd;
            cos_t1 = -dot(nd,beam.v); %angle of incidence with right normal vector        

            sin_t1 = sin(acos(cos_t1));
            vm = zeros(3,length(interaction.m)*size(beam.v,2));
            for k=1:length(interaction.m)
                sin_t2 = interaction.m(k).*beam.wavelength./interaction.d - sin_t1; %grating equation
                tan_t2 = sin_t2./sqrt(1-sin_t2.^2);
                % v' = v - dot(nl) * n - g * tan(t2)
                v = nd - repmat(interaction.g,1,size(vp,2)) .* repmat(tan_t2,3,1);
            
                vm(:,((k-1)*size(beam.v,2)+1) : (k*size(beam.v,2))) = v./repmat(MV3norm(v),3,1);
            end           
        end
        
        function v = reflect(n,beam)
            %returns direction of specular reflected beam on obj
            %3x1 or 3xs n (s=number of beams) vector
            %Beamclass beam
            
            %e.g. see
            %https://en.wikipedia.org/wiki/Reflection_%28mathematics%29 
            %"Reflection through a hyperplane in n dimensions"
            if(size(n,2) == 1)
                n = repmat(n,1,size(beam.v,2));
            end
            % v' = v - 2* nl n
            v = beam.v-2*repmat(dot(beam.v,n ),3,1).*n./repmat(MV3norm(n).^2,3,1);
            v = v./repmat(MV3norm(v),3,1);
        end
        
        function [v,eta,c12,n_beam] = refract(n,n_normal,n_back,beam)
            %returns direction of refraction following snells law at locally plane
            %surfaces
            %3x1 or 3xs n (s=number of beams) vector
            %double n_normal, refractive index on side of normal vector
            %double n_backside, refractive index on backside of normal
            %vector
            %Beamclass beam
            
            %e.g. see
            %https://en.wikipedia.org/wiki/Snell's_law
            
            %Vector form, not valid for sin_t2_sq > 1 - returns imag, "critical angle"    
            if(size(n,2) == 1)
                nd = repmat(n,1,size(beam.v,2));
            else
                nd = n;
            end
            cos_t1 = -dot(nd,beam.v); 
            
            inv0 = cos_t1 > 0;
            invS = (inv0 - 0.5)*2 ; %angle of incidence > pi/2 -> beam arrives from backward direction
            n_beam = n_back.*inv0 + n_normal.*~inv0; %refractive index for propagation distance calculation
            
            eta = (n_normal./n_back).^(invS);
            
            cos_t1 = invS.*cos_t1;
            nd = repmat(invS,3,1).*nd;
           
            sin_t2_sq = eta.^2 .* (1 - cos_t1.^2);
            v = repmat(eta,3,1).*beam.v + repmat(eta.*cos_t1 - sqrt(1 - sin_t2_sq),3,1 ) .* nd;
            
            c12 = cos_t1./sqrt(1-sin_t2_sq);
            v = v./repmat(MV3norm(v),3,1);
        end
        
        function beam_interact = returnBeam(intersection,beam,v,Eff,ID,t,nm)
             %returns Beamclass after input specs
             sel = intersection.isIntersect;
             if(sum(sel) > 0)
                 beam_interact = Beamclass(intersection.p,...
                     v,...
                     zeros(1,size(intersection.p,2)),...
                     beam.wavelength,...
                     beam.polarization,...
                     beam.E0.*abs(Eff),...
                     ID*ones(1,size(intersection.p,2)),...
                     beam.dist + MV3norm(beam.v .* repmat(t.*beam.nm,3,1)) + double(Eff < 0).*beam.wavelength/2,... %pi/2 phase jump at dense surface (reflection)
                     nm);
             else
                 beam_interact = Beamclass();
             end
        end
                        
        function beam_interact = interactOpticalSurf(obj,beam,intersection)
            %set normal vector
            switch obj.geometry.type
                case OpticalSurf.geometryTypes{1}
                    n = obj.geometry.n;
                case OpticalSurf.geometryTypes{2}
                    n = obj.geometry.n;
                case OpticalSurf.geometryTypes{3}
                    %point on quad surface is defined through n = div F(x,y,z) =
                    %[dF/dx,dF/dy,dF/dy] at position p
                    p = intersection.p - repmat( obj.geometry.center,1,size(intersection.p,2) ); %point on surf
                    A = obj.geometry.A;
                    n = [A(1,1)*p(1,:) + A(1,2)*p(2,:) + A(1,3)*p(3,:);...
                         A(2,2)*p(2,:) + A(2,1)*p(1,:) + A(2,3)*p(3,:);...
                         A(3,3)*p(3,:) + A(3,1)*p(1,:) + A(3,2)*p(2,:)] ;
                    n = 2*(n + repmat(A(1:3,4),1,size(p,2)));
                    n = n./repmat(MV3norm(n),3,1);
                otherwise
                    fprintf('ERROR: geometry type not implemented\n');
            end
            
            %calculate new beam
            switch obj.interactionType
                case OpticalSurf.interactionTypes{1} %reflection
                    v = OpticalSurf.reflect(n,beam);
                    beam_interact = OpticalSurf.returnBeam(intersection,beam,v,1,obj.ID,intersection.t,beam.nm);
                case OpticalSurf.interactionTypes{2} %refraction
                    [n_normal,n_back] = OpticalSurf.refractiveIndex(obj,beam);
                    [v,eta,c12,n_beam] = OpticalSurf.refract(n,n_normal,n_back,beam);  
                    [E_refr,~] = OpticalSurf.efficiency(obj,beam,eta,c12);
                    beam_interact = OpticalSurf.returnBeam(intersection,beam,v,E_refr,obj.ID,intersection.t,n_beam);
                case OpticalSurf.interactionTypes{3} %refraction and reflection
                    vr = OpticalSurf.reflect(n,beam);
                    [n_normal,n_back] = OpticalSurf.refractiveIndex(obj,beam);
                    [vt,eta,c12,n_beam] = OpticalSurf.refract(n,n_normal,n_back,beam);
                    [E_refr,E_refl] = OpticalSurf.efficiency(obj,beam,eta,c12);
                    beam_r = OpticalSurf.returnBeam(intersection,beam,vr,E_refl,obj.ID,intersection.t,beam.nm);
                    beam_t = OpticalSurf.returnBeam(intersection,beam,vt,E_refr,obj.ID,intersection.t,n_beam);
                    beam_interact = mergeBeams(beam_r,beam_t);
                case OpticalSurf.interactionTypes{4} %diffraction
                    v = OpticalSurf.diffract(obj.interaction,n,beam);
                    beam_interact(length(obj.interaction.m)) = Beamclass();
                    for k=1:length(obj.interaction.m)
                        beam_interact(k) = OpticalSurf.returnBeam(intersection,beam,v(:,((k-1)*size(beam.v,2)+1):(k*size(beam.v,2))),1,obj.ID,intersection.t,beam.nm);
                    end
                    beam_interact = mergeBeams(beam_interact);
                    S = sum(abs(imag(v))) == 0;
                    beam_interact = mergeBeams(selectBeamByIndex(beam_interact,S));
                otherwise
                    fprintf('ERROR: interaction type not implemented!\n');
            end
        end
    end
end
