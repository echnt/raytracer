function beam = beam_template(v,r,beamDiameter,N,E,pol,wavelength,type,intensity_function)
%rudimentary options for quick generation of equispaced beams, see
%Beamclass for description of parameters

switch type
    case 'circular' %circular disc shaped beam profile
        phi = linspace(0,2*pi,N(2));
        x = linspace(0,beamDiameter/2,N(1));
        if(N(1) == 1)            
            if(N(2) == 1)
                x = 0;
                phi = 0;
            else
                x = beamDiameter/2; %constant distance from center for different phi
            end
        elseif(N(2) == 2)
            phi = [0,pi];
        end
        
        [v_sph(1),v_sph(2),v_sph(3)] = cart2sph(v(1),v(2),v(3));
        [v1(1),v1(2),v1(3)] = sph2cart(v_sph(1),v_sph(2)+pi/2,1);
        v2 = cross(v1,v/norm(v));
        p = zeros(3,N(2));
        for k=1:3
            p(k,:) = (sin(phi)*v1(k)+cos(phi)*v2(k));
        end
        beam = Beamclass();
        switch intensity_function
            case 'exp'
                Er = E.*exp(-(x/beamDiameter).^2);
            case 'gauss'
                Er = E.*exp(-(x/beamDiameter).^2).^2;
            case 'flat'
                Er = E*ones(size(x));
        end
        for k=1:length(x)
            for j=1:length(phi)
                for lI=1:length(wavelength)
                    beam = mergeBeams(Beamclass(p(:,j)*x(k)+r,v,0,wavelength(lI),pol,Er(k),0,0,1),beam); %char(java.util.UUID.randomUUID))
                end
            end
        end        
end
