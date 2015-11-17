function shape = getDisc(r,N,normal_vector)
    %approximate disc shaped structure by polygon plane
    %r center position
    %N number of edges
    %normal_vector plane orientation
    [v_sph(1),v_sph(2),v_sph(3)] = cart2sph(normal_vector(1),normal_vector(2),normal_vector(3));
    [v1(1),v1(2),v1(3)] = sph2cart(v_sph(1),v_sph(2)+pi/2,1);
    v2 = cross(v1,normal_vector/norm(normal_vector));
    phi = linspace(0,2*pi,N);
    shape = zeros(3,N);
    for k=1:3
        shape(k,:) = r*(sin(phi)*v1(k)+cos(phi)*v2(k));
    end
