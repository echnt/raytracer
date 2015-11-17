function h = plotRaytracing(OO,beam,I_t0_col,t_t0_len,t_tmax_len)
    %Plot beam and Optical Objects
    %I_t0_col color of beams with t==0
    %t_t0_len minimum length of beam if t==0
    %t_tmax_le maximum length of beam

    
%% Optical Object
    globalAlpha = 0.4;
    
    for k=1:length(OO)
        for j=1:length(OO(k).optElements)
            switch OO(k).optElements(j).geometry.type
                case OpticalSurf.geometryTypes{1} %disc plane
                    disc = getDisc(OO(k).optElements(j).geometry.radius,20,OO(k).optElements(j).geometry.n);
                    fillHandle = fill3(disc(1,:)+OO(k).optElements(j).geometry.center(1),disc(2,:)+OO(k).optElements(j).geometry.center(2),disc(3,:)+OO(k).optElements(j).geometry.center(3),disc(3,:),'Edgecolor',[1,1,1]); hold on;
                    alpha(fillHandle,globalAlpha*2)
                    hold on;
                case OpticalSurf.geometryTypes{2} %rect plane
                    r = OO(k).optElements(j).geometry.center;
                    n = OO(k).optElements(j).geometry.n;
                    h = OO(k).optElements(j).geometry.height;
                    w = OO(k).optElements(j).geometry.width*cross(n,h/norm(h));
                    s(1,:) = r+h/2+w/2;
                    s(2,:) = r+h/2-w/2;
                    s(3,:) = r-h/2-w/2;
                    s(4,:) = r-h/2+w/2;
                    fillHandle = fill3(s(:,1),...
                                       s(:,2),...
                                       s(:,3),...
                                       s(:,3),'Edgecolor',[1,1,1]); 
                    hold on;
                    alpha(fillHandle,globalAlpha)
                case OpticalSurf.geometryTypes{3} %plot quad surf
                    D = OO(k).optElements(j).geometry.D;
                    r = OO(k).optElements(j).geometry.center;
                    
                    x = linspace(-D/2,D/2,100);
                    y = x;
                    [x,y] = meshgrid(x,y);
                    A = OO(k).optElements(j).geometry.A;
                    %F = Ax^2+By^2+Cz^2 + 2Dxy + 2Exz + 2Fyz + 2Gx + 2Hy + 2Iz +J = 0
                    a = A(3,3);
                    b = 2*(A(3,1)*x + A(3,2)*y + A(3,4));
                    c = A(1,1)*x.^2 + A(2,2)*y.^2 + 2*(A(1,2)*x.*y + A(1,4)*x + A(2,4)*y + A(4,4));
                    if(a ~= 0)
                        F1 = -b./(2*a) - sqrt(b.^2./(4*a.^2) - c./a);
                        F2 = -b./(2*a) + sqrt(b.^2./(4*a.^2) - c./a);
                        S1 = abs(imag(F1)) == 0 & ((x.^2 + y.^2 + F1.^2) < (D/2)^2);
                        S2 = abs(imag(F2)) == 0 & ((x.^2 + y.^2 + F2.^2) < (D/2)^2);
                        
                        x = x+r(1);
                        y = y+r(2);
                        F1 = F1+r(3);
                        F2 = F2+r(3);
                        
                        if(sum(S1 == 1) == 0)
                            surf(x,y,F1,'edgecolor','none');
                        else
                            plot3(x(S1),y(S1),F1(S1),'b.');
                        end
                            hold on;
                        if(sum(S2 == 1) == 0)
                            surf(x,y,F2,'edgecolor','none');                            
                        else
                            plot3(x(S2),y(S2),F2(S2),'b.');
                        end
                    else
                        F = -c./b;
                        surf(x,y,F,'edgecolor','none');
                        hold on;
                    end   
            end
        end
    end
%% Beam    
    I = abs(beam.E0).^2;
    if(nargin == 5)
        I(beam.t == 0) = I_t0_col;
        beam.t(beam.t == 0) = t_t0_len;
        beam.t(beam.t > t_tmax_len) = t_tmax_len;
    else
        I(beam.t == 0) = 0.4;
        beam.t(beam.t == 0) = 0.05;
        beam.t(beam.t > 1) = 1;
    end
    I(I > 1) = 1;
    I = I.^0.4;
    I = I-min(I);
    I = I/max(I);
    
    X = [beam.r(1,:);beam.r(1,:)+beam.v(1,:).*beam.t];
    Y = [beam.r(2,:);beam.r(2,:)+beam.v(2,:).*beam.t];
    Z = [beam.r(3,:);beam.r(3,:)+beam.v(3,:).*beam.t];
    Ncol = 1024;
    cm = hot(Ncol);
    cm = flipud(cm);    
    
    for k=1:size(beam.r,2)     
        beamColor = cm(round(I(k)*(Ncol-1))+1,:).^2;
        if(I(k) > 1e-8)
            h = line(X(:,k),Y(:,k),Z(:,k),'LineWidth',1,'Color',beamColor);
            hold on;
        end
    end

    xlabel('X /m');
    ylabel('Y /m');
    zlabel('Z /m');
    grid on;
