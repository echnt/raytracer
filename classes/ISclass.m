classdef ISclass
    %finds and saves 3D intersections of Beamclass and OpticalSurf class.
	properties
        p; %intersect point
        t; %multiplier from r+v*t
        angle_rad; %angle of incidence to surf normal vector in radiants
        isIntersect; %array of intersections in beam
        index; %min index in array of ISclass
        ID; %Optical Surf ID
    end
    methods
        function obj = ISclass(p,t,angle_rad,isIntersect,index)
            switch nargin
                case 4
                    obj.p = p;
                    obj.t = t;
                    obj.angle_rad = angle_rad;
                    obj.isIntersect = isIntersect;
                    obj.index = [];
                    obj.ID = [];
                case 5
                    obj.p = p;
                    obj.t = t;
                    obj.angle_rad = angle_rad;
                    obj.isIntersect = isIntersect;
                    obj.index = index;
                    obj.ID = [];
                otherwise
                    obj.p = [];
                    obj.t = [];
                    obj.angle_rad = [];       
                    obj.isIntersect = [];
                    obj.index = [];
                    obj.ID = [];
            end
        end
        
        function obj = getMinTfromArray(in)
            %compiles ISclass of shortest intersections from ISclass array, expects same sized
            %array lengths in ISclass objects
            obj = in(1);
            obj.index = ones(size(obj.t));
            for k=2:length(in)
                %test if beam is smallest, nonzero, noncomplex, valid solution
                %on surf
                v = ((in(k).t < obj.t | ~obj.isIntersect) & (in(k).t > 0) & in(k).isIntersect & isreal(in(k).t));     
                obj.p(:,v) = in(k).p(:,v);
                obj.t(v) = in(k).t(v);
                obj.angle_rad(v) = in(k).angle_rad(v);
                obj.isIntersect(v) = in(k).isIntersect(v);
                obj.index(v) = k;
                obj.ID(v) = in(k).ID(v);
            end
        end
        
        function obj = selectISbyIndex(in,S)
            %array style access of ISclass in with index S
            obj.p = in.p(:,S);
            obj.t = in.t(S);
            obj.angle_rad = in.angle_rad(S);
            obj.isIntersect = in.isIntersect(S);
            obj.index = in.index(S);
            obj.ID = in.ID(S);
        end
        
        function obj1 = mergeISclass(obj1,obj2)
            %merge obj1 in array style with obj2
            switch nargin
                case 1 %array of ISclass
                    for k=2:length(obj1)
                        obj1(1) = mergeISclass(obj1(1),obj1(k));
                    end
                    obj1 = obj1(1);
                case 2 %two ISclass
                    N1 = size(obj1.p,2);
                    N2 = size(obj2.p,2);
                    obj1.p(:,N1+(1:N2)) = obj2.p;
                    obj1.t(N1+(1:N2)) = obj2.t;
                    obj1.angle_rad(:,N1+(1:N2)) = obj2.angle_rad;
                    obj1.isIntersect(:,N1+(1:N2)) = obj2.isIntersect;
                    if(~isempty(obj2.index))
                        obj1.index(:,N1+(1:N2)) = obj2.index;
                    end
                    if(~isempty(obj2.ID))
                        obj1.ID(:,N1+(1:N2)) = obj2.ID;
                    end
            end
        end
    end
    
    methods(Static)
        
        function obj = intersectPlane(optic,beam,mode)
            %finds the intersection of a line with a plane fullfilling the
            %equations:
            %Plane: (point on plane-support point of plane)*n = 0 (Hesse
            %Normal Form)
            %Line:  r+v*t = point on Line;
            dn = dot(beam.v,repmat(optic.geometry.n,1,size(beam.v,2)));
            tDP = dot((repmat(optic.geometry.center,1,size(beam.r,2)) - beam.r),repmat(optic.geometry.n,1,size(beam.r,2)) )./dn;
            pDP = beam.r+repmat(tDP,3,1).*beam.v;
            angle_radDP = acos( dn./MV3norm(beam.v)/norm(optic.geometry.n) );
            switch mode
                case 'disc'
                    isIntersectDP = (tDP > eps*10) .* (optic.geometry.radius >= MV3norm(repmat(optic.geometry.center,1,size(beam.r,2)) - pDP)); %eps*10: exclude beams starting on the surface
                case 'rect'
                    dist = repmat(optic.geometry.center,1,size(beam.r,2)) - pDP;
                    h = optic.geometry.height/norm(optic.geometry.height);
                    ph = dot(dist,repmat(h,1,size(beam.v,2)));
                    pw = dot(dist,repmat( cross(h,optic.geometry.n),1,size(beam.v,2)));
                    isIntersectDP = (tDP > eps*10) & abs(ph) <= norm(optic.geometry.height)/2 & abs(pw) < optic.geometry.width/2; %eps*10: exclude beams starting on the surface
            end
            isIntersectDP = boolean(isIntersectDP);
            obj = ISclass(pDP,tDP,angle_radDP,isIntersectDP);
        end
        
        function obj = intersectQuadSurf(optic,beam)
            %returns if an quad surf structure fullfilling the condition
            %(vt+r-c)' A (vt+r-c) = 0 where A consists of eigenvalues
            %[1/ai^2,-1] for ellispoids/spheres representing the radial extend of the object. Or in general:
            %A = [ADEG;
            %     DBFH;
            %     EFCI;
            %     GHIJ]
            %fullfilling the second order equation:
            %Ax^2+By^2+Cz^2+2Dxy+2Exy+2Fxy+2Gx+2Hy+2Jz+I = 0;
            %A can be transformed afterwards to account for object rotations. v is
            %the beam vector and c (the center of the object) - r shift vector of beam.
            %Note: define normal vector and radius to limit intersection to
            %parts of spheres/ellipsoids, otherwise set inf
            %use plane intersects for quicker computation
            
            %for an overview of second order surfaces see e.g.
            %http://www.cs.mtu.edu/~shene/COURSES/cs3621/NOTES/geometry/simple.html
            %Ax^2+By^2+Cz^2 + 2Dxy + 2Exz + 2Fyz + 2Gx + 2Hy + 2Iz +J = 0
            
            %Problem to solve: get t from (vt+r-C)' A (vt+r-C) = 1
            %0 = a*t^2 + b*t + c with
            %a = v' A v
            %b = v' A r + r'T A v
            %c = r' A r
            %t = -b/2a +- sqrt( (b/2a)^2 - c/a ) if a ~= 0!
            
            if(size(optic.geometry.A,1) == 3 && size(optic.geometry.A,2) == 3)
                A = [ [optic.geometry.A,zeros(3,1)] ; [zeros(1,3),-1]];
            else
                A = optic.geometry.A;
            end
            r = [beam.r;ones(1,size(beam.r,2))];
            v = [beam.v;zeros(1,size(beam.v,2))];
            
            m = size(r,2);
            R = r - repmat([optic.geometry.center;0],1,m);
            a = v' * A * v;
            a = a( sub2ind(size(a),1:size(a,1),1:size(a,2)) );
            b = v' * A * R + R' * A * v; %yes, A is symmetric
            b = b( sub2ind(size(b),1:size(b,1),1:size(b,2)) );
            c = R' * A * R;
            c = c( sub2ind(size(c),1:size(c,1),1:size(c,2)) );
            
            S = a~=0;
            t(1,S) = (-b(S)-sqrt(b(S).^2-4*a(S).*c(S)))./(2*a(S)); % solve parametric quadratic equation with 2 possible solutions
            t(2,S) = (-b(S)+sqrt(b(S).^2-4*a(S).*c(S)))./(2*a(S));
            
            t(1,~S) = -c(~S)./b(~S);
            t(2,~S) = 0;

            
            %select if hit and process first hit
            isIntersect = isreal(t) & t > 10*eps;  %eps*10: exclude beams starting on the surface
            
            p = repmat(beam.v,1,2).*repmat([t(1,:),t(2,:)],3,1)+repmat(beam.r,1,2);
            pc = p-repmat(optic.geometry.center,1,2*m);
            pd = dot(pc,repmat(optic.geometry.n,1,2*m));
            
            %intersection condition: p must be within radius/cylinder of n
            %with radius
            dist = sqrt(MV3norm(pc).^2 - pd.^2); %distance from optical axis            
            dist = [dist(1:size(beam.v,2));dist( (size(beam.v,2)+1) : (size(beam.v,2)*2))];
            pd = [pd(1:size(beam.v,2));pd( (size(beam.v,2)+1) : (size(beam.v,2)*2))];
            isIntersect = isIntersect & (dist <= optic.geometry.D/2) & (pd >= -10*eps);
                        
            t(~isIntersect) = NaN;
            [ts,I] = min(t);    
            Ip = (I-1)*size(beam.v,2)+(1:size(beam.v,2));
            p = p(:,Ip);
            isIntersect = isIntersect(sub2ind(size(isIntersect),I,1:size(isIntersect,2)));
            
            obj = ISclass(p,ts,zeros(size(ts)),isIntersect); %angle is not returned as the interaction plane would have to be computed
        end
        
        function obj = intersectOpticalSurf(optic,beam)
            %function to select correct intersection function
            switch optic.geometry.type
                case OpticalSurf.geometryTypes{1}
                    obj = ISclass.intersectPlane(optic,beam,'disc');
                case OpticalSurf.geometryTypes{2}
                    obj = ISclass.intersectPlane(optic,beam,'rect');
                case OpticalSurf.geometryTypes{3}
                    obj = ISclass.intersectQuadSurf(optic,beam);
                otherwise
                    fprintf('Error: geometry type unkown!\n');
            end
        end
    end

end
