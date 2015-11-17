classdef Beamclass
    %provides a class for beam handling including some basic handling
    %functions like merge and array like access, expects vectors as 3xn for
    %raytrace handling
	properties
        r; %start point, 3xn vec
        v; %direction of beam, 3xn vec
        t; %multiplier r+t*v, scalar 1xn
        
        wavelength; %scalar 1xn
        polarization; %0=undef,1=S,2=P; should be a jones vector 4x1, 1xn scalar
        E0; %E field Amplitude, scalar 1xn
        ID; %random double ID, scalar 1xn,
            %should be unique 128bit ID by java.util.UUID.randomUUID
        dist; % propagated distance, scalar 1xn
        nm; %refractive index of medium, scalar 1xn
    end
    methods
        function obj = Beamclass(r,v,t,wavelength,polarization,E0,ID,dist,nm)
            switch nargin
                case 0
                    obj.r = [];
                    obj.v = [];
                    obj.t = [];
                    obj.wavelength = [];
                    obj.polarization = [];
                    obj.E0 = [];
                    obj.ID = [];
                    obj.dist = 0;
                    obj.nm = 0;
                otherwise
                    obj.r = r;
                    obj.v = v;
                    obj.t = t;
                    obj.wavelength = wavelength;
                    obj.polarization = polarization;
                    obj.E0 = E0;
                    obj.ID = ID;
                    obj.dist = dist;
                    obj.nm = nm;
            end
        end
        
        function obj1 = mergeBeams(obj1, obj2)
            %merge beams obj1,obj2 or array of obj1
            switch nargin
                case 1
                    for k=2:length(obj1)
                        obj1(1) = mergeBeams(obj1(1),obj1(k));
                    end
                    obj1 = obj1(1);
                case 2
                    N1 = size(obj1.v,2);
                    N2 = size(obj2.v,2);
                    if(N2 > 0)
                        obj1.r(:,N1+(1:N2)) = obj2.r;
                        obj1.v(:,N1+(1:N2)) = obj2.v;
                        obj1.t(N1+(1:N2)) = obj2.t;

                        obj1.wavelength(N1+(1:N2)) = obj2.wavelength;
                        obj1.polarization(:,N1+(1:N2)) = obj2.polarization;
                        obj1.E0(N1+(1:N2)) = obj2.E0;
                        obj1.ID(N1+(1:N2)) = obj2.ID;
                        obj1.dist(N1+(1:N2)) = obj2.dist;
                        obj1.nm(N1+(1:N2)) = obj2.nm;
                    end
            end
        end
        
        function obj = selectBeamByIndex(obj2,S)
            %select (S) beams in obj2 in array like style
            obj = Beamclass(obj2.r(:,S),obj2.v(:,S),obj2.t(S),...
                            obj2.wavelength(S),obj2.polarization(:,S),...
                            obj2.E0(S),obj2.ID(S),obj2.dist(S),obj2.nm(S));
        end
        
        function obj = updateTvals(obj,intersection)
            %updates beam.t with values from intersection
            obj.t = intersection.t;
        end
    end
end
