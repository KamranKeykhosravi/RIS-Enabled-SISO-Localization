classdef SupPoint
    % This object represents a super point which have X, Y, and Z
    % dimensions, where the three dimensions are matrices of same size 
    %(Z can be a scalar which indicates all the points have the same Z)
    % this class is made to deal with multipe locations points at once which makes
    
    properties
        X                           %matrix of X positions
        Y                           %matrix of Y positions
        Z                           %z positions
        sinAz                       %sin of Azimuth 
        cosAz                       %cos of Azimuth 
        sinEl                       %sin of Elevation
        cosEl                       %cos of Elevation
        XY2                         %X.^2+Y.^2
        XY                          %sqrt(X.^2+Y.^2)
        abs2                        %X.^2+Y.^2+z.^2
        abs                         %sqrt(X.^2+Y.^2+z.^2)
    end
    methods
        function obj=SupPoint(X,Y,Z)
            % Constructor gets either three matrices or one 1X3 vector
            
            if nargin ==1 % in case a 1X3 vector is passed
                Y = X(2);
                Z = X(3);
                X = X(1);
            end
            obj.X = X;
            obj.Y = Y;
            obj.Z = Z;
            obj = obj.calcProperties;
        end
        function obj = calcProperties(obj)
            % calculates the attributes of superpoint
            
           obj.XY2 = obj.X.^2+obj.Y.^2;
           obj.XY = sqrt(obj.XY2);
           obj.abs2 = obj.XY2+obj.Z.^2;
           obj.abs = sqrt(obj.abs2);
           obj.cosAz = obj.X./obj.XY;
           obj.sinAz = obj.Y./obj.XY;
           obj.cosEl = obj.Z./obj.abs;
           obj.sinEl = obj.XY./obj.abs;
        end

        %%
        
        function k = getKv(obj,lambda)
            % calculates the wavenumber vector 
            k = -1*[(2*pi/lambda)*obj.sinEl(:).*obj.cosAz(:),(2*pi/lambda)*obj.sinEl(:).*obj.sinAz(:),(2*pi/lambda)*obj.cosEl(:)];
        end
        function k = getKDazv(obj,lambda)
            % calculates the derivative of the wavenumber vector with
            % respect to Azimuth angle
            k = -1*[-(2*pi/lambda)*obj.sinEl(:).*obj.sinAz(:),(2*pi/lambda)*obj.sinEl(:).*obj.cosAz(:),zeros(length(obj.cosAz(:)),1)];
        end
        function k = getKDelv(obj,lambda)
            % calculates the derivative of the wavenumber vector with
            % respect to Elevation angle
            k = -1*[(2*pi/lambda)*obj.cosEl(:).*obj.cosAz(:),(2*pi/lambda)*obj.cosEl(:).*obj.sinAz(:),-(2*pi/lambda)*obj.sinEl(:)];
        end
        
        %%
        function v = getV(obj,i,j)
            % outputs the 1X3 vector of dimensions of the point i,j
            
            if length(obj.Z) == 1
                z=obj.Z;
            else
                z=obj.Z(i,j);
            end
            v=[obj.X(i,j), obj.Y(i,j), z];
        end
    end
end