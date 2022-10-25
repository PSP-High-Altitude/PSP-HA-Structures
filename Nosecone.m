classdef Nosecone < handle
    %NOSECONE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        % NOSCONE SPECIFIC GEOMETRY PARAMETERS
        SHAPE;
        PARAMETER;
        PROFILE;
        EF;
        
        % GEOMETRIC PROPERTIES
        t;
        AF_t;
        L;
        tip_L;
        AF_ID;
        AF_OD;
       
        % MATERIAL PROPERTIES
        rho1;
        rho2;
        E;
        v;
        y;
        
        % OPERATING PROPERTIES
        NET_F;
        
        % DESIGN CRITERIA
        MASS;
        NCCL;
        NCBS;
        CoM;
        MoIy_tip;
        MoIx_tip;
        MoIz_tip;
        MoIx_nc;
        MoIy_nc;
        MoIz_nc;
    end
    
    methods
        function obj = Nosecone(input1,input2,input3,input4,input5,input6,input7,input8,input9)
            %NOSECONE Construct an instance of this class
            %   Detailed explanation goes here
            obj.SHAPE = input1;
            obj.PARAMETER = input2;
            obj.t = input3;
            obj.AF_t = input4;
            obj.tip_L = input5;
            obj.AF_ID = input6;
            obj.E = input7(1);
            obj.y = input7(2);
            obj.v = input7(3);
            obj.rho2 = input7(4);
            obj.rho1 = input8;
            obj.NET_F = input9;
            
            obj.AF_OD = obj.AF_ID + 2*obj.AF_t;
            obj.L = obj.AF_OD*5;
            obj.Nosecone_Mass(obj.SHAPE, obj.PARAMETER, obj.t, obj.L, obj.AF_t, obj.AF_ID, obj.rho1, obj.rho2, obj.tip_L);
            obj.Center_of_Mass(0, obj.L+obj.tip_L, 0, obj.PROFILE, obj.L);
        end
        
        function Nosecone_Mass(obj, shape, shapeParameter, wallThickness, length, AF_thickness, AF_inner_diameter, tipDensity, bodyDensity, tipLength)

            %calc tip volume
            %calc frustum volume
            %multiply by respective material densities
            %add

            baseDiameter = 2 * AF_thickness + AF_inner_diameter;

            tipVol = 0; % reassigned if there is a solid tip
            baseR = baseDiameter / 2;

            % assign shape equation
            if shape == 1 % power series
                f = @(x) baseR * (x / length) .^ shapeParameter;
            elseif shape == 2 % haak series
                angle = @(x) acos(1 - ((2 * x) / length));
                %obj.EF = angle;
                f = @(x) (baseR / sqrt(pi)) * sqrt(angle(x) - (sin(2 * angle(x)) / 2) + shapeParameter * (sin(angle(x))) .^ 3);
            elseif shape == 3 % conical
                f = @(x) (x * baseR) / length;
            elseif shape == 4 % tangent ogive
                rho = (baseR .^ 2 + length .^ 2) / (2 * baseR);
                f = @(x) sqrt( (rho .^ 2) - (length - x) .^ 2) + baseR - rho;
            elseif shape == 5 %elliptical
                f = @(x) baseR * sqrt(1 - (x .^ 2 / length .^ 2));
            end

            if tipLength > 0
                tipVol = calcTipVol(f, tipLength);
            end

            hollowVol = calcHollowVol(f, wallThickness, tipLength, length);

            tipMass = tipVol * tipDensity;
            hollowMass = hollowVol * bodyDensity;

            obj.MASS = tipMass + hollowMass;
            obj.PROFILE = f;
        end
        
        function Center_of_Mass(obj, xmin, xmax, ymin, ymax, length)
            Moment_fun = @(x, y) x;
            Moment = integral2(Moment_fun, xmin, xmax, ymin, ymax);
            Mass = integral(ymax, 0, length);
            obj.CoM = Moment/Mass;
            
        end
        
    end
end

% CALC VOLUME OF HOLLOW PORTION (FRUSTUM OR CONE)
function hollowVol = calcHollowVol(shapeFunction, wallThickness, tipLength, length)
    % outer f squared - inner f squared inside the integral

    outerF = @(x) shapeFunction(x) .^ 2;
    innerF = @(x) (shapeFunction(x) - wallThickness) .^ 2;

    integrand = @(x) outerF(x) - innerF(x); 

    hollowVol = pi * integral(integrand, tipLength, length);
end


% CALC VOLUME OF TIP
function tipVol = calcTipVol(shapeFunction, tipLength)
    integrand = @(x) shapeFunction(x) .^2;
    tipVol = pi * integral(integrand, 0, tipLength);
end

