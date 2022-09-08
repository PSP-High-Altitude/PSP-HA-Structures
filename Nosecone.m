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
        function obj = Nosecone(input1,input2, input3, input4, input5, input6, input7, input8, input9, input10, input11, input12, input13)
            %NOSECONE Construct an instance of this class
            %   Detailed explanation goes here
            obj.SHAPE = input1;
            obj.PARAMETER = input2;
            obj.t = input3;
            obj.AF_t = input4;
            obj.L = input5;
            obj.tip_L = input6;
            obj.AF_ID = input7;
            obj.rho1 = input8;
            obj.rho2 = input9;
            obj.E = input10;
            obj.v = input11;
            obj.y = input12;
            obj.NET_F = input13;
            obj.Nosecone_Mass(obj.SHAPE, obj.PARAMETER, obj.t, obj.L, obj.AF_t, obj.AF_ID, obj.rho1, obj.rho2, obj.tip_L);
            obj.NC_Buckling(obj.E, obj.t, obj.v);
            obj.Center_of_Mass(0, obj.L+obj.tip_L, 0, obj.PROFILE, obj.L);
            Moment_of_Inertia_nc();
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

        function NC_Buckling(obj, E, t, v)
            
            Psi = .25/cosd(5.74);
            obj.NCCL = 1000 * 0.33 * (2*pi*E*t^2*cosd(5.74)^2)/(sqrt(3*(1-v^2)));
            
            obj.NCBS = obj.NET_F/(2*pi*Psi*t*cosd(5.74)^2);
           
        end
        
        function Center_of_Mass(obj, xmin, xmax, ymin, ymax, length)
            Moment_fun = @(x, y) x;
            Moment = integral2(Moment_fun, xmin, xmax, ymin, ymax);
            Mass = integral(ymax, 0, length);
            obj.CoM = Moment/Mass;
            
        end
        
        %Moment of Inertia 
        function Moment_of_Inertia_nc(obj, MASS,tip_L, AF_ID, AF_t)
            AF_OD = AF_ID + 2 * AF_t;
            c = 0.7; %constant for assuming the moment for the nosecone is equivalent to that of a cylinder

            MoI_z = (3/10) * Mass * ((AF_OD / 2) ^ 2 - (AF_ID / 2) ^ 2);
            MoI_x = c * 1/12 * MASS * ((3*((AF_OD / 2)^2+(AF_IF / 2)^2)) + tip_L^2); 
            MoI_y = c * 1/12 * MASS * ((3*((AF_OD / 2)^2+(AF_IF / 2)^2)) + tip_L^2); 

            obj.MoIz_nc = MoI_z;
            obj.MoIy_nc = MoI_y;
            obj.MoIx_nc = MoI_x;

        end


        function Moment_of_Inertia_tip(obj, MASS_tip, tip_L)

            OR = tip_L / 5; %5 to 1 ratio
            c = 0.7; %constant for assuming the moment for the nosecone tip is equivalent to that of a cylinder
            MoI_z = (3/10) * Mass_tip * (OR ^ 2);
            MoI_x = c * ((1/4) * Mass_tip * (OR) ^ 2 + 1/12 * Mass_tip * tip_L ^ 2);
            MoI_y = c * ((1/4) * Mass_tip * (OR) ^ 2 + 1/12 * Mass_tip * tip_L ^ 2);

            obj.MoIz_tip = MoI_z;
            obj.MoIx_tip = MoI_x;
            obj.MoIy_tip = MoI_y;

        end 



%Combine both moments of inertia to get a final one for the nosecone and
%the nosecone tip combined
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



