classdef Airframe < handle
    %AIRFRAME: Given material properties and operating inputs, this class
    %   will automatically calculate geometrical properties for a sub 
    %   minimum diameter aiframe using built in safety factors. Along with 
    %   those dimensions, certain strucutral and physical properties are
    %   calculated to be used in the analysis and evaluation of the design.
    % 
    
    properties
        % GEOMETRIC PROPERTIES
        %inputs
        OD;     % Outer Diameter, in
        
        %constants
        NZ_t;   % Nozzle thickness, in
        IT_L;   % Internals length, in
        
        %calculated
        OR;     % Outer Radius, in
        ID;     % Inner Diameter, in
        IR;     % Inner Radius, in
        t;      % Wall Thickness, in
        L;      % Length, in
        CSA;    % Cross Sectional Area, in^2
        FBH_t;  % Forward Bulkhead thickness, in
        PL;     % Prop Mass length, in
        Rs;         % Slenderness Ratio, unitless
        Rt;         % Transition Ratio, unitless
        
        % PROPULSION/OPERATING PROPERTIES
        %inputs
        MEOP;   % Maximum Expected Operating Pressure of Propellant, ksi
        PM;     % Propellant mass, lbs
        RhoProp % Propellant Density, lbs/in^3
        Max_dynamic_F;
        
        % MATERIAL PROPERTIES
        %inputs
        E;      % Young's Modulus, ksi
        y;      % Yield Stress, ksi
        v;      % Poisson's ratio, unitless
        rho;    % Density, lbs/in^3
       
        %DESIGN CRITERIA
        MASS;       % Airframe mass, lbs
        FBH_MASS;   % Forward Bulkhead mass, lbs
        NZ_MASS;    % Mass of the nozzle, lbs
        LBS;        % Local Buckling Stress, ksi
        JBS;        % Johnson Buckling Stress, ksi
        EBS;        % Euler Buckling Stress, ksi
        AS;         % Applied stress due to force at max Q, ksi
        CT;         % Column Type, unitless
        CoM;        % Center of Mass, in
        MoIx;       % Mass moment of inertia about x-axis, lb-ft*s^2   
        MoIy;       % Mass moment of inertia about y-axis, lb-ft*s^2
        MoIz;       % Mass moment of inertia about z-axis, lb-ft*s^2
        
    end
    
    methods
        function obj = Airframe(input1,input2, input3, input4, input5, input6, input7, input8, input9)
            %AIRFRAME Construct an instance of this class
            %   Detailed explanation goes here
            % INPUTS
            obj.OD = input1;        % inputs
            obj.MEOP = input2;      % "    "
            obj.PM = input3;        % "    "
            obj.RhoProp = input4;   % "    "
            obj.E = input5;         % "    "
            obj.y = input6;         % "    "
            obj.v = input7;         % "    "
            obj.rho = input8;       % "    "
            obj.Max_dynamic_F = input9;
            
            % CONSTANTS
            obj.NZ_t = 5;
            obj.IT_L = 60;
            
            % GEOMETRIC PROPERTIES
            obj.Inner_Diameter(obj.OD, obj.y, obj.MEOP); % ID func call
            obj.t = (obj.OD - obj.ID)/2; % Wall thickness calc
            obj.OR = obj.OD/2;      % Outer radius calc
            obj.IR = obj.ID/2;      % Inner radius calc
            
            % FUNCTION CALLS(misc.)
            obj.Prop_Length(obj.PM, obj.RhoProp, obj.IR);
            obj.Local_Buckling(obj.E, obj.t, obj.OR, obj.IR, obj.v); 
            obj.CSA = (pi .* obj.OR.^2) - (pi .* obj.IR.^2);
            obj.Moment_of_Inertia(obj.OR, obj.IR, obj.L, obj.MASS);
            obj.Bulkhead_Thickness(obj.MEOP, obj.IR, obj.y);
            obj.Length(obj.PL, obj.FBH_t, obj.NZ_t, obj.IT_L);
            obj.Column_Ratios(obj.OR, obj.IR, obj.L, obj.E, obj.y)
            obj.Airframe_Mass(obj.OR, obj.IR, obj.L, obj.rho);
            obj.Bulkhead_Masses(obj.FBH_t, obj.IR, obj.rho, obj.NZ_t);
            obj.Center_of_Mass();
            if obj.Rs < obj.Rt
                obj.Johnson_Buckling(obj.y, obj.L, obj.E, obj.OR, obj.IR);
                obj.CT = 0;
            elseif obj.Rt < obj.Rs
                obj.Euler_Buckling(obj.L, obj.E, obj.OR, obj.IR);
                obj.CT = 1;
            end
            obj.AS = obj.Max_dynamic_F/obj.CSA;
            
        end
        
        function Inner_Diameter(obj, OD, Y, MEOP)
            % This function takes the outer diameter, yield strength, and
            % the maximum expected operating pressure of the
            % propellant/motor. It should be noted that, for the formula
            % being used to calculate inner diameter in this fucntion,
            % standards say that 1.5x or 2x the actual MEOP be used as
            % input, which is what has been assumed MEOP is here.
            %  
            % Inputs:
            % Property          Variable Name           Units
            % Outer Diameter    OD                      in
            % Yield strength    Y                       ksi
            % Maximum expected  MEOP                    ksi
            %  - operating pressure
            % 
            % Outputs(assingment):
            % Property          Variable Name           Units
            % Inner diameter    ID                      in
            
            HS_max = Y/5;
            obj.ID = (2 * HS_max * OD - MEOP * OD)/(MEOP + 2 * HS_max);

        end
        
        function Prop_Length(obj, PropMass, PropDensity, IR)
            % This function takes the propellant mass, propellant density,
            % and inner diameter to calculate the length of the Airframe
            % that is necessary to contain the propellant mass.
            %  
            % Inputs:
            % Property          Variable Name           Units
            % Propellant mass   PropMass                lbs
            % Prop. density     PropDensity             lbs/in^3
            % Inner Diameter    ID                      in
            % 
            % Outputs(assingment):
            % Property          Variable Name           Units
            % Prop. length      PL                      in
            
            X_sec_A = pi*(IR^2 - 0.25^2);
            PropVol = PropMass/PropDensity;
            obj.PL = PropVol/X_sec_A ;
        end
        
        function Bulkhead_Thickness(obj, MEOP, R, Y)
            FBCS = Y/2;
            obj.FBH_t = sqrt((3 * MEOP * R^2)/(4 * FBCS));
        end
        
        function Length(obj, PL, BH_L, NZ_L, IT)
            % This function takes the propellant mass, propellant density,
            % and inner diameter to calculate the length of the Airframe
            % that is necessary to contain the propellant mass.
            %  
            % Inputs:
            % Property          Variable Name           Units
            % Prop. length      PL                      in
            % Bulkhead length   BH_L                    in
            % Nozzle length     NZ_L                    in
            % 
            % Outputs(assingment):
            % Property          Variable Name           Units
            % Prop. length      PL                      in
            
            obj.L = PL + BH_L + NZ_L + IT;
        end
        
        function Bulkhead_Masses(obj, t, r, rho, NZ_t)
            X_sec_A = pi*r^2;
            Volume = t * X_sec_A;
            obj.FBH_MASS = Volume * rho;
            %obj.NZ_MASS = rho * X_sec_A * NZ_t;
            obj.NZ_MASS = rho * Volume;
        end
        
        function Airframe_Mass(obj, OR, IR, L, density)
            % This function uses cross-sectional area, and length to
            % calculate volume, then applies density to the volume to
            % calculate mass.
            %
            % Inputs:
            % Property          Variable Name           Units
            % Outer Radius      OR                      in
            % Inner Radius      IR                      in
            % Length            L                       in
            % Density           density                 lbs/in^3
            % 
            % Outputs(assingment):
            % Property          Variable Name           Units
            % Mass              mass                    lbs
            
            X_sec_A = (pi .* OR.^2) - (pi .* IR.^2);  % Cross Sectional Area, in^2
            Volume = X_sec_A .* L;  % Volume, in^3
            mass = Volume * density;
            obj.MASS = mass;  % Mass, lbs
            
        end
        
        function Center_of_Mass(obj)
            
            obj.CoM = obj.L/2;
            
        end
        
        function Local_Buckling(obj, E, t, OR, IR, v)
            % This function uses Young's Modulus, wall thickness, outer
            % diameter, and poisson's ratio to calculate the local buckling
            % stress of the airframe.
            %
            % Inputs:
            % Property          Variable Name           Units
            % Young's Modulus   E                       ksi
            % Wall Thickness    t                       in
            % Outer Diameter    OD                      in
            % Poisson's ratio   v                       unitless
            %
            % Outputs(assignment):
            % Property          Variable                Units
            % Local Buckling    lbs                     ksi
            %   Stress
            
            r = (OR + IR)/2;  % Mean radius, in
            phi = 1/16 * sqrt(r/t);  % 
            gamma = 1 - 0.901 * (1 - exp(-1 * phi));  % Buckling knockdown factor
            
            obj.LBS = (gamma * E)/sqrt(3 * (1 - v^2)) * t/r;  % Critical Stress, ksi
            
        end
        
        function Column_Ratios(obj, OR, IR, L, E, Y)
            % This function uses outer radius, inner radius, Young's
            % modulus, yield strength, and length to calculate the
            % slenderness ratio and the transition ratio for the airframe
            % to determine the buckling stress case.
            % Inputs:
            % Property          Variable Name           Units
            % Outer radius      OR                      in
            % Inner radius      IR                      in
            % Length            L                       in
            % Young's Modulus   E                       ksi
            % Yield Strength    Y                       ksi
            % 
            % Outputs(assingment):
            % Property          Variable Name           Units
            % Slenderness ratio Rs                      unitless
            % Transition ratio  Rg                      unitless
            
            S = Y/2;
            r = sqrt((OR^2+IR^2)/2); % Radius of Gyration along axis of rotational symmetry
            obj.Rs = L/r;
            obj.Rt = sqrt((2*pi^2*E)/(4*S));
            
        end
        
        function Johnson_Buckling(obj, Y, L, E, OR, IR)
            % This function uses the material yield strength, length,
            % Young's Modulus, outer radius, and inner radius to calculate
            % the Johnson buckling stress.
            % 
            % Inputs:
            % Property          Variable Name           Units
            % Yield Strength    Y                       ksi
            % Length            L                       in
            % Young's Modulus   E                       ksi
            % Outer radius      OR                      in
            % Inner radius      IR                      in
            % 
            % Outputs(assingment):
            % Property          Variable Name           Units
            % Johnson Buckling  JBS                     ksi
            
            r = sqrt((OR^2+IR^2)/2); % Radius of Gyration along Y axis(axis of symmetry)
            half_Y = Y/2;  % Half of the yield Strength
            
            obj.JBS = half_Y - ((half_Y*2*L)/(2*pi*r))^2 * (1/E);
            
        end
        
        function Euler_Buckling(obj, L, E, OR, IR)
            % This function uses the material yield strength, length,
            % Young's Modulus, outer radius, and inner radius to calculate
            % the Johnson buckling stress.
            % 
            % Inputs:
            % Property          Variable Name           Units
            % Length            L                       in
            % Young's Modulus   E                       ksi
            % Outer radius      OR                      in
            % Inner radius      IR                      in
            % 
            % Outputs(assingment):
            % Property          Variable Name           Units
            % Euler Buckling    EBS                     ksi
            
            r = sqrt((OR^2+IR^2)/2);
            obj.EBS = (pi^2 * E)/((2*(L/r))^2); 
            
        end
        
        function Moment_of_Inertia(obj, OR, IR, L, Mass)
            % This function uses the inner radius and outer radius to
            % calculate the 2nd moment of area(MoI) of the airframe about the 3
            % principle axes(X, Y, Z).
            % 
            % Inputs:
            % Property          Variable Name           Units
            % Inner radius      IR                      in
            % Outer radius      OR                      in
            % 
            % Outputs(assingments):
            % Property          Variable Name           Units
            % MoI about x       moix                    lb-ft*s^2
            % MoI about y       moiy                    lb-ft*s^2
            % MoI about z       moiz                    lb-ft*s^2
            
            moix = 1/12 * Mass * (3*(OR^2+IR^2) + L^2);
            moiy = 1/12 * Mass * (3*(OR^2+IR^2) + L^2);
            moiz = 1/2 * Mass * (OR^2+IR^2);
            
            obj.MoIx = moix;
            obj.MoIy = moiy;
            obj.MoIz = moiz;
            
        end
    end
end

