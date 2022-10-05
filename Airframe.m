classdef Airframe < handle
    %AIRFRAME Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        % GEOMETRIC PROPERTIES
        %inputs
        OD;     % Outer Diameter, m
        
        %constants
        NZ_t;   % Nozzle thickness, m
        IT_L;   % Internals length, m
        
        %calculated
        OR;     % Outer Radius, m
        ID;     % Inner Diameter, m
        IR;     % Inner Radius, m
        t;      % Wall Thickness, m
        L;      % Length, 
        CSA;    % Cross Sectional Area, m^2
        FBH_t;  % Forward Bulkhead thickness, m
        PL;     % Prop Mass length, m
        Rs;         % Slenderness Ratio, unitless
        Rt;         % Transition Ratio, unitless
        
        % PROPULSION/OPERATING PROPERTIES
        %inputs
        MEOP;   % Maximum Expected Operating Pressure of Propellant, Mpa
        PM;     % Propellant mass, Kg
        PD % Propellant Density, g/cm^3
        Max_F;
        
        % MATERIAL PROPERTIES
        %inputs
        E;      % Young's Modulus, Gpa
        y;      % Yield Stress, Mpa
        v;      % Poisson's ratio, unitless
        rho;    % Density, Kg/m^3
       
        %DESIGN CRITERIA
        MASS;       % Airframe mass, Kg
        FBH_MASS;   % Forward Bulkhead mass, Kg
        NZ_MASS;    % Mass of the nozzle, Kg
        LBS;        % Local Buckling Stress, Kg
        JBS;        % Johnson Buckling Stress, Mpa
        EBS;        % Euler Buckling Stress, Mpa
        AS;         % Applied stress due to force at max Q, Mpa
        CT;         % Column Type, unitless
        CoM;        % Center of Mass, cm
        MoIx_af;    % Mass moment of inertia about x-axis, Kg-m^2   
        MoIy_af;    % Mass moment of inertia about y-axis, Kg-m^2
        MoIz_af;    % Mass moment of inertia about z-axis, Kg-m^2
        MoIx_prop;  % Mass moment of inertia about x-axis, Kg-m^2
        MoIy_prop;  % Mass moment of inertia about y-axis, Kg-m^2
        MoIz_prop;  % Mass moment of inertia about z-axis, Kg-m^2
    end
    
    methods
        function obj = Airframe(input1,input2,input3,input4,input5,input6)
            %AIRFRAME Construct an instance of this class
            %   Detailed explanation goes here
            obj.OD = input1;
            obj.MEOP = input2;
            obj.PM = input3;
            obj.PD = input4;
            obj.E = input5(1);
            obj.y = input5(2);
            obj.v = input5(3);
            obj.rho = input5(4);
            obj.Max_F = input6;
            
            obj.NZ_t = 0.127;
            obj.IT_L = 1.016;
            
            obj.Inner_Diameter(obj.OD, obj.y, obj.MEOP);
%             if obj.ID < obj.OD - 0.00254
%                 obj.ID = obj.OD - 0.00254;
%             end
            obj.OR = obj.OD/2;
            obj.IR = obj.ID/2;
            obj.t = obj.OR - obj.IR;
            if obj.t < 0.00127
                obj.t = 0.00127;
                obj.ID = obj.OD - 2*obj.t;
                obj.OR = obj.OD/2;
                obj.IR = obj.ID/2;
            end
            
            obj.Prop_Length(obj.PM, obj.PD, obj.IR)
            obj.Bulkhead_Thickness(obj.MEOP, obj.IR, obj.y);
            obj.Length(obj.PL, obj.FBH_t, obj.NZ_t, obj.IT_L);
            obj.Bulkhead_Masses(obj.FBH_t, obj.IR, obj.rho, obj.NZ_t);
            obj.Airframe_Mass(obj.OR, obj.IR, obj.L, obj.rho);
            obj.Center_of_Mass();
            obj.Local_Buckling(obj.E, obj.t, obj.OR, obj.IR, obj.v);
            obj.Column_Ratios(obj.OR, obj.IR, obj.L, obj.E, obj.y)
            if obj.Rs < obj.Rt
                obj.Johnson_Buckling(obj.y, obj.L, obj.E, obj.OR, obj.IR);
                obj.CT = 0;
            elseif obj.Rt < obj.Rs
                obj.Euler_Buckling(obj.L, obj.E, obj.OR, obj.IR);
                obj.CT = 1;
            end
            obj.AS = obj.Max_F/obj.CSA;
            
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
            % Outer Diameter    OD                      m
            % Yield strength    Y                       Mpa
            % Maximum expected  MEOP                    Mpa
            %  - operating pressure
            % 
            % Outputs(assingment):
            % Property          Variable Name           Units
            % Inner diameter    ID                      m
            
            HS_max = Y/2;
            obj.ID = (2 * HS_max * OD - MEOP * OD)/(MEOP + 2 * HS_max);

        end
        
        function Prop_Length(obj, PropMass, PropDensity, IR)
            % This function takes the propellant mass, propellant density,
            % and inner diameter to calculate the length of the Airframe
            % that is necessary to contain the propellant mass.
            %  
            % Inputs:
            % Property          Variable Name           Units
            % Propellant mass   PropMass                Kg
            % Prop. density     PropDensity             Kg/m^3
            % Inner Diameter    ID                      m
            % 
            % Outputs(assingment):
            % Property          Variable Name           Units
            % Prop. length      PL                      m
            
            X_sec_A = pi*(IR^2 - 0.00635^2);
            PropVol = PropMass/PropDensity;
            obj.PL = PropVol/X_sec_A ;
        end
        
        function Bulkhead_Thickness(obj, MEOP, R, Y)
            FBCS = Y/2;
            obj.FBH_t = sqrt((3 * MEOP * R^2)/(4 * FBCS));
        end
        
        function Length(obj, PL, BH_L, NZ_L, IT_L)
            % This function takes the propellant mass, propellant density,
            % and inner diameter to calculate the length of the Airframe
            % that is necessary to contain the propellant mass.
            %  
            % Inputs:
            % Property          Variable Name           Units
            % Prop. length      PL                      m
            % Bulkhead length   BH_L                    m
            % Nozzle length     NZ_L                    m
            % 
            % Outputs(assingment):
            % Property          Variable Name           Units
            % Prop. length      PL                      m
            
            obj.L = PL + BH_L + NZ_L + IT_L;
        end
        
        function Bulkhead_Masses(obj, t, r, rho, NZ_t)
            % This function takes the wall thickness, airframe inner
            % radius, inner radius, and density to calculate the mass of
            % the bulkheads that contain the propellant mass.
            %  
            % Inputs:
            % Property          Variable Name           Units
            % thickness         t                       m
            % Inner radius      r                       m
            % Density           rho                     m
            % 
            % Outputs(assingment):
            % Property          Variable Name           Units
            % Bulkhead masses   FBH_t                   m
            
            X_sec_A = pi*r^2;
            Volume = t * X_sec_A;
            obj.FBH_MASS = Volume * rho;
            %obj.NZ_MASS = rho * X_sec_A * NZ_t;
            obj.NZ_MASS = Volume * rho;
        end
        
        function Airframe_Mass(obj, OR, IR, L, density)
            % This function uses cross-sectional area, and length to
            % calculate volume, then applies density to the volume to
            % calculate mass.
            %
            % Inputs:
            % Property          Variable Name           Units
            % Outer Radius      OR                      m
            % Inner Radius      IR                      m
            % Length            L                       m
            % Density           density                 kg/m^3
            % 
            % Outputs(assingment):
            % Property          Variable Name           Units
            % Mass              mass                    kg
            
            X_sec_A = (pi .* OR.^2) - (pi .* IR.^2);  % Cross Sectional Area, in^2
            Volume = X_sec_A .* L;  % Volume, in^3
            mass = Volume * density;
            obj.MASS = mass;  % Mass, lbs
            obj.CSA = X_sec_A;
            
        end
        
        function Center_of_Mass(obj)
            % This function uses the AF instance length to calculate the
            % center of mass of the component.
            %  
            % Inputs: (none)
            % Property          Variable Name           Units
            %
            % 
            % Outputs(assingment):
            % Property          Variable Name           Units
            % Center of Mass    CoM                     m (from the top of
            %                                             the airframe)
            
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
            % Outer radius      OR                      m
            % Inner radius      IR                      m
            % Length            L                       m
            % Young's Modulus   E                       Mpa
            % Yield Strength    Y                       Mpa
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
            % This function uses the length, material yield strength,
            % Young's Modulus, outer radius, and inner radius to calculate
            % the Euler buckling stress.
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
        
    end
end

