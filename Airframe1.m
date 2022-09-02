classdef Airframe1 < handle
    %AIRFRAME2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        %inputs
        OD_a;
        
        %constants
        NZ_t;   % Nozzle thickness, in
        IT_L;   % Internals length, in
        
        %calculated parameters
        ID_a;
        OD_c;
        ID_c;
        t_a;
        t_c;
        FBH_t;
        Rs;
        Rt;
        CSA;
        
        
        %casing material props
        rho_c;
        E_c;
        v_c;
        y_c;
        
        %aerostructure material props
        rho_a;
        E_a;
        v_a;
        y_a;
        
        %operating parameters
        MEOP;
        PM;
        PD;
        Max_dynamic_F;
        
        %propellant properties
        rhoProp
        
        
        %design criteria
        MASS_a;
        MASS_c;
        FBH_MASS;
        NZ_MASS;
        JBS;
        EBS;
        AS;
        SF;
        L_a;
        L_c;
        CoM;
    end
    
    methods
        function obj = Airframe1(input1, input2, input3, input4, input5, input6, input7)
            %AIRFRAME2 Construct an instance of this class
            %   Detailed explanation goes here
            
            % INPUTS
            obj.OD_a = input1;
            
            obj.rho_c = input2(1);
            obj.E_c = input2(2);
            obj.v_c = input2(3);
            obj.y_c = input2(4);
            
            obj.rho_a = input3(1);
            obj.E_a = input3(2);
            obj.v_a = input3(3);
            obj.y_a = input3(4);
            
            obj.MEOP = input4;
            obj.rhoProp = input5;
            obj.PM = input6;
            obj.Max_dynamic_F = input7;
            
            % CONSTANTS
            obj.NZ_t = 5;
            obj.IT_L = 60;
            
            
            obj.t_a = .125;
            obj.ID_a = obj.OD_a - 2*obj.t_a;
            obj.OD_c = obj.ID_a;
            obj.Inner_Diameter_c(obj.OD_c, obj.y_c, obj.MEOP);
            obj.t_c = (obj.OD_c - obj.ID_c)/2;
            obj.Prop_Length(obj.ID_c, obj.PM, obj.rhoProp);
            obj.Bulkhead_Thickness(obj.MEOP, obj.ID_c/2, obj.y_c);
            obj.Length(obj.IT_L, obj.L_c, obj.FBH_t, obj.NZ_t);
            obj.MASS_a = obj.Airframe_Mass(obj.OD_a/2, obj.ID_a/2, obj.L_a, obj.rho_a);
            obj.MASS_c = obj.Airframe_Mass(obj.OD_c/2, obj.ID_c/2, obj.L_c, obj.rho_c);
            obj.Bulkhead_Masses(obj.FBH_t, obj.OD_c/2, obj.rho_c, obj.NZ_t);
            obj.Center_of_Mass();
            obj.Column_Ratios(obj.OD_a, obj.ID_a, obj.L_a, obj.E_a, obj.y_a);
            obj.AS = obj.Max_dynamic_F/obj.CSA;
            if obj.Rs < obj.Rt
                obj.Johnson_Buckling(obj.y_a, obj.L_a, obj.E_a, obj.OD_a/2, obj.ID_a/2);
                obj.SF = obj.EBS/obj.AS;
            elseif obj.Rt < obj.Rs
                obj.Euler_Buckling(obj.L_a, obj.E_a, obj.OD_a/2, obj.ID_a/2);
                obj.SF = obj.EBS/obj.AS;
            end
            
            
            
        end
        
        function Inner_Diameter_c(obj, OD, Y, MEOP)
            % This function takes the outer diameter of the casing, yield 
            % strength, and the maximum expected operating pressure of the
            % propellant/motor. It will output the inner diameter of the 
            % casing. It should be noted that, for the formula
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
            
            HS_max = Y/2;
            obj.ID_c = (2 * HS_max * OD - MEOP * OD)/(MEOP + 2 * HS_max);

        end
        
        function Prop_Length(obj, ID, PropMass, PropDensity)
            X_sec_A =  pi*((ID/2)^2 - 0.25^2);
            PropVol = PropMass/PropDensity;
            obj.L_c = PropVol/X_sec_A;
        end
        
        function Bulkhead_Thickness(obj, MEOP, R, Y)
            FBCS = Y/2;
            obj.FBH_t = sqrt((3 * MEOP * R^2)/(4 * FBCS));
        end
        
        function Length(obj, IT_L, PL, BH_L, NZ_L)
            obj.L_a = IT_L + PL + BH_L + NZ_L;
        end
        
        function Column_Ratios(obj, OD, ID, L, E, Y)
            % This function uses outer radius, inner radius, Young's
            % modulus, yield strength, and length to calculate the
            % slenderness ratio and the transition ratio for the airframe
            % to determine the buckling stress case.
            % Inputs:
            % Property          Variable Name           Units
            % Outer diameter    OD                      in
            % Inner diameter    ID                      in
            % Length            L                       in
            % Young's Modulus   E                       ksi
            % Yield Strength    Y                       ksi
            % 
            % Outputs(assingment):
            % Property          Variable Name           Units
            % Slenderness ratio Rs                      unitless
            % Transition ratio  Rg                      unitless
            
            OR = OD/2;
            IR = ID/2;
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
        
        function [mass] = Airframe_Mass(obj, OR, IR, L, density)
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
            mass = Volume * density;  % Mass, lbs
            obj.CSA = X_sec_A;
            
        end
        
        function Bulkhead_Masses(obj, t, r, rho, NZ_t)
            X_sec_A = pi*r^2;
            Volume = t * X_sec_A;
            obj.FBH_MASS = Volume * rho;
            %obj.NZ_MASS = rho * X_sec_A * NZ_t;
            obj.NZ_MASS = rho * Volume;
        end
        
        function Center_of_Mass(obj)
            obj.CoM = obj.L_a/2;
        end
        
    end
end

