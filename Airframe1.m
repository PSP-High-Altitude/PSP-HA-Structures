classdef Airframe1 < handle
    %AIRFRAME1 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        % GEOMETRIC PROPERTIES
        %inputs
        OD;
        OR;
        %constants
        IT_L;
        NZ_t;
        %calculated
        ID_m;
        IR_m;
        t_m;
        ID_t;
        IR_t;
        t_t;
        L;
        
        CSA;
        
        
        
        PL;
        FBH_th;
        
        % MOTOR CASING MATERIAL PROPERTIES
        E_m;
        y_m;
        v_m;
        rho_m;
        
        % AIRFRAME TUBE MATERIAL PROPERTIES
        E_t;
        y_t;
        v_t;
        rho_t;
        
        
        % PROPULSION/OPERATING PROPERTIES
        MEOP;
        PD;
        PM;
        Max_F;
        
        % DESIGN CRITERIA
        MASS_m;
        MASS_t;
        MASS;
        FBH_MASS;
        NZ_MASS;
    end
    
    methods
        function obj = Airframe1(input1,input2, input3, input4, input5, input6, input7)
            %AIRFRAME1 Construct an instance of this class
            %   Detailed explanation goes here
            obj.OD = input1;
            
            obj.E_m = input2(1);
            obj.y_m = input2(2);
            obj.v_m = input2(3);
            obj.rho_m = input2(4);
            
            obj.E_t = input3(1);
            obj.y_t = input3(2);
            obj.v_t = input3(3);
            obj.rho_t = input3(4);
            
            obj.MEOP = input4;
            obj.PD = input5;
            obj.PM = input6;
            obj.Max_F = input7;
            
            obj.IT_L = 0.9906;
            obj.NZ_t = 0.127;
            
            obj.Inner_Diameter_m(obj.OD, obj.y_m, obj.MEOP);
            obj.OR = obj.OD/2;
            obj.IR_m = obj.ID_m/2;
            obj.t_m = obj.OR-obj.ID_m;
            if obj.t_m < 0.00127
                obj.t_m = 0.00127;
                obj.ID_m = obj.OD - 2*obj.t_m;
                obj.OR = obj.OD/2;
                obj.IR_m = obj.ID_m/2;
            end
            obj.ID_t = obj.ID_m;
            obj.IR_t = obj.ID_t/2;
            obj.t_t = 0.003175;
            
            
            obj.Prop_Length(obj.ID_m, obj.PM, obj.PD);
            obj.Bulkhead_Thickness(obj.MEOP, obj.IR_m, obj.y_m);
            obj.Length(obj.IT_L, obj.PL, obj.FBH_th, obj.NZ_t);
            obj.MASS_m = obj.Airframe_Mass(obj.OR, obj.IR_m, obj.PL+obj.NZ_t+obj.FBH_th, obj.rho_m);
            obj.MASS_t = obj.Airframe_Mass(obj.OR, obj.IR_t, obj.IT_L, obj.rho_t);
            obj.Bulkhead_Masses(obj.FBH_th, obj.IR_m, obj.rho_m, obj.NZ_t);
        end
        
        function Inner_Diameter_m(obj, OD, Y, MEOP)
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
            % Outer Diameter    OD                      m
            % Yield strength    Y                       Mpa
            % Maximum expected  MEOP                    Mpa
            %  - operating pressure
            % 
            % Outputs(assingment):
            % Property          Variable Name           Units
            % Inner diameter    ID                      m
            
            HS_max = Y/2;
            obj.ID_m = (2 * HS_max * OD - MEOP * OD)/(MEOP + 2 * HS_max);

        end
        
        function Prop_Length(obj, ID, PropMass, PropDensity)
            X_sec_A =  pi*((ID/2)^2 - 0.00635^2);
            PropVol = PropMass/PropDensity;
            obj.PL = PropVol/X_sec_A;
        end
        
        function Bulkhead_Thickness(obj, MEOP, R, Y)
            FBCS = Y/2;
            obj.FBH_th = sqrt((3 * MEOP * R^2)/(4 * FBCS));
        end
        
        function Length(obj, IT_L, PL, BH_L, NZ_L)
            obj.L = IT_L + PL + BH_L + NZ_L;
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
        
    end
end

