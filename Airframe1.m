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
        E;
        y;
        v;
        rho;
        
        
        % PROPULSION/OPERATING PROPERTIES
        MEOP;
        PD;
        PM;
        Max_F;
        
        % DESIGN CRITERIA
        MASS;
        FBH_MASS;
        NZ_MASS;
    end
    
    methods
        function obj = Airframe1(input1,input2, input3, input4, input5, input6)
            %AIRFRAME1 Construct an instance of this class
            %   Detailed explanation goes here
            obj.OD = input1;
            
            obj.E = input2(1);
            obj.y = input2(2);
            obj.v = input2(3);
            obj.rho = input2(4);
            
            obj.MEOP = input3;
            obj.PD = input4;
            obj.PM = input5;
            obj.Max_F = input6;
            
            obj.NZ_t = 0.127;
            
            %obj.Inner_Diameter_m(obj.OD, obj.y_m, obj.MEOP);
            obj.t_m = 0.003175;
            obj.OR = obj.OD/2;
            obj.ID_m = obj.OD - obj.t_m*2;
            obj.IR_m = obj.ID_m/2;
            %obj.t_m = obj.OR-obj.ID_m;
%             if obj.t_m < 0.00127
%                 obj.t_m = 0.00127;
%                 obj.ID_m = obj.OD - 2*obj.t_m;
%                 obj.OR = obj.OD/2;
%                 obj.IR_m = obj.ID_m/2;
%             end
            
            obj.Prop_Length(obj.ID_m, obj.PM, obj.PD);
            obj.Bulkhead_Thickness(obj.MEOP, obj.IR_m, obj.y);
            obj.Length(obj.PL, obj.FBH_th, obj.NZ_t);
            obj.MASS = obj.Airframe_Mass(obj.OR, obj.IR_m, obj.PL+obj.NZ_t+obj.FBH_th, obj.rho);
            obj.Bulkhead_Masses(obj.FBH_th, obj.IR_m, obj.rho, obj.NZ_t);
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
        
        function Length(obj, PL, BH_L, NZ_L)
            obj.L = PL + BH_L + NZ_L;
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

