classdef Interstage < handle
    %INTERSTAGE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        OD1;
        OD2;
        L;
        t;
        rho;
        PROFILE;
        MASS;
        CoM;
        MoIx_IS;
        MoIy_IS;
        MoIz_IS;
    end
    
    methods
        function obj = Interstage(input1,input2, input3, input4, input5)
            %INTERSTAGE Construct an instance of this class
            %   Detailed explanation goes here
            obj.OD1 = input1;
            obj.OD2 = input2;
            obj.L = input3;
            obj.t = input4;
            obj.rho = input5;
            
            obj.Interstage_Mass(obj.OD1, obj.OD2, obj.L, obj.t, obj.rho);
            obj.Center_of_Mass(0, obj.L, 0, obj.PROFILE, obj.L);
            obj.MoI_IS(obj.OD1, obj.L, obj.t, obj.MASS);
        end
        
        function Interstage_Mass(obj, D1, D2, L, t, density)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            R1 = D1/2;
            R2 = D2/2;
            m = (R2-R1)/L;
            y_int = R1;
            Outer_Profile = @(x) (m*x + y_int);
            Inner_Profile = @(x) (m*x + y_int - t);
            Integrand = @(x) pi.*(Outer_Profile(x).^2 - Inner_Profile(x).^2);
            Volume = integral(Integrand, 0, L);
            obj.PROFILE = Outer_Profile;
            obj.MASS = Volume * density;
            
        end
        
        function Center_of_Mass(obj, xmin, xmax, ymin, ymax, length)
            Moment_fun = @(x, y) x;
            Moment = integral2(Moment_fun, xmin, xmax, ymin, ymax);
            Mass = integral(ymax, 0, length);
            obj.CoM = Moment/Mass;
            
        end
       
        function MoI_IS(obj, OD1, L, t, MASS)
            c = 0.9; %Constant for assuming that the interstage is a cylinder
            OR = OD1 / 2;
            IR = OD1 / 2 - t;
            
            moix = c *1/12 * MASS * (3*(OR^2+IR^2) + L^2);
            moiy = c * 1/12 * MASS * (3*(OR^2+IR^2) + L^2);
            moiz = c * 1/2 * MASS * (OR^2+IR^2);
        
            obj.MoIx_IS = moix;
            obj.MoIy_IS = moiy;
            obj.MoIz_IS = moiz;
            
        end
        
    end
    
end

