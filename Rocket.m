classdef Rocket < handle
    %ROCKET Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        % INPUTS
        NC2;                    % Second stage NC object
        AF2;                    % Second stage AF object
        FN2;                    % Second stage FN object
        IS1;                     % Interstage object
        AF1;                    % First stage AF object
        FN1;                    % First stage FN object
        Max_dynamic_F;          % Drag force at Max_Q, N
        Max_Q;                  % Max dynamic pressure, Mpa
        PL_MASS;                % Payload Mass, kg
        PM;
        PM_NEW;
        PM_Hist_2nd_Stage;      % Propellant mass history for 2nd stage, kg(t)
        PM_Hist_1st_Stage;      % Propellant mass history for 2nd stage, kg(t)
        
        % DESIGN CRITERIA
        
        L_2nd_Stage;            % Length of the 2nd stage, m
        L_1st_Stage;            % Length of the 1st stage, m
        
        MASS_2nd_Stage;         % Pad mass of the 2nd stage, kg
        MASS_Hist_2nd_Stage;    % Mass history of the 2nd stage, kg(t)
        MASS_1st_Stage;         % Pad mass of the 1st stage, kg
        MASS_Hist_1st_Stage;    % Pad mass of the 1st stage, kg(t)
        Rocket_CoM_Hist;        % CoM hist of the rocket, m(t)
        
        CoM_2nd_Stage;          % Pad CoM of the 2nd stage, m
        CoM_Hist_2nd_Stage;     % CoM hist of the 2nd stage, m(t)
        CoM_1st_Stage;          % Pad CoM of the 1st stage, m
        CoM_Hist_1st_Stage;     % CoM hist of the 1st stage, m(t)

        MoIx_2nd_Stage;         
        MoIx_Hist_2nd_Stage;
        MoIx_1st_Stage;
        MoIx_Hist_1st_Stage;

        MoIz_2nd_Stage;
        MoIz_Hist_2nd_Stage;
        MoIz_1st_Stage;
        MoIz_Hist_1st_Stage;
    end
    
    methods
        function obj = Rocket(input1,input2,input3,input4,input5,input6,input7,input8,input9,input10,input11)
            %ROCKET Construct an instance of this class
            %   Detailed explanation goes here
            obj.NC2 = input1;
            obj.AF2 = input2;
            obj.FN2 = input3;
            obj.IS1 = input4;
            obj.AF1 = input5;
            obj.FN1 = input6;
            obj.Max_dynamic_F = input7;
            obj.Max_Q = input8;
            obj.PL_MASS = input9;
            obj.PM_Hist_2nd_Stage = input10;
            obj.PM_Hist_1st_Stage = input11;
            
            obj.Second_Stage_CoM_History(obj.PM_Hist_2nd_Stage);
            obj.CoM_1st_Stage = obj.First_Stage_CoM(obj.PM_Hist_1st_Stage(1));
            obj.First_Stage_CoM_History(obj.PM_Hist_1st_Stage);
        end
        
        function [CoM, Mass_out] = Second_Stage_CoM(obj, PM)
            NC = [obj.NC2.MASS, obj.NC2.CoM];
            AF = [obj.AF2.MASS, obj.NC2.L+obj.AF2.CoM];
            RC = [1.814, obj.NC2.L + 0.6096/2];
            SP = [0.453, obj.NC2.L + 0.6096 + 0.0381];
            AV = [0.800, obj.NC2.L + 0.6096 + 0.0762 + 0.127];
            FBH = [obj.AF2.FBH_MASS, obj.NC2.L+obj.AF2.L-obj.AF2.NZ_t-obj.AF2.PL-obj.AF2.FBH_t/2];
            PROP = [PM, obj.NC2.L+obj.AF2.L-obj.AF2.NZ_t-obj.AF2.PL/2];
            LBH = [obj.AF2.FBH_MASS, obj.NC2.L+obj.AF2.L-obj.AF2.NZ_t+obj.AF2.FBH_t/2];
            FN = [obj.FN2.MASS, obj.NC2.L+obj.AF2.L-obj.FN2.CoM_y];
            
            Rocket = [NC; AF; RC; SP; AV; FBH; PROP; LBH; FN];
            
            Moment = 0;
            Mass = 0;
            for i = 1:size(Rocket, 1)
                Moment = Moment + Rocket(i,1)*Rocket(i,2);
                Mass = Mass + Rocket(i,1);
            end
            Mass_out = Mass;
            CoM = Moment/Mass;
        end
        
        function Second_Stage_CoM_History(obj, PM_Hist)
            CoM_hist = [];
            Mass_hist = [];
            steps = length(PM_Hist);
            for i = 1:steps
                [CoM_i, Mass_i] = obj.Second_Stage_CoM(PM_Hist(i));
                Mass_hist = [Mass_hist, Mass_i];
                CoM_hist = [CoM_hist, CoM_i];
            end
            
            obj.L_2nd_Stage = obj.NC2.L+obj.NC2.tip_L+obj.AF2.L;
            obj.MASS_Hist_2nd_Stage = Mass_hist;
            obj.CoM_Hist_2nd_Stage = CoM_hist;
        end
        
        function [CoM, Mass_out] = First_Stage_CoM(obj, PM)
            IS = [obj.IS1.MASS, obj.IS1.CoM];
            FBH = [obj.AF1.FBH_MASS, obj.IS1.L+obj.AF1.L-obj.AF1.NZ_t-obj.AF1.PL-obj.AF1.FBH_th/2];
            AF = [obj.AF1.MASS, obj.IS1.L+obj.AF1.L/2];
            PROP = [PM, obj.IS1.L+obj.AF1.L-obj.AF1.NZ_t-obj.AF1.PL/2];
            LBH = [obj.AF1.FBH_MASS, obj.IS1.L+obj.AF1.L-obj.AF1.NZ_t+obj.AF1.FBH_th/2];
            FN = [obj.FN1.MASS, obj.IS1.L+obj.AF1.L-obj.FN1.CoM_y];
            
            Rocket = [IS; FBH; AF; PROP; LBH; FN];
            
            Moment = 0;
            Mass = 0;
            for i = 1:size(Rocket, 1)
                Moment = Moment + Rocket(i,1)*Rocket(i,2);
                Mass = Mass + Rocket(i,1);
            end
            Mass_out = Mass;
            CoM = Moment/Mass;
        end
        
        function First_Stage_CoM_History(obj, PM_Hist)
            CoM_hist = [];
            Mass_hist = [];
            steps = length(PM_Hist);
            for i = 1:steps
                [CoM_i, Mass_i] = obj.First_Stage_CoM(PM_Hist(i));
                Mass_hist = [Mass_hist, Mass_i];
                CoM_hist = [CoM_hist, CoM_i];
            end
            
            obj.L_1st_Stage = obj.IS1.L+obj.AF1.L;
            obj.MASS_Hist_1st_Stage = Mass_hist;
            obj.CoM_Hist_1st_Stage = CoM_hist;
        end
        
    end
end

