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
        Max_dynamic_F;          % Drag force at Max_Q, lbf
        Max_Q;                  % Max dynamic pressure, ksi
        PL_MASS;                % Payload Mass, lbs
        PM;
        PM_NEW;
        PM_Hist_2nd_Stage;      % Propellant mass history for 2nd stage
        PM_Hist_1st_Stage;      % Propellant mass history for 2nd stage
        
        % DESIGN CRITERIA
        
        L_2nd_Stage;
        L_1st_Stage;
        
        MASS_2nd_Stage;
        MASS_Hist_2nd_Stage;
        MASS_1st_Stage;
        MASS_Hist_1st_Stage;
        Rocket_CoM_Hist;
        
        CoM_2nd_Stage;
        CoM_Hist_2nd_Stage;
        CoM_1st_Stage;
        CoM_Hist_1st_Stage;

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
        function obj = Rocket(input1,input2, input3, input4, input5, input6, input7, input8, input9, input10, input11)
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
            
            
            obj.Second_Stage_CoM();
            obj.Second_Stage_CoM_History(obj.PM_Hist_2nd_Stage);
            obj.First_Stage_CoM();
            obj.First_Stage_CoM_History(obj.PM_Hist_1st_Stage);
            %obj.Rocket_CoM_Hist(obj.CoM_Hist_1st_Stage, obj.MASS_Hist_1st_Stage);
            
            obj.Second_Stage_MoI();
            %obj.Second_Stage_MoI_History(obj.PM_Hist_2nd_Stage);
            obj.First_Stage_MoI();
            %obj.First_Stage_MoI_History(obj.PM_Hist_1st_Stage);
            
        end
        
        function Second_Stage_CoM(obj)
            NC = [obj.NC2.MASS, obj.NC2.CoM];               % Nosecone
            AF = [obj.AF2.MASS, obj.NC2.L + obj.AF2.CoM];   % Airframe
            MC = [0.535, obj.NC2.L + 4.724 + 5];            % Main Chute(MC)
            MCSC = [0.048, obj.NC2.L + 14.961 + 1];         % MC Shock cord 
            BH1 = [1.355, obj.NC2.L + 17 + 1];              % First Bulkhead
            PL = [obj.PL_MASS, obj.NC2.L + 20 + 7];         % Payload
            BH2 = [1.355, obj.NC2.L + 35 + 1];              % Second bulkhead
            DCSC = [0.013, obj.NC2.L + 37 + 1];             % DC shock cord
            DC = [0.104, obj.NC2.L + 47 + 5];               % Drogue Chute(DC)
            FBH = [obj.AF2.FBH_MASS, obj.NC2.L+(obj.AF2.L-obj.AF2.NZ_t-obj.AF2.PL-obj.AF2.FBH_t/2)];    % Forward Bulkhead
            PROP = [obj.AF2.PM, obj.NC2.L + obj.AF2.L - obj.AF2.NZ_t - obj.AF2.PL/2];   % Propellant Mass
            %NZ = [obj.AF2.NZ_MASS, obj.NC2.L + obj.AF2.L - obj.AF2.NZ_t/2];             % Nozzle
            NZ = [obj.AF2.NZ_MASS, obj.NC2.L+obj.AF2.L-obj.AF2.NZ_t-obj.AF2.FBH_t/2];     % Nozzle
            FN = [obj.FN2.MASS, obj.NC2.L + obj.AF2.L - obj.FN2.RC + obj.FN2.CoM_y];       % Fins
           
            Rocket = [NC; AF; MC; MCSC; BH1; PL; BH2; DCSC; DC; FBH; PROP; NZ; FN];
            %Rocket = [NC; AF; MC; MCSC; BH1; PL; BH2; DCSC; DC; FBH; NZ; FN];
            
            Moment = 0;
            Mass = 0;
            for i = 1:size(Rocket, 1)
                Moment = Moment + Rocket(i,1)*Rocket(i,2);
                Mass = Mass + Rocket(i,1);
            end
            
            obj.MASS_2nd_Stage = Mass;
            obj.CoM_2nd_Stage = Moment/Mass;
        end
        
        function [CoM, Mass] = Second_Stage_CoM_Hist_Version(obj, PM)
            NC = [obj.NC2.MASS, obj.NC2.CoM];               % Nosecone
            AF = [obj.AF2.MASS, obj.NC2.L + obj.AF2.CoM];   % Airframe
            MC = [0.535, obj.NC2.L + 4.724 + 5];            % Main Chute(MC)
            MCSC = [0.048, obj.NC2.L + 14.961 + 1];         % MC Shock cord 
            BH1 = [1.355, obj.NC2.L + 17 + 1];              % First Bulkhead
            PL = [obj.PL_MASS, obj.NC2.L + 20 + 7];         % Payload
            BH2 = [1.355, obj.NC2.L + 35 + 1];              % Second bulkhead
            DCSC = [0.013, obj.NC2.L + 37 + 1];             % DC shock cord
            DC = [0.104, obj.NC2.L + 47 + 5];               % Drogue Chute(DC)
            FBH = [obj.AF2.FBH_MASS, obj.NC2.L+(obj.AF2.L-obj.AF2.NZ_t-obj.AF2.PL-obj.AF2.FBH_t/2)];    % Forward Bulkhead
            PROP = [PM, obj.NC2.L + obj.AF2.L - obj.AF2.NZ_t - obj.AF2.PL/2];   % Propellant Mass
            %NZ = [obj.AF2.NZ_MASS, obj.NC2.L + obj.AF2.L - obj.AF2.NZ_t/2];             % Nozzle
            NZ = [obj.AF2.NZ_MASS, obj.NC2.L+obj.AF2.L-obj.AF2.NZ_t-obj.AF2.FBH_t/2];     % Nozzle
            FN = [obj.FN2.MASS, obj.NC2.L + obj.AF2.L - obj.FN2.RC + obj.FN2.CoM_y];       % Fins
           
            Rocket = [NC; AF; MC; MCSC; BH1; PL; BH2; DCSC; DC; FBH; PROP; NZ; FN];
            %Rocket = [NC; AF; MC; MCSC; BH1; PL; BH2; DCSC; DC; FBH; NZ; FN];
            
            Moment = 0;
            Mass = 0;
            for i = 1:size(Rocket, 1)
                Moment = Moment + Rocket(i,1)*Rocket(i,2);
                Mass = Mass + Rocket(i,1);
            end
            
            
            CoM = Moment/Mass;
        end
        
        function Second_Stage_CoM_History(obj, PM_Hist)
            CoM_hist = [];
            Mass_hist = [];
            steps = length(PM_Hist);
            for i = 1:steps
                [CoM_i, Mass_i] = obj.Second_Stage_CoM_Hist_Version(PM_Hist(i));
                Mass_hist = [Mass_hist, Mass_i];
                CoM_hist = [CoM_hist, CoM_i];
            end
            
            obj.MASS_Hist_2nd_Stage = Mass_hist;
            obj.CoM_Hist_2nd_Stage = CoM_hist;
        end
        
        function First_Stage_CoM(obj)
            IS = [obj.IS1.MASS, obj.IS1.CoM];               % Interstage
            AF = [obj.AF1.MASS_a, obj.IS1.L+obj.AF1.CoM];   % Airframe(aero)
            CA = [obj.AF1.MASS_c, obj.IS1.L+obj.AF1.L_a-obj.AF1.NZ_t-obj.AF1.L_c/2];    % Motor Casing
            FBH = [obj.AF1.FBH_MASS, obj.IS1.L+obj.AF1.L_a-obj.AF1.NZ_t-obj.AF1.L_c-obj.AF1.FBH_t/2];   % Forward Bulkhead
            NZ = [obj.AF1.NZ_MASS, obj.IS1.L+obj.AF1.L_a-obj.AF1.NZ_t-obj.AF1.FBH_t/2];   % Nozzle
            MC = [0.535, obj.IS1.L + 4.724 + 5];            % Main Chute(MC)
            MCSC = [0.048, obj.IS1.L + 14.961 + 1];         % MC Shock cord 
            BH1 = [1.355, obj.IS1.L + 17 + 1];              % First Bulkhead
            PL = [obj.PL_MASS, obj.IS1.L + 20 + 7];         % Payload
            BH2 = [1.355, obj.IS1.L + 35 + 1];              % Second bulkhead
            DCSC = [0.013, obj.IS1.L + 37 + 1];             % DC shock cord
            DC = [0.104, obj.IS1.L + 47 + 5];               % Drogue Chute(DC)
            PROP = [obj.AF1.PM, CA(1,2)];                   % Propellant mass
            FN = [obj.FN1.MASS, obj.IS1.L+obj.AF1.L_a-obj.FN1.RC+obj.FN1.CoM_y];
            
            Rocket = [IS; AF; CA; FBH; NZ; MC; MCSC; BH1; PL; BH2; DCSC; DC; PROP; FN];
            Moment = 0;
            Mass = 0;
            for i = 1:size(Rocket, 1)
                Moment = Moment + Rocket(i,1)*Rocket(i,2);
                Mass = Mass + Rocket(i,1);
            end
            
            obj.L_2nd_Stage = obj.NC2.L+obj.NC2.tip_L+obj.AF2.L;
            obj.L_1st_Stage = obj.IS1.L+obj.AF1.L_a;
            obj.MASS_1st_Stage = Mass;
            obj.CoM_1st_Stage = Moment/Mass;
            
        end
        
        function [CoM, Mass] = First_Stage_CoM_Hist_Version(obj, PM)
            IS = [obj.IS1.MASS, obj.IS1.CoM];               % Interstage
            AF = [obj.AF1.MASS_a, obj.IS1.L+obj.AF1.CoM];   % Airframe(aero)
            CA = [obj.AF1.MASS_c, obj.IS1.L+obj.AF1.L_a-obj.AF1.NZ_t-obj.AF1.L_c/2];    % Motor Casing
            FBH = [obj.AF1.FBH_MASS, obj.IS1.L+obj.AF1.L_a-obj.AF1.NZ_t-obj.AF1.L_c-obj.AF1.FBH_t/2];   % Forward Bulkhead
            NZ = [obj.AF1.NZ_MASS, obj.IS1.L+obj.AF1.L_a-obj.AF1.NZ_t-obj.AF1.FBH_t/2];   % Nozzle
            MC = [0.535, obj.IS1.L + 4.724 + 5];            % Main Chute(MC)
            MCSC = [0.048, obj.IS1.L + 14.961 + 1];         % MC Shock cord 
            BH1 = [1.355, obj.IS1.L + 17 + 1];              % First Bulkhead
            PL = [obj.PL_MASS, obj.IS1.L + 20 + 7];         % Payload
            BH2 = [1.355, obj.IS1.L + 35 + 1];              % Second bulkhead
            DCSC = [0.013, obj.IS1.L + 37 + 1];             % DC shock cord
            DC = [0.104, obj.IS1.L + 47 + 5];               % Drogue Chute(DC)
            PROP = [PM, CA(1,2)];                   % Propellant mass
            FN = [obj.FN1.MASS, obj.IS1.L+obj.AF1.L_a-obj.FN1.RC+obj.FN1.CoM_y];
            
            Rocket = [IS; AF; CA; FBH; NZ; MC; MCSC; BH1; PL; BH2; DCSC; DC; PROP; FN];
            Moment = 0;
            Mass = 0;
            for i = 1:size(Rocket, 1)
                Moment = Moment + Rocket(i,1)*Rocket(i,2);
                Mass = Mass + Rocket(i,1);
            end
            
            CoM = Moment/Mass;
        end
        
        function First_Stage_CoM_History(obj, PM_Hist)
            CoM_hist = [];
            Mass_hist = [];
            steps = length(PM_Hist);
            for i = 1:steps
                [CoM_i, Mass_i] = obj.First_Stage_CoM_Hist_Version(PM_Hist(i));
                Mass_hist = [Mass_hist, Mass_i];
                CoM_hist = [CoM_hist, CoM_i];
            end
            
            obj.MASS_Hist_1st_Stage = Mass_hist;
            obj.CoM_Hist_1st_Stage = CoM_hist;
        end


        %this function can calculate the MoI of all objects using
        %parallel axis theorm. IF you don't need Parallax, then just input
        %the RokCOM as a zero instead of calling the CoM function. 
        function [MoI] = Parallax(obj, RokCoM, Mass, MoIs, CoMs)
            MoI = 0;
            %RokCoM
            for i = 1:size(MoIs)
                
                MoI = (MoI + MoIs(i)) + (Mass(i) * (RokCoM - CoMs(i))^2);
                %here = CoMs(i)
            end
            
        end

        function First_Stage_MoI(obj)
            %3D Objects
            IS = [obj.IS1.MASS, obj.IS1.MoIx_IS, obj.IS1.MoIy_IS, obj.IS1.MoIz_IS, obj.IS1.CoM];               % Interstage
            AF = [obj.AF1.MASS_a, obj.AF1.MoIx_a, obj.AF1.MoIy_a, obj.AF1.MoIz_a, obj.IS1.L+obj.AF1.CoM];   % Airframe(aero)
            CA = [obj.AF1.MASS_c, obj.AF1.MoIx_c, obj.AF1.MoIy_c, obj.AF1.MoIz_c, obj.IS1.L+obj.AF1.L_a-obj.AF1.NZ_t-obj.AF1.L_c/2];    % Motor Casing
            PROP = [obj.AF1.PM, obj.AF1.MoIx_prop, obj.AF1.MoIy_prop, obj.AF1.MoIz_prop, CA(1,5)];                   % Propellant
            
            
            
            %Point Masses
            FBH = [obj.AF1.FBH_MASS, obj.IS1.L+obj.AF1.L_a-obj.AF1.NZ_t-obj.AF1.L_c-obj.AF1.FBH_t/2];   % Forward Bulkhead
            NZ = [obj.AF1.NZ_MASS, obj.IS1.L+obj.AF1.L_a-obj.AF1.NZ_t-obj.AF1.FBH_t/2];   % Nozzle
            MC = [0.535, obj.IS1.L + 4.724 + 5];            % Main Chute(MC)
            MCSC = [0.048, obj.IS1.L + 14.961 + 1];         % MC Shock cord 
            BH1 = [1.355, obj.IS1.L + 17 + 1];              % First Bulkhead
            PL = [obj.PL_MASS, obj.IS1.L + 20 + 7];         % Payload
            BH2 = [1.355, obj.IS1.L + 35 + 1];              % Second bulkhead
            DCSC = [0.013, obj.IS1.L + 37 + 1];             % DC shock cord
            DC = [0.104, obj.IS1.L + 47 + 5];               % Drogue Chute(DC)
            FN = [obj.FN1.MASS, obj.IS1.L+obj.AF1.L_a-obj.FN1.RC+obj.FN1.CoM_y];
            
            S1_CoM = obj.CoM_1st_Stage;           %First_Stage_CoM();

            %split up the physical and point masses because they have different vector sizes.
            Rocket3D = [IS; AF; CA; PROP];
            RocketPt = [FBH; NZ; MC; MCSC; BH1; PL; BH2; DCSC; DC; FN];
            %disp("HELLO THERE");            
            MoIx3Da = obj.Parallax(S1_CoM, Rocket3D(:,1), Rocket3D(:,2), Rocket3D(:,5))
            MoIxPt = obj.Parallax(S1_CoM, RocketPt(:,1), zeros(1, 2), RocketPt(:,2));
            obj.MoIx_1st_Stage = MoIx3Da + MoIxPt;
            
                    
            %point masses are centered around the centerline still, need to
            %write a function that generates cylinders for them if time
            %permits. 
            obj.MoIz_1st_Stage = obj.Parallax(0, zeros(1, 4), Rocket3D(:,4), Rocket3D(:,5));
        end
        
        function Second_Stage_MoI(obj)
            %3D Objects
            NC = [obj.NC2.MASS, obj.NC2.MoIx_nc, obj.NC2.MoIy_nc, obj.NC2.MoIz_nc, obj.NC2.CoM];               % Nosecone
            AF = [obj.AF2.MASS, obj.AF2.MoIx_af, obj.AF2.MoIy_af, obj.AF2.MoIz_af, obj.NC2.L + obj.AF2.CoM];   % Airframe
            PROP = [obj.AF2.PM, obj.AF2.MoIx_prop, obj.AF2.MoIy_prop, obj.AF2.MoIz_prop, obj.NC2.L + obj.AF2.L - obj.AF2.NZ_t - obj.AF2.PL/2];   %Propellant
            
            %Point Masses
            MC = [0.535, obj.NC2.L + 4.724 + 5];            % Main Chute(MC)
            MCSC = [0.048, obj.NC2.L + 14.961 + 1];         % MC Shock cord 
            BH1 = [1.355, obj.NC2.L + 17 + 1];              % First Bulkhead
            PL = [obj.PL_MASS, obj.NC2.L + 20 + 7];         % Payload
            BH2 = [1.355, obj.NC2.L + 35 + 1];              % Second bulkhead
            DCSC = [0.013, obj.NC2.L + 37 + 1];             % DC shock cord
            DC = [0.104, obj.NC2.L + 47 + 5];               % Drogue Chute(DC)
            FBH = [obj.AF2.FBH_MASS, obj.NC2.L+(obj.AF2.L-obj.AF2.NZ_t-obj.AF2.PL-obj.AF2.FBH_t/2)];    % Forward Bulkhead
            %NZ = [obj.AF2.NZ_MASS, obj.NC2.L + obj.AF2.L - obj.AF2.NZ_t/2];             % Nozzle
            NZ = [obj.AF2.NZ_MASS, obj.NC2.L+obj.AF2.L-obj.AF2.NZ_t-obj.AF2.FBH_t/2];     % Nozzle
            FN = [obj.FN2.MASS, obj.NC2.L + obj.AF2.L - obj.FN2.RC + obj.FN2.CoM_y];       % Fins
            
            S2_CoM = obj.CoM_2nd_Stage;

            %split up the physical and point masses because they have different vector sizes.
            Rocket3D = [NC; AF; PROP];
            RocketPt = [FBH; NZ; MC; MCSC; BH1; PL; BH2; DCSC; DC; FN];
            
            MoIx3D = obj.Parallax(S2_CoM, Rocket3D(:,1), Rocket3D(:,2), Rocket3D(:,5)); 
            MoIxPt = obj.Parallax(S2_CoM, RocketPt(:,1), zeros(1, 2), RocketPt(:,2));
            obj.MoIx_2nd_Stage = MoIx3D + MoIxPt;
            
            
            %point masses are centered around the centerline still, need to
            %write a function that generates cylinders for them if time
            %permits. 
            obj.MoIz_2nd_Stage = obj.Parallax(0, zeros(1,4), Rocket3D(:,4), Rocket3D(:,5));
        end
        
        
%         function Rocket_CoM_History(obj, FS_CoM_Hist, FS_Mass_Hist)
%             CoM_hist = [];
%             steps = size(FS_Mass_Hist);
%             for i = 1:steps(2)
%                 moment = obj.CoM_2nd_Stage*obj.MASS_2nd_Stage + FS_CoM_Hist(i)*FS_Mass_Hist(i);
%                 mass = obj.MASS_2nd_Stage + FS_Mass_Hist(i);
%                 CoM_i = moment/mass;
%                 CoM_hist = [CoM_hist, CoM_i];
%             end
%             
%             obj.Rocket_CoM_Hist = CoM_Hist;
%         end
    end
end

