
clc;
point_design_ID = 010407111317;
data = read_Data("test") * 2.20462;
data2 = read_Data("test2") * 2.20462;
data3 = read_Data("test3");
stg2_PM = data(1);
stg1_PM = data2(1);
stg2_diameter = data3(1);
stg1_diameter = data3(2);
stg2_MEOP = data3(3) * 1.5/1000;
stg1_MEOP = data3(4) * 1.5/1000;
max_a = data3(5);
max_drag = data3(6);
payload = data3(7) * 2.20462;
Net_Force = 1.5 * ((198/2.20462)* max_a + max_drag)/(4.448*1000);  % 1.5x (M_rocket * Max_acc + Max_F_drag), kilopound-force
AF2 = Airframe(stg2_diameter, stg2_MEOP, stg2_PM, 0.0686, 14500, 72.5, 0.3, 0.162, Net_Force);  % 2nd stage airframe obj
NC2 = Nosecone(2, 0, AF2.t, AF2.t, 5*AF2.OD, 0, AF2.ID, 0, 0.162, 14500, 0.3, 72.5, Net_Force);  % 2nd stage nosecone obj
FN2 = Fins(15, 5, 5, .125, 4, .162, 8);  % 2nd stage fins obj

AF1 = Airframe1(stg1_diameter, [0.162,14500,0.3,72.5], [0.0813,14489,0.286,106.7],stg1_MEOP, 0.0686, stg1_PM, Net_Force);
IS1 = Interstage(AF2.OD, AF1.OD_a, 6, 0.125, 0.0813);
FN1 = Fins(15, 5, 5, .125, 4, 0.0813, 8);
RK = Rocket(NC2, AF2, FN2, IS1, AF1, FN1, 0, 0, 11.02, data, data2);


disp("POINT DESIGN ID#: " + point_design_ID);
%% AIRFRAME
disp("-----------------2ND STAGE AIRFRAME STATS-----------------")
disp("PHYSICAL PROPERTIES:")
disp("Thickness: " + AF2.t + " in");
disp("Foward Bulkhead minimum thickness: " + AF2.FBH_t + " in");
disp("Length: " + AF2.L + " in");
disp("Mass: " + AF2.MASS + " lbs");
%disp("Forward Bulkhead minimum mass: " + AF2.FBH_MASS + " lbs");
disp("STRUCTURAL STABILITY:")
disp("Local Buckling SF: " + AF2.LBS/AF2.AS);
if AF2.CT == 0
    disp("Johnson Buckling(Intermediate) SF: " + AF2.JBS/AF2.AS);
elseif AF2.CT == 1
    disp("Euler Buckling(Long) SF: " + AF2.EBS/AF2.AS);
end


%% NOSECONE
disp("-----------------2ND STAGE NOSECONE STATS-----------------");
disp("PHYSICAL PROPERTIES:");
disp("Mass: " + NC2.MASS + " lbs");
disp("STRUCTURAL STABILITY:");
disp("Critical Axial Load: " + NC2.NCCL + " lbs");
disp("Buckling SF: " + NC2.E/NC2.NCBS);

%% FINS
disp("-------------------2ND STAGE FINS STATS-------------------")
disp("PHYSICAL PROPERTIES:");
disp("Mass: " + FN2.MASS + " lbs");

%% 2ND STAGE ROCKET
disp("---------------------2ND STAGE ROCKET---------------------")
disp("Rocket Mass: " + RK.MASS_2nd_Stage + " lbs");
disp("Center of Mass: " + RK.CoM_2nd_Stage + " in");

%% ROCKET
disp("-----------------------FULL ROCKET------------------------")
mass = RK.MASS_2nd_Stage + RK.MASS_1st_Stage;
moment = RK.MASS_2nd_Stage*RK.CoM_2nd_Stage+RK.MASS_1st_Stage*(RK.L_2nd_Stage+RK.CoM_1st_Stage);
CoM = moment/mass;
length = RK.L_2nd_Stage + RK.L_1st_Stage;
disp("Length: " + length + " in")
disp("Mass: " + mass + " lbs")
disp("Center of Mass: " + CoM + " in from NC tip")

CoM_hist = [];
steps = size(RK.CoM_Hist_1st_Stage);
for i = 1:steps(2)
    moment_i = RK.CoM_2nd_Stage*RK.MASS_2nd_Stage + (RK.L_2nd_Stage+RK.CoM_Hist_1st_Stage(i))*RK.MASS_Hist_1st_Stage(i);
    mass_i = RK.MASS_2nd_Stage + RK.MASS_Hist_1st_Stage(i);
    CoM_i = moment_i/mass_i;
    CoM_hist = [CoM_hist, CoM_i];    
end

% MoIx = RK.Parallax(CoM_hist(1), RK.MASS_1st_Stage, RK.MoIx_1st_Stage, RK.CoM_1st_Stage);
% MoIx = MoIx + RK.Parallax(CoM_hist(1), RK.MASS_2nd_Stage, RK.MoIx_2nd_Stage, RK.CoM_2nd_Stage);
% RK.Parallax(CoM_hist(1), RK.MASS_2nd_Stage, RK.MoIx_2nd_Stage, RK.CoM_2nd_Stage);
% 
% RK.MoIx_1st_Stage
% RK.MoIx_2nd_Stage

%AF1.MoI to check values
%FN2.CoM_x

% disp(" ")
% disp(" ")
% disp("---------2nd Stage NC----------")
% disp("MoIx: " + NC2.MoIx_nc)
% disp("MoIy: " + NC2.MoIy_nc)
% disp("MoIz: " + NC2.MoIz_nc)
% disp("---------2nd Stage AF----------")
% disp("MoIx: " + AF2.MoIx_af)
% disp("MoIy: " + AF2.MoIy_af)
% disp("MoIz: " + AF2.MoIz_af)
% disp("---------2nd Stage PM----------")
% disp("MoIx: " + AF2.MoIx_prop)
% disp("MoIy: " + AF2.MoIy_prop)
% disp("MoIz: " + AF2.MoIz_prop)
% disp("---------1st Stage IS----------")
% disp("MoIx: " + IS1.MoIx_IS)
% disp("MoIy: " + IS1.MoIy_IS)
% disp("MoIz: " + IS1.MoIz_IS)
% disp("---------1st Stage AF----------")
% disp("MoIx: " + AF1.MoIx_a)
% disp("MoIy: " + AF1.MoIy_a)
% disp("MoIz: " + AF1.MoIz_a)
% disp("---------1st Stage CA----------")
% disp("MoIx: " + AF1.MoIx_c)
% disp("MoIy: " + AF1.MoIy_c)
% disp("MoIz: " + AF1.MoIz_c)
% disp("---------1st Stage PM----------")
% disp("MoIx: " + AF1.MoIx_prop)
% disp("MoIy: " + AF1.MoIy_prop)
% disp("MoIz: " + AF1.MoIz_prop)
% 
% 
% 
% 
% 
