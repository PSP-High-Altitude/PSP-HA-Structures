
clc;
point_design_ID = 10509101316;
data = read_Data("test") * 2.20462;
data2 = read_Data("test2") * 2.20462;
stg2_PM = data(1);
stg1_PM = data2(1);
Net_Force = 1.5 * ((118.5014/2.20462)* 347.8159 + 3037.4)/(4.448*1000);  % 1.5x (M_rocket * Max_acc + Max_F_drag), kilopound-force
AF2 = Airframe(4.5, 1.2, stg2_PM, 0.0686, 14500, 72.5, 0.3, 0.162, Net_Force);  % 2nd stage airframe obj
NC2 = Nosecone(2, 0, AF2.t, AF2.t, 22.5, 0, AF2.ID, 0, 0.162, 14500, 0.3, 72.5, Net_Force);  % 2nd stage nosecone obj
FN2 = Fins(15, 5, 5, .125, 4, .162, 8);  % 2nd stage fins obj

AF1 = Airframe1(4.5, [0.162,14500,0.3,72.5], [0.0813,14489,0.286,106.7],2.205, 0.0686, stg1_PM, Net_Force);
IS1 = Interstage(AF2.OD, AF1.OD_a, 6, 0.125, 0.162);
FN1 = Fins(15, 5, 5, .125, 4, 0.0813, 8);
RK = Rocket(NC2, AF2, FN2, IS1, AF1, FN1, 0, 0, 2.204, data, data2);


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

%AF1.MoI to check values
%FN2.CoM_x

