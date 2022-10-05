
data = read_Data("test");       % Mass hist of the 2nd stage prop mass, Kg
data2 = read_Data("test2");     % Mass hist of the 1st stage prop mass, Kg
data3 = read_Data("test3");     % Various inputs(see below)
stg2_PM = data(1);              % On-pad 2nd stage prop mass, Kg
stg1_PM = data2(1);             % On-pad 1st stage prop mass, Kg
stg2_diameter = data3(1)/39.37; % 2nd stage diameter, in -> m
stg1_diameter = data3(2)/39.37; % 1st stage diameter, in -> m
stg2_MEOP = data3(3)/145.038;   % 2nd stage max chamber pressure, psi -> Mpa
stg1_MEOP = data3(4)/145.038;   % 1st stage max chamber pressure, psi -> Mpa
max_a = data3(5);         % Max acceleration, m/s^2
max_drag = data3(6);            % Max drag force, N
prop_mass_density = 1890;       % Prop mass desnsity, Kg/m^3

mass = 25.301;                  % Mass of the 2nd stage rocket
%mass = 26.1723;
Max_C_Force = (mass*max_a + max_drag)/1000000;  % Maximum compressive Force

AF2 = Airframe(stg2_diameter, stg2_MEOP, stg2_PM, prop_mass_density, [110316.12, 786, 0.37, 4456.4647], Max_C_Force);
%AF2 = Airframe(stg2_diameter, stg2_MEOP, stg2_PM, prop_mass_density, [75842.33, 372.317, 0.343, 2795.6704], Max_C_Force);
NC2 = Nosecone(2, 0, AF2.t, AF2.t, 0, AF2.ID, [110316.12, 786, 0.37, 4456.4647], 0, 0);
FN2 = Fins(0.381, 0.127, 0.127, 0.003175, 4, 4456.4647, 0.2032);

IS1 = Interstage(AF2.OD, stg1_diameter, 0.0762, 0.00127, [66603, 239.9, 0.325, 2693.25473]);
AF1 = Airframe1(stg1_diameter,[110316.12,786,0.37,4456.4647],[21029.01,206.843,0.315,1970],stg1_MEOP,prop_mass_density,stg1_PM, Max_C_Force);
FN1 = Fins(0.381, 0.127, 0.127, 0.003175, 4, 4456.4647, 0.2032);
RK = Rocket(NC2, AF2, FN2, IS1, AF1, FN1, 0, 0, 0, data, data2);

CoM_hist = [];
Mass_hist = [];
steps = size(RK.CoM_Hist_1st_Stage);
for i = 1:steps(2)
    moment_i = RK.CoM_Hist_2nd_Stage(1)*RK.MASS_Hist_2nd_Stage(1) + (RK.L_2nd_Stage+RK.CoM_Hist_1st_Stage(i))*RK.MASS_Hist_1st_Stage(i);
    mass_i = RK.MASS_Hist_2nd_Stage(1) + RK.MASS_Hist_1st_Stage(i);
    CoM_i = moment_i/mass_i;
    Mass_hist = [Mass_hist, mass_i];
    CoM_hist = [CoM_hist, CoM_i];    
end

vector1 = RK.MASS_Hist_2nd_Stage';
vector2 = RK.CoM_Hist_2nd_Stage';
SS = [vector1, vector2];
writeOutputs("Second_Stage", SS);

vector3 = RK.MASS_Hist_1st_Stage';
vector4 = RK.CoM_Hist_1st_Stage';
FS = [vector3, vector4];
writeOutputs("First_Stage", FS);

vector5 = Mass_hist';
vector6 = CoM_hist';
FR = [vector5, vector6];
writeOutputs("Rocket", FR);

RK.MASS_Hist_2nd_Stage(end)
RK.MASS_Hist_1st_Stage(end)
RK.MASS_Hist_2nd_Stage(end)+RK.MASS_Hist_1st_Stage(end) - 18.39896054

AF2.JBS/AF2.AS
AF2.LBS/AF2.AS
AF2.y/AF2.AS
