%%Boost-Glide Model Baseline Main Script
% ASEN 2804
% Author: John Mah
% Modifications by: Preston Tee
% Date: Initiated - 15 Dec. 2023
% Last modified - 15 Jan. 2024

%% Clean Workspace and Housekeeping
clear; clc; close all;

% removes warnings for table variable names for a cleaner output
warning('OFF', 'MATLAB:table:ModifiedAndSavedVarnames')

%% Import and Read Aircraft Design File
Design_Input = readtable("Design Input File.xlsx",'Sheet','Input','ReadRowNames',true); %Read in Aircraft Geometry File
Count = 9; % Change This

%height(Design_Input); %Number of different aircraft configurations in geometry file

% Import Airfoil Data File
Airfoil = readtable("Design Input File.xlsx",'Sheet','Airfoil_Data'); %Read in Airfoil Data

% Import Benchmark Aircraft Truth Data
Benchmark = readtable("Design Input File.xlsx",'Sheet','Benchmark_Truth'); %Read in Benchmark "Truth" Data for model validation only

% Import Material Properties Data
Material_Data = readtable("Design Input File.xlsx",'Sheet','Materials'); %Read in prototyp material densities for weight model

%% Caluations - Conditions and Sizing
% US Standard Atmophere - uses provided MATLAB File Exchange function
[rho,a,T,P,nu,z]= atmos(Design_Input.altitude_o(:,:)); 

ATMOS = table(rho,a,T,P,nu,z); % Reorganize atmopheric conditions into a table for ease of passing into functions
clearvars rho a T P nu z % Clear original variables now that they are in a table

% Call Wing Geometry Calcuation Function
WingGeo_Data = WingGeo(Design_Input,Count); %Calculate specific wing geometry from wing configuration parameters

%% Quick Explainer - Tables
% This code heavily utilizes tables for data organization. You should think
% of a table as a spreadsheet. Tables are also very similar to a 2D array
% of data exept that the columns can be named so that it is clear what data
% is in those columns. There are multiple ways to get the data out of
% columns, through standard indexing, and through dot indexing. 
%
% Standard indexing:
%
% Like when indexing into an array you can get data out of a table by using
% parenthasis, (), which will return another table, which is often a 
% problem for calulations and plotting
%   Example:
%       NewTable = OriginalTable(:,:)
%
% Alternativly, if you index in the same way but with curly braces, {}, a
% standard array will be returned
%   Example:
%       NewArray = OriginalTable{:,:}
%
% Dot indexing:
%
% Finally, if you would like to access just one column, tables support dot
% indexing using the name of the column header. This takes the form of the 
% name of the variable, then a dot, then the name of the column, This will
% return a standard 1D array of data
%   Example: 
%       NewArray = OriginalTable.ColumnName_1
%
% The tables in this code have been purposly organized such that the rows
% will ALWAYS correspond to the different configuration inputs in the input
% file. In other words, the first input row of the input file (1st row not
% including the header row) will match up with row 1 of the tables in this
% code, the second input will be the 2nd row of tables here, etc. Columns
% will always be variables of interest and will be named appropriately.
%
% More MATLAB documentation on getting data out of tables:
% https://www.mathworks.com/help/matlab/matlab_prog/access-data-in-a-table.html
%
% We have provided the necessary code for packaging the data into tables to
% output from and input to functions in order to keep the size of the
% function headers reasonable. It will be your responsibility to unpack and
% use data passed into functions in tables correctly. Please ensure you
% are using the preallocated varaible names and do not modify the code that
% creates the tables. We want to help you with the math, not with general
% coding

%% Calculations - Lift and Drag
% Call Wing Lift & Drag Model Function
[WingLiftModel,AoA,AoA_Count,AirfoilLiftCurve,WingLiftCurve,WingDragCurve] =...
    WingLiftDrag(Design_Input,Airfoil,Count); 

% Call Parasite Drag Buildup Model Function
[Parasite_Drag_Data,FF_Table] = ...
    ParasiteDrag(Design_Input,Airfoil,WingGeo_Data,ATMOS,Count);

% Call Induced Drag Model Function
InducedDrag_Data = ...
    InducedDrag(Design_Input,WingLiftModel,WingLiftCurve,WingDragCurve,WingGeo_Data,Count);

% Call Complete Drag Polar Function
[DragPolar_mod1,DragPolar_mod2,DragPolar_mod3] = ...
    DragPolar(Parasite_Drag_Data,InducedDrag_Data,Design_Input,AoA_Count,WingLiftCurve,Count);

% Call L/D Analysis Function
[LD_mod1,LD_mod2,LD_mod3,LD_benchmark] = ...
    LD(Benchmark,DragPolar_mod1,DragPolar_mod2,DragPolar_mod3,WingLiftCurve,AoA_Count,Count);

% Call Weight Model
[Weight_Data,CG_Data] = ...
    Weight(Design_Input,Count,WingGeo_Data,Airfoil,Material_Data);

%% Calculations - Dynamic Models
% Call Thrust Model
[ThrustCurves, Time, statsStruct] = Thrust();
% Call Boost-Ascent Flight Dynamics Model
[apogee, hApogee, stateStruct] = BoostAscent(Design_Input, ATMOS, Parasite_Drag_Data, Weight_Data, ThrustCurves, Time, Count);
% Call Glide Flight Dynamics Model
[GlideRange] = GlideDescent(LD_mod1, apogee, Design_Input, ATMOS, Weight_Data, WingLiftModel, WingLiftCurve, Count);


%% Plotting

%Lift Curves
figure
hold on
plot(AoA,Airfoil{1,(5:22)});
plot(Benchmark.AoA(:),Benchmark.CL(:));
plot(AoA,AirfoilLiftCurve{1,:});
plot(AoA,WingLiftCurve{1,:});
xlabel('Angle of Attack (deg)');
ylabel('Coefficient of Lift (CL)');
title('Lift Curve Slope Modeling Comparison');
legend('Airfoil Data','Benchmark Aircraft','Airfoil Model','Wing Model','Location','southeast');
hold off

%Drag Polar Curves
figure
hold on
plot(AirfoilLiftCurve{1,:},Airfoil{1,(24:41)});
plot(WingLiftCurve{1,:},WingDragCurve{1,:});
plot(Benchmark.CL(:),Benchmark.CD(:));
plot(WingLiftCurve{1,:},DragPolar_mod1{1,:});
plot(WingLiftCurve{1,:},DragPolar_mod2{1,:});
plot(WingLiftCurve{1,:},DragPolar_mod3{1,:});
xlabel('Coefficient of Lift (CL)');
ylabel('Coefficient of Drag (CD)');
title('Drag Polar Model Comparison');
legend('Airfoil Drag Polar','3D Wing Drag Polar','Benchmark Drag Polar','Drag Polar Cavallo','Drag Polar Obert', 'Drag Polar Kroo', 'Location','northwest');
hold off

figure
hold on
for n = 1:Count
    plot(WingLiftCurve{1,:},DragPolar_mod1{n,:}); % brace indexing for plotting tables
end
xlabel('Coefficient of Lift (CL)');
ylabel('Coefficient of Drag (CD)');
title('Drag Polar Configuration Comparison - Cavallo Oswald Model 1');
legend(Design_Input.Config(:),'Location','northwest');
hold off

%Lift over Drag Analysis Plots (Configuration 1)
figure
hold on
plot(WingLiftCurve{1,:},LD_mod1{1,:},'--');
plot(WingLiftCurve{1,:},LD_mod2{1,:},'--');
plot(WingLiftCurve{1,:},LD_mod3{1,:},'--');
plot(WingLiftCurve{1,:},LD_benchmark{1,:});
%plot(WingLiftCurve{1,:},(WingLiftCurve{1,:}./WingDragCurve{1,:}));

xlabel('Coefficient of Lift (CL)');
ylabel('L/D Ratio (-)');
title('Lift over Drag Comparisons - Configuration 1.01');
legend('L/D Cavallo','L/D Obert','L/D Kroo','L/D Benchmark', 'L/D  3D Wing', 'Location','southeast');
hold off
%% New Graph
figure
hold on
induced1=(InducedDrag_Data.k1_mod1(n) * (WingLiftCurve{n,:}).^2);
induced2=(InducedDrag_Data.k1_mod2(n) * (WingLiftCurve{n,:}).^2);
induced3=(InducedDrag_Data.k1_mod3(n) * (WingLiftCurve{n,:}).^2);

Benchmark_matrix=table2array(Benchmark);
CD_benchmark_matrix = Benchmark_matrix(:, 3);
CD_3dcd_matrix=transpose(table2array(WingDragCurve));
Parasite_Matrix = 2*.01323215.*ones([18 1]);
induced_benchmark = transpose((CD_benchmark_matrix - Parasite_Matrix));
Parasite_Matrix_w = 2*.0052711.*ones([18 1]);
induced_benchmark_w = transpose((CD_3dcd_matrix - Parasite_Matrix_w));

plot(WingLiftCurve{1,:},induced1(1,:),'--');
plot(WingLiftCurve{1,:},induced2(1,:),'--');
plot(WingLiftCurve{1,:},induced3(1,:),'--');
WingLiftCurve=table2array(WingLiftCurve);
plot(WingLiftCurve(1,:),induced_benchmark(1,:));
plot(WingLiftCurve(1,:),induced_benchmark_w(1,:));

xlabel('Coefficient of Lift (C_L)');
ylabel('Induced Drag (C_D_i)');
title('Induced Drag vs. Coefficient of Lift');
legend('Cavallo','Obert','Kroo','Benchmark', '3D Wing', 'Location','southeast');
hold off

%% Ryan Testing
e_benchmark.CL_B = Benchmark_matrix(:,2);
e_benchmark.CD_B = Benchmark_matrix(:,3);
e_benchmark.CL_zero = interp1(e_benchmark.CL_B,e_benchmark.CD_B,0);
e_benchmark.CDo = e_benchmark.CL_zero;
e_benchmark.AR = Design_Input.AR_w(n);

e_benchmark.e_bench = zeros(size(e_benchmark.CL_B));
for i = 1:length(e_benchmark.CL_B)
e_benchmark.e_bench(i) = (e_benchmark.CL_B(i)^2)/(pi*e_benchmark.AR*(e_benchmark.CD_B(i)-e_benchmark.CDo));
end

e_benchmark.k_B = 1/(pi*e_benchmark.e_bench(2)*Design_Input.AR_w(n));

figure()
plot(e_benchmark.CD_B,e_benchmark.CL_B)

% Calculate error for each angle of attack
e_benchmark.test1 = abs(((induced1 - induced_benchmark)./induced_benchmark).*100);
e_benchmark.test2 = abs(((induced2 - induced_benchmark)./induced_benchmark).*100);
e_benchmark.test3 = abs(((induced3 - induced_benchmark)./induced_benchmark).*100);

% Calculate error for CL = 0
e_benchmark.value = interp1(e_benchmark.CL_B,AoA,0);
e_benchmark.B1 = interp1(AoA,induced1,e_benchmark.value);
e_benchmark.B2 = interp1(AoA,induced2,e_benchmark.value);
e_benchmark.B3 = interp1(AoA,induced3,e_benchmark.value);

%Calculate oswalds for AoA at L=0
e_benchmark.I = interp1(AoA,induced_benchmark,e_benchmark.value);
e_benchmark.ZL1 = abs(((e_benchmark.B1 - e_benchmark.I)./e_benchmark.I).*100);
e_benchmark.ZL2 = abs(((e_benchmark.B2 - e_benchmark.I)./e_benchmark.I).*100);
e_benchmark.ZL3 = abs(((e_benchmark.B3 - e_benchmark.I)./e_benchmark.I).*100);

figure
hold on
plot(AoA,induced_benchmark)
plot(AoA,induced1)
plot(AoA,induced2)
plot(AoA,induced3)
legend('benchmark','i1','i2','i3');
hold off

%% Static Test Thrust Profle Model Plots
num2Ltests = 8; % Just hard code the number of 2L tests in the ThrustCurves table
cmap = colormap(parula(num2Ltests));
set(0,'DefaultAxesColorOrder',cmap)
set(gca(),'ColorOrder',cmap);
figure(10)
for i = 1:num2Ltests % Looping over just the 2L data in ThrustCurves
    plot(Time, ThrustCurves{:, i}, DisplayName = [ThrustCurves.Properties.VariableNames{i}(6:end), ' mL'])
    if i == 1
        hold on
    end
end
xlabel('Test Time [s]');
ylabel('Fitted Thrust [N]');
title('2L Fitted Thrust Profiles');
legend();
grid on
hold off

cmap = colormap(parula(width(ThrustCurves)-num2Ltests));
set(0,'DefaultAxesColorOrder',cmap)
set(gca(),'ColorOrder',cmap);
figure(110)
for i = num2Ltests+1:width(ThrustCurves) % Looping over the rest of the data in ThrustCurves
    plot(Time, ThrustCurves{:, i}, DisplayName = [ThrustCurves.Properties.VariableNames{i}(6:end), ' mL'])
    if i == num2Ltests+1
        hold on
    end
end
xlabel('Test Time [s]');
ylabel('Fitted Thrust [N]');
title('1.25L Fitted Thrust Profiles');
legend();
grid on
hold off

%% Boost_Ascent Flight Profile Plots
% Some setup to make plots more readable in color, look up the
% documentation for 'cmap' for other color map options
figure(20)
cmap = colormap(parula(Count));
set(0,'DefaultAxesColorOrder',cmap)
set(gca(),'ColorOrder',cmap);

fields = fieldnames(stateStruct);
for n = 1:Count
    distBoost = vecnorm([stateStruct.(fields{n}).data(:, 4), stateStruct.(fields{n}).data(:, 5)], 2, 2);
    plot(distBoost,...
            -stateStruct.(fields{n}).data(:, 6), ...
            DisplayName=Design_Input.Properties.RowNames{n}, Color=cmap(n, :))
    if n == 1
        hold on
    end
end
xlabel('Total Distance Traveled [m]');
ylabel('Height Achieved [m]');
title('Boost 2D Total Distance Traveled');
legend();
grid on
hold off

figure(21)
for n = 1:Count
    plot(stateStruct.(fields{n}).data(:, 4),...
            stateStruct.(fields{n}).data(:, 5), ...
            DisplayName=Design_Input.Properties.RowNames{n}, Color=cmap(n, :))
    if n == 1
        hold on
    end
end
xlabel('y [m] - Positive = East');
ylabel('x [m] - Positive = North');
title('Boost Ground Track');
legend();
grid on
hold off

figure(22)
for n = 1:Count    
    plot3(stateStruct.(fields{n}).data(:, 4),...
            stateStruct.(fields{n}).data(:, 5),...
            stateStruct.(fields{n}).data(:, 6), ...
            DisplayName=Design_Input.Properties.RowNames{n})
    if n == 1
        hold on
    end
end
xlabel('x [m] - Positive = North');
ylabel('y [m] - Positive = East');
zlabel('z [m] - Positive = Down');
title('Boost Trajectory Plots');
set(gca, 'ZDir','reverse')
set(gca, 'YDir','reverse')
legend();
grid on
axis equal
hold off

figure(23)
for n = 1:Count
    Wx = -Design_Input.V_wind(n)*cosd(Design_Input.Wind_Az(n)); 
    Wy = -Design_Input.V_wind(n)*sind(Design_Input.Wind_Az(n));
    quiver3(stateStruct.(fields{n}).data(:, 4), ... % x
        stateStruct.(fields{n}).data(:, 5), ... % y
        stateStruct.(fields{n}).data(:, 6), ... % z
        stateStruct.(fields{n}).data(:, 1)-Wx, ... % Vax
        stateStruct.(fields{n}).data(:, 2)-Wy, ... % Vay
        stateStruct.(fields{n}).data(:, 3), ... % Vz
        DisplayName=Design_Input.Properties.RowNames{n}, ...
        LineWidth=2)
    if n == 1
        hold on
    end
end
xlabel('x [m] - Positive = North');
ylabel('y [m] - Positive = East');
zlabel('z [m] - Positive = Down');
title('Boost Trajectory Plots with Heading Vetors');
set(gca, 'ZDir','reverse')
set(gca, 'YDir','reverse')
legend();
grid on
axis equal
hold off

%% Reset default color order
set(0,'DefaultAxesColorOrder','default')
