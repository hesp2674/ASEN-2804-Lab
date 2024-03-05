% clear
% clc
% close all
% 
% [ThrustCurves, Time,raw_data,time_relevant,milestone1_data] = Thrustm();
% end_point = size(ThrustCurves);
% thrust = table2array(ThrustCurves);
% names = ThrustCurves.Properties.VariableNames;
%   t = Time';
% 
% for i = 1:end_point(2)
%     string_from_cell = char(ThrustCurves.Properties.VariableNames{i});
%     name_splitter(i,:)=strsplit(string_from_cell,'_');
%     milestone1_data(i,1:2) = [str2double(name_splitter(i,1)),str2double(name_splitter(i,2))];
%     name_actual = strcat(char(name_splitter(i,1)), ' mL Bottle-' , char(name_splitter(i,2)) , ' mL Water');
%     figure(i)
%     plotThrust = (thrust(:,i));
%     plot(t,plotThrust,LineWidth=2)
%     title(name_actual)
%     hold on
%     plot(time_relevant,raw_data(:,i))
%     legend('Fourier Approximation','Raw Data')
%     xlabel("Time(s)");
%     ylabel("Thrust(lbs)");
% end


function [ThrustCurves, Time,Thrust_Raw_Data,time_relevant,milestone1_data] = Thrust()
%% Thrust Summary
% This funciton will take in the file location of the two test setups and
% using the file names in those directories, will pull out all of the
% avialable tests, cendition their data, and fit that data into a standard
% formatting for output. Note that statistics are also requested for the
% student deliverable, but how to pass those out will be left up to the
% students as they are not needed to be passed into any later functions.
% Despite this, the first two outputs of the funciton are not permitted to
% have their form modified.

%% Outputs:
% ThrustCurves:
%   A table containing 0.5 seconds of thrust data for each of the cases
%   available, this data will have formatting such that there are 501
%   evenly spaced thrust data points (rows) for each test (columns). The
%   ordering of the columns will go from max to min water volume in the 2L
%   bottle and then max to min in the 1.25L bottle
%
% Time:
%   A 1D array corresponding to the times of the thrust data points in the
%   ThrustCurves table
%
% <User defined variable(s) for statistics>
%

%% Define data locations
% This is hard coded!!!
fileLoc_2L = 'Thrust_Test_Data/2000mL 60 psi/'; % path to the data files, be sure to include a trailing slash
fileLoc_1pt25L = 'Thrust_Test_Data/1250mL 60 psi/'; % path to the data files, be sure to include a trailing slash

%% Read in all of the avilable data and find what data there is
testInfo_2L = getThrustTestNames(fileLoc_2L);
    configs_2L = unique(testInfo_2L.waterVol);
numConfigs_2L = length(configs_2L);

testInfo_1pt25L = getThrustTestNames(fileLoc_1pt25L);
configs_1pt25L = unique(testInfo_1pt25L.waterVol);
numConfigs_1pt25L = length(configs_1pt25L);

numConfigs = numConfigs_2L + numConfigs_1pt25L;

% Set known sampling frequency
f= 1652; % [Hz]

%% Preallocate variables of interest
Time = 0:0.001:0.5; % just go ahead and define this, note that it will be 501 long
ThrustCurves = zeros(length(Time),numConfigs);
Thrust_Raw_Data = zeros(f/2+1,numConfigs);

ThrustCurvesNames = {};

%% Loop over all of the configurations
for N = 1:numConfigs % use upper case N to distiguish that it is counting something different from the aerodynamic modeling loops
    %% Dertemine what configuration to use for this iteration in the loop
    if N <=  numConfigs_2L % determine if we should be reading 2L or 1.25L data
        bottleSize = '2000'; % [ml]
        waterSize = configs_2L(N);
        testIndexes = find(testInfo_2L.waterVol == waterSize); % finds the index of the relavant tests
        numTests = length(testIndexes); % finds the number of tests performed
        testNames = testInfo_2L.fileNames(testIndexes, :); % pulls all of the test names of interest, weird indexing is due to string arrays
    else
        bottleSize = '1250'; % [ml]
        waterSize = configs_1pt25L(N-numConfigs_2L);
        testIndexes = find(testInfo_1pt25L.waterVol == waterSize); % finds the index of the relavant tests
        numTests = length(testIndexes); % finds the number of tests performed
        testNames = testInfo_1pt25L.fileNames(testIndexes, :); % pulls all of the test names of interest, weird indexing is due to string arrays
    end

    % /////////////////////////////////////////////////////////////////////////
    % MODIFY THIS SECTION
    % /////////////////////////////////////////////////////////////////////////
    % Notice that there is little to no guidance in place for this
    % function. This is on purpose as there are many different and equally
    % valid ways to process data (not to say that any way is valid though).
    % The lack of guidance is therefore to encourage you to think about,
    % discuss, and potentially debate as a group the best set of steps to
    % extract just the meaningful part of the thrust profile
    data_relevant = zeros(f/2+1,numTests);
    for i = 1:numTests
    %% Load data
        % The folloowing three lines will pull all of the files in each
        % test setup for you and give the array "data" which should be
        % conditioned. You should not need to modify any of this section of
        % code
        fileName = testNames(i, :); % again weird indexing is due to string arrays, we have to ask for all the characters in a row
        data = readmatrix(fileName); % load the data
        data = data(:,3)*4.448; % take only the third column and converting from lbf to N

    %% Data Conditioning
    %start point (manual)/ end point (start point +.5s)
    %get correct time stamp and match with num samples
    time_actual = 1/f:1/f:(length(data)/f);
    time_actual = time_actual';
    found =0;
    found2= 0;


    for j = 1:(length(data)-100)
        if found ==0 && data(j) < (mean(data(j+1:j+5)-10))
            start_index =j-10;
            end_index = start_index+1652/2;
        end
        if data(j) < (mean(data(j+1:j+30)-50)) && found == 0 &&data(j) > -15
            start_index =j-10;
            end_index = start_index+1652/2;
            found =1;
        end
        if found ==1 && found2 ==0 && data(j)<= mean(data(j+1:j+30)) +.05 && data(j)>= mean(data(j+1:j+30)) -.05
            flat_index(i) = j-start_index;
            found2 = 1;
        end
    end
    time_relevant  = time_actual(start_index:end_index) - time_actual(start_index);
    data = data';
    data_relevant(:,i)  = data(start_index:end_index)';

    
    data_relevant(:,i) = data_relevant(:,i) - mean(data_relevant(1:5,i));
    if mean(data_relevant(:,i)) < 0
        data_relevant(:,i)  = -data_relevant(:,i);
    end
    time_flat(:,i) = time_relevant(flat_index(i));
    
    %Adjustment calculations: first declare zeroes
    adjustment = zeros(end_index-start_index+1,1);
    %find slope between 0 and thrust value when it flattens
    avg_weight_after = mean(data_relevant(700:827,i));
    adj_slope = (avg_weight_after)/(flat_index(i));
    %create scalar of indeces for adj_slope (it is in thrust/index)
    adjuster = 1:1:flat_index(i);
    %find linear part of adjustment
    adjustment(1:flat_index(i)) = adj_slope.*adjuster;
    %for second half of values just use the final value
    adjustment(flat_index(i)+1:827)= avg_weight_after;
    %subtract the adjustment
    data_relevant(:,i) = data_relevant(:,i) - adjustment;

    % %create line of best fit
    % data_fit_help = polyfit(time_relevant,data_relevant,8);
    % 
    % data_fit = polyval(data_fit_help,Time);
        
    %% Averaging
    %adjust for water weight offset





    end



    % Note that averaging should accour before data fitting. Technically
    % either can be done, but the output will be much more smooth if the
    % fit is applied to an average
    %% Data Fitting
    final_data_fit = mean(data_relevant,2);
    peak_thrusts = max(data_relevant,[],1);
    mean_max = mean(peak_thrusts);
    std_max = std(peak_thrusts);
    time_mean = mean(time_flat);
    std_time = std(time_flat);
    for i = 1:length(peak_thrusts)
        avg_thrust(i) = mean(data_relevant(10:flat_index(i), i) );
        time_thrust(i) = time_relevant(flat_index(i),1)- time_relevant(10,1);
    end
    impulse = mean(avg_thrust.*time_thrust);

    Fs = 1652;
    L = 827;
    xaxis = Fs/L*(0:L-1);


    freq = fft(final_data_fit);
    freq_norm = 1/L * freq;
    freq_recon = freq;
    k = 30;
    for i = (k+2):(L-k)
        freq_recon(i) = 0;
    end

    transformed_data_recon = ifft(freq_recon);
    transformed_data_orig = ifft(freq);

    for i = length(transformed_data_recon)-170:length(transformed_data_recon)
        transformed_data_recon(i) = 0;
    end
    interped = interp1(time_relevant, transformed_data_recon,Time);
    q = interped <0;
    interped(q) = 0;
    
    if(interped(6)<20)
        w = find(q == 1, 1);
        interped(1:w) = 0;
        interped(1:10) =0;
    end

    thrustOut = interped;
    


    %% Sample onto the standard output array format
    milestone1_data(N,:) =[0 0 mean_max std_max time_mean std_time impulse];




% /////////////////////////////////////////////////////////////////////////
% END OF SECTION TO MODIFY
% /////////////////////////////////////////////////////////////////////////
    %% Convert to table for output
    % It is very important that the data is 501 elements long corresponding
    % to 0-0.5 seconds of time at this point!!!
    ThrustCurves(:, N) = thrustOut;
    Thrust_Raw_Data(:,N)= final_data_fit;
    % Header naming convention of <bottle size (in ml)>_<water volume (in ml)>
    ThrustCurvesNames{N} = [bottleSize, '_', num2str(waterSize)];
end
ThrustCurves = array2table(ThrustCurves);
ThrustCurves.Properties.VariableNames = ThrustCurvesNames;
end
