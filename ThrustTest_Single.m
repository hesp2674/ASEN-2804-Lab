clear; clc; close all;

%% ThrustTest_Single Summary
% First, some of the summary from Thrust:
% This funciton will take in the file location of the two test setups and
% using the file names in those directories, will pull out all of the
% avialable tests, process them, and fit the data into a standard
% formatting for output
%
% ThrustTest_Single is meant to be a testing function where you can test
% the conditioning steps that your group develops. This version of the code
% loads in only one data set at a time so that the workspace is less
% cluttered. It is suggested that students make plots of their thrust data
% often throughout their development in this code section to visually check
% that their conditioning is working as desired. Once conditioning is
% working on one set of data, students are engcouraged to try other single
% data sets. Once groups are satisfied, the full funciton has the same form
% as this, so the conditioning code can simply be dragged and dropped into
% the full "Thrust.m" funciton.
% 
% Note that there is one major step missing from this testing funciton.
% That is averaging over multiple tests from the same setup. This will need
% to be added in the full "Thrust.m" function


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
% <User defined variables for statistics>

% Set known sampling frequency
f= 1652; % [Hz]

%% Preallocate variables of interest
Time = 0:0.001:0.5; % just go ahead and define this, note that it will be 501 long
ThrustCurves = zeros(length(Time),1);

%% List what configuration is being used
bottleSize = 1250; % [ml]
waterSize = 400;% [ml]
% Be sure that the above matches up with the file name specified below. In
% the full "Thrust" script, all of this data will be loaded automatically
% for you
testName = 'Thrust_Test_Data/2000mL 60 psi/Group10Test05_W0600_B2000'; % Should be a string of the path to the data (including the data file name)

% /////////////////////////////////////////////////////////////////////////
% MODIFY THIS SECTION
% /////////////////////////////////////////////////////////////////////////

%% Load data
% This should not have to be modified
fileName = testName; % again weird indexing is due to string arrays, we have to ask for all the characters in a row
data = readmatrix(fileName); % load the data
data = data(:,3)*4.448; % take only the third column and converting from lbf to N


%% Data Conditioning
i =1;
N = 1;
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




% freq_one_side = freq(1:length(freq)/2);
% freq_norm_one_side = 1/L * freq_one_side;
% xaxis2 = Fs/(L/2)*(0:(L/2)-1);
% plot(xaxis2,abs(freq_one_side))
% 
% figure();
% stem(xaxis2,abs(freq_norm_one_side),'x','LineWidth',0.2)
% xlabel('frequency')
% ylabel('amplitude')
% grid on




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

    for i = 1:length(peak_thrusts)
        avg_thrust2(i) = mean(transformed_data_recon(10:flat_index(i), i) );
        time_thrust2(i) = time_relevant(flat_index(i),1)- time_relevant(10,1);
    end
    impulse2 = mean(avg_thrust2.*time_thrust2);
    interped = interp1(time_relevant, transformed_data_recon,Time);
     
    thrustOut = interped;

    figure()

    plot(Time,interped)
    hold on
   % plot(time_relevant,data_relevant)


    %% Sample onto the standard output array format
    milestone1_data(N,:) =[0 0 mean_max std_max time_mean std_time impulse];


% hold on
% plot(time_relevant,data_relevant)
% plot(Time,data_fit);
% scatter(time_flat,data(flat_index), LineWidth=10)
% hold off

%% Sample onto the standard output array format

% /////////////////////////////////////////////////////////////////////////
% END OF SECTION TO MODIFY
% /////////////////////////////////////////////////////////////////////////


