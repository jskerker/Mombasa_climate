%% January 2021
% This file includes multiple tests for a wrapper function to run different
% static or flex dam sizing options, or using a genetic algorithm to pick
% dam sizes

%% Runs the scp_climate or scp_climate_rcp45 function with different static storage sizes
clear all;
close all;
storage = [60 100];

while storage(1) < 105
    addpath('/Users/jenniferskerker/Documents/GradSchool/Research/Project1/Code/Models/Fletcher_2019_Learning_Climate/SDP');
    addpath('/Users/jenniferskerker/Documents/GradSchool/Research/Project1/Code/Models/Fletcher_2019_Learning_Climate/Streamflow_yield_analysis');
    run('sdp_climate.m');
    storage = storage + 5;
end

%% run sizing for all possible combinations for RCP8.5 
clear all;
close all;
%x = [70 90];
addpath('/Users/jenniferskerker/Documents/GradSchool/Research/Project1/Code/Models/Fletcher_2019_Learning_Climate/SDP');
addpath('/Users/jenniferskerker/Documents/GradSchool/Research/Project1/Code/Models/Fletcher_2019_Learning_Climate/Streamflow_yield_analysis');
addpath('/Users/jenniferskerker/Documents/GradSchool/Research/Project1/Code/Models/Fletcher_2019_Learning_Climate'); 
for j= 65:5:105
    x(1) = j;
    for k= 95:5:135
        x(2) = k;
        run('sdp_climate_opttest.m'); % this is a version of the sdp_climate.m 
        %file where I modified the storage variables to take in x, and I
        %commented out the display message line in the sdp step to make the
        %model run raster
        disp(x);
        meanCostFlex = mean(sum(totalCostTime(:,:,3), 2))
    end
end
%% run sizing for all possible combinations for RCP4.5 
clear all;
close all;
x = [70 90];
addpath('/Users/jenniferskerker/Documents/GradSchool/Research/Project1/Code/Models/Mombasa_climate/SDP');
addpath('/Users/jenniferskerker/Documents/GradSchool/Research/Project1/Code/Models/Mombasa_climate/Streamflow_yield_analysis');
addpath('/Users/jenniferskerker/Documents/GradSchool/Research/Project1/Code/Models/Mombasa_climate'); 
for j= 65:5:105
    x(1) = j;
    for k= 95:5:135
        x(2) = k;
        run('sdp_climate_RCP45_opttest.m');
        disp(x);
        meanCostFlex = mean(sum(totalCostTime(:,:,3), 2))
    end
end


%% create function to optimize flex sizing using ga
% NOTE: This piece is not up to date, but is essentially the same as the
% updated piece of code below.
clear all;
close all;

nvars = 2;
ObjectiveFunction = @run_model;
%ConstraintFunction = @simple_constraint;

% int means that the input variables 1 and 2 are integers
Int = [1 2];
% lower bound, upper bound, and starting values
LB2 = [1 1];
UB2 = [9 10];
X02 = [5 7];

%options.InitialPopulationMatrix = X0;
options = optimoptions('ga', 'Display', 'iter', ...
    'InitialPopulationMatrix', X02, 'MutationFcn', @mutationadaptfeasible);

[x, fval] = ga(ObjectiveFunction, nvars,[],[],[],[],LB2,UB2, [],Int, options);

%% create function to optimize flex sizing- RCP4.5

clear all;
close all;

nvars = 2;
addpath('/Users/jenniferskerker/Documents/GradSchool/Research/Project1/Code/Models/Mombasa_climate/SDP');
ObjectiveFunction = @run_model_RCP45;
%ConstraintFunction = @simple_constraint;

% int means that the input variables 1 and 2 must be integers
Int = [1 2];

% lower bound, upper bound, and starting values
% in the run_model_RCP45 script, I multiply the values to get the actual
% dam sizes I want, but having the input variables be low integer numbers
% was easier
LB2 = [1 1]; % numbers correspond to counting values of the above
UB2 = [9 10];
X02 = [5 7];

%options.InitialPopulationMatrix = X0;
options = optimoptions('ga', 'InitialPopulationMatrix', X02,...
    'Display', 'iter', 'PlotFcn', {@gaplotbestf, @gaplotstopping}, ...
    'MaxGenerations', 100);
options.Display = 'iter';
[x, fval, exitFlag, Output] = ga(ObjectiveFunction, nvars,[],[],[],[],LB2,UB2,[], Int, options);
