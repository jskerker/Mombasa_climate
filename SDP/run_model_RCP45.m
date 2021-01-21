function act = run_model(x)

addpath('/Users/jenniferskerker/Documents/GradSchool/Research/Project1/Code/Models/Mombasa_climate/SDP');
addpath('/Users/jenniferskerker/Documents/GradSchool/Research/Project1/Code/Models/Mombasa_climate/Streamflow_yield_analysis');
x
x(1) = 60 + (x(1)-1)*5;
x(2) = 90 + (x(2)-1)*5;
run('sdp_climate_RCP45_opttest.m');

% looked at different ways of measuring the optimal dam size- currently
% thinking we could measure it based on the action chosen- the number of
% times the flex dam is chosen

%medCostFlex = median(sum(totalCostTime(:,:,3), 2));
%meanCostFlex = mean(sum(totalCostTime(:,:,3), 2));
%perc75 = prctile(sum(totalCostTime(:,:,3), 2),75);
%perc90 = prctile(sum(totalCostTime(:,:,3), 2),90);

act = 1-sum(action(:,1,end)==3)/10000;
disp(act)