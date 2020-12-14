%Setup
clear all; close all
%addpath('/Users/jenniferskerker/Documents/GradSchool/Research/Project1/Code/Models/Mombasa_climate/BMA_code')
%load('BMA_results_RCP45_2020-11-11.mat')

% RCP8.5
addpath('/Users/jenniferskerker/Documents/GradSchool/Research/Project1/Code/Models/Fletcher_2019_Learning_Climate/BMA_code')
load('BMA_results_RCP85_2020-11-11.mat')

N = 5;

% Percent change in precip from one time period to next
climParam.P_min = -.3;
climParam.P_max = .3;
climParam.P_delta = .02; 
s_P = climParam.P_min : climParam.P_delta : climParam.P_max;
climParam.P0 = s_P(15);
climParam.P0_abs = 77; %mm/month
M_P = length(s_P);

% Change in temperature from one time period to next
climParam.T_min = 0;
climParam.T_max = 1.5;
climParam.T_delta = 0.05; % deg C
s_T = climParam.T_min: climParam.T_delta: climParam.T_max;
climParam.T0 = s_T(1);
climParam.T0_abs = 26;
M_T = length(s_T);

% Absolute temperature values
T_abs_max = max(s_T) * N;
s_T_abs = climParam.T0_abs : climParam.T_delta : climParam.T0_abs+ T_abs_max;
M_T_abs = length(s_T_abs);
T_bins = [s_T_abs-climParam.T_delta/2 s_T_abs(end)+climParam.T_delta/2];

% Absolute percip values
P_abs_max = max(s_P) * N;
s_P_abs = 66:1:97;
M_P_abs = length(s_P_abs);
P_bins = [s_P_abs-climParam.P_delta/2 s_P_abs(end)+climParam.P_delta/2];

climParam = struct;
climParam.numSamp_delta2abs = 100000;
climParam.numSampTS = 100;
climParam.checkBins = true;

% Percent change in precip from one time period to next
climParam.P_min = -.3;
climParam.P_max = .3;
climParam.P_delta = .02; 
s_P = climParam.P_min : climParam.P_delta : climParam.P_max;
climParam.P0 = s_P(15);
climParam.P0_abs = 77; %mm/month
M_P = length(s_P);

% Change in temperature from one time period to next
climParam.T_min = 0;
climParam.T_max = 1.5;
climParam.T_delta = 0.05; % deg C
s_T = climParam.T_min: climParam.T_delta: climParam.T_max;
climParam.T0 = s_T(1);
climParam.T0_abs = 26;
M_T = length(s_T);

[T_Temp, T_Precip, ~, ~, ~, ~] = bma2TransMat( NUT, NUP, s_T, s_P, N, climParam);

numSamp = 25000;
decades = { '1990', '2010', '2030', '2050', '2070', '2090'};

% Starting point
T0 = s_T(1);
T0_abs = 26;
P0 = s_P(15);
P0_abs = 75;

%% Updating over time with runoff
addpath('/Users/jenniferskerker/Documents/GradSchool/Research/Project1/Code/Models/Mombasa_climate/SDP')
load('runoff_by_state_Mar16_knnboot_1t')
load('shortage_costs_01_Dec_2020_19_31_36', 'yield', 'unmet_dom', 'shortageCost')

% Run for loop to create sample sets
numRuns = 1000;

for M=1:numRuns
    % Set time series
    state_ind_P = zeros(1,N);
    state_ind_T = zeros(1,N);
    state_ind_P(1) =  find(P0_abs==s_P_abs);
    state_ind_T(1) = find(T0_abs==s_T_abs);
    randGen = true;
    %state_ind_P(2:N) = [12 17 19 22];
    state_ind_P(2:N) = [10 11 14 17];
    %state_ind_T(2:N) = [9 17 26 33];
    state_ind_T(2:N) = [10 19 26 33];

    MAR = cellfun(@(x) mean(mean(x)), runoff);
    p = randi(numSamp,N-1);
    T_over_time = cell(1,N);
    P_over_time = cell(1,N);
    MAR_over_time = cell(1,N);

    for t = 1:N
        % Sample forward distribution given current state
        T_current = s_T_abs(state_ind_T(t));
        P_current = s_P_abs(state_ind_P(t));
        [T_over_time{t}] = T2forwardSimTemp(T_Temp, s_T_abs, N, t, T_current, numSamp, false);
        [P_over_time{t}] = T2forwardSimTemp(T_Precip, s_P_abs, N, t, P_current, numSamp, false);

        % Lookup MAR and yield for forward distribution
        T_ind = arrayfun(@(x) find(x == s_T_abs), T_over_time{t});
        P_ind = arrayfun(@(x) find(x == s_P_abs), P_over_time{t});
        [~,t_steps] = size(T_ind);
        MAR_over_time{t} = zeros(size(T_ind));
        yield_over_time{t} = zeros(size(T_ind));
        for i = 1:numSamp
            for j = 1:t_steps
                MAR_over_time{t}(i,j) = MAR(T_ind(i,j), P_ind(i,j), 1);
                yield_over_time{t}(i,j) = unmet_dom(T_ind(i,j), P_ind(i,j),1, 1) ;   % 80 MCM storage
            end
        end

        % Sample next time period
        if randGen
            state_ind_T(t+1) = find(T_over_time{t}(p(t),2)==s_T_abs);
            state_ind_P(t+1) = find(P_over_time{t}(p(t),2)==s_P_abs);
        end
    end

    % Calculate CI for first time period and t2-1 time period
    t1 = 1; %1990
    t2 = 4; %2050
    CIFcn = @(x,p)prctile(x,abs([0,100]-(100-p)/2)); %Source: https://www.mathworks.com/matlabcentral/answers/159417-how-to-calculate-the-confidence-interval
    CI1 = CIFcn(yield_over_time{t1},80);
    CI2 = CIFcn(yield_over_time{t2-1},80);
    CI3 = CIFcn(yield_over_time{end},80);
    
    % Find the percent change between the two CIs for 2050
    diff1 = CI1(2,t2)-CI1(1,t2); % looking at the differences in CIs in 2050
    diff2 = CI2(2,6-t2)-CI2(1,6-t2);
    percChange2050(M) = (diff2 - diff1)/diff1 * 100; % negative value means CI is getting smaller

    % Find the percent change between the two CIs for 2090
    diff1 = CI1(2,end)-CI1(1,end); % looking at the differences in CIs in 2090
    diff2 = CI3(2,end)-CI3(1,end);
    percChange2090(M) = (diff2 - diff1)/diff1 * 100; % negative value means CI is getting smaller
    
    if N==5
        state_ind_P_runs(M,:) = state_ind_P;
        state_ind_T_runs(M,:) = state_ind_T;
    end
    stateMsg = strcat('M= ', num2str(M));
    disp(stateMsg);
end

%% Create boxplot of 2050 data for Precip

f = figure('Position', [173 213 716 457]);

% load RCP85 data for certain parameters
load('bayesian_learning_RCP85_1000_CI80.mat', 's_P_abs', 't2', 'state_ind_P_runs')
state_ind_P_abs_RCP85 = s_P_abs(state_ind_P_runs(:,t2-1));

% RCP4.5 for perc Change in 2050
load('bayesian_learning_RCP45_1000_CI80.mat');

% put data in matrix for use in boxplot
% make matrix wide enough to encompass minP to maxP (we do not know which
% dataset the min and max are in)
state_ind_P_abs_RCP45 = s_P_abs(state_ind_P_runs(:,t2-1));
minP = min([state_ind_P_abs_RCP45, state_ind_P_abs_RCP85]);
maxP = max([state_ind_P_abs_RCP45, state_ind_P_abs_RCP85]);
rng = maxP - minP;
xlab = minP:1:maxP;
boxData = NaN(500,rng+1);
sum = 0;
for i=1:(rng+1)
    % if there are no values for a certain P value, skip it
    if(isempty(find(state_ind_P_abs_RCP45==minP+i-1)) == 0)
        ind = find(state_ind_P_abs_RCP45 == minP+i-1);
        boxData(1:numel(ind),i) = percChange2050(ind);
        sum = sum + numel(ind);
        disp({'hello', num2str(i), num2str(sum)});
    end
    
end

% create boxplot
pos_1=1:1:(rng+1);
b1 = boxplot(boxData, 'colors', 'r', 'width', 0.3,'symbol', '',...
    'position', pos_1);
hold on

% RCP8.5 for perc Change in 2050
load('bayesian_learning_RCP85_1000_CI80.mat');

% put data in matrix for use in boxplot

%rng = range(state_ind_P_abs);
boxData = NaN(500,rng+1);
for i=1:(rng+1)
    % if there are no values for a certain P value, skip it
    if(isempty(find(state_ind_P_abs_RCP85==minP+i-1)) == 0)
        ind = find(state_ind_P_abs_RCP85 == minP+i-1);
        boxData(1:numel(ind),i) = percChange2050(ind);
        disp({'hello', num2str(i)})
    end
    
end

pos_2 = 1.:1:(rng+1);
b2 = boxplot(boxData, 'colors', 'b', 'width', 0.3,...
    'symbol', '', 'position', pos_2);
xticklabels(xlab);
xlabel('2030 Precip (mm/mo)');
ylabel('% Change in CI MAS (MCM/y)');
legend([b1(4,1), b2(4,1)], {'RCP4.5', 'RCP8.5'}, 'Location', 'northeast');
title('%\Delta in 2050 CI (80% CI) of Mean Annual Shortage (beyond 10%) from 1990 to 2030 Compared to Precip Data');

%% Create boxplot of 2090 data for Precip

f = figure('Position', [173 213 716 457]);

% load RCP85 data for certain parameters
load('bayesian_learning_RCP85_1000_CI80.mat', 's_P_abs', 't2', 'state_ind_P_runs')
sz = size(state_ind_P_runs, 2);
state_ind_P_abs_RCP85 = s_P_abs(state_ind_P_runs(:,sz-1));

% RCP4.5 for perc Change in 2050
load('bayesian_learning_RCP45_1000_CI80.mat');

% put data in matrix for use in boxplot
% make matrix wide enough to encompass minP to maxP (we do not know which
% dataset the min and max are in)
state_ind_P_abs_RCP45 = s_P_abs(state_ind_P_runs(:,sz-1));
minP = min([state_ind_P_abs_RCP45, state_ind_P_abs_RCP85]);
maxP = max([state_ind_P_abs_RCP45, state_ind_P_abs_RCP85]);
rng = maxP - minP;
xlab = minP:1:maxP;
boxData = NaN(500,rng+1);
sum = 0;
for i=1:(rng+1)
    % if there are no values for a certain P value, skip it
    if(isempty(find(state_ind_P_abs_RCP45==minP+i-1)) == 0)
        ind = find(state_ind_P_abs_RCP45 == minP+i-1);
        boxData(1:numel(ind),i) = percChange2090(ind);
        sum = sum + numel(ind);
        disp({'hello', num2str(i), num2str(sum)});
    end
    
end

% create boxplot
pos_1=1:1:(rng+1);
b1 = boxplot(boxData, 'colors', 'r', 'width', 0.3,'symbol', '',...
    'position', pos_1);
hold on

% RCP8.5 for perc Change in 2050
load('bayesian_learning_RCP85_1000_CI80.mat');

% put data in matrix for use in boxplot

%rng = range(state_ind_P_abs);
boxData = NaN(500,rng+1);
for i=1:(rng+1)
    % if there are no values for a certain P value, skip it
    if(isempty(find(state_ind_P_abs_RCP85==minP+i-1)) == 0)
        ind = find(state_ind_P_abs_RCP85 == minP+i-1);
        boxData(1:numel(ind),i) = percChange2090(ind);
        disp({'hello', num2str(i)})
    end
    
end

pos_2 = 1.:1:(rng+1);
b2 = boxplot(boxData, 'colors', 'b', 'width', 0.3,...
    'symbol', '', 'position', pos_2);
xticklabels(xlab);
xlabel('2070 Precip (mm/mo)');
ylabel('% Change in CI MAS (MCM/y)');
legend([b1(4,1), b2(4,1)], {'RCP4.5', 'RCP8.5'}, 'Location', 'northeast');
title('%\Delta in 2090 CI (80% CI) of Mean Annual Shortage (beyond 10%) from 1990 to 2070 Compared to Precip Data');

%% Create boxplot of 2050 data for Temp

f = figure('Position', [173 213 716 457]);
load('bayesian_learning_RCP85_1000_CI80.mat', 'state_ind_T_runs', 't2', 's_T_abs');
state_ind_T_abs_RCP85 = s_T_abs(state_ind_T_runs(:,t2-1));

% RCP4.5 for perc Change in 2050
load('bayesian_learning_RCP45_1000_CI80.mat');

% put data in matrix for use in boxplot
state_ind_T_abs_RCP45 = s_T_abs(state_ind_T_runs(:,t2-1));
rng = round((max(state_ind_T_abs_RCP85) - min(state_ind_T_abs_RCP45))*20);
xlab = min(state_ind_T_abs_RCP45):0.05:max(state_ind_T_abs_RCP85);
boxData = NaN(1000,rng+1);
sum = 0;
for i=1:(rng+1)
    % if there are no values for a certain T value, skip it
    if(isempty(find(round(state_ind_T_abs_RCP45,2)==round(min(state_ind_T_abs_RCP45)+(i-1)/20,2))) == 0)
        ind = find(round(state_ind_T_abs_RCP45,2) == (round(min(state_ind_T_abs_RCP45)+(i-1)/20,2)));
        boxData(1:numel(ind),i) = percChange2050(ind);
        sum = sum + numel(ind);
        disp({'hello', num2str(i), num2str(sum)});
    end
    
end
   
pos_1=1:1:(rng+1);
b1 = boxplot(boxData, 'colors', 'r', 'width', 0.3,'symbol', '',...
    'position', pos_1);
hold on

% RCP8.5 for perc Change in 2050
load('bayesian_learning_RCP85_1000_CI80.mat');

% put data in matrix for use in boxplot
%state_ind_T_abs_RCP85 = s_P_abs(state_ind_T_runs(:,t2-1));
%rng = range(state_ind_P_abs);
boxData = NaN(1000,rng+1);
for i=1:(rng+1)
    % if there are no values for a certain P value, skip it
    if(isempty(find(round(state_ind_T_abs_RCP85,2)==round(min(state_ind_T_abs_RCP45)+(i-1)/20,2))) == 0)
        ind = find(round(state_ind_T_abs_RCP85,2) == (round(min(state_ind_T_abs_RCP45)+(i-1)/20,2)));
        boxData(1:numel(ind),i) = percChange2050(ind);
        disp({'hello', num2str(i)})
    end
    
end

pos_2 = 1.4:1:(rng + 1.4);
b2 = boxplot(boxData, 'colors', 'b', 'width', 0.3,...
    'symbol', '', 'position', pos_2);
ax=gca;
ax.XTick = ax.XTick - 0.2;
ax.XLim = ax.XLim - 0.2;
xticklabels(xlab);
xlabel('2030 Temp (^oC)');
ylabel('% Change in CI MAS (MCM/y)');
legend([b1(4,1), b2(4,1)], {'RCP4.5', 'RCP8.5'}, 'Location', 'northeast');
title('%\Delta in 2050 CI (80%) of Mean Annual Shortage (beyond 10%) from 1990 to 2030 Compared to Temp Data');

%% Create boxplot of 2090 data for Temp

f = figure('Position', [114 157 1100 420]);
load('bayesian_learning_RCP85_1000_CI80.mat', 'state_ind_T_runs', 's_T_abs');
sz = size(state_ind_T_runs, 2);
% get abs and relative temperature states
state_ind_T_abs_RCP85 = s_T_abs(state_ind_T_runs(:,sz-1));
minT_RCP85 = min(state_ind_T_runs(:,sz-1));
state_ind_T_rel_RCP85 = state_ind_T_runs(:,sz-1) - minT_RCP85;

% RCP4.5 for perc Change in 2090
load('bayesian_learning_RCP45_1000_CI80.mat');

% put data in matrix for use in boxplot
minT_RCP45 = min(state_ind_T_runs(:,sz-1));
state_ind_T_rel_RCP45 = state_ind_T_runs(:,sz-1) - minT_RCP45;
state_ind_T_abs_RCP45 = s_T_abs(state_ind_T_runs(:,sz-1));
% use relative temp states so that can put RCP4.5 and RCP8.5 on same plot
minT = min(min([state_ind_T_rel_RCP45, state_ind_T_rel_RCP85]));
maxT = max(max([state_ind_T_rel_RCP45, state_ind_T_rel_RCP85]));
rng = maxT - minT;
boxData = NaN(1000,rng+1);
sum = 0;
for i=1:(rng+1)
    % if there are no values for a certain T value, skip it
    if(isempty(find(state_ind_T_rel_RCP45==(minT+i-1))) == 0)
        ind = find(state_ind_T_rel_RCP45 == (minT+i-1));
        boxData(1:numel(ind),i) = percChange2090(ind);
        sum = sum + numel(ind);
        disp({'hello', num2str(i), num2str(sum)});
    end
    
end
   
pos_1=1:1:(rng+1);
b1 = boxplot(boxData, 'colors', 'r', 'width', 0.3,'symbol', '',...
    'position', pos_1);
hold on

% RCP8.5 for perc Change in 2050
load('bayesian_learning_RCP85_1000_CI80.mat');

% put data in matrix for use in boxplot
boxData = NaN(1000,rng+1);
for i=1:(rng+1)
    % if there are no values for a certain P value, skip it
    if(isempty(find(state_ind_T_rel_RCP85==(minT+i-1))) == 0)
        ind = find(state_ind_T_rel_RCP85 == (minT+i-1));
        boxData(1:numel(ind),i) = percChange2090(ind);
        disp({'hello', num2str(i)})
    end
    
end

pos_2 = 1.4:1:(rng+1.4);
b2 = boxplot(boxData, 'colors', 'b', 'width', 0.3,...
    'symbol', '', 'position', pos_2);
% format x axis
ax=gca;
ax.XTick = ax.XTick - 0.2;
ax.XLim = ax.XLim - 0.2;
orig = strcat(num2str(s_T_abs(minT_RCP45)), '/', num2str(s_T_abs(minT_RCP85)));
rep = strcat('+', string((1:1:rng) * 0.05));
xlab = [orig, rep];
xticklabels(xlab);
ylim([-110 100]);
xlabel('2070 Temp (^oC)');
ylabel('% Change in CI MAS (MCM/y)');
legend([b1(4,1), b2(4,1)], {'RCP4.5', 'RCP8.5'}, 'Location', 'northeast');
sgtitle('%\Delta in 2090 CI (80%) of Mean Annual Shortage (beyond 10%) from 1990 to 2070 Compared to Temp Data');

%% 2 x 2 boxplots looking at 2050 and 2090 data for P and T
f1 = figure('Position', [126 59 1019 634])
% load RCP4.5 data
load('bayesian_learning_RCP45.mat');

subplot(2,2,1)
state_ind_P_abs = s_P_abs(state_ind_P_runs(:,t2-1));
b1 = boxplot(percChange2050, state_ind_P_abs, 'colors', 'r', 'width', 0.15,...
    'symbol', '');

subplot(2,2,2)
state_ind_P_abs = s_P_abs(state_ind_P_runs(:,end-1));
b2 = boxplot(percChange2090, state_ind_P_abs, 'colors', 'r', 'width', 0.15,...
    'symbol', '');

% load RCP8.5 data
load('bayesian_learning_RCP85.mat');
subplot(2,2,3)
state_ind_P_abs = s_P_abs(state_ind_P_runs(:,t2-1));
b3 = boxplot(percChange2050, state_ind_P_abs, 'colors', 'b', 'width', 0.15,...
    'symbol', '');

subplot(2,2,4)
state_ind_P_abs = s_P_abs(state_ind_P_runs(:,end-1));
b4 = boxplot(percChange2090, state_ind_P_abs, 'colors', 'b', 'width', 0.15,...
    'symbol', '');

f2 = figure('Position', [126 59 1019 634])
% load RCP4.5 data
load('bayesian_learning_RCP45.mat');

subplot(2,2,1)
state_ind_T_abs = s_T_abs(state_ind_T_runs(:,t2-1));
b1 = boxplot(percChange2050, state_ind_T_abs, 'colors', 'r', 'width', 0.15,...
    'symbol', '');

subplot(2,2,2)
state_ind_T_abs = s_T_abs(state_ind_T_runs(:,end-1));
b2 = boxplot(percChange2090, state_ind_T_abs, 'colors', 'r', 'width', 0.15,...
    'symbol', '');

% load RCP8.5 data
load('bayesian_learning_RCP85.mat');
subplot(2,2,3)
state_ind_T_abs = s_T_abs(state_ind_T_runs(:,t2-1));
b3 = boxplot(percChange2050, state_ind_T_abs, 'colors', 'b', 'width', 0.15,...
    'symbol', '');

subplot(2,2,4)
state_ind_T_abs = s_T_abs(state_ind_T_runs(:,end-1));
b4 = boxplot(percChange2090, state_ind_T_abs, 'colors', 'b', 'width', 0.15,...
    'symbol', '');

%% 

figure;
load('bayesian_learning_RCP45.mat');
subplot(2,1,1)
h1 = histogram(percChange2050, 'FaceColor', 'r', 'BinWidth', 5,...
    'Normalization', 'probability');
xlabel('% Change in 2050 MAR CI');
ylabel('Probability');
title('% Change in 2050 MAR CI from 1990 to 2030')
hold on

subplot(2,1,2)
h2 = histogram(percChange2090, 'FaceColor', 'r', 'BinWidth', 10,...
    'Normalization', 'probability');
xlabel('% Change in 2090 MAR CI');
ylabel('Probability');
title('% Change in 2090 MAR CI from 1990 to 2070')
hold on

load('bayesian_learning_RCP85.mat');
subplot(2,1,1)
h3 = histogram(percChange2050, 'FaceColor', 'b', 'BinWidth', 5,...
    'Normalization', 'probability');
legend([h1, h3], {'RCP4.5', 'RCP8.5'}, 'Location', 'northeast');

subplot(2,1,2)
h4 = histogram(percChange2090, 'FaceColor', 'b', 'BinWidth', 10,...
    'Normalization', 'probability');

%% Create figure

fig = figure;
addpath('/Users/jenniferskerker/Documents/GradSchool/Research/Project1/Code/Models/Fletcher_2019_Learning_Climate/Plots/plotting utilities/cbrewer/cbrewer/cbrewer');
[clrmp1]=cbrewer('seq', 'Reds', N);
[clrmp2]=cbrewer('seq', 'Blues', N);
[clrmp3]=cbrewer('seq', 'Greens', N);
[clrmp4]=cbrewer('seq', 'Purples', N);
%set(fig,'Position', [680 558 1400 750])
set(fig, 'Position', [357 38 648 660]);

for t =1:N
    x = t:N+1;
    X=[x,fliplr(x)];
    T_p01 = prctile(T_over_time{t},.01);
    T_p995 = prctile(T_over_time{t},99.9);
    P_p01 = prctile(P_over_time{t},.01);
    P_p995 = prctile(P_over_time{t},99.9);
    MAR_p01 = prctile(MAR_over_time{t},.01);
    MAR_p995 = prctile(MAR_over_time{t},99.9);
    yield_p01 = prctile(yield_over_time{t},.01);
    yield_p995 = prctile(yield_over_time{t},99.9);
    
    subplot(2,2,1)
    Y=[T_p01,fliplr(T_p995)];
    Y = Y - s_T_abs(state_ind_T(1));
    hold on
    fill(X,Y,clrmp1(t,:), 'LineWidth', 1);
    scatter(t,s_T_abs(state_ind_T(t)) - s_T_abs(state_ind_T(1)), 'k', 'MarkerFaceColor', 'k') 
    xticks(1:6)
    xticklabels(decades)
    ylabel('Degrees C')
    title('Cumulative T Change')
    %ylim([-.1 3.75])
    ylim([-.1 2.5])
    ytickformat('%.1f')
    %set(gca,'Units','normalized')
    if t ==1
    yLabelHandle = get( gca ,'YLabel' );
    pos  = get( yLabelHandle , 'position' )
    pos1 = pos - [0.15 0 0]; 
    set( yLabelHandle , 'position' , pos1 );
    end
    
    subplot(2,2,2)
    Y=[P_p01,fliplr(P_p995)];
    hold on
    fill(X,Y-P0_abs,clrmp2(t,:), 'LineWidth', 1);
    scatter(t,s_P_abs(state_ind_P(t))-P0_abs, 'k', 'MarkerFaceColor', 'k') 
    xticks(1:6)
    xticklabels(decades)
    ylabel('mm/month')
    title('Cumulative P Change')
    ylim([-18 20])
    set(gca,'Units','normalized')
    if t ==1
    yLabelHandle = get( gca ,'YLabel' );
    pos  = get( yLabelHandle , 'position' )
    pos1 = pos - [0.15 0 0]; 
    set( yLabelHandle , 'position' , pos1 );
    end
    
    subplot(2,2,3)
    Y=[MAR_p01,fliplr(MAR_p995)];
    hold on
    fill(X,Y,clrmp3(t,:), 'LineWidth', 1);
    scatter(t,MAR_over_time{t}(1,1), 'k', 'MarkerFaceColor', 'k') 
    xticks(1:6)
    xticklabels(decades)
    ylabel('MCM/y')
    title('Mean Annual Runoff')
    set(gca,'Units','normalized')
    if t ==1
    yLabelHandle = get( gca ,'YLabel' );
    pos  = get( yLabelHandle , 'position' )
    pos1 = pos - [0.15 0 0]; 
    set( yLabelHandle , 'position' , pos1 );
    end
    
    subplot(2,2,4)
    Y=[yield_p01,fliplr(yield_p995)];
    hold on
    fill(X,Y,clrmp4(t,:), 'LineWidth', 1);
    scatter(t,yield_over_time{t}(1,1), 'k', 'MarkerFaceColor', 'k') 
    xticks(1:6)
    xticklabels(decades)
    ylabel('MCM/y')
    title('Mean Annual Shortage (beyond 10%) ')
    ylim([0 20])
    set(gca,'Units','normalized')
    if t ==1
    yLabelHandle = get( gca ,'YLabel' );
    pos  = get( yLabelHandle , 'position' )
    pos1 = pos - [0.15 0 0]; 
    set( yLabelHandle , 'position' , pos1 );
    end
        
    frames(t) = getframe(gcf);
    
end

%% 


% Example time series 2: dry path

% Set time series
state_ind_P = zeros(1,N);
state_ind_T = zeros(1,N);
state_ind_P(1) =  find(P0_abs==s_P_abs);
state_ind_T(1) = find(T0_abs==s_T_abs);
randGen = false;
state_ind_P(2:N) = [9 9 8 7];
%state_ind_T(2:N) = [9 17 26 33];
state_ind_T(2:N) = [10 19 26 33];


MAR = cellfun(@(x) mean(mean(x)), runoff);
p = randi(numSamp,N-1);
T_over_time = cell(1,N);
P_over_time = cell(1,N);
MAR_over_time = cell(1,N);

for t = 1:N
    % Sample forward distribution given current state
    T_current = s_T_abs(state_ind_T(t));
    P_current = s_P_abs(state_ind_P(t));
    [T_over_time{t}] = T2forwardSimTemp(T_Temp, s_T_abs, N, t, T_current, numSamp, false);
    [P_over_time{t}] = T2forwardSimTemp(T_Precip, s_P_abs, N, t, P_current, numSamp, false);
    
    % Lookup MAR and yield for forward distribution
    T_ind = arrayfun(@(x) find(x == s_T_abs), T_over_time{t});
    P_ind = arrayfun(@(x) find(x == s_P_abs), P_over_time{t});
    [~,t_steps] = size(T_ind);
    MAR_over_time{t} = zeros(size(T_ind));
    yield_over_time{t} = zeros(size(T_ind));
    for i = 1:numSamp
        for j = 1:t_steps   
            MAR_over_time{t}(i,j) = MAR(T_ind(i,j), P_ind(i,j), 1);
            yield_over_time{t}(i,j) = unmet_dom(T_ind(i,j), P_ind(i,j),1, 1) ;   % 80 MCM storage
        end
    end
    
    % Sample next time period
    if randGen
        state_ind_T(t+1) = find(T_over_time{t}(p(t),2)==s_T_abs);
        state_ind_P(t+1) = find(P_over_time{t}(p(t),2)==s_P_abs);
    end
end



for t =1:N
    x = t:N+1;
    X=[x,fliplr(x)];
    T_p01 = prctile(T_over_time{t},.01);
    T_p995 = prctile(T_over_time{t},99.9);
    P_p01 = prctile(P_over_time{t},.01);
    P_p995 = prctile(P_over_time{t},99.9);
    MAR_p01 = prctile(MAR_over_time{t},.01);
    MAR_p995 = prctile(MAR_over_time{t},99.9);
    yield_p01 = prctile(yield_over_time{t},.01);
    yield_p995 = prctile(yield_over_time{t},99.9);
    
    subplot(4,2,5)
    Y=[T_p01,fliplr(T_p995)];
    Y = Y - s_T_abs(state_ind_T(1));
    hold on
    fill(X,Y,clrmp1(t,:), 'LineWidth', 1);
    scatter(t,s_T_abs(state_ind_T(t)) - s_T_abs(state_ind_T(1)), 'k', 'MarkerFaceColor', 'k') 
    xticks(1:6)
    xticklabels(decades)
    ylabel('Degrees C')
    title('Cumulative T Change')
    %ylim([-.1 3.75])
    ylim([-.1 2.5])
    ytickformat('%.1f')
    set(gca,'Units','normalized')
    if t ==1
    yLabelHandle = get( gca ,'YLabel' );
    pos  = get( yLabelHandle , 'position' )
    pos1 = pos - [0.15 0 0]; 
    set( yLabelHandle , 'position' , pos1 );
    end
    
    subplot(4,2,6)
    Y=[P_p01,fliplr(P_p995)];
    hold on
    fill(X,Y-P0_abs,clrmp2(t,:), 'LineWidth', 1);
    scatter(t,s_P_abs(state_ind_P(t))-P0_abs, 'k', 'MarkerFaceColor', 'k') 
    xticks(1:6)
    xticklabels(decades)
    ylabel('mm/month')
    title('Cumulative P Change')
    ylim([-18 20])
    set(gca,'Units','normalized')
    if t ==1
    yLabelHandle = get( gca ,'YLabel' );
    pos  = get( yLabelHandle , 'position' )
    pos1 = pos - [0.15 0 0]; 
    set( yLabelHandle , 'position' , pos1 );
    end
    
    subplot(4,2,7)
    Y=[MAR_p01,fliplr(MAR_p995)];
    hold on
    fill(X,Y,clrmp3(t,:), 'LineWidth', 1);
    scatter(t,MAR_over_time{t}(1,1), 'k', 'MarkerFaceColor', 'k') 
    xticks(1:6)
    xticklabels(decades)
    ylabel('MCM/y')
    title('Mean Annual Runoff')
    set(gca,'Units','normalized')
    if t ==1
    yLabelHandle = get( gca ,'YLabel' );
    pos  = get( yLabelHandle , 'position' )
    pos1 = pos - [0.15 0 0]; 
    set( yLabelHandle , 'position' , pos1 );
    end
    
    subplot(4,2,8)
    Y=[yield_p01,fliplr(yield_p995)];
    hold on
    fill(X,Y,clrmp4(t,:), 'LineWidth', 1);
    scatter(t,yield_over_time{t}(1,1), 'k', 'MarkerFaceColor', 'k') 
    xticks(1:6)
    xticklabels(decades)
    ylabel('MCM/y')
    title('Mean Annual Shortage (beyond 10%) ')
    ylim([0 20])
    set(gca,'Units','normalized')
    if t ==1
    yLabelHandle = get( gca ,'YLabel' );
    pos  = get( yLabelHandle , 'position' )
    pos1 = pos - [0.15 0 0]; 
    set( yLabelHandle , 'position' , pos1 );
    end
        
    frames(t) = getframe(gcf);
    
end




set(findall(fig.Children,'Type', 'Scatter'), 'SizeData', 10)

