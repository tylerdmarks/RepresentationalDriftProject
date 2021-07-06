% select stimulus to study
stim = 'NatMov';
switch stim
    case 'NatMov'
        F0 = 31:50;
        offidx = logical([ones(1, 50) zeros(1, 300)]);
    case 'PDG'
        F0 = logical(repmat([zeros(1, 20) ones(1, 20) zeros(1, 20)], 1, 12));
        offidx = logical(repmat([ones(1, 40) zeros(1, 20)], 1, 12));
end

all = StabilityAnalyzer;            % initialize object
all.importData();                   % import first field
all.addData();                      % add each subsequent field
all.setQualityThreshold(3);         % set quality threshold
all.cellSelection;                  % select inclusion criteria
num_sessions = all.num_sessions;

qual = all.getUse_cells;
sessions = all.getUse_sessions;
sessions(isnan(sessions)) = 1;              % NaN values to 1s for logical vector purposes
present = sum(sessions, 2) == num_sessions;
good_cells = qual' & present;

Resp = all.RespData.(stim).RespMat_Full(:, :, good_cells, :);

catResp = [];
for kk = 1:all.num_sessions
    catResp = cat(1, catResp, Resp(:, :, :, kk));
end

% initialize another object that will pool spike data from the above fields 
% (run deconvolution of DFF data to generate spike data for each field, 
% then run the stability data preprocessing pipeline on it to generate datafiles for use here (generate 'RespData' data structure using spike data)
TDM077an_spdata = StabilityAnalyzer;     
TDM077an_spdata.importData();               % import datafiles generated from spike data for first mouse
TDM077an_spdata.addData();                  % add datafiles from subsequent mice

all.setSpikedata(TDM077an_spdata);       % feed spike data StabilityAnalyzer object into our first object for use in eventDetector

eventdata = all.extractWaveforms(stim, 5, 0.01, 2);             % input args: stimulus, min event length, sig threshold, min sessions
% save eventdata_pdg eventdata                          % save data

event_vectors = eventdata.event_vectors(good_cells, :);         % extract data based on inclusion criteria
ts = eventdata.waveform_ts(good_cells);
waveforms = eventdata.waveforms(good_cells);

% Get rid of neurons whose max responses occur during the offtime
% mean_catResp = squeeze(nanmean(catResp, 1));            % averaged across trials
% [~, max_idx] = max(mean_catResp, [], 1);
% 
% remove_idx = false(1, sum(good_cells));
% for ii = 1:sum(good_cells)
%     if any(max_idx(ii) == find(offidx))
%         remove_idx(ii) = true;
%     end
% end
% 
% event_vectors(remove_idx, :) = [];
% ts(remove_idx) = [];
% waveforms(remove_idx) = [];
% Resp(:, :, remove_idx, :) = [];
% catResp(:, :, remove_idx) = [];

%% Events per neuron for both stimuli

all.cellSelection;               % select single stimulus responsive (for desired stimulus), + 'quality' + 'presence'
num_sessions = all.num_sessions;
qual = all.getUse_cells;
sessions = all.getUse_sessions;
sessions(isnan(sessions)) = 1;
present = sum(sessions, 2) == num_sessions;
good_cells = qual' & present;
waveforms = eventdata.waveforms(good_cells);
num_waves = cellfun(@length, waveforms);

figure
hold on
histogram(num_waves(num_waves > 0), 'FaceAlpha', 1)

all.cellSelection;           % select dual responsive (+ 'quality' + 'presence')
num_sessions = all.num_sessions;
qual = all.getUse_cells;
sessions = all.getUse_sessions;
sessions(isnan(sessions)) = 1;
present = sum(sessions, 2) == num_sessions;
good_cells = qual' & present;
waveforms = eventdata.waveforms(good_cells);
num_waves = cellfun(@length, waveforms);

histogram(num_waves(num_waves > 0), 'FaceAlpha', 1)
xlim([0 13])


%% Evidence of events appearing, disappearing, staying the same, or time-shifting

% For each neurons's event waveforms, test whether first and last session frame-avg DF/F distributions (trials) are significantly different.
% If first session significantly > last assign -1, if last session significantly > first assign 1, if no difference assign 0
num_cells = size(event_vectors, 1);
zScores = cell(1, num_cells);
deltas = cell(1, num_cells);
offResp = Resp(:, F0, :, :);


%zscore trial-by-trial waveforms using off response
for ii = 1:num_cells
    num_waves = length(waveforms{ii});
    for ww = 1:num_waves
        curr_wave = waveforms{ii}{ww};
        curr_off = squeeze(offResp(:, :, ii, 1:size(curr_wave, 3)));
        zScores{ii}{ww} = squeeze(nanmean(curr_wave, 2))./squeeze(nanstd(curr_off, [], 2));
    end
end

%determine the 'delta' of the waveform across sessions
for ii = 1:num_cells
    num_waves = length(waveforms{ii});
    for ww = 1:num_waves
        curr_wave = zScores{ii}{ww};
        ses1_vals = cat(1, curr_wave(:, 1), curr_wave(:, 2));
        sesf_vals = cat(1, curr_wave(:, end-1), curr_wave(:, end));
        p = ranksum(ses1_vals, sesf_vals);
        if p < 0.01
            if mean(ses1_vals) > mean(sesf_vals)
                deltas{ii}(ww) = -1;
            else
                deltas{ii}(ww) = 1;
            end
        else
            deltas{ii}(ww) = 0;
        end
    end
end

%visualizing deltas (coloring events)
visualizer = zeros(length(deltas), size(event_vectors, 2), 3);          % create RGB image matrix, fill in colors based on deltas of each waveform
for ii = 1:length(deltas)
    deltalist = deltas{ii};
    for ww = 1:length(deltalist)
        curr_ts = ts{ii}(ww, :);
        curr_delta = deltalist(ww);
        switch curr_delta
            case -1
                visualizer(ii, curr_ts(1):curr_ts(2), 1) = 1;
            case 0 
                visualizer(ii, curr_ts(1):curr_ts(2), 2) = 1;
            case 1
                visualizer(ii, curr_ts(1):curr_ts(2), 3) = 1;
                visualizer(ii, curr_ts(1):curr_ts(2), 2) = 0.5;
        end
    end
end

mean_catResp = squeeze(nanmean(catResp, 1));            % averaged across trials
[maxima, max_idx] = max(mean_catResp, [], 1);
meanmean_catResp = mean(mean_catResp, 1);           % averaged across trials and time

% sorted by max response time
[times, I] = sort(max_idx);                    
visualizer = visualizer(I, :, :);
visualizer(times < 51, :, :) = [];

% or sorted by max response magnitude
[~, I] = sort(maxima, 'descend');    
% or sorted by average response magnitude
[~, I] = sort(meanmean_catResp, 'descend');      


visualizer = visualizer(I, :, :);


figure
imagesc(visualizer);
box off
xlim([50 350])


% proportions of event types surrounding 1s of a neuron's peak response
mean_catResp = squeeze(nanmean(catResp, 1));            % averaged across trials
[maxima, max_idx] = max(mean_catResp, [], 1);

%iterate through every event, find its response time, if it's sufficiently close to the neuron's max respone add this to one of three vectors based on the event type
ups = 0;
downs = 0;
statics = 0;
num_waves = cellfun(@length, waveforms);

for ii = 1:num_cells
    for ww = 1:num_waves(ii)
        if max_idx(ii) > 50
            curr_ts = ts{ii}(ww, :);
            dis = floor(abs(mean(curr_ts) - max_idx(ii)));          % distance from event's center to neuron's max response (in frames)
            if dis <= 5             %frames
                curr_delta = deltas{ii}(ww);
                switch curr_delta
                    case 1
                        ups = ups + 1;
                    case -1
                        downs = downs + 1;
                    case 0
                        statics = statics + 1;
                end
            end
        end
    end
end

figure
pie([ups downs statics])
title('Strongest event type proportions')


%% for every field separately, determine proportion of events with positive, negative, and 0 delta
total_events = sum(cellfun(@length, deltas));
posevents = sum(cellfun(@length, cellfun(@(x) find(x == 1), deltas, 'UniformOutput', 0)))/total_events;      %fraction of positive deltas
negevents = sum(cellfun(@length, cellfun(@(x) find(x == -1), deltas, 'UniformOutput', 0)))/total_events;     %fraction of negative deltas
zeroevents = sum(cellfun(@length, cellfun(@(x) find(x == 0), deltas, 'UniformOutput', 0)))/total_events;        %fraction of no deltas

%pool data

% cat_posevents = [];  % initialize once
% cat_negevents = [];
% cat_zeroevents = [];

% add each field
cat_posevents = [cat_posevents posevents];
cat_negevents = [cat_negevents negevents];
cat_zeroevents = [cat_zeroevents zeroevents];

save cat_posevents_natmov cat_posevents             % save separate files for pdg and natmov
save cat_negevents_natmov cat_negevents
save cat_zeroevents_natmov cat_zeroevents

%plotting 
% cat_posevents_pdg = importdata('cat_posevents_pdg.mat');
% cat_negevents_pdg = importdata('cat_negevents_pdg.mat');
% cat_zeroevents_pdg = importdata('cat_zeroevents_pdg.mat');
% cat_posevents_natmov = importdata('cat_posevents_natmov.mat');
% cat_negevents_natmov = importdata('cat_negevents_natmov.mat');
% cat_zeroevents_natmov = importdata('cat_zeroevents_natmov.mat');

mean_posevents_pdg = mean(cat_posevents_pdg);
mean_negevents_pdg = mean(cat_negevents_pdg);
mean_zeroevents_pdg = mean(cat_zeroevents_pdg);
SEM_posevents_pdg = std(cat_posevents_pdg)/sqrt(length(cat_posevents_pdg));
SEM_negevents_pdg = std(cat_negevents_pdg)/sqrt(length(cat_negevents_pdg));
SEM_zeroevents_pdg = std(cat_zeroevents_pdg)/sqrt(length(cat_zeroevents_pdg));

mean_posevents_natmov = mean(cat_posevents_natmov);
mean_negevents_natmov = mean(cat_negevents_natmov);
mean_zeroevents_natmov = mean(cat_zeroevents_natmov);
SEM_posevents_natmov = std(cat_posevents_natmov)/sqrt(length(cat_posevents_natmov));
SEM_negevents_natmov = std(cat_negevents_natmov)/sqrt(length(cat_negevents_natmov));
SEM_zeroevents_natmov = std(cat_zeroevents_natmov)/sqrt(length(cat_zeroevents_natmov));

bar_vector = [mean_posevents_pdg mean_posevents_natmov; mean_negevents_pdg mean_negevents_natmov; mean_zeroevents_pdg mean_zeroevents_natmov];
error_vector = [SEM_posevents_pdg SEM_posevents_natmov; SEM_negevents_pdg SEM_negevents_natmov; SEM_zeroevents_pdg SEM_zeroevents_natmov];

figure
bar(bar_vector);
hold on

ngroups = size(bar_vector, 1);
nbars = size(bar_vector, 2);
% Calculating the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
for nn = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*nn-1) * groupwidth / (2*nbars);
    errorbar(x, bar_vector(:,nn), error_vector(:,nn), '.');
end

xlabel('1 is up, 2 is down, 3 is no change')
ylabel('Proportion of events')
[~, p1, ~, stats] = ttest(cat_posevents_pdg, cat_posevents_natmov);
tstat1 = stats.tstat;
[~, p2, ~, stats] = ttest(cat_negevents_pdg, cat_negevents_natmov);
tstat2 = stats.tstat;
[~, p3, ~, stats] = ttest(cat_zeroevents_pdg, cat_zeroevents_natmov);
tstat3 = stats.tstat;


%% What influences the way in which an event changes over time?
% Delta zScore over time as a function of 
% 1) Ses-average event zscore (or 1st session event zscore)
% 2) Length of the event (unused)
% 3) Time of event (unused)
% 4) 'Redundancy' of an event

num_waves = cellfun(@length, waveforms);
num_cells = length(waveforms);
% 1) Calculate by-session zscore for every event
offResp = Resp(:, F0, :, :);
sesZScores = cell(1, num_cells);
meanZScores = cell(1, num_cells);
firstZScores = cell(1, num_cells);

for ii = 1:num_cells
    for ww = 1:num_waves(ii)
        curr_wave = squeeze(nanmean(waveforms{ii}{ww}, 1));
        curr_off = squeeze(nanmean(squeeze(offResp(:, :, ii, 1:size(waveforms{ii}{ww}, 3))), 1));
        for kk = 1:size(waveforms{ii}{ww}, 3)
            sesZScores{ii}{ww}(kk) = nanmean(curr_wave(:, kk), 1)/nanstd(curr_off(:, kk), [], 1);                   % rectify negative values to 0
        end
        meanZScores{ii}(ww) = mean(sesZScores{ii}{ww});
        firstZScores{ii}(ww) = sesZScores{ii}{ww}(1);
    end
end


% % 2) Calculate length of each event
% event_lengths = cell(1, num_cells);
% 
% for ii = 1:num_cells
%     for ww = 1:num_waves(ii)
%         event_lengths{ii}(ww) = size(waveforms{ii}{ww}, 2);
%     end
% end

% 3) Log start time of each event
% start_times = cell(1, num_cells);
% 
% for ii = 1:num_cells
%     for ww = 1:num_waves(ii)
%         start_times{ii}(ww) = ts{ii}(ww, 1);
%     end
% end

% 4) Calculate 'redundancy' of each event, defined by [number of neurons with an event that overlaps at least n percent with a given event / total number of neurons]

%find number of good cells per mouse
total_cells_per_mouse = all.num_cells_bymouse;
good_cells_bymouse = zeros(1, length(total_cells_per_mouse));
ct = 0;
for mm = 1:length(total_cells_per_mouse)
    curr_idx = ct+1:ct+total_cells_per_mouse(mm);
    curr_goodcells = good_cells(curr_idx);
    good_cells_bymouse(mm) = sum(curr_goodcells);
    ct = ct + total_cells_per_mouse(mm);
end 
    
redundancy = cell(1, num_cells);
pct = 0.75;
ct = 0;
bymouse_idx = 1;
for ii = 1:num_cells                % iterate through every neuron
    ii
    if ii > ct + good_cells_bymouse(bymouse_idx)                 % when ii exceeds the count for the current population, jump the counter that iterates jj and the population idx
        ct = ct + good_cells_bymouse(bymouse_idx);
        bymouse_idx = bymouse_idx + 1;
    end
    for ww = 1:num_waves(ii)        % iterate through each of its events
        curr_overlapping = 0;
        curr_vector = zeros(1, size(Resp, 2));
        curr_vector(ts{ii}(ww, 1):ts{ii}(ww, 2)) = 1;
        for jj = ct+1:ct+good_cells_bymouse(bymouse_idx)         % iterate only through all the neurons that belong to neuron ii's population of good_cells
            if ii ~= jj
                stop_flag = 0;
                for qq = 1:num_waves(jj)             % iterate through this neuron's events
                    if ~stop_flag
                        temp_vector = zeros(1, size(Resp, 2));
                        temp_vector(ts{jj}(qq, 1):ts{jj}(qq, 2)) = 1;
                        comp = curr_vector & temp_vector;
                        if sum(comp)/sum(curr_vector) >= pct
                            curr_overlapping = curr_overlapping + 1;
                            stop_flag = 1;
                        end
                    end
                end
            end
        end
        redundancy{ii}(ww) = curr_overlapping/good_cells_bymouse(bymouse_idx);
    end
end

% X) Calculate difference between early and later zScores for every event
ZScorediff = cell(1, num_cells);
pseudoRDI = cell(1, num_cells);
for ii = 1:num_cells
    for ww = 1:num_waves(ii)
        late = mean(sesZScores{ii}{ww}(end-1:end));
        early = mean(sesZScores{ii}{ww}(1:2));
        if late < 0 || early < 0
            ZScorediff{ii}(ww) = NaN;
            pseudoRDI{ii}(ww) = NaN;
        else
            ZScorediff{ii}(ww) = abs(late-early);         % Comparing average zScore of first two sessions to average zScore of last two sessions, 
            pseudoRDI{ii}(ww) = abs(late-early)/(late+early);
        end
    end
end


% Generate y vector for plotting
y = pseudoRDI;        
cat_y = [];
for ii = 1:num_cells
    for ww = 1:num_waves(ii)
        cat_y = [cat_y y{ii}(ww)];                 % use either delta zScore or std zScore here
    end
end

%Generate x vector for plotting
x = redundancy;        % choose: meanZScores firstZScores start_times redundancy
xlab = 'Event average zScore';
cat_x = [];     
for ii = 1:num_cells
    for ww = 1:num_waves(ii)
        cat_x = [cat_x x{ii}(ww)];
    end
end


figure
axis square
hold on
edge_vector = 0:10:100;                         
edges = prctile(cat_x, edge_vector);
binned_x = discretize(cat_x, edges);
CI = zeros(1, length(edges)-1);
binned_y = zeros(1, length(edges)-1);
midpoints = zeros(1, length(edges)-1);
std_binned_y = zeros(1, length(edges)-1);
SEM_binned_y = zeros(1, length(edges)-1);
t = zeros(length(edges)-1, 2);
for bb = 1:length(edges)-1
    binned_y(bb) = nanmean(cat_y(binned_x == bb));
    SEM_binned_y(bb) = nanstd(cat_y(binned_x == bb))/sqrt(sum(binned_x == bb));
    std_binned_y(bb) = nanstd(cat_y(binned_x == bb));
    midpoints(bb) = (edges(bb+1) + edges(bb))/2;
    t(bb, :) = tinv([0.025  0.975],sum(binned_x == bb)-1);      % T-Score
    temp = t(bb, :)*SEM_binned_y(bb);          % not true CI, just want the +-
    CI(bb) = temp(:, 2);
%     temp = bootci(100, {@(x)nanmean(x), cat_y(binned_x == bb)}, 'alpha', 0.05); 
%     CI(bb) = temp(2) - temp(1);
end
% %OPTIONAL: run this if you want to force axes limits
% xlimits = [0 20];
% ylimits = [];
% if ~isempty(xlimits) 
%     cat_x(cat_x < xlimits(1)) = xlimits(1);
%     cat_x(cat_x > xlimits(2)) = xlimits(2);
%     midpoints(midpoints < xlimits(1)) = xlimits(1);
%     midpoints(midpoints > xlimits(2)) = xlimits(2);
% end
% if ~isempty(ylimits)
%     cat_y(cat_y < ylimits(1)) = ylimits(1);
%     cat_y(cat_y > ylimits(2)) = ylimits(2);
%     binned_y(binned_y < ylimits(1)) = ylimits(1);
%     binned_y(binned_y > ylimits(2)) = ylimits(2);
%     SEM_binned_y(SEM_binned_y < ylimits(1)) = ylimits(1);
%     SEM_binned_y(SEM_binned_y > ylimits(2)) = ylimits(2);
%     std_binned_y(std_binned_y < ylimits(1)) = ylimits(1);
%     std_binned_y(std_binned_y > ylimits(2)) = ylimits(2);
% end

scatter(cat_x, cat_y, 'filled', 'MarkerFaceAlpha', 0.2);
errorbar(midpoints, binned_y, CI, 'LineStyle', 'none');
scatter(midpoints, binned_y, 'filled');
xlabel(xlab)
ylabel('Norm. |Delta zScore|')
xlim([0 0.4])
ylim([0 1])
set(gca, 'xscale', 'log');

%% ------------mean zScore vs redundancy-------------------%              
y = redundancy;          

cat_y = [];
for ii = 1:num_cells
    for ww = 1:num_waves(ii)
        cat_y = [cat_y y{ii}(ww)];               
    end
end

x = meanZScores;       
xlab = 'Event mean zScore';
cat_x = [];     
for ii = 1:num_cells
    for ww = 1:num_waves(ii)
        cat_x = [cat_x x{ii}(ww)];
    end
end

figure
axis square
hold on
edges = prctile(cat_x, 0:10:100);
binned_x = discretize(cat_x, edges);
CI = zeros(length(edges)-1, 2);
binned_y = zeros(1, length(edges)-1);
for bb = 1:length(edges)-1
    binned_y(bb) = nanmean(cat_y(binned_x == bb));
    SEM_binned_y(bb) = nanstd(cat_y(binned_x == bb))/sqrt(sum(binned_x == bb));
    std_binned_y(bb) = nanstd(cat_y(binned_x == bb));
    midpoints(bb) = (edges(bb+1) + edges(bb))/2;
    t(bb, :) = tinv([0.025  0.975],sum(binned_x == bb)-1);      % T-Score
    CI(bb, :) = t(bb, :)*SEM_binned_y(bb);          % not true CI, just want the +-
end

scatter(cat_x, cat_y, 'filled', 'MarkerFaceAlpha', 0.2);
errorbar(midpoints, binned_y, CI(:, 2), 'LineStyle', 'none');
scatter(midpoints, binned_y, 'filled');
xlabel(xlab)
ylabel('Event redundancy')
xlim([1 100])
ylim([0 0.35])
set(gca, 'xscale', 'log');

%% ------------------------------------------------------%

%basic ranksum test on top and bottom quartiles
top = cat_y(cat_x >= prctile(cat_x, 75));
bottom = cat_y(cat_x <= prctile(cat_x, 25));
[p, ~, stats] = ranksum(bottom, top);
figure
boxplot(cat(2, bottom, top), [zeros(1, length(bottom)) ones(1, length(top))]);
title(sprintf('p = %.2e', p))
xlabel(xlab);
ylabel('|delta zscore|');
ylim([-0.05 1])
set(gcf, 'Position', [400 400 200 500]);

%binning data by 1 percentiles, ranksum testing top and bottom quartiles of this binned data
% edges = prctile(cat_x, 0:1:100);
% binned_x = discretize(cat_x, edges);
% for bb = 1:length(edges)-1
%     binned_y_100(bb) = nanmean(cat_y(binned_x == bb));
% end
% bottom = binned_y_100(1:25);       %first 'quartile'
% top = binned_y_100(76:100);
% p = ranksum(top, bottom);
% figure
% boxplot(cat(2, bottom, top), [zeros(1, length(bottom)) ones(1, length(top))]);
% title(sprintf('p = %.2e', p))
% xlabel('Binned event avg zscore quartiles');
% ylabel('Norm. |delta zscore|');
% set(gcf, 'Position', [400 400 200 500]);
% ylim([0 0.2])

%1-way anova on the start time stuff for natmov
% deciles = [];
% groups = [];
% for bb = 1:length(edges)-1
%     curr_decile = cat_y(binned_x == bb);
%     groupsize = length(curr_decile);
%     deciles = [deciles curr_decile];
%     groups = [groups bb*ones(1, groupsize)]; 
% end
% [p, ~, stats] = anova1(deciles, groups);
% c = multcompare(stats, 'CType', 'bonferroni');

%comparing each bin to average with one sample ttest 
% chance = nanmean(cat_y);
% p = [];
% for bb = 1:max(binned_x)
%     curr_group = cat_y(binned_x == bb);
%     [~, p(bb)] = ttest(curr_group, chance);
% end

%save the 10th prctile binned averages from both stimuli to show on same plot for last subpanel       (mean zscore vs instability)

% for PDG
pdgplotdata.cat_x = cat_x;
pdgplotdata.cat_y = cat_y;
pdgplotdata.midpoints = midpoints;
pdgplotdata.binned_y = binned_y;
pdgplotdata.CI = CI;
save pdgplotdata pdgplotdata

% for MOV
natmovplotdata.cat_x = cat_x;
natmovplotdata.cat_y = cat_y;
natmovplotdata.midpoints = midpoints;
natmovplotdata.binned_y = binned_y;
natmovplotdata.CI = CI;
save natmovplotdata natmovplotdata

load pdgplotdata
load natmovplotdata

% plotting original binned data on same plot
figure

fa = polyfit(log(pdgplotdata.midpoints), pdgplotdata.binned_y, 1);
ypts = polyval(fa,log(pdgplotdata.midpoints));
semilogx(pdgplotdata.midpoints, ypts,'-')

hold on

fb = polyfit(log(natmovplotdata.midpoints), natmovplotdata.binned_y, 1);
ypts = polyval(fb,log(natmovplotdata.midpoints));
semilogx(natmovplotdata.midpoints, ypts,'-')

errorbar(pdgplotdata.midpoints, pdgplotdata.binned_y, pdgplotdata.CI(:, 2), 'LineStyle', 'none');
scatter(pdgplotdata.midpoints, pdgplotdata.binned_y, 'filled');
errorbar(natmovplotdata.midpoints, natmovplotdata.binned_y, natmovplotdata.CI(:, 2), 'LineStyle', 'none');
scatter(natmovplotdata.midpoints, natmovplotdata.binned_y, 'filled');

axis square
xlim([1 100])
ylim([0 0.5])
xlabel('Session-average event zScore (binned)')
ylabel('Norm |delta zScore|');

% plotting bins from combined event data from both stimuli (to make bins with equivalent boundaries between stimuli)
cat_x = [pdgplotdata.cat_x natmovplotdata.cat_x];
edges = prctile(cat_x, 0:10:100);                   % bin edges calculated from all data

figure

%PDG
cat_x = pdgplotdata.cat_x; cat_y = pdgplotdata.cat_y;
binned_x = discretize(cat_x, edges);            % binning done separately for each stimulus based on pooled edges
CI = zeros(1, length(edges)-1);
binned_y = zeros(1, length(edges)-1);
for bb = 1:length(edges)-1
    binned_y(bb) = nanmean(cat_y(binned_x == bb));
    SEM_binned_y(bb) = nanstd(cat_y(binned_x == bb))/sqrt(sum(binned_x == bb));
    std_binned_y(bb) = nanstd(cat_y(binned_x == bb));
    midpoints(bb) = (edges(bb+1) + edges(bb))/2;
    t(bb, :) = tinv([0.025  0.975],sum(binned_x == bb)-1);      % T-Score
    temp = t(bb, :)*SEM_binned_y(bb);          % not true CI, just want the +-
    CI(bb) = temp(2);
%     temp = bootci(2000, @(x)nanmean(x), cat_y(binned_x == bb)); 
%     CI(bb) = temp(2) - temp(1);
    eles(bb) = sum(binned_x == bb);
end

f = polyfit(log(midpoints), binned_y, 1);
ypts = polyval(f,log(midpoints));
semilogx(midpoints, ypts,'-')
hold on
errorbar(midpoints, binned_y, CI, 'LineStyle', 'none');
scatter(midpoints, binned_y, 'filled');

%MOV
cat_x = natmovplotdata.cat_x; cat_y = natmovplotdata.cat_y;
binned_x = discretize(cat_x, edges);            % binning done separately for each stimulus based on pooled edges
CI = zeros(1, length(edges)-1);
binned_y = zeros(1, length(edges)-1);
for bb = 1:length(edges)-1
    binned_y(bb) = nanmean(cat_y(binned_x == bb));
    SEM_binned_y(bb) = nanstd(cat_y(binned_x == bb))/sqrt(sum(binned_x == bb));
    std_binned_y(bb) = nanstd(cat_y(binned_x == bb));
    midpoints(bb) = (edges(bb+1) + edges(bb))/2;
    t(bb, :) = tinv([0.025  0.975],sum(binned_x == bb)-1);      % T-Score
    temp = t(bb, :)*SEM_binned_y(bb);          % not true CI, just want the +-
    CI(bb) = temp(2);
%     temp = bootci(2000, @(x)nanmean(x), cat_y(binned_x == bb)); 
%     CI(bb) = temp(2) - temp(1);
    eles(bb) = sum(binned_x == bb);
end

f = polyfit(log(midpoints), binned_y, 1);
ypts = polyval(f,log(midpoints));
semilogx(midpoints, ypts,'-')
errorbar(midpoints, binned_y, CI, 'LineStyle', 'none');
scatter(midpoints, binned_y, 'filled');

axis square
xlim([1 100])
ylim([0 0.5])
xlabel('Session-average event zScore (binned)')
ylabel('Norm |delta zScore|');

%stat comparisons of each bin
cat_x = pdgplotdata.cat_x; cat_y = pdgplotdata.cat_y;
binned_x = discretize(cat_x, edges);            % binning done separately for each stimulus based on pooled edges
bins_y_pdg = cell(1, length(edges)-1);
for bb = 1:length(edges)-1
    bins_y_pdg{bb} = cat_y(binned_x == bb);
end

cat_x = natmovplotdata.cat_x; cat_y = natmovplotdata.cat_y;
binned_x = discretize(cat_x, edges);            % binning done separately for each stimulus based on pooled edges
bins_y_natmov = cell(1, length(edges)-1);
for bb = 1:length(edges)-1
    bins_y_natmov{bb} = cat_y(binned_x == bb);
end

for bb = 1:length(bins_y_pdg)
    p(bb) = ranksum(bins_y_pdg{bb}, bins_y_natmov{bb});
end

%% bootstrap analysis of natmov vs pdg waveform mean zscore distributions
% natmov_wavemeans = cat_x;            % compute each of these first
% pdg_wavemeans = cat_x;
% 
% sample_size = length(pdg_wavemeans);
% iterations = 1000;
% ps = [];
% for tt = 1:iterations
%     tt
%     bootsample = datasample(natmov_wavemeans, sample_size);
%     [~, ps(tt)] = ttest2(bootsample, pdg_wavemeans);
% end
% p = mean(ps);
% 
% figure; hold on
% histogram(pdg_wavemeans, 0:2:100)
% histogram(natmov_wavemeans, 0:2:100)
% title(sprintf('p = %.2e', p));
% xlim([0 60])

%% Event independency testing (unused)
       
% First: run delta calculations above 

% Calcualte a "diversity score" for every neuron, within groups of neurons with different # of events (only consider neurons with 2+ events)
% diversity = zeros(1, num_cells);
% max_same = zeros(1, num_cells);
% num_events = cellfun(@length, deltas);
% 
% for ii = 1:num_cells    
%     num_up = sum(deltas{ii} == 1);
%     num_down = sum(deltas{ii} == -1);
%     num_static = sum(deltas{ii} == 0);
%     
%     max_same(ii) = max([num_up num_down num_static]);
% %     diversity(ii) = (1-max_same/num_events(ii))/(2/3);          % sets 0.333 (best score) to 1 and sets 1 (worst score) to 0
%     diversity(ii) = max_same(ii)/num_events(ii);
% end
% 
% % plotting
% cat_div = [];
% g = [];
% for yy = 3:8
%     if yy < 8
%         cat_div = cat(2, cat_div, diversity(num_events == yy));
%         g = cat(2, g, yy*ones(1, sum(num_events == yy)));
%     else
%         cat_div = cat(2, cat_div, diversity(num_events >= yy));
%         g = cat(2, g, yy*ones(1, sum(num_events >= yy)));
%     end
% end
% 
% histogram(diversity(num_events >= 3), 10)
% 
% figure
% violinplot(cat_div, g, 'ShowBox', false, 'ShowNotches', false, 'ShowMedian', false, 'ShowWhisker', false, 'ShowViolin', false);
% xlabel('Number of events')
% ylabel('Diversity score')
% 
% %visualizing proportions of diversity scores within each 'event class'
% figure                          % using diversity score
% for yy = 3:8
%     if yy < 9
%         curr_scores = diversity(num_events == yy);
%     else
%         curr_scores = diversity(num_events >= yy);
%     end
%     x = categorical(curr_scores);
%     
%     subplot(3, 2, yy-2)
%     pie(x)
%     if yy < 9
%         title(sprintf('Neurons with %d events', yy));
%     else
%         title(sprintf('Neurons with %d+ events', yy));
%     end
% end
% 
% figure                          % using max_same
% for yy = 3:8
%     if yy < 9
%         curr_scores = max_same(num_events == yy);
%     else
%         curr_scores = max_same(num_events >= yy);
%     end
%     x = categorical(curr_scores);
%     
%     subplot(2, 3, yy-2)
%     pie(x)
%     if yy < 9
%         title(sprintf('Neurons with %d events', yy));
%     else
%         title(sprintf('Neurons with %d+ events', yy));
%     end
% end
% 
% 
% figure
% histogram(num_events, 0.5:1:9.5)
% xlabel('Number of events')
% ylabel('Neurons')
% 
% 
% % Is diversity score dependent on anything?
% 
% %Average DFF zScore
% frame_catResp = [];
% for rr = 1:size(Resp, 1)
%     frame_catResp = cat(1, frame_catResp, squeeze(Resp(rr, :, :, :)));
% end
%     
% fullzScores = zeros(size(frame_catResp, 2), size(frame_catResp, 3));
% for kk = 1:all.num_sessions
%     fullzScores(:, kk) = nanmean(frame_catResp(:, :, kk), 1)./nanstd(frame_catResp(:, :, kk), [], 1);
% end
% fullzScores = nanmean(fullzScores, 2);
% 
% figure
% scatter(fullzScores(num_events >= 6), diversity(num_events >= 6), 'filled');
% [~, p] = corr(fullzScores(num_events >= 6), diversity(num_events >= 6)', 'Type', 'Spearman');
% title(sprintf('All neurons >= 6 events, Spearman p = %.2d', p));
% xlabel('Ses-average zScore of response')
% ylabel('Diversity score')
% axis square
% 
% %Mean spike rate
% natmov_misr = importdata('cat_natmov_misr.mat');
% natmov_misr = natmov_misr(1:length(good_cells));
% natmov_misr = natmov_misr(good_cells);
% 
% figure
% scatter(natmov_misr(num_events >= 6), diversity(num_events >= 6), 'filled');
% [~, p] = corr(natmov_misr(num_events >= 6), diversity(num_events >= 6)', 'Type', 'Spearman');
% title(sprintf('All neurons >= 6 events, Spearman p = %d', p));
% xlabel('Ses-average MISR')
% ylabel('Diversity score')
% axis square
% 
% %% Independency continued
% 
% % Find average pairwise correlations between event zscore curves within neurons (in neurons with n or more events)
% % compare to correlations resulting from using session-shuffled zscore data
% 
% %shuffling zscores
% shuffled_sesZScores = cell(1, num_cells);
% for ii = 1:num_cells
%     curr_waves = sesZScores{ii};
%     if ~isempty(curr_waves)
%         curr_sessions = length(curr_waves{1});
%         for ww = 1:num_waves(ii)
%             neworder = randperm(curr_sessions);
%             shuffled_sesZScores{ii}{ww} = curr_waves{ww}(neworder);
%         end
%     end
% end
% 
% 
% avg_cc = zeros(1, num_cells);
% avg_cc_shuffled = zeros(1, num_cells);
% for ii = 1:num_cells
%     curr_waves = sesZScores{ii};
%     curr_waves_shuffled = shuffled_sesZScores{ii};
%     pairwise_ccs = [];
%     pairwise_ccs_shuffled = [];
%     for ww = 1:num_waves(ii)
%         for mm = ww:num_waves(ii)
%             if ww ~= mm
%                 pairwise_ccs = [pairwise_ccs corr(curr_waves{ww}', curr_waves{mm}')];
%                 pairwise_ccs_shuffled = [pairwise_ccs_shuffled corr(curr_waves_shuffled{ww}', curr_waves_shuffled{mm}')];
%             end
%         end
%     end
%     avg_cc(ii) = nanmean(pairwise_ccs);
%     avg_cc_shuffled(ii) = nanmean(pairwise_ccs_shuffled);
% end
%                 
% figure
% histogram(avg_cc(num_waves >= 4))
% hold on
% histogram(avg_cc_shuffled(num_waves >= 4))
% title('neurons with 4 or more events')
% xlabel('Average pairwise cc of event zscore curves')
% ylabel('neurons');
% 
% 
% % average pairwise zscore cc versus a neuron's responsivity (ses-avg zscore or first session zscore)
% 
% bycell_sesZScores = zeros(num_cells, all.num_sessions);
% for kk = 1:all.num_sessions
%     curr_Resp = permute(Resp(:, :, :, kk), [2 3 4 1]);
%     framecat_Resp = squeeze(all.catTrials(curr_Resp));            % [total frames x cells]
%     if strcmp(stim, 'PDG')
%         exidx = logical(repmat([ones(1, 40) zeros(1, 20)], 1, 96));
%         framecat_Resp(exidx, :) = [];                        %excise offtime for PDG
%     end
%     bycell_sesZScores(:, kk) = nanmean(framecat_Resp, 1)./std(framecat_Resp, [], 1);
% end
% bycell_meanZScores = nanmean(bycell_sesZScores, 2);
% bycell_firstZScores = bycell_sesZScores(:, 1);
% 
% %
% figure
% subplot(1, 2, 1)
% x = bycell_meanZScores(num_waves >= 4);
% y = avg_cc(num_waves >= 4);
% scatter(x, y, 'filled');
% axis square
% xlabel('Ses-avg zscore')
% ylabel('avg pairwise cc of event zscores')
% [r, p] = corr(x, y', 'Type', 'Spearman');
% title(sprintf('neurons with 4 or more events,\n r = %.2f p = %.2e', r, p))
% hold on
% edges = prctile(x, 0:10:100);
% binned_x = discretize(x, edges);
% for bb = 1:length(edges)-1
%     binned_y(bb) = nanmean(y(binned_x == bb));
%     SEM_binned_y(bb) = nanstd(y(binned_x == bb))/sqrt(sum(binned_x == bb));
%     midpoints(bb) = (edges(bb+1) + edges(bb))/2;
% end
% scatter(midpoints, binned_y, 'filled');
% errorbar(midpoints, binned_y, SEM_binned_y, 'LineStyle', 'none');
% 
% subplot(1, 2, 2)
% bot_thresh = prctile(x, 10);
% top_thresh = prctile(x, 90);
% bot_vals = y(x <= bot_thresh);
% top_vals = y(x >= top_thresh);
% [~, p] = ttest2(bot_vals, top_vals);
% boxplot([bot_vals; top_vals]')
% title(sprintf('p = %.2e', p))
% xlabel('Bottom vs top 10th prctile')
% ylabel('avg ccs')
% axis square
% 
% %
% figure
% subplot(1, 2, 1)
% x = bycell_firstZScores(num_waves >= 4);
% y = avg_cc(num_waves >= 4);
% scatter(x, y, 'filled');
% axis square
% xlabel('First session zscore')
% ylabel('avg pairwise cc of event zscores')
% [r, p] = corr(x, y', 'Type', 'Spearman');
% title(sprintf('neurons with 4 or more events,\n r = %.2f p = %.2e', r, p))
% hold on
% edges = prctile(x, 0:10:100);
% binned_x = discretize(x, edges);
% for bb = 1:length(edges)-1
%     binned_y(bb) = nanmean(y(binned_x == bb));
%     SEM_binned_y(bb) = nanstd(y(binned_x == bb))/sqrt(sum(binned_x == bb));
%     midpoints(bb) = (edges(bb+1) + edges(bb))/2;
% end
% scatter(midpoints, binned_y, 'filled');
% errorbar(midpoints, binned_y, SEM_binned_y, 'LineStyle', 'none');
% 
% subplot(1, 2, 2)
% bot_thresh = prctile(x, 10);
% top_thresh = prctile(x, 90);
% bot_vals = y(x <= bot_thresh);
% top_vals = y(x >= top_thresh);
% [~, p] = ttest2(bot_vals, top_vals);
% boxplot([bot_vals; top_vals]')
% title(sprintf('p = %.2e', p))
% xlabel('Bottom vs top 10th prctile')
% ylabel('avg ccs')
% axis square

%% plotting events

for ii = 1:num_cells
    figure
    subplot(2, 1, 1)
    imagesc(catResp(:, :, ii));
    
    subplot(2, 1, 2)
    plot(mean(catResp(:, :, ii), 1));
    hold on
    
    for ww = 1:size(ts{ii}, 1)
        p = patch([ts{ii}(ww, 1) ts{ii}(ww, 1) ts{ii}(ww, 2) ts{ii}(ww, 2)], [min(mean(catResp(:, :, ii))) max(mean(catResp(:, :, ii))) max(mean(catResp(:, :, ii))) min(mean(catResp(:, :, ii)))], 'c');
        p.EdgeAlpha = 0;
        p.FaceAlpha = 0.4;
    end
    
%     subplot(3, 1, 3)
%     plot(event_vectors(ii, :));
%     ylim([-0.1 1.1])
    pause
    close
end

%plot change of activity over time for each event for a given neuron
neuron = 4;                 %with respect to good_cells
figure
trace = zeros(210, num_waves(neuron));
for ww = 1:num_waves(neuron)
    subplot(num_waves(neuron), 1, ww)
    curr_wave = waveforms{neuron}{ww};
    curr_wave_cattrials = [];
    offResp_cattrials = [];
    for kk = 1:size(curr_wave, 3)
        curr_wave_cattrials = cat(1, curr_wave_cattrials, curr_wave(:, :, kk));
        offResp_cattrials = cat(1, offResp_cattrials, offResp(:, :, neuron, kk));
    end
    curr_trace = mean(curr_wave_cattrials, 2)./std(offResp_cattrials, [], 2);
%     trace = [trace' repmat(trace(end), 1, 30)];          % zero padding on the edges
    plot(movmean(curr_trace, 30));
    trace(:, ww) = movmean(curr_trace, 30);
    ylim([0 11])
    xlim([0 210])
    axis square
    box off
end

ii = 4;
figure
for ii = 1:num_cells
    for ww = 1:num_waves(ii)
        curr_wave = squeeze(mean(waveforms{ii}{ww}, 1));
        for kk = 1:size(curr_wave, 2)
            subplot(size(curr_wave, 2), 1, kk)
            plot(curr_wave(:, kk))
            ylim([min(min(curr_wave)) max(max(curr_wave))])
        end
        pause
    end
end
            


