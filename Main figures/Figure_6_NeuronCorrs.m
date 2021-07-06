
% create Analyzer (data import is in constructor in this one)
SKKS066 = EnsembleAnalyzer();
% set quality threshold and assign inclusion criteria
SKKS066.setQualityThreshold(3);
SKKS066.cellSelection;
num_sessions = SKKS066.num_sessions;
num_cells = SKKS066.num_cells;

qual = SKKS066.getUse_cells;            % reflects included neurons based on above selections
sessions = SKKS066.getUse_sessions;     % presence data for each neuron
% consider neurons present on all sessions (some have NaN values where data was not collected, set to 1 for logical purposes)
sessions(isnan(sessions)) = 1;
present = sum(sessions, 2) == num_sessions;
good_cells = qual' & present;           % included neurons that are present on all sessions

% if using same neurons for both stimuli, these will be identical
% if not, run lines 5-15 separately for each stim, and save the good_cells vector separately for each stimulus
save good_cells_natmov good_cells
save good_cells_pdg good_cells


%% ------- Signal corrs -----------%

% Full matrices
Scorr_natmov = zeros(num_cells, num_cells, num_sessions);
Scorr_pdg = zeros(num_cells, num_cells, num_sessions);
for kk = 1:num_sessions
    fprintf('Session %d of %d\n', kk, num_sessions);
    Scorr_natmov(:, :, kk) = SKKS066.SignalCorr(SKKS066.RespData.NatMov.RespMat_Full(:, :, :, kk));
    Scorr_pdg(:, :, kk) = SKKS066.SignalCorr(SKKS066.RespData.PDG.RespMat_Full(:, :, :, kk));
end
save Scorr_natmov Scorr_natmov
save Scorr_pdg Scorr_pdg

% Within session control matrices
% Scorr_natmov_WS = zeros(num_cells, num_cells, num_sessions, 2);
% Scorr_pdg_WS = zeros(num_cells, num_cells, num_sessions, 2);
% 
% clear NatMov_RespMat
% clear PDG_RespMat
% for kk = 1:num_sessions
%     fprintf('Session %d of %d\n', kk, num_sessions);
%     NatMov_RespMat(:, :, :, 1) = SKKS066.RespData.NatMov.RespMat_Full(2:2:end, :, :, kk);
%     NatMov_RespMat(:, :, :, 2) = SKKS066.RespData.NatMov.RespMat_Full(1:2:end, :, :, kk);
%     PDG_RespMat(:, :, :, 1) = SKKS066.RespData.PDG.RespMat_Full(2:2:end, :, :, kk);
%     PDG_RespMat(:, :, :, 2) = SKKS066.RespData.PDG.RespMat_Full(1:2:end, :, :, kk);
%     
%     for eo = 1:2
%         Scorr_natmov_WS(:, :, kk, eo) = SKKS066.SignalCorr(NatMov_RespMat(:, :, :, eo));
%         Scorr_pdg_WS(:, :, kk, eo) = SKKS066.SignalCorr(PDG_RespMat(:, :, :, eo));
%     end
% end
% 
% save Scorr_natmov_WS Scorr_natmov_WS
% save Scorr_pdg_WS Scorr_pdg_WS

%% ------- Noise corrs 201210 ---------%

% Full matrices
Ncorr_natmov = zeros(sum(good_cells), sum(good_cells), num_sessions);           % need to limit to only good_cells to save time
Ncorr_pdg = zeros(sum(good_cells), sum(good_cells), num_sessions);

% matrices for within-session comparisons (D0 only)
% Ncorr_natmov_WS = zeros(sum(good_cells), sum(good_cells), num_sessions, 2);               % going to yield two different versions, one using even bins, other using odd bins
% Ncorr_pdg_WS = zeros(sum(good_cells), sum(good_cells), num_sessions, 2);

frames = 5760;      % to match amount of data used for corr calculation between stimuli
for kk = 1:num_sessions
    fprintf('Session %d out of %d\n', kk, num_sessions);
    PDG = SKKS066.RespData.PDG.RespMat_Full(:, :, :, kk);               
    NatMov = SKKS066.RespData.NatMov.RespMat_Full(:, :, :, kk);
    
    fprintf('PDG...\n');
    Ncorr_pdg(:, :, kk) = SKKS066.NoiseCorr(PDG, frames);       % the NoiseCorr method only performs calculations on good_cells, as established above.
%     Ncorr_pdg_WS(:, :, kk, 1) = MJG013w1.NoiseCorr(PDG(2:2:end, :, :), frames);   %  squeeze(mean(PDG, 1))
%     Ncorr_pdg_WS(:, :, kk, 2) = MJG013w1.NoiseCorr(PDG(1:2:end, :, :), frames);

    
    fprintf('NatMov...\n');
    Ncorr_natmov(:, :, kk) = SKKS066.NoiseCorr(NatMov, frames);
%     Ncorr_natmov_WS(:, :, kk, 1) = MJG013w1.NoiseCorr(NatMov(2:2:end, :, :), frames);   %squeeze(mean(NatMov, 1))
%     Ncorr_natmov_WS(:, :, kk, 2) = MJG013w1.NoiseCorr(NatMov(1:2:end, :, :), frames); 
end

save Ncorr_natmov Ncorr_natmov
save Ncorr_pdg Ncorr_pdg
% save Ncorr_natmov_WS Ncorr_natmov_WS
% save Ncorr_pdg_WS Ncorr_pdg_WS

%% Visualization

% load either Scorrs or Ncorrs
corr_natmov = importdata('Ncorr_natmov.mat');
corr_pdg = importdata('Ncorr_pdg.mat');

% rectify identity line to 0 for visualization
num_cells = size(corr_natmov, 1);
num_sessions = size(corr_natmov, 3);
I = logical(eye(num_cells));
for kk = 1:num_sessions                         
    natmov_curr = corr_natmov(:, :, kk);
    natmov_curr(I) = 0;
    corr_natmov(:, :, kk) = natmov_curr;
    
    pdg_curr = corr_pdg(:, :, kk);
    pdg_curr(I) = 0;
    corr_pdg(:, :, kk) = pdg_curr;
end

% extract only good_cells -- only necessary for Scorrs, since Ncorrs are output with good_cells only already
load good_cells_natmov
corr_natmov = corr_natmov(good_cells, good_cells, :);
load good_cells_pdg
corr_pdg = corr_pdg(good_cells, good_cells, :);

% sort
PDGResp = SKKS066.RespData.PDG.RespMat_Full(:, :, good_cells, :);
NatResp = SKKS066.RespData.NatMov.RespMat_Full(:, :, good_cells, :);
PDGcatResp = [];
NatcatResp = [];
for kk = 1:SKKS066.num_sessions
    PDGcatResp = cat(1, PDGcatResp, PDGResp(:, :, :, kk));
    NatcatResp = cat(1, NatcatResp, NatResp(:, :, :, kk));
end
[~, pdgmaxidcs] = max(squeeze(mean(PDGcatResp, 1)), [], 1);
[~, natmaxidcs] = max(squeeze(mean(NatcatResp, 1)), [], 1);

[~, pdgorder] = sort(pdgmaxidcs);
[~, natorder] = sort(natmaxidcs);

corr_pdg_sorted = corr_pdg(pdgorder, pdgorder, :);
corr_natmov_sorted = corr_natmov(natorder, natorder, :);


%%plotting

% A: plot corr matrices D0 and D42 for both stimuli
coloraxis = ([-0.2 0.2]);
figure
colormap(redblue);
subplot(2, 2, 1)
imagesc(corr_pdg_sorted(:, :, 1))
axis square
caxis(coloraxis)
colorbar
xlabel('D0')
ylabel('D0')
subplot(2, 2, 2)
imagesc(corr_pdg_sorted(:, :, num_sessions))
axis square
caxis(coloraxis)
colorbar
xlabel('D final')
ylabel('D final')
subplot(2, 2, 3)
imagesc(corr_natmov_sorted(:, :, 1))
axis square
caxis(coloraxis)
colorbar
xlabel('D0')
ylabel('D0')
subplot(2, 2, 4)
imagesc(corr_natmov_sorted(:, :, num_sessions))
axis square
caxis(coloraxis)
colorbar
xlabel('D final')
ylabel('D final')

% difference matrix between D0 and D42 for both stimuli
coloraxis = ([-0.2 0.2]);
figure
colormap(redblue);
subplot(1, 2, 1)
imagesc(corr_pdg_sorted(:, :, num_sessions) - corr_pdg_sorted(:, :, 1));
axis square
caxis(coloraxis)
colorbar
subplot(1, 2, 2)
imagesc(corr_natmov_sorted(:, :, num_sessions) - corr_natmov_sorted(:, :, 1));
axis square
caxis(coloraxis)
colorbar


%% changes in connectivity magnitude

%Mean connectivity
pdgdiff = zeros(1, sum(good_cells));
natdiff = zeros(1, sum(good_cells));
for ii = 1:sum(good_cells)
    pdgdiff(ii) = mean(corr_pdg(1:end~=ii, ii, num_sessions), 1) - mean(corr_pdg(1:end~=ii, ii, 1), 1);         % D42 mean conn. - D0 mean conn.
    natdiff(ii) = mean(corr_natmov(1:end~=ii, ii, num_sessions), 1) - mean(corr_natmov(1:end~=ii, ii, 1), 1);
end

% % B: Example field - distribution of neurons' D42 mean conn. - D0 mean conn., PDG vs NatMov
figure
histogram(abs(pdgdiff), 0:0.005:0.08);            % bins here are for Scorrs
hold on
histogram(abs(natdiff), 0:0.005:0.08);
xline(mean(abs(pdgdiff)));
xline(mean(abs(natdiff)));
xlabel('|Delta mean conn.|')
ylabel('Neurons')
[p, ~, stats] = ranksum(abs(pdgdiff), abs(natdiff));
title(sprintf('|D42 mean conn. - D0 mean conn.|, PDG v NatMov, p = %.2e', p))

% % C: pooling data for all fields - distribution of fields' mean [D42 mean conn. - D0 mean conn.] for PDG vs NatMov
% initialize these vectors once
pdg_meanconn_absdiff_cat = [];
nat_meanconn_absdiff_cat = [];

% % add each mouse 
pdg_meanconn_absdiff_cat = [pdg_meanconn_absdiff_cat mean(abs(pdgdiff))];          %mean of the distribution of connectivity difference abs values for each of the stimuli
nat_meanconn_absdiff_cat = [nat_meanconn_absdiff_cat mean(abs(natdiff))];

save pdg_meanconn_absdiff_cat pdg_meanconn_absdiff_cat
save nat_meanconn_absdiff_cat nat_meanconn_absdiff_cat


%% Correlating corr matrices

% CCbs

%import either Scorr or Ncorr
corr_natmov = importdata('Ncorr_natmov.mat');
corr_pdg = importdata('Ncorr_pdg.mat');

%rectify identity line to 0 
num_cells = size(corr_natmov, 1);
num_sessions = size(corr_natmov, 3);
I = logical(eye(num_cells));
for kk = 1:num_sessions                         
    natmov_curr = corr_natmov(:, :, kk);
    natmov_curr(I) = 0;
    corr_natmov(:, :, kk) = natmov_curr;
    
    pdg_curr = corr_pdg(:, :, kk);
    pdg_curr(I) = 0;
    corr_pdg(:, :, kk) = pdg_curr;
end

%extract only good_cells if necessary
load good_cells_natmov
corr_natmov = corr_natmov(good_cells, good_cells, :);
load good_cells_pdg
corr_pdg = corr_pdg(good_cells, good_cells, :);

nat_btwses_cc = zeros(1, num_sessions-1);
pdg_btwses_cc = zeros(1, num_sessions-1);

%2D correlation
for kk = 2:num_sessions
    nat_btwses_cc(kk-1) = corr2(corr_natmov(:, :, kk), corr_natmov(:, :, 1));
    pdg_btwses_cc(kk-1) = corr2(corr_pdg(:, :, kk), corr_pdg(:, :, 1));
end
% gives the same result
% for kk = 2:num_sessions
%     nat_btwses_cc(kk-1) = max(max(normxcorr2(corr_natmov(:, :, kk), corr_natmov(:, :, 1))));
%     pdg_btwses_cc(kk-1) = max(max(normxcorr2(corr_pdg(:, :, kk), corr_pdg(:, :, 1))));
% end


% CCws

% %import either Scorr or Ncorr
% corr_natmov = importdata('Ncorr_natmov_WS.mat');
% corr_pdg = importdata('Ncorr_pdg_WS.mat');
% 
% %rectify identity line to 0 for visualization
% num_cells = size(corr_natmov, 1);
% num_sessions = size(corr_natmov, 3);
% I = logical(eye(num_cells));
% for eo = 1:2
%     for kk = 1:num_sessions                         
%         natmov_curr = corr_natmov(:, :, kk, eo);
%         natmov_curr(I) = 0;
%         corr_natmov(:, :, kk, eo) = natmov_curr;
% 
%         pdg_curr = corr_pdg(:, :, kk, eo);
%         pdg_curr(I) = 0;
%         corr_pdg(:, :, kk, eo) = pdg_curr;
%     end
% end
% 
% %extract only good_cells -- IF NECESSARY
% % load good_cells_natmov
% % corr_natmov = corr_natmov(good_cells, good_cells, :, :);
% % load good_cells_pdg
% % corr_pdg = corr_pdg(good_cells, good_cells, :, :);
% 
% nat_winses_cc = zeros(1, num_sessions);
% pdg_winses_cc = zeros(1, num_sessions);
% 
% %2D correlation
% for kk = 1:num_sessions
%     nat_winses_cc(kk) = corr2(corr_natmov(:, :, kk, 1), corr_natmov(:, :, kk, 2));
%     pdg_winses_cc(kk) = corr2(corr_pdg(:, :, kk, 1), corr_pdg(:, :, kk, 2));
% end


%% pooling data across fields

% CC bs
% initialize these cell arrays once
nat_btwses_cc_cat = cell(1, 6);
pdg_btwses_cc_cat = cell(1, 6);

% run this for each mouse after calculating its values above
weeks = [1 1 1 1 1 1];                  % starting on week 2, indicates weeks for which this mouse has data. e.g. [1 1 1 1 0 0] means mouse was recorded for only 5 weeks.
if ~isempty('weeks')
    num_sessions = length(weeks);
    ct = 1;
    for kk = 1:num_sessions
        if weeks(kk)
            nat_btwses_cc_cat{kk} = [nat_btwses_cc_cat{kk} nat_btwses_cc(ct)];
            pdg_btwses_cc_cat{kk} = [pdg_btwses_cc_cat{kk} pdg_btwses_cc(ct)];
            ct = ct + 1;
        end
    end
end

% save data once all mice have been added
save nat_btwses_cc_cat nat_btwses_cc_cat
save pdg_btwses_cc_cat pdg_btwses_cc_cat

% % CC ws
% % 
% % nat_winses_cc_cat = cell(1, 7);
% % pdg_winses_cc_cat = cell(1, 7);
% 
% weeks = [1 1 1 1 1 1 1];                  % starting on week 1
% if ~isempty('weeks')
%     num_sessions = length(weeks);
%     ct = 1;
%     for kk = 1:num_sessions
%         if weeks(kk)
%             nat_winses_cc_cat{kk} = [nat_winses_cc_cat{kk} nat_winses_cc(ct)];
%             pdg_winses_cc_cat{kk} = [pdg_winses_cc_cat{kk} pdg_winses_cc(ct)];
%             ct = ct + 1;
%         end
%     end
% end
% % 
% save nat_winses_cc_cat nat_winses_cc_cat
% save pdg_winses_cc_cat pdg_winses_cc_cat



%% plotting pooled data

% C
pdgabs = importdata('pdg_meanconn_absdiff_cat.mat');
natabs = importdata('nat_meanconn_absdiff_cat.mat');

figure
hold on
for kk = 1:length(pdgabs)
    plot(1:2, [pdgabs(kk) natabs(kk)], 'k');
    scatter(1, pdgabs(kk), 'filled', 'MarkerFaceColor', 'b');
    scatter(2, natabs(kk), 'filled', 'MarkerFaceColor', 'r');
end
[~, p, ~, stats] = ttest(pdgabs, natabs);
xlabel('Stimulus')
ylabel('Mean |D42 mean conn. - D0 mean conn.|')
title(sprintf('p = %.2e', p))
ylim([0 0.12])
xlim([0 3])
set(gcf, 'Position', [200 200 250 500])

% D
metric = 'cc';              % cc or mse
switch metric
    case 'cc'
        nat = importdata('nat_btwses_cc_cat.mat');
        pdg = importdata('pdg_btwses_cc_cat.mat');
%     case 'mse'     
%         nat = importdata('nat_btwses_mse_cat.mat');
%         pdg = importdatat('pdg_btwses_mse_cat.mat');
end

nat = cellfun(@(x)1-x, nat, 'UniformOutput', false);
pdg = cellfun(@(x)1-x, pdg, 'UniformOutput', false);
mean_nat = cellfun(@mean, nat);
mean_pdg = cellfun(@mean, pdg);
SEM_nat = cellfun(@std, nat)./sqrt(cellfun(@length, nat));
SEM_pdg = cellfun(@std, pdg)./sqrt(cellfun(@length, pdg));

for kk = 1:length(nat)
    [~, p(kk), ~, stats(kk)] = ttest(nat{kk}, pdg{kk});
end

figure
hold on
errorbar(1:6, mean_pdg, SEM_pdg, 'b');
errorbar(1:6, mean_nat, SEM_nat, 'r');
xlim([0 7])
ylim([0.5 1])
axis square
title(sprintf('p = %.2e ', p))






