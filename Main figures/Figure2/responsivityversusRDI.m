all = StabilityAnalyzer;
all.importData();               % import first field
all.addData();                  % add subsequent fields
all.setQualityThreshold(3);
all.cellSelection();


corr_type = 'Pearson';
stimtype = 'NatMov';
timespan = 'average';

binsize = 10;        % percent


%% responsivity vs RDI
close all
[plotdata, stats] = all.responsivity('responsivity', stimtype, timespan, binsize);   

ylim([-0.12 0.8])
xlim([0 1])

% %% reliability vs RDI
% 
% %load a StabilityAnalyzer with StabilityData files containing CCs from subsampled NatMov
% close all
% an.responsivity('reliability', stimtype, timespan, 10);
% 
% %% pooling spike data from first sessions (might pool average spike rate over all sessions later)
% natmov_filename = uigetfile('.mat');
% pdg_filename = uigetfile('.mat');
% natmovdata = importdata(natmov_filename);
% pdgdata = importdata(pdg_filename);
% 
% natmov_spikes = natmovdata.spikes(:, 1:10500);
% pdg_spikes = pdgdata.spikes(:, 1:5760);        % get rid of extra frames at the end
% 
% % natmov_exidx = logical(repmat([ones(1, 50) zeros(1, 300)], 1, 30));       % for excising out offtimes (probably not necessary)
% % pdg_exidx = logical(repmat([ones(1, 40) zeros(1, 20)], 1, 96));
% % 
% % natmov_spikes(:, natmov_exidx) = [];
% % pdg_spikes(:, pdg_exidx) = [];
% 
% mean_natmov_isr = mean(natmov_spikes, 2);
% mean_pdg_isr = mean(pdg_spikes, 2);
% 
% 
% % cat_natmov_misr = [];
% % cat_pdg_misr = [];
% 
% cat_natmov_misr = cat(1, cat_natmov_misr, mean_natmov_isr);
% cat_pdg_misr = cat(1, cat_pdg_misr, mean_pdg_isr);
% 
% save cat_natmov_misr cat_natmov_misr
% save cat_pdg_misr cat_pdg_misr
% 
% %% pooling spike data from all sessions (average)
% natmov_filenames = uigetfile('.mat', 'MultiSelect', 'on');
% pdg_filenames = uigetfile('.mat', 'MultiSelect', 'on');
% 
% clear natmov_spikes pdg_spikes
% for kk = 1:length(natmov_filenames)
%     kk
%     curr_data = importdata(natmov_filenames{kk});
%     natmov_spikes(:, :, kk) = curr_data.spikes(:, 1:10500);
%     curr_data = importdata(pdg_filenames{kk});
%     pdg_spikes(:, :, kk) = curr_data.spikes(:, 1:5760);
% end
% 
% mean_natmov_isr = mean(mean(natmov_spikes, 2), 3);
% mean_pdg_isr = mean(mean(pdg_spikes, 2), 3);
% 
% % 
% % cat_natmov_misr = [];
% % cat_pdg_misr = [];
% 
% cat_natmov_misr = cat(1, cat_natmov_misr, mean_natmov_isr);
% cat_pdg_misr = cat(1, cat_pdg_misr, mean_pdg_isr);
% 
% save cat_natmov_misr cat_natmov_misr
% save cat_pdg_misr cat_pdg_misr
% 
% 
% 
% %% cell selection
% an.cellSelection;                      %change this if you want to include only responsive cells, all cells, etc.
% qual = an.getUse_cells;            
% sessions = an.getUse_sessions;
% sessions(isnan(sessions)) = 1;
% present = sum(sessions, 2) == 7;
% good_cells = qual' & present;
% 
% %% spike rate distributions
% figure
% histogram(cat_pdg_misr(good_cells), 0:0.15:4.5);
% hold on
% histogram(cat_natmov_misr(good_cells), 0:0.15:4.5);
% % axis square
% [~, p] = ttest2(cat_natmov_misr(good_cells), cat_pdg_misr(good_cells));
% title(sprintf('natmov mean: %.2f, pdg mean: %.2f, p = %.2e', mean(cat_natmov_misr(good_cells)), mean(cat_pdg_misr(good_cells)), p));
% 
% % figure
% % boxplot([cat_natmov_misr(good_cells) cat_pdg_misr(good_cells)], [zeros(1, sum(good_cells)) ones(1, sum(good_cells))]);
% 
% figure
% scatter(cat_pdg_misr(good_cells), cat_natmov_misr(good_cells), 'filled', 'MarkerFaceAlpha', 0.4);
% [r, p] = corr(cat_pdg_misr(good_cells), cat_natmov_misr(good_cells));
% title(sprintf('r = %.2f, p = %.2e', r, p))
% hold on
% scatter(mean(cat_pdg_misr(good_cells)), mean(cat_natmov_misr(good_cells)), 'filled', 'MarkerFaceColor', 'r');
% axis square
% xlabel('PDG Mean ISR (spikes/s)');
% ylabel('NatMov Mean ISR (spikes/s)');
% set(gca, 'xscale', 'log', 'yscale', 'log');
% refline(1, 0)
% 
% 
% % plot RDI curves
% an.RDIplotter('PDG', 'Average');
% hold on
% an.RDIplotter('NatMov', 'Average');
% ylim([-0.2 1])
% 
% %% ISR vs reliability 
% % pdg_ccs = all.StabilityData.PDG.CCs(1, :);          % first session's values
% % natmov_ccs = all.StabilityData.NatMov.CCs(1, :);
% 
% pdg_ccs = nanmean(an.StabilityData.PDG.CCs(:, good_cells), 1);          % average values
% natmov_ccs = nanmean(an.StabilityData.NatMov.CCs(:, good_cells), 1);
% 
% %
% figure
% scatter(cat_pdg_misr(good_cells), pdg_ccs, 'filled');
% axis square
% hold on
% edges = prctile(cat_pdg_misr(good_cells), 0:binsize:100);
% binned_misr = discretize(cat_pdg_misr(good_cells), edges);
% for ee = 1:length(edges)-1
%     binned_ccs(ee) = mean(pdg_ccs(binned_misr == ee));
%     binned_cc_sem(ee) = std(pdg_ccs(binned_misr == ee))/sqrt(sum(binned_misr == ee));
%     midpoints(ee) = (edges(ee)+edges(ee+1))/2;
% end
% errorbar(midpoints, binned_ccs, binned_cc_sem, 'vertical', 'LineStyle', 'none');
% scatter(midpoints, binned_ccs, 'filled');
% % f = polyfit(midpoints, binned_ccs, 1);
% % refline(f(1), f(2));
% [r_binned, p_binned] = corr(midpoints', binned_ccs', 'Type', corr_type);
% % [r, p] = corr(cat_pdg_misr(good_cells), pdg_ccs(good_cells)', 'Type', corr_type);
% xlabel('PDG Mean ISR (spikes/s)')
% ylabel('PDG Reliability')
% title(sprintf('r = %.2f, p = %.2e', r_binned, p_binned))
% 
% %
% figure
% scatter(cat_natmov_misr(good_cells), natmov_ccs, 'filled');
% axis square
% hold on
% edges = prctile(cat_natmov_misr(good_cells), 0:binsize:100);
% binned_misr = discretize(cat_natmov_misr(good_cells), edges);
% for ee = 1:length(edges)-1
%     binned_ccs(ee) = mean(natmov_ccs(binned_misr == ee));
%     binned_cc_sem(ee) = std(natmov_ccs(binned_misr == ee))/sqrt(sum(binned_misr == ee));
%     midpoints(ee) = (edges(ee)+edges(ee+1))/2;
% end
% errorbar(midpoints, binned_ccs, binned_cc_sem, 'vertical', 'LineStyle', 'none');
% scatter(midpoints, binned_ccs, 'filled');
% % f = polyfit(midpoints, binned_ccs, 1);
% % refline(f(1), f(2));
% [r_binned, p_binned] = corr(midpoints', binned_ccs', 'Type', corr_type);
% % [r, p] = corr(cat_natmov_misr(good_cells), natmov_ccs(good_cells)', 'Type', corr_type);
% xlabel('NatMov Mean ISR (spikes/s)')
% ylabel('NatMov Reliability')
% title(sprintf('r = %.2f, p = %.2e', r_binned, p_binned))
% 
% %% ISR vs RDI
% % pdg_mean_RDI = zeros(1, length(good_cells));
% % natmov_mean_RDI = zeros(1, length(good_cells));
% % 
% % pdg_RDI = all.StabilityData.PDG.RDI;
% % natmov_RDI = all.StabilityData.NatMov.RDI;
% % 
% % for ii = 1:length(good_cells)
% %     pdg_mean_RDI(ii) = nanmean(pdg_RDI(logical(sessions(ii, 2:end)), ii), 1);
% %     natmov_mean_RDI(ii) = nanmean(natmov_RDI(logical(sessions(ii, 2:end)), ii), 1);
% % end
% 
% pdg_mean_RDI = nanmean(an.StabilityData.PDG.RDI(:, good_cells), 1);
% natmov_mean_RDI = nanmean(an.StabilityData.NatMov.RDI(:, good_cells), 1);
% 
% %
% figure
% x = cat_pdg_misr(good_cells);
% scatter(x, pdg_mean_RDI, 'filled', 'MarkerFaceAlpha', 0.4);
% axis square
% [r, p] = corr(x, pdg_mean_RDI', 'Type', corr_type);
% xlabel('PDG mean ISR (spikes/s)')
% ylabel('PDG mean RDI')
% hold on
% 
% edges = prctile(x, 0:binsize:100);
% binned_misr = discretize(x, edges);
% binned_pdg_RDI = zeros(1, length(edges)-1);
% binned_pdg_RDI_SEM = zeros(1, length(edges)-1);
% pdg_midpoints = zeros(1, length(edges)-1);
% for ee = 1:length(edges)-1
%     binned_pdg_RDI(ee) = mean(pdg_mean_RDI(binned_misr == ee));
%     binned_pdg_RDI_SEM(ee) = std(pdg_mean_RDI(binned_misr == ee))/sqrt(sum(binned_misr == ee));
%     pdg_midpoints(ee) = (edges(ee)+edges(ee+1))/2;
% end
% 
% mdl = fitlm(pdg_midpoints', binned_pdg_RDI');            % get confidence intervals for plotting later
% Xnew = linspace(min(pdg_midpoints)-1, max(pdg_midpoints)+1, 1000)';
% [~, CI] = predict(mdl, Xnew);
% 
% f = polyfit(pdg_midpoints, binned_pdg_RDI, 1);
% 
% %basic test
% bot_thresh = prctile(x, 25);
% top_thresh = prctile(x, 75);
% bot_vals = pdg_mean_RDI(x <= bot_thresh);
% top_vals = pdg_mean_RDI(x >= top_thresh);
% [~, p] = ttest2(bot_vals, top_vals);
% 
% %bootstrapped test
% bot_thresh = prctile(x, 25);
% top_thresh = prctile(x, 75);
% bot_vals = pdg_mean_RDI(x <= bot_thresh);
% top_vals = pdg_mean_RDI(x >= top_thresh);
% 
% sample_size = length(bot_vals);
% iterations = 1000;
% ps = [];
% for tt = 1:iterations
%     bootsample = datasample(bot_vals, sample_size);
%     ps(tt) = ranksum(bootsample, top_vals);
% end
% p = mean(ps);
% 
% %
% errorbar(pdg_midpoints, binned_pdg_RDI, binned_pdg_RDI_SEM, 'vertical', 'LineStyle', 'none');
% scatter(pdg_midpoints, binned_pdg_RDI, 'filled');
% [r_binned, p_binned] = corr(pdg_midpoints', binned_pdg_RDI', 'Type', corr_type);
% title(sprintf('r = %.2f, p = %.2e', r_binned, p_binned))
% % plot(Xnew, CI(:, 1));
% % plot(Xnew, CI(:, 2));
% xlim([min(cat_pdg_misr(good_cells)) max(cat_pdg_misr(good_cells))])
% set(gca, 'xscale', 'log')
% refline(f(1), f(2))
% ylim([-0.1 0.8])
% % f = @(b,a) exp(b(1)) .* exp(b(2).*a) + b(3);                % Exponential Fit With Y-Offset
% % B = fminsearch(@(b) norm(binned_pdg_RDI - f(b, pdg_midpoints)), [1; 1; 1]);         % Estimate Parameters
% % plot(pdg_midpoints, f(B, binned_pdg_RDI), '-r')
% 
% 
% 
% 
% figure
% boxplot([bot_vals; top_vals]')
% title(sprintf('p = %.2e', p))
% xlabel('Bottom vs top 25th prctile')
% ylabel('avg RDI')
% set(gcf, 'Position', [400 400 200 500]);
% ylim([-0.12 0.8])
% 
% % figure
% % hold on
% % errorbar(pdg_midpoints, binned_pdg_RDI, binned_pdg_RDI_SEM, 'vertical', 'LineStyle', 'none');
% % scatter(pdg_midpoints, binned_pdg_RDI, 'filled');
% % f = polyfit(pdg_midpoints, binned_pdg_RDI, 1);
% % refline(f(1), f(2));
% % 
% % mdl = fitlm(pdg_midpoints', binned_pdg_RDI');
% % Xnew = linspace(-0.2, 1, 1000)';
% % [~, CI] = predict(mdl, Xnew);
% % plot(Xnew, CI(:, 1));
% % plot(Xnew, CI(:, 2));
% % xlabel('PDG mean ISR, binned (spikes/s)')
% % ylabel('PDG mean RDI')
% 
% %
% figure
% x = cat_natmov_misr(good_cells);
% scatter(x, natmov_mean_RDI, 'filled', 'MarkerFaceAlpha', 0.4);
% axis square
% [r, p] = corr(x, natmov_mean_RDI', 'Type', corr_type);
% xlabel('NatMov mean ISR (spikes/s)')
% ylabel('NatMov mean RDI')
% hold on
% 
% edges = prctile(x, 0:binsize:100);
% binned_misr = discretize(x, edges);
% binned_natmov_RDI = zeros(1, length(edges)-1);
% binned_natmov_RDI_SEM = zeros(1, length(edges)-1);
% natmov_midpoints = zeros(1, length(edges)-1);
% for ee = 1:length(edges)-1
%     binned_natmov_RDI(ee) = mean(natmov_mean_RDI(binned_misr == ee));
%     binned_natmov_RDI_SEM(ee) = std(natmov_mean_RDI(binned_misr == ee))/sqrt(sum(binned_misr == ee));
%     natmov_midpoints(ee) = (edges(ee)+edges(ee+1))/2;
% end
% 
% errorbar(natmov_midpoints, binned_natmov_RDI, binned_natmov_RDI_SEM, 'vertical', 'LineStyle', 'none');
% scatter(natmov_midpoints, binned_natmov_RDI, 'filled');
% f = polyfit(natmov_midpoints, binned_natmov_RDI, 1);
% [r_binned, p_binned] = corr(natmov_midpoints', binned_natmov_RDI', 'Type', corr_type);
% title(sprintf('r = %.2f, p = %.2e', r_binned, p_binned))
% set(gca, 'xscale', 'log')
% refline(f(1), f(2))
% xlim([min(x) max(x)])
% ylim([-0.1 0.8])
% 
% %bootstrapped test
% bot_thresh = prctile(x, 25);
% top_thresh = prctile(x, 75);
% bot_vals = natmov_mean_RDI(x <= bot_thresh);
% top_vals = natmov_mean_RDI(x >= top_thresh);
% 
% sample_size = length(bot_vals);
% iterations = 1000;
% ps = [];
% for tt = 1:iterations
%     bootsample = datasample(bot_vals, sample_size);
%     ps(tt) = ranksum(bootsample, top_vals);
% end
% p = mean(ps);
% 
% 
% %basic test
% figure
% bot_thresh = prctile(x, 25);
% top_thresh = prctile(x, 75);
% bot_vals = natmov_mean_RDI(x <= bot_thresh);
% top_vals = natmov_mean_RDI(x >= top_thresh);
% [~, p] = ttest2(bot_vals, top_vals);
% boxplot([bot_vals; top_vals]')
% title(sprintf('p = %.2e', p))
% xlabel('Bottom vs top 25th prctile')
% ylabel('avg RDI')
% set(gcf, 'Position', [400 400 200 500]);
% ylim([-0.12 0.8])




% figure
% hold on
% errorbar(natmov_midpoints, binned_natmov_RDI, binned_natmov_RDI_SEM, 'vertical', 'LineStyle', 'none');
% scatter(natmov_midpoints, binned_natmov_RDI, 'filled');
% f = polyfit(pdg_midpoints, binned_natmov_RDI, 1);
% 
% mdl = fitlm(natmov_midpoints', binned_natmov_RDI');
% Xnew = linspace(-0.5, 4, 1000)';
% [~, CI] = predict(mdl, Xnew);
% plot(Xnew, CI(:, 1));
% plot(Xnew, CI(:, 2));
% xlabel('NatMov mean ISR, binned (spikes/s)')
% ylabel('NatMov mean RDI')
% refline(f(1), f(2));



