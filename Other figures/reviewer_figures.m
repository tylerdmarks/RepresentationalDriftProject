% REVIEWER FIGURES


%% age at first recording vs. final RDI
% set age, create analysis and import data for each mouse
age = 250;
an = StabilityAnalyzer();
an.importData();
% an.addData();                 % for mice with mulitple fields, combine into same object by adding subsequent fields
an.setQualityThreshold(3);


figure
an.RDIplotter('PDG', 'Average');         % running RDIplotter produces the RDIdata structure, which we can then extract
hold on
an.RDIplotter('NatMov', 'Average');

RDIdata = an.getRDIdata;
PDG_RDIavg = RDIdata.PDG.RDI_avg(end);
NatMov_RDIavg = RDIdata.NatMov.RDI_avg(end);

% initialize once
ages = [];
PDG_RDIs = [];
MOV_RDIs = [];

% find values for each mouse (above), then concatenate
ages = [ages age];
PDG_RDIs = [PDG_RDIs PDG_RDIavg];
MOV_RDIs = [MOV_RDIs NatMov_RDIavg];

% plot pooled data
figure; 
scatter(ages, PDG_RDIs, 60, 'filled');
hold on
scatter(ages, MOV_RDIs, 60, 'filled');
xlabel('Age (days)');
ylabel('Final RDI value')
[r_pdg, p_pdg] = corr(ages', PDG_RDIs');
[r_mov, p_mov] = corr(ages', MOV_RDIs'); 
title(sprintf('PDG r = %.2f, p = %.2f; MOV r = %.2f, p = %.2f', r_pdg, p_pdg, r_mov, p_mov));
axis square




%% Inhibitory neuron categorization

an = StabilityAnalyzer();
an.importData();
an.addData();
an.setQualityThreshold(3);
an.cellSelection;
qual = an.getUse_cells;
sessions = an.getUse_sessions;
num_cells = an.num_cells;
num_sessions = an.num_sessions;

sessions(isnan(sessions)) = 1;            % nan values converted to zeros for logical purposes
present = sum(sessions, 2) == an.num_sessions;
good_cells = qual' & present;


figure
an.RDIplotter('PDG', 'Average');
hold on
an.RDIplotter('NatMov', 'Average');
oridata = an.getOrientationData;
RDIdata = an.getRDIdata;

thresh = 0.4;
osi = oridata.osiMat;
avgosi = zeros(1, num_cells);
stdosi = zeros(1, num_cells);
for ii = 1:size(osi, 2)
    avgosi(ii) = mean(osi(logical(sessions(ii, :)), ii));
    stdosi(ii) = std(osi(logical(sessions(ii, :)), ii));
end
sharp = avgosi >= thresh;
broad = avgosi < thresh;

pdg_sharp = zeros(2, num_sessions-1);
pdg_broad = zeros(2, num_sessions-1);
mov_sharp = zeros(2, num_sessions-1);
mov_broad = zeros(2, num_sessions-1);
n_sharp = zeros(1, num_sessions-1);
n_broad = zeros(2, num_sessions-1);
for kk = 1:num_sessions-1
    pdg_sharp(1, kk) = nanmean(RDIdata.PDG.RDI_included_scatter(kk, sharp));
    pdg_sharp(2, kk) = nanstd(RDIdata.PDG.RDI_included_scatter(kk, sharp))/sqrt(sum(sharp & qual & sessions(:, kk+1)'));
    pdg_broad(1, kk) = nanmean(RDIdata.PDG.RDI_included_scatter(kk, broad));
    pdg_broad(2, kk) = nanstd(RDIdata.PDG.RDI_included_scatter(kk, broad))/sqrt(sum(broad & qual & sessions(:, kk+1)'));
    
    mov_sharp(1, kk) = nanmean(RDIdata.NatMov.RDI_included_scatter(kk, sharp));
    mov_sharp(2, kk) = nanstd(RDIdata.NatMov.RDI_included_scatter(kk, sharp))/sqrt(sum(sharp & qual & sessions(:, kk+1)'));
    mov_broad(1, kk) = nanmean(RDIdata.NatMov.RDI_included_scatter(kk, broad));
    mov_broad(2, kk) = nanstd(RDIdata.NatMov.RDI_included_scatter(kk, broad))/sqrt(sum(broad & qual & sessions(:, kk+1)'));
    n_sharp(kk) = sum(sharp & qual & sessions(:, kk+1)');
    n_broad(kk) = sum(broad & qual & sessions(:, kk+1)');
    
    curr_pdg_sharp = RDIdata.PDG.RDI_included_scatter(kk, sharp);
    curr_pdg_sharp(isnan(curr_pdg_sharp)) = [];
    curr_natmov_sharp = RDIdata.NatMov.RDI_included_scatter(kk, sharp);
    curr_natmov_sharp(isnan(curr_natmov_sharp)) = [];
    curr_pdg_broad = RDIdata.PDG.RDI_included_scatter(kk, broad);
    curr_pdg_broad(isnan(curr_pdg_broad)) = [];
    curr_natmov_broad = RDIdata.NatMov.RDI_included_scatter(kk, broad);
    curr_natmov_broad(isnan(curr_natmov_broad)) = [];
    p_pdg(kk) = ranksum(curr_pdg_sharp, curr_pdg_broad, 'tail', 'both');
    p_mov(kk) = ranksum(curr_natmov_sharp, curr_natmov_broad, 'tail', 'both');
end

figure
hold on
errorbar(1:num_sessions-1, pdg_sharp(1, :), pdg_sharp(2, :), '-b');
errorbar(1:num_sessions-1, pdg_broad(1, :), pdg_broad(2, :), '--c');
errorbar(1:num_sessions-1, mov_sharp(1, :), mov_sharp(2, :), '-r');
errorbar(1:num_sessions-1, mov_broad(1, :), mov_broad(2, :), '--m');
xlim([0 7])
ylim([-0.1 0.5])
axis square






