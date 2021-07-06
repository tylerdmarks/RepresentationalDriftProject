% initialize Analyzer
all = StabilityAnalyzer();
% First mouse: import data (RespData, RoiINFO, StabilityData)
all.importData();
% For subsequent mice, add data
all.addData();
% if the mouse has non-consecutive sessions, indicate weeks for which the mouse has data in the following way
all.addData([1 1 0 1 1 1]);             % indicates data present on weeks 1, 2, 4, 5, 6


% distribution of quality values
qual = all.RoiINFO.quality;
figure
histogram(qual);
xline(2.5)      % threshold indicator (2.5 == 3)
xlim([0 6])
axis square


%% Generate RDI curves using different quality thresholds


% below is old, all pooled-data RDI curves should be generated using the MixedEffectsModels_RDI function

% figure
% hold on
% all.visualize_quality = true;
% all.setQualityThreshold(5);
% all.RDIplotter('PDG', 'Average');        
% all.RDIplotter('NatMov', 'Average');
% RDIdata = all.getRDIdata;
% Fdata.qual5.PDG.curve = RDIdata.PDG.RDI_avg;
% Fdata.qual5.PDG.errorbars = RDIdata.PDG.RDI_SEM;
% Fdata.qual5.NatMov.curve = RDIdata.NatMov.RDI_avg;
% Fdata.qual5.NatMov.errorbars = RDIdata.NatMov.RDI_SEM;
% 
% all.setQualityThreshold(4);
% all.RDIplotter('PDG', 'Average');        
% all.RDIplotter('NatMov', 'Average');
% RDIdata = all.getRDIdata;
% Fdata.qual4.PDG.curve = RDIdata.PDG.RDI_avg;
% Fdata.qual4.PDG.errorbars = RDIdata.PDG.RDI_SEM;
% Fdata.qual4.NatMov.curve = RDIdata.NatMov.RDI_avg;
% Fdata.qual4.NatMov.errorbars = RDIdata.NatMov.RDI_SEM;
% 
% all.setQualityThreshold(3);
% all.RDIplotter('PDG', 'Average');       
% all.RDIplotter('NatMov', 'Average');
% RDIdata = all.getRDIdata;
% Fdata.qual3.PDG.curve = RDIdata.PDG.RDI_avg;
% Fdata.qual3.PDG.errorbars = RDIdata.PDG.RDI_SEM;
% Fdata.qual3.NatMov.curve = RDIdata.NatMov.RDI_avg;
% Fdata.qual3.NatMov.errorbars = RDIdata.NatMov.RDI_SEM;
% 
% all.setQualityThreshold(2);
% all.RDIplotter('PDG', 'Average');        
% all.RDIplotter('NatMov', 'Average');
% ylim([-0.1 0.5])
% RDIdata = all.getRDIdata;
% Fdata.qual2.PDG.curve = RDIdata.PDG.RDI_avg;
% Fdata.qual2.PDG.errorbars = RDIdata.PDG.RDI_SEM;
% Fdata.qual2.NatMov.curve = RDIdata.NatMov.RDI_avg;
% Fdata.qual2.NatMov.errorbars = RDIdata.NatMov.RDI_SEM;

% these control lines are still used though
all.setQualityThreshold(3);
all.cellSelection;
good_cells = all.getUse_cells;
refline(0, nanmean(all.StabilityData.PDG.RDI_control(good_cells)));
refline(0, nanmean(all.StabilityData.NatMov.RDI_control(good_cells)));


