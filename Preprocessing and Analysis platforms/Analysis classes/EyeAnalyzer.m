classdef EyeAnalyzer < StabilityAnalyzer
    
    properties
    	pupilinfo
    	pupilraw
    end
    
    methods

    	function obj = EyeAnalyzer()
    		%empty (use importData)
    	end

    	function density = drawCentroidMap(obj)
    		ysize = obj.pupilinfo.movie_size(1);
    		xsize = obj.pupilinfo.movie_size(2);
			ypts = 1:ysize;
			xpts = 1:xsize;
    		num_plots = ceil(obj.num_sessions/3);
    		%plot bounds
    		top = 18;
    		bottom = 82;
    		left = 40;
    		right = 120;

    		dmaps = figure;
    		frames = figure;

    		for ii = 1:obj.num_sessions
    			currpoints = obj.pupilinfo.centroid(ii).resampled;

    			density(:, :, ii) = histcounts2(currpoints(:, 2), currpoints(:, 1), ypts, xpts);
    			smoothmap = imgaussfilt(density(:, :, ii), 2);

				% outline = obj.getEyeOutline(ii);
				% smoothmap(outline) = max(max(smoothmap));

    			figure(dmaps)
				dm(ii) = subplot(num_plots, 3, ii);
				imagesc(xpts, ypts, smoothmap);
				% set(gca, 'XLim', xpts([left right]), 'YLim', ypts([top bottom]));
				set(gca, 'XLim', xpts([1 end]), 'YLim', ypts([1 end]));
				colormap parula
                pbaspect(dm(ii), [xsize ysize 1])
				hold on
				obj.drawAveragePupilBoundary(ii);

				figure(frames)
				fr(ii) = subplot(num_plots, 3, ii);
				% imagesc(obj.pupilinfo.mean_frame(top:bottom, left:right, ii));
				imagesc(obj.pupilinfo.mean_frame(:, :, ii));	%should normalize all of these		
				colormap gray
                pbaspect(fr(ii), [xsize ysize 1])
    		end

    		avgmap = mean(density, 3);
    		smoothmap = imgaussfilt(avgmap, 2);
    		figure
    		imagesc(smoothmap);
			% set(gca, 'XLim', xpts([left right]), 'YLim', ypts([top bottom]));
    		hold on
    		obj.drawAveragePupilBoundary(1);
    		mean_centroids = obj.getMeanCentroids();
            scatter(mean_centroids(:, 1), mean_centroids(:, 2));
            pbaspect(gca, [xsize ysize 1])
    	end

        %% All of the following is old, unused code
        
    	% function [RDIavgs, BS_distances, WS_distances, std_centroids, pupil_diff] = deltaPupil(obj)
    	% 	%BS_distances: distance between each session's mean centroid and the reference session's mean centroid
    	% 	%WS_distances: average distance between every frame's centroid and the mean centroid for given session
    	% 	%std_centroids: standard deviation of x and y deflections on every session
    	% 	mean_centroids = obj.getMeanCentroids();
    	% 	std_centroids = obj.getStdCentroids();
     %        pupilsize = obj.getPupilSize();

     %        BS_distances = pdist2(mean_centroids(2:end, :), mean_centroids(1, :));

    	% 	WS_distances = zeros(obj.num_sessions, 1);
    	% 	for jj = 1:obj.num_sessions
    	% 		WS_distances_curr = pdist2(obj.pupilinfo.centroid(jj).raw, mean_centroids(jj, :));
    	% 		WS_distances(jj) = mean(WS_distances_curr);
    	% 	end

     %        mean_pupilsize = cellfun(@mean, pupilsize);
     %        pupil_diff.mean_pupilsize = mean_pupilsize;
     %        % pupil_diff.BS_diff = mean_pupilsize(2:end) - mean_pupilsize(1);
     %        % pupil_diff.BS_diff_norm = pupil_diff.BS_diff/mean_pupilsize(1);
     %        pupil_diff.std_pupilsize = cellfun(@std, pupilsize);


     %        obj.cellSelection;      %select cells to use from criteria list
     %        [curr_RDIdata, ~] = obj.extractRDI(obj.pupilinfo.stim_type);
     %        obj.setRDIdata(structPacker(obj.getRDIdata, curr_RDIdata, obj.pupilinfo.stim_type));
     %        RDIavgs = curr_RDIdata.RDI_avg;
    	% end

    	% function [BS_distance_CCs, p_bsdist, WS_distance_CCs, p_wsdist, y_CCs, p_y, x_CCs, p_x] = bycellCC(obj)
    	% 	[~, BS_distances, WS_distances, std_centroids] = obj.deltaCentroid;
    	% 	use_sessions = obj.getUse_sessions;

    	% 	y_stds = std_centroids(2:end, 2);
    	% 	x_stds = std_centroids(2:end, 1);

    	% 	BS_distance_CCs = zeros(1, obj.num_cells);			% CC between a cell's RDI and distance from mean centroid on reference (avg across sessions)
    	% 	WS_distance_CCs = zeros(1, obj.num_cells);
    	% 	y_CCs = zeros(1, obj.num_cells);				% CC between a cell's RDI and std of y deflection (avg across sessions)
    	% 	x_CCs = zeros(1, obj.num_cells);				% '' x deflection
    	% 	RDI = obj.StabilityData.(obj.pupilinfo.stim_type).RDI;
    	% 	for ii = 1:obj.num_cells
    	% 		[BS_distance_CCs(ii), p_bsdist(ii)] = corr(RDI(:, ii), BS_distances);
    	% 		[WS_distance_CCs(ii), p_wsdist(ii)] = corr(RDI(:, ii), WS_distances(2:end));
    	% 		[y_CCs(ii), p_y(ii)] = corr(RDI(:, ii), y_stds);
    	% 		[x_CCs(ii), p_x(ii)] = corr(RDI(:, ii), x_stds);
    	% 	end

    	% 	% figure
    	% 	% scatter(ones(1, sum(obj.use_cells & ~isnan(BS_distance_CCs) & p_bsdist < 0.05)),...
    	% 	% 	BS_distance_CCs(obj.use_cells & ~isnan(BS_distance_CCs) & p_bsdist < 0.05), 'filled');
    	% 	% ylim([-1 1])
    	% 	% percent_sig_bsdist = sum(p_bsdist < 0.05 & obj.use_cells)/sum(obj.use_cells)

    	% 	% figure
    	% 	% scatter(ones(1, sum(obj.use_cells & ~isnan(WS_distance_CCs) & p_wsdist < 0.05)),...
    	% 	% 	WS_distance_CCs(obj.use_cells & ~isnan(WS_distance_CCs) & p_wsdist < 0.05), 'filled');
    	% 	% ylim([-1 1])
    	% 	% percent_sig_wsdist = sum(p_wsdist < 0.05 & obj.use_cells)/sum(obj.use_cells)

    	% 	% figure
    	% 	% scatter(ones(1, sum(obj.use_cells & ~isnan(y_CCs) & p_y < 0.05)),...
    	% 	% 	y_CCs(obj.use_cells & ~isnan(y_CCs) & p_y < 0.05), 'filled');
    	% 	% ylim([-1 1])
    	% 	% percent_sig_y = sum(p_y < 0.05 & obj.use_cells)/sum(obj.use_cells)

    	% 	% figure
    	% 	% scatter(ones(1, sum(obj.use_cells & ~isnan(x_CCs) & p_x < 0.05)),...
    	% 	% 	x_CCs(obj.use_cells & ~isnan(x_CCs) & p_x < 0.05), 'filled');
    	% 	% ylim([-1 1])
    	% 	% percent_sig_x = sum(p_x < 0.05 & obj.use_cells)/sum(obj.use_cells)
    	% end

        % function delta_area = percentfromMeanArea(obj, comp_type)

        %     switch comp_type
        %     case 'WS'
        %         ref_idx = 1:obj.num_sessions;
        %     case 'BS'
        %         ref_idx = ones(1, obj.num_sessions);
        %     end

        %     for kk = 1:obj.num_sessions
        %         mean_area(kk) = mean(obj.pupilinfo.area(kk).raw, 2);
        %     end

        %     for kk = 1:obj.num_sessions
        %         delta_area{kk} = ((obj.pupilinfo.area(kk).raw - mean_area(ref_idx(kk)))/mean_area(ref_idx(kk)))*100;
        %     end
        % end

        % function out = weeksApartArea(obj)
        %     corrected_area = obj.alignWeeks(obj.pupilinfo.area, 0, 0, 'pupilinfo');

        %     sessions = obj.getUse_sessions;
        %     corrected_sessions = obj.alignWeeks(sessions, 2, 1, 'RoiINFO');

        %     max_distance = size(corrected_area, 2)-1;
        %     out = cell(1, max_distance);        
        %     weeks = obj.getWeeks;
        %     for kk = 1:size(corrected_area, 2)
        %         kk
        %         for qq = 1:size(corrected_area, 2)
        %             try
        %             if (qq > kk) && weeks(qq) && weeks(kk)
        %                 curr_dist = qq - kk;
        %                 area_diff = ((mean(corrected_area(qq).raw) - mean(corrected_area(kk).raw))/mean(corrected_area(kk).raw))*100;
        %                 out{curr_dist} = [out{curr_dist} area_diff];
        %             end
        %             catch
        %                 keyboard
        %             end
        %         end
        %     end
        % end

        % function out = weeksApartdeltacentroid(obj, type)
        %     corrected_centroid = obj.alignWeeks(obj.pupilinfo.centroid, 0, 0, 'pupilinfo');

        %     sessions = obj.getUse_sessions;
        %     corrected_sessions = obj.alignWeeks(sessions, 2, 1, 'RoiINFO');

        %     max_distance = size(corrected_centroid, 2)-1;
        %     out = cell(1, max_distance);        % output is cell array of vectors containing avg RDI values calculated for every combination of 'weeks apart', i.e. out{1} contains all avg RDI values for sessions 1 week apart
        %     weeks = obj.getWeeks;
        %     mean_centroid = obj.getMeanCentroids(corrected_centroid);

        %     for kk = 1:size(corrected_centroid, 2)
        %         kk
        %         for qq = 1:size(corrected_centroid, 2)
        %             if (qq > kk) && weeks(qq) && weeks(kk)
        %                 curr_dist = qq - kk;
        %                 switch type
        %                     case 'WS'
        %                         ref_meandelta = mean(pdist2(corrected_centroid(kk).raw, mean_centroid(kk, :)));
        %                         currtest_meandelta = mean(pdist2(corrected_centroid(qq).raw, mean_centroid(qq, :)));
        %                         delta_diff = ((currtest_meandelta - ref_meandelta)/ref_meandelta)*100;
        %                         out{curr_dist} = [out{curr_dist} delta_diff];
        %                     case 'BS'
        %                         delta = pdist2(mean_centroid(qq, :), mean_centroid(kk, :));
        %                         out{curr_dist} = [out{curr_dist} delta];
        %                 end
        %             end
        %         end
        %     end
        % end





    	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  %   	function out = sdovertime(obj)
  %   		for ii = 1:length(obj.pupilinfo.area)
  %   			Area(:, ii) = obj.pupilinfo.area(ii).resampled;
  %   			centroid_x(:, ii) = obj.pupilinfo.centroid(ii).resampled(:, 1);
  %   			centroid_y(:, ii) = obj.pupilinfo.centroid(ii).resampled(:, 2);
  %   		end

  %   		Area_std = std(Area, [], 1)./mean(Area, 1);
  %   		centroidx_std = std(centroid_x, [], 1)./mean(centroid_x, 1);
  %   		centroidy_std = std(centroid_y, [], 1)./mean(centroid_y, 1);
  %   		out.Area_std = Area_std;
  %   		out.centroidx_std = centroidx_std;
  %   		out.centroidy_std = centroidy_std;

  %   		for kk = 1:obj.num_cells
	 %    		out.PDG_area_corrs(kk) = corr(obj.StabilityData.PDG.RDI(:, kk), Area_std(2:end)');
	 %    		out.PDG_centroidx_corrs(kk) = corr(obj.StabilityData.PDG.RDI(:, kk), centroidx_std(2:end)');
	 %    		out.PDG_centroidy_corrs(kk) = corr(obj.StabilityData.PDG.RDI(:, kk), centroidy_std(2:end)');

		% 		out.NatMov_area_corrs(kk) = corr(obj.StabilityData.NatMov.RDI(:, kk), Area_std(2:end)');
	 %    		out.NatMov_centroidx_corrs(kk) = corr(obj.StabilityData.NatMov.RDI(:, kk), centroidx_std(2:end)');
	 %    		out.NatMov_centroidy_corrs(kk) = corr(obj.StabilityData.NatMov.RDI(:, kk), centroidy_std(2:end)');    		
  %   		end
  %   	end
    	
  %   	function pupilAreaBS(obj)
 	% 		%Extracting cells to be used in analysis 
  %   		obj.cellSelection();
  %   		cell_indicator = obj.getUse_cells;

  %   		meanTrial_temp = zeros(obj.num_sessions-1, obj.num_cells);
  %   		trialAveraged_temp = zeros(obj.num_sessions-1, obj.num_cells) ;
  %   		for ii = 2:obj.num_sessions
  %   			[meanTrial_temp(ii, :), trialAveraged_temp(ii, :)] = obj.pupilAreaWS(ii, cell_indicator, 0);
  %   		end
  %   		for jj = 1:obj.num_cells
	 %    		meanTrial_corr(jj) = nanmean(meanTrial_temp(obj.use_sessions_full(jj, :), jj), 1);
	 %    		trialAveraged_corr(jj) = nanmean(trialAveraged_temp(obj.use_sessions_full(jj, :), jj), 1);
	 %    	end

  %   		stim_type = obj.eyemovie_stimdata(1).stim_type;

		% 	Resp = obj.RespData.(stim_type).RespMat_Full;
  %   		if session ~= 1		%if we're not evaluateing the ref session, we can look at RDI
  %   			RDI = obj.StabilityData.(stim_type).RDI(session-1, :);
  %   		else
  %   			fprintf('Warning: using the reference session means you are unable to evaluate RDI.\nPlease choose another session if you wish to do so.\n');
  %   		end
  %   		CCs = obj.StabilityData.(stim_type).CCs(session, :);
  %   		keyboard
  %   		%%%%%%
  %   		figure
 	% 		histogram(meanTrial_corr(cell_indicator), 10);	% 4

 	% 		[pos_corr, highpos_corr, neg_corr, highneg_corr] = obj.splithist(meanTrial_corr, 0.2);
		% 	A = RDI(highneg_corr)';
		% 	B = RDI(neg_corr)';
		% 	C = RDI(pos_corr)';
		% 	D = RDI(highpos_corr)';
		% 	tog = [A; B; C; D];
		% 	groups = [zeros(length(A), 1); ones(length(B), 1); 2*ones(length(C), 1); 3*ones(length(D), 1)];
		% 	fprintf('--corr = %d cells\n-corr = %d cells\n+corr = %d cells\n++corr = %d cells\n',...
		% 	 sum(highneg_corr), sum(neg_corr), sum(pos_corr), sum(highpos_corr));

		% 	figure
		% 	boxplot(tog, groups);
		% 	set(gca, 'XTickLabel', {'--corr', '-corr', '+corr', '++corr'});
		% 	%%%%%

		% 	figure
		% 	histogram(trialAveraged_corr(cell_indicator), 10);
 	% 		[pos_corr, highpos_corr, neg_corr, highneg_corr] = obj.splithist(trialAveraged_corr', 0.1);
		% 	E = RDI(highneg_corr)';
		% 	F = RDI(neg_corr)';
		% 	G = RDI(pos_corr)';
		% 	H = RDI(highpos_corr)';
		% 	tog = [E; F; G; H];
		% 	groups = [zeros(length(E), 1); ones(length(F), 1); 2*ones(length(G), 1); 3*ones(length(H), 1)];
		% 	fprintf('--corr = %d cells\n-corr = %d cells\n+corr = %d cells\n++corr = %d cells\n',...
		% 	 sum(highneg_corr), sum(neg_corr), sum(pos_corr), sum(highpos_corr));

		% 	figure
		% 	boxplot(tog, groups);
		% 	set(gca, 'XTickLabel', {'--corr', '-corr', '+corr', '++corr'});
		% end


  %   	function [meanTrial_corr, trialAveraged_corr] = pupilAreaWS(obj, session, plot_flag, cell_indicator)		
  %   		%This method is dedicated to analysis done using data only WITHIN a given session (including questions regarding RDI on a given session)
  %   		if nargin < 4
	 %    		obj.cellSelection();
	 %    		cell_indicator = obj.getUse_cells;
	 %    	end

  %   		stim_type = obj.eyemovie_stimdata(session).stim_type;
    		
  %   		movieframes = size(obj.eyeinfo(session).cropped_movie, 3);
    		
  %   		for ff = 1:movieframes
  %   			pupil_area(ff) = obj.eyeinfo(session).pupil(ff).Area;
  %   		end
  %   		currmovie_stimdata = obj.eyemovie_stimdata(session);

  %   		Resp = obj.RespData.(stim_type).RespMat_Full(:, :, :, session);
  %   		if session ~= 1		%if we're not evaluateing the ref session, we can look at RDI
  %   			RDI = obj.StabilityData.(stim_type).RDI(session-1, :);
  %   		else
  %   			fprintf('Warning: using the reference session means you are unable to evaluate RDI.\nPlease choose another session if you wish to do so.\n');
  %   		end
  %   		CCs = obj.StabilityData.(stim_type).CCs(session, :);

  %   		pupil_area_resampled = obj.downsampleEye(pupil_area, size(Resp, 1)*size(Resp, 2));
  %   		pupil_area_sorted = obj.sortFrames(pupil_area_resampled, currmovie_stimdata, size(Resp, 2));
  %   		%%


		% 	%average pupil size vs average response on trial by trial basis
 	% 		%for every cell:
 	% 		%1. find frame-average neural response on every trial
 	% 		%2. find frame-average pupil size on every trial
 	% 		%3. find correlation across trials between these two
 	% 		%Then:
 	% 		%4. look at distribution of these correlations
 	% 		%5. compare RDI of cells with high correlations to that of cells with low correlations

 	% 		meanTrial_Resp = squeeze(mean(Resp, 2));	% 1
 	% 		meanTrial_pupilsize = mean(pupil_area_sorted, 2); % 2
 	% 		meanTrial_corr = corr(meanTrial_Resp, meanTrial_pupilsize);	% 3
 	% 		if plot_flag
	 % 			figure
	 % 			histogram(meanTrial_corr(cell_indicator), 10);	% 4

	 % 			[pos_corr, highpos_corr, neg_corr, highneg_corr] = obj.splithist(meanTrial_corr, 0.2);
	 % 			if session ~= 1
	 % 				A = RDI(highneg_corr & cell_indicator)';
	 % 				B = RDI(neg_corr & cell_indicator)';
	 % 				C = RDI(pos_corr & cell_indicator)';
	 % 				D = RDI(highpos_corr & cell_indicator)';
	 % 				tog = [A; B; C; D];
	 % 				groups = [zeros(length(A), 1); ones(length(B), 1); 2*ones(length(C), 1); 3*ones(length(D), 1)];
	 % 				fprintf('--corr = %d cells\n-corr = %d cells\n+corr = %d cells\n++corr = %d cells\n',...
	 % 				 length(A), length(B), length(C), length(D));

	 % 				figure
	 % 				boxplot(tog, groups);
	 % 				set(gca, 'XTickLabel', {'--corr', '-corr', '+corr', '++corr'});
	 % 			end
	 % 		end

 	% 		%%%%%%%%%

 	% 		%correlations between pupil size and individual cell response
 	% 		%for every cell:
 	% 		%1. find correlation between neural response and pupil size on every trial
 	% 		%2. find average correlation across trials
 	% 		%Then:
 	% 		%3. look at distribution of these correlations 
 	% 		%4. compare RDI of cells with high correlations to that of cells with low corrlations 

 	% 		for ii = 1:obj.num_cells
 	% 			for jj = 1:size(Resp, 1)
	 % 				curr_corrs(jj) = corr(Resp(jj, :, ii)', pupil_area_sorted(jj, :)');
	 % 			end
 	% 			trialAveraged_corr(ii) = mean(curr_corrs);
 	% 		end

 	% 		if plot_flag
 	% 			figure
 	% 			histogram(trialAveraged_corr(cell_indicator), 10);
	 % 			[pos_corr, highpos_corr, neg_corr, highneg_corr] = obj.splithist(trialAveraged_corr', 0.1);
	 % 			if session ~= 1
	 % 				E = RDI(highneg_corr & cell_indicator)';
	 % 				F = RDI(neg_corr & cell_indicator)';
	 % 				G = RDI(pos_corr & cell_indicator)';
	 % 				H = RDI(highpos_corr & cell_indicator)';
	 % 				tog = [E; F; G; H];
	 % 				groups = [zeros(length(E), 1); ones(length(F), 1); 2*ones(length(G), 1); 3*ones(length(H), 1)];
	 % 				fprintf('--corr = %d cells\n-corr = %d cells\n+corr = %d cells\n++corr = %d cells\n',...
	 % 				 length(A), length(B), length(C), length(D));

	 % 				figure
	 % 				boxplot(tog, groups);
	 % 				set(gca, 'XTickLabel', {'--corr', '-corr', '+corr', '++corr'});
	 % 			end
	 % 		end
  %   	end

   
    end

	methods (Access = private)

		function drawAveragePupilBoundary(obj, session)
			center = mean(obj.pupilinfo.centroid(session).resampled, 1);		%average centroid location 
			major = mean(obj.pupilinfo.majoraxis(session).resampled);
			minor = mean(obj.pupilinfo.minoraxis(session).resampled);
			ori = mean(obj.pupilinfo.orientation(session).resampled);

			theta = 0 : 0.01 : 2*pi;
            x = minor/2 * cos(theta);
            y = major/2 * sin(theta);
            
            % rotation?
            R = [cosd(ori), -sind(ori);... % create rotation matrix
            sind(ori), cosd(ori)];
            
            rCoords = R * [x; y]; % apply transform
            
            xr = rCoords(1, :)';
            yr = rCoords(2, :)';
            
            plot(yr + center(1), xr + center(2), 'LineWidth', 2, 'Color', [1, 0, 0]);
        end

        function outline = getEyeOutline(obj, session)
        	currimg = obj.pupilinfo.mean_frame(:, :, session);
        	dark = currimg < min(min(currimg)) + 3.5*std(std(currimg, [], 1), [], 2);
        	outline = bwmorph(dark, 'remove');
        end

        function centroids = getMeanCentroids(obj, input)
            if nargin < 2
                input = obj.pupilinfo.centroid;
            end
            for jj = 1:length(input)
                curr_mean = mean(input(jj).resampled, 1);
                if isempty(curr_mean)
                    centroids(jj, :) = NaN;
                else
                    centroids(jj, :) = curr_mean;
                end
            end
        end

        function stds = getStdCentroids(obj)
        	for jj = 1:obj.num_sessions
        		stds(jj, 1) = std(obj.pupilinfo.centroid(jj).resampled(:, 1), [], 1);
        		stds(jj, 2) = std(obj.pupilinfo.centroid(jj).resampled(:, 2), [], 1);
        	end
        end

        function size = getPupilSize(obj)
            for kk = 1:obj.num_sessions
                size{kk, :} = obj.pupilinfo.area(kk).raw;
            end
        end
    end
end