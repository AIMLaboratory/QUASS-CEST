
% clc; clear; close all;
clear

OutputDir = ".\Out";

filefolder = "..\invivo_code\Outprocess_multi_batch";
load(filefolder+"\MeanArr.mat","QUASS_HC","QUASS_MS","Raw_HC","Raw_MS")
    % structure, [1,nsubject]
    %   Raw_MS(idxsub).[GM_mean/WM_mean/lesion_mean]: [nSatpara, nCEST, nsubregion]
data_HC = cat(1,Raw_HC,QUASS_HC); % [2,nsubject]
data_MS = cat(1,Raw_MS,QUASS_MS);

%% 1. data
ytagList = ["Amp.","","","",""]; 
% ytagList = ["Amide Amp.","CEST@2ppm Amp.","NOE Amp.","MT Amp."]; 
ymaxList = [0.035, 0.035, 0.09, 0.12;...
            0.035, 0.035, 0.09, 0.12];

for idxSat = 1:3

    Fig = figure();set(gcf,'position',[50,50,1500,800])
    tiledlayout(2,4,'TileSpacing','loose','Padding','compact')
    
    for idxQUASS = 1:2 % 1 for raw, 2 for QUASS
        for idxCEST = 1:4
            % GM region
            GM_HC = arrayfun(@(x) x.GM_mean(idxSat,idxCEST), data_HC(idxQUASS,:));
            GM_MS = arrayfun(@(x) x.GM_mean(idxSat,idxCEST), data_MS(idxQUASS,:));
    
            % WM region
            WM_HC = arrayfun(@(x) x.WM_mean(idxSat,idxCEST), data_HC(idxQUASS,:));
            WM_MS = arrayfun(@(x) x.WM_mean(idxSat,idxCEST), data_MS(idxQUASS,:));
            WM_lesion = cell2mat(arrayfun(@(x) squeeze(x.lesion_mean(idxSat,idxCEST,:))', data_MS(idxQUASS,:), 'UniformOutput', false));
            WM_lesion = WM_lesion(~isnan(WM_lesion));
    
            if idxQUASS==1 && idxCEST ==1
                legend_flag = 1;
            else
                legend_flag = 0;
            end
    
            nexttile
            plotRegionHist(GM_HC(:), GM_MS(:), WM_HC(:), WM_MS(:), WM_lesion(:), ytagList(idxCEST), legend_flag, ymaxList(idxQUASS,idxCEST));
        end
    end
    
    exportgraphics(Fig, OutputDir+"\Sat#"+num2str(idxSat)+"_sublesion.png", 'BackgroundColor', 'white', 'Resolution', 600);
end

function plotRegionHist(GM_HC, GM_MS, WM_HC, WM_MS, WM_lesion, ytag, legend_flag,ymax_range)
% INPUT
%   [GM/WM]_[HC/MS/lesion]: Nx1 cell double array

    %% 2. Merge data
    data_GM = {GM_HC, GM_MS}; % 1x2 cell array
    data_WM = {WM_HC, WM_MS, WM_lesion}; % 1x3 cell array
    
    %% 3. Statistical tests (ANOVA per group)
    % Two GM measurements (HC, MS)
    % (1) Independent-samples t-test: WM_HC vs WM_MS
    [~, p_GM, ~, ~] = ttest2(GM_HC, GM_MS, 'Vartype', 'unequal');
    c_GM = [1,2,p_GM];

    % Three WM measurements (WM, NAWM, lesion)
    % (1) Independent-samples t-test: WM_HC vs WM_MS
    [~, p_WM_NAWM] = ttest2(WM_HC, WM_MS, 'Vartype', 'unequal');

    
    % (2) Independent-samples t-test: WM_HC vs WM_lesion
    [~, p_WM_lesion] = ttest2(WM_HC, WM_lesion, 'Vartype', 'unequal');

    
    % (3) Independent-samples t-test: WM_MS vs WM_lesion
    [~, p_NAWM_lesion] = ttest2(WM_MS, WM_lesion, 'Vartype', 'unequal');
    
    % (4) Simple Bonferroni correction (optional)
    % p_ind1 = p_ind1 * 3; % Number of comparisons = 3
    % p_ind2 = p_ind2 * 3;
    % p_paired = p_paired * 3;

    c_WM = zeros(3,3); % 3 group for [index1, index2, pvalue]
    c_WM(1,:) = [1,2,p_WM_NAWM]; % WM - NAWM
    c_WM(2,:) = [1,3,p_WM_lesion]; % WM - Lesion
    c_WM(3,:) = [2,3,p_NAWM_lesion]; % NAWM - Lesion
    
    %% 4. Compute mean and standard error
    means_GM = cellfun(@mean, data_GM);
    sems_GM = cellfun(@(x) std(x)/sqrt(length(x)), data_GM);
    
    means_WM = cellfun(@mean, data_WM);
    sems_WM = cellfun(@(x) std(x)/sqrt(length(x)), data_WM);
    
    % Combine all means/SEMs for plotting
    all_means = [means_GM, means_WM];
    all_sems = [sems_GM, sems_WM];
    
    %% 5. Plot bar chart
    x_center_GM = 1;
    x_center_WM = 3;
    bar_width = 0.4;     % width of bar
    bar_gap = 0.1;       % gap within bar pair
    x_pos_GM = x_center_GM + [-1,1]*(bar_width/2+bar_gap/2); % Positions for two GM measurements
    x_pos_WM = x_center_WM + [-1,0,1]*(bar_width+bar_gap); % Positions for three WM measurements
    x_all = [x_pos_GM, x_pos_WM];
    
    hold on;
    
    % Draw each bar separately and set DisplayName
    bar_colors = [0.929 0.694 0.125;  % GM_HC, yellow
                  0.850 0.325 0.098;   % GM_MS, orange
                  0.929 0.694 0.125;  % WM_HC, yellow
                  0.850 0.325 0.098;   % WM_MS, orange
                  0.000 0.447 0.741];   % WM_MS lesion, blue
    
    % Draw each bar
    for i = 1:length(x_all)
        h_bar(i) = bar(x_all(i), all_means(i), bar_width, ...
            'FaceColor', bar_colors(i,:),'EdgeColor','k','LineWidth',1);
        hold on;
    end
    
    % Add error bars
    errorbar(x_all, all_means, all_sems, 'k.', 'LineWidth', 1.5, 'CapSize', 10);
    
    %% 6. Set x-axis labels and styles
    ylabel(ytag, 'FontSize', 12);
    set(gca, 'XTick', x_all, 'XTickLabel', ''); % Clear existing labels
    
    % Define group label positions and text
    metric_groups = {'GM', 'WM'};
    group_positions = [x_center_GM, x_center_WM]; % Group centers
    
    % Add group labels below x-axis
    text(group_positions, ...
         zeros(1,2) - 0.05 * max(ylim), ... % Y position (below axis)
         metric_groups, ...
         'HorizontalAlignment', 'center', ...
         'VerticalAlignment', 'top', ...
         'FontSize', 12, 'FontWeight', 'bold');
    set(gca, 'FontName', 'Times New Roman', 'FontSize', 14,'FontWeight','bold');
    
    %% 7. Add significance markers
    if ~exist('ymax_range','var') | isempty(ymax_range)
        y_max = max(all_means + all_sems);
    else
        y_max = ymax_range;
    end
    y_offset = 0.3 * y_max;
    
    % (1) Significance for two GM bars
    for i = 1:size(c_GM, 1)
        group1 = c_GM(i, 1);
        group2 = c_GM(i, 2);
        p_val = c_GM(i, end);
        
        % Determine significance symbol
        if p_val < 0.001
            sig_symbol = '***';
        elseif p_val < 0.01
            sig_symbol = '**';
        elseif p_val < 0.05
            sig_symbol = '*';
        else
            sig_symbol = 'ns';
        end
        
        % Draw connecting line and marker
        x1 = x_pos_GM(group1);
        x2 = x_pos_GM(group2);
        y_line = y_max + y_offset * (i * 0.5);
    
        plot([x1, x2], [y_line, y_line], '-k', 'LineWidth', 2);
        text((x1+x2)/2, y_line + 0.05*y_max, sig_symbol, ...
             'HorizontalAlignment', 'center', 'FontSize', 12);
    end
    
    % (2) Significance for three WM bars
    for i = 1:size(c_WM, 1)
        group1 = c_WM(i, 1);
        group2 = c_WM(i, 2);
        p_val = c_WM(i, end);
        
        if p_val < 0.001
            sig_symbol = '***';
        elseif p_val < 0.01
            sig_symbol = '**';
        elseif p_val < 0.05
            sig_symbol = '*';
        else
            sig_symbol = 'ns';
        end
        
        % Draw connecting line and marker
        x1 = x_pos_WM(group1);
        x2 = x_pos_WM(group2);
        y_line = y_max + y_offset * (i * 0.5);
    
        plot([x1, x2], [y_line, y_line], '-k', 'LineWidth', 2);
        text((x1+x2)/2, y_line + 0.05*y_max, sig_symbol, ...
             'HorizontalAlignment', 'center', 'FontSize', 12);
    end
    
    % Adjust Y-axis range
    ylim([0, y_max + y_offset * 2]);

    hold off;
    
    %% 8. Add legend (by groups)
    if legend_flag == 1
        legend_labels = {'HC', 'MS', 'MS lesions'};
        legend(h_bar(3:5), legend_labels(1:3),'Location', 'northwest', 'NumColumns', 1, 'FontSize', 10);
    end
end
