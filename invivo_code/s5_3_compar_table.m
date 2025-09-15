clear

filefolder = ".\Outprocess_multi_batch";
load(filefolder+"\MeanArr.mat","QUASS_HC","QUASS_MS","Raw_HC","Raw_MS")
    % structure, [1,nsubject]
    %   Raw_MS(idxsub).[GM_mean/WM_mean/lesion_mean]: [nSatpara, nCEST, nsubregion]
data_HC = cat(1,Raw_HC,QUASS_HC); % [2,nsubject]
data_MS = cat(1,Raw_MS,QUASS_MS);

%% 
CEST_names = {"Amide", "CEST@2ppm", "NOE", "MT"};
postfix_names = {"_raw", "_QUASS"};

% table size: [2(raw/QUASS) * 4(CESTs) * 6(Sat#), 8(Columns)]
result_table = table(...
    'Size', [2*4*6, 8], ...
    'VariableTypes', {'string', 'string', 'string', 'string', 'string', 'string', 'string', 'string'}, ...
    'VariableNames', {'SatGroup','CEST contrast', 'Control', 'NAWM', 'Lesion', '% Change (Control/NAWM)', '% Change (Control/Lesion)', '% Change (NAWM/Lesion)'} ...
);

i = 1;
for idxSat = 1:6
    
    for idxCEST = 1:4
        for idxQUASS = 1:2 % 1 for raw, 2 for QUASS
        
            % WM region
            WM_HC = arrayfun(@(x) x.WM_mean(idxSat,idxCEST), data_HC(idxQUASS,:));
            WM_MS = arrayfun(@(x) x.WM_mean(idxSat,idxCEST), data_MS(idxQUASS,:));
            WM_lesion = cell2mat(arrayfun(@(x) squeeze(x.lesion_mean(idxSat,idxCEST,:))', data_MS(idxQUASS,:), 'UniformOutput', false));
            WM_lesion = WM_lesion(~isnan(WM_lesion));
    
    
            % Call helper to get formatted strings
            [mean_std_A, mean_std_B, mean_std_C, change_BA, change_CA, change_CB] = format_group_stats(WM_HC, WM_MS, WM_lesion);
            
            % Write into table
            result_table{i,1} = "Sat#"+num2str(idxSat);
            result_table{i,2} = CEST_names{idxCEST} + postfix_names{idxQUASS};
            result_table{i,3} = mean_std_A;
            result_table{i,4} = mean_std_B;
            result_table{i,5} = mean_std_C;
            result_table{i,6} = change_BA;
            result_table{i,7} = change_CA;
            result_table{i,8} = change_CB;
            i = i+1;
        end
    end
end
% Display table
disp(result_table);

% Export to Excel (optional)
writetable(result_table, filefolder+"\stat_results.xlsx");

%%
function [mean_std_A, mean_std_B, mean_std_C, change_BA, change_CA, change_CB] = format_group_stats(A, B, C)
    % Compute mean and standard deviation
    mean_A = mean(A); std_A = std(A);
    mean_B = mean(B); std_B = std(B);
    mean_C = mean(C); std_C = std(C);
    
    % Compute significance (t-test)
    [~, p_BA] = ttest2(B, A, 'Vartype', 'unequal');
    [~, p_CA] = ttest2(C, A, 'Vartype', 'unequal');
    [~, p_CB] = ttest2(C, B, 'Vartype', 'unequal');
    
    % Format strings
    mean_std_A = sprintf("%.3f ± %.3f", mean_A, std_A);
    mean_std_B = sprintf("%.3f ± %.3f", mean_B, std_B);
    mean_std_C = sprintf("%.3f ± %.3f", mean_C, std_C);
    
    % Add significance asterisk *
    change_BA = sprintf("%.1f%%%s", (mean_B - mean_A)/mean_A*100, ifelse(p_BA < 0.05, '*', ''));
    change_CA = sprintf("%.1f%%%s", (mean_C - mean_A)/mean_A*100, ifelse(p_CA < 0.05, '*', ''));
    change_CB = sprintf("%.1f%%%s", (mean_C - mean_B)/mean_B*100, ifelse(p_CB < 0.05, '*', ''));
end

% --- Helper function ---
function out = ifelse(condition, true_val, false_val)
    if condition
        out = true_val;
    else
        out = false_val;
    end
end