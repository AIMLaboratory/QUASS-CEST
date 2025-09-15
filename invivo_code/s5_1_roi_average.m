
% clc; clear; close all;
clear

fittingDataDir = ".\Outprocess_multi_batch";

%% load all fitting results
dateDir = dir(fittingDataDir);
dateDir = dateDir([dateDir.isdir]);
dateDir = dateDir(~ismember({dateDir.name}, {'.', '..'}));

hc_folders = {};
ms_folders = {};
for i = 1:length(dateDir)
    folderName = dateDir(i).name;
    if startsWith(folderName, 'QUASS_preproc_hc_')
        hc_folders{end+1} = folderName;
    elseif startsWith(folderName, 'QUASS_preproc_ms_')
        ms_folders{end+1} = folderName;
    end
end
hc_folderList = fullfile(fittingDataDir,hc_folders);
ms_folderList = fullfile(fittingDataDir,ms_folders);

Raw_MS = calGroupMeanValue(ms_folderList,'ms','raw');
Raw_HC = calGroupMeanValue(hc_folderList,'hc','raw');

QUASS_MS = calGroupMeanValue(ms_folderList,'ms','QUASS');
QUASS_HC = calGroupMeanValue(hc_folderList,'hc','QUASS');

save(fittingDataDir+"\MeanArr.mat","QUASS_MS","QUASS_HC","Raw_MS","Raw_HC");

%%
function [meanResults] = calGroupMeanValue(folderList,class,preptag)
% INPUT
%   output_DMPLF.mat
%       CESTfit_data: [1,6] structure
%           TR
%           Ts
%           B1
%           raw_fitpara: [256,256,16] fitting parameters
%           QUASS_fitpara: [256,256,16] fitting parameters
%   class: 'hc', 'ms'
%   preptag: 'raw', 'QUASS'
% OUTPUT:
%   [contrast]MeanArr: [nsubject,nSatpara,nregion]
    groupPathList = folderList;
    
    nsubject = length(groupPathList);
    switch class
        case 'hc'
            nregion = 2; % GM/WM
        case 'ms'
            nregion = 3; % GM/WM/lesion
    end
    nSatpara = 6; % saturation settings 
    meanResults = struct();
    
    calMeanfromMap = @(map_2D, roi_2D) mean(map_2D(roi_2D~=0),'all');
    subjectIDList = {};
    for folderidx = 1:nsubject
        % (0) subject ID
        subjectID = strsplit(groupPathList{folderidx},'\');
        subjectID = strrep(subjectID{end},'QUASS_preproc_ms_','');
        subjectIDList{folderidx} = subjectID;

        % (1) mask
        load(groupPathList{folderidx}+"\CESTinput.mat",'roi_gm','roi_wm','roi_csf');
        roi_wm(roi_csf==1) = 0;
        roi_wm(roi_gm==1) = 0;
        roi_gm(roi_csf==1) = 0;
        if strcmp(class,'ms') == 1
            load(".\MSlesionROI_8MS_20250702\roi_lesion_"+subjectID+".mat",'roi_lesion');
            roi_csf(roi_lesion==1) = 0;
            roi_gm(roi_lesion==1) = 0;
            roi_wm(roi_lesion==1) = 0;
        end

        % (2) CEST contrast mean value
        meanResults(folderidx).GM_mean = zeros(nSatpara,4,1); % [nSat, nCEST, nsubROI]
        meanResults(folderidx).WM_mean = zeros(nSatpara,4,1); 

        if strcmp(class,'ms') == 1
            if 1 % seperate roi_lesion base on connectivity
                [roi_sublesion, numLesions] = bwlabeln(roi_lesion); % roi value from 1 to numLesions
                meanResults(folderidx).lesion_mean = zeros(nSatpara,4,max(numLesions,1)); % avoid empty
            else
                roi_sublesion = roi_lesion;
                numLesions = 1;
            end
        end
        

        load(groupPathList{folderidx}+"\output_DMPLF.mat",'CESTfit_data');
        for idxSatpara = 1:nSatpara
            dataTemp = CESTfit_data(idxSatpara);
            switch preptag
                case 'raw'
                    fitpara = dataTemp.raw_fitpara;
                case 'QUASS'
                    fitpara = dataTemp.QUASS_fitpara;
            end
            
            % Water (2-4), Amide (5-7), NOE (8-10), MT (11-13), Guan (14-16)
            amide_map = fitpara(:,:,5);
            noe_map = fitpara(:,:,8);
            mt_map = fitpara(:,:,11);
            guan_map = fitpara(:,:,14);
    
            meanResults(folderidx).GM_mean(idxSatpara,1,1) = calMeanfromMap(amide_map, roi_gm);
            meanResults(folderidx).GM_mean(idxSatpara,2,1) = calMeanfromMap(guan_map, roi_gm);
            meanResults(folderidx).GM_mean(idxSatpara,3,1) = calMeanfromMap(noe_map, roi_gm);
            meanResults(folderidx).GM_mean(idxSatpara,4,1) = calMeanfromMap(mt_map, roi_gm);

            meanResults(folderidx).WM_mean(idxSatpara,1,1) = calMeanfromMap(amide_map, roi_wm);
            meanResults(folderidx).WM_mean(idxSatpara,2,1) = calMeanfromMap(guan_map, roi_wm);
            meanResults(folderidx).WM_mean(idxSatpara,3,1) = calMeanfromMap(noe_map, roi_wm);
            meanResults(folderidx).WM_mean(idxSatpara,4,1) = calMeanfromMap(mt_map, roi_wm);

            if strcmp(class,'ms') == 1
                if numLesions == 0
                    meanResults(folderidx).lesion_mean(idxSatpara,1,1) = NaN;
                    meanResults(folderidx).lesion_mean(idxSatpara,2,1) = NaN;
                    meanResults(folderidx).lesion_mean(idxSatpara,3,1) = NaN;
                    meanResults(folderidx).lesion_mean(idxSatpara,4,1) = NaN;
                else
                    for idxsubles = 1:numLesions
                        calMeanfromLes = @(map_2D, roi_2D) mean(map_2D(roi_2D == idxsubles),'all');
                        meanResults(folderidx).lesion_mean(idxSatpara,1,idxsubles) = calMeanfromLes(amide_map, roi_sublesion);
                        meanResults(folderidx).lesion_mean(idxSatpara,2,idxsubles) = calMeanfromLes(guan_map, roi_sublesion);
                        meanResults(folderidx).lesion_mean(idxSatpara,3,idxsubles) = calMeanfromLes(noe_map, roi_sublesion);
                        meanResults(folderidx).lesion_mean(idxSatpara,4,idxsubles) = calMeanfromLes(mt_map, roi_sublesion);
                    end
                end
            end
            
        end
    end

end
