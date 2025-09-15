% convert dicom to mat:
% (1) coregister T1 to CEST@m0
% (2) segment brain on registered T1
% (3) TBD: realign/motion-correction on CEST
% (4) convert nii to mat

% Input:
%   \[date]\A\*_3D Ax T1 MPRAGE, 72 slices
%   \[date]\A\*_2D-T1 FA7 TR 30ms, 1 slice
%   \[date]\A\*_2D-T1 FA20 TR 30ms, 1 slice
%   \[date]\A\*_2D WASABI 3.7uT 5ms, 1 slice

%   \[date]\A\*_CEST 23s-0.5uT, 1 slice
%   \[date]\A\*_CEST 1.52.3s-0.5uT, 1 slice
%   \[date]\A\*_CEST 11.6s-0.5uT, 1 slice
%   \[date]\A\*_CEST 11.6s-1.0uT, 1 slice
%   \[date]\A\*_CEST 1.52.3s-1.0uT, 1 slice
%   \[date]\A\*_CEST 11.6s-1.0uT, 1 slice


% Huabin Zhang - huabinz1202@connect.hku.hk, 2025-0908

clear
close all
clc

seqDcmList = ["*_3D Ax T1 MPRAGE\";...
           "*_3D Ax T2 Cube\";...
           "*_3D-T1 FA7 TR30ms\";...
           "*_3D-T1 FA20 TR30ms\";...
           "*_2D WASABI 3.7uT 5ms\";...
           "*_CEST 23s-0.5uT\";...
           "*_CEST1.52.3s-0.5uT\";...
           "*_CEST 11.6s-0.5uT\";...
           "*_CEST 23s-1.0uT\";...
           "*_CEST1.52.3s-1.0uT\";...
           "*_CEST 11.6s-1.0uT\"]; 

dateDirList = ["\ms\20250206\";...
               "\ms\20250213\";...
               "\ms\20250220\";...
               "\ms\20250227\";...
               "\ms\20250306\";...
               "\ms\20250502\";...
               "\ms\20250520\";...
               "\ms\20250522\";...
               "\hc\20250214_1\";...
               "\hc\20250214_2\";...
               "\hc\20250228\";...
               "\hc\20250307_1\";...
               "\hc\20250307_2\";...
               "\hc\20250313\";...
               "\hc\20250509_1\";...
               "\hc\20250509_2\"];
dateDirList = "D:\InvivoData\GE\MS\MotionCorr(8MS8HC)" + dateDirList;

% dateDirList = [".\sampleData\hc\20250228\"];
for folderidx = 1:length(dateDirList)
dataDir = dateDirList(folderidx); 
dataID = strsplit(dataDir,'\');
dataID = dataID(end-2)+"_"+dataID(end-1);
fprintf("preprocessing "+dataID+"\n");

    %% (0) dcm to nii
    seqNiiList = [];
    % [fileNameList,dicomPath] = uigetfile({'*.*','All Files (*.*)'},'Select multiple CEST DICOM files','.\sampleData\20250228\A\7_3D CEST, CW 0.8uT2s\', 'MultiSelect', 'on');
    % dicomFileDir = reshape(fullfile(dicomPath,fileNameList),[],1); % cell, full dir of all CEST DICOM files
    
    tic
    fprintf("\tconvert dicom to nii: ");
    for idx = 1:length(seqDcmList)
        temp = dir(dataDir+"A\"+seqDcmList(idx));
        temp = temp(~ismember({temp.name}, {'.', '..'})); % eliminate ".\" and "..\"
        dicomFolder = temp(1).folder;
        dicomFileDir = dicomFolder + "\" + {temp.name};
        dicomFileDir = cellstr(reshape(dicomFileDir,[],1));
        
        strTemp = strsplit(dicomFolder,'\');
        strTemp{end-1} = char(strTemp{end-1}+"pro");
        AproPath = string(join(strTemp(1:end-1),'\')); % store preprocessed .nii data in folder Apro
        if ~exist(AproPath,'dir')
            mkdir(AproPath);
        end
        
        OutputDir = string(join(strTemp,'\')) + "\niiTemp"; % temporarily store nii data
        seqNiiList = [seqNiiList;OutputDir];
        if exist(OutputDir,'dir')
            rmdir(OutputDir,'s'); % remove nonempty folder
        end
        mkdir(OutputDir);
        
        matlabbatch_dcm2nii{1}.spm.util.import.dicom.data = dicomFileDir;
        matlabbatch_dcm2nii{1}.spm.util.import.dicom.root = 'flat';
        matlabbatch_dcm2nii{1}.spm.util.import.dicom.outdir = {char(OutputDir)};
        matlabbatch_dcm2nii{1}.spm.util.import.dicom.protfilter = '.*';
        matlabbatch_dcm2nii{1}.spm.util.import.dicom.convopts.format = 'nii';
        matlabbatch_dcm2nii{1}.spm.util.import.dicom.convopts.meta = 0;
        matlabbatch_dcm2nii{1}.spm.util.import.dicom.convopts.icedims = 0;
        evalc("spm_jobman('run', matlabbatch_dcm2nii)");
    end
    fprintf('%.4f s\n', toc);

    %% (1) segment on HR T1w
    fprintf("\tsegment on HR T1w: ");
    tic
    HRT1Name = dir(seqNiiList(1)+"\*-000001-01.nii"); % first one is T1w
    
    spm_tpmDir = fullfile(spm('Dir'),'tpm','TPM.nii');
    matlabbatch_segment{1}.spm.spatial.preproc.channel.vols = {char(HRT1Name.folder+"\"+HRT1Name.name+",1")};
    matlabbatch_segment{1}.spm.spatial.preproc.channel.biasreg = 0.001;
    matlabbatch_segment{1}.spm.spatial.preproc.channel.biasfwhm = 60;
    matlabbatch_segment{1}.spm.spatial.preproc.channel.write = [0 0];
    matlabbatch_segment{1}.spm.spatial.preproc.tissue(1).tpm = {char(spm_tpmDir+",1")};
    matlabbatch_segment{1}.spm.spatial.preproc.tissue(1).ngaus = 1;
    matlabbatch_segment{1}.spm.spatial.preproc.tissue(1).native = [1 0];
    matlabbatch_segment{1}.spm.spatial.preproc.tissue(1).warped = [0 0];
    matlabbatch_segment{1}.spm.spatial.preproc.tissue(2).tpm = {char(spm_tpmDir+",2")};
    matlabbatch_segment{1}.spm.spatial.preproc.tissue(2).ngaus = 1;
    matlabbatch_segment{1}.spm.spatial.preproc.tissue(2).native = [1 0];
    matlabbatch_segment{1}.spm.spatial.preproc.tissue(2).warped = [0 0];
    matlabbatch_segment{1}.spm.spatial.preproc.tissue(3).tpm = {char(spm_tpmDir+",3")};
    matlabbatch_segment{1}.spm.spatial.preproc.tissue(3).ngaus = 2;
    matlabbatch_segment{1}.spm.spatial.preproc.tissue(3).native = [1 0];
    matlabbatch_segment{1}.spm.spatial.preproc.tissue(3).warped = [0 0];
    matlabbatch_segment{1}.spm.spatial.preproc.tissue(4).tpm = {char(spm_tpmDir+",4")};
    matlabbatch_segment{1}.spm.spatial.preproc.tissue(4).ngaus = 3;
    matlabbatch_segment{1}.spm.spatial.preproc.tissue(4).native = [1 0];
    matlabbatch_segment{1}.spm.spatial.preproc.tissue(4).warped = [0 0];
    matlabbatch_segment{1}.spm.spatial.preproc.tissue(5).tpm = {char(spm_tpmDir+",5")};
    matlabbatch_segment{1}.spm.spatial.preproc.tissue(5).ngaus = 4;
    matlabbatch_segment{1}.spm.spatial.preproc.tissue(5).native = [1 0];
    matlabbatch_segment{1}.spm.spatial.preproc.tissue(5).warped = [0 0];
    matlabbatch_segment{1}.spm.spatial.preproc.tissue(6).tpm = {char(spm_tpmDir+",6")};
    matlabbatch_segment{1}.spm.spatial.preproc.tissue(6).ngaus = 2;
    matlabbatch_segment{1}.spm.spatial.preproc.tissue(6).native = [0 0];
    matlabbatch_segment{1}.spm.spatial.preproc.tissue(6).warped = [0 0];
    matlabbatch_segment{1}.spm.spatial.preproc.warp.mrf = 1;
    matlabbatch_segment{1}.spm.spatial.preproc.warp.cleanup = 1;
    matlabbatch_segment{1}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
    matlabbatch_segment{1}.spm.spatial.preproc.warp.affreg = 'eastern';
    matlabbatch_segment{1}.spm.spatial.preproc.warp.fwhm = 0;
    matlabbatch_segment{1}.spm.spatial.preproc.warp.samp = 3;
    matlabbatch_segment{1}.spm.spatial.preproc.warp.write = [0 0];
    evalc("spm_jobman('run', matlabbatch_segment)");
    fprintf('%.4f s\n', toc);

    %% (2) reslice mask/T1w/T2w/DFA to CEST@m0
    fprintf("\treslice mask/T1w/T2w/DFA to CEST@m0: ");
    tic 
    CESTm0 = dir(seqNiiList(end)+"\*-000001-01.nii");

    T1wName = dir(seqNiiList(1)+"\*-000001-01.nii"); % "HRT1wName.nii" and "c[1-5]HRT1wName.nii"
    T2wName = dir(seqNiiList(2)+"\*-000001-01.nii");
    FA7Name = dir(seqNiiList(3)+"\*-000001-01.nii");
    FA20Name = dir(seqNiiList(4)+"\*-000001-01.nii");

    targetName = [T1wName;FA7Name;FA20Name]; % first one is T1w
    matlabbatch{1}.spm.spatial.coreg.write.ref = {char(CESTm0.folder+"\"+CESTm0.name+",1")};
    matlabbatch{1}.spm.spatial.coreg.write.source = reshape(cellstr({targetName.folder}+"\"+{targetName.name}+",1"),[],1);
    matlabbatch{1}.spm.spatial.coreg.write.roptions.interp = 4; % n-spline interpolation
    matlabbatch{1}.spm.spatial.coreg.write.roptions.wrap = [0 0 0];
    matlabbatch{1}.spm.spatial.coreg.write.roptions.mask = 0;
    matlabbatch{1}.spm.spatial.coreg.write.roptions.prefix = 'CESTspace_';
    evalc("spm_jobman('run', matlabbatch)");
    fprintf('%.4f s\n', toc);

    %% (3) convert nii to mat
    fprintf("\tsave preprocessed to .mat: ");
    tic

    % T1w, in uint16
    T1w = read_nii(HRT1Name.folder+"\CESTspace_"+HRT1Name.name);
    [nx,ny,nz] = size(T1w);
    
    % masks for GM, WM, and CSF, in uint8
    dir_roi = dir(seqNiiList(1)+"\CESTspace_c1*.nii");
    roi_gm_uint = read_nii(dir_roi.folder+"\"+dir_roi.name);
    dir_roi = dir(seqNiiList(1)+"\CESTspace_c2*.nii");
    roi_wm_uint = read_nii(dir_roi.folder+"\"+dir_roi.name);
    dir_roi = dir(seqNiiList(1)+"\CESTspace_c3*.nii");
    roi_csf_uint = read_nii(dir_roi.folder+"\"+dir_roi.name);

    % T2w
    dir_T2w = dir(seqNiiList(2)+"\CESTspace_*-000001-01.nii");
    T2w = read_nii(dir_T2w.folder+"\"+dir_T2w.name);

    % DFA at FA = 7, 20
    dir_FA7 = dir(seqNiiList(3)+"\CESTspace_*-000001-01.nii");
    GRE_FA7 = read_nii(dir_FA7.folder+"\"+dir_FA7.name);
    
    dir_FA20 = dir(seqNiiList(4)+"\CESTspace_*-000001-01.nii");
    GRE_FA20 = read_nii(dir_FA20.folder+"\"+dir_FA20.name);
    
    % WASABI
    WASABIimg = loadCESTfromNifti(seqNiiList(5),nx,ny,nz);

    % CEST images
    CEST_data(1).raw = loadCESTfromNifti(seqNiiList(6),nx,ny,nz);
    CEST_data(1).TR = 3; % s
    CEST_data(1).Ts = 2; % s
    CEST_data(1).B1 = 0.5; % uT

    CEST_data(2).raw = loadCESTfromNifti(seqNiiList(7),nx,ny,nz);
    CEST_data(2).TR = 2.3; % s
    CEST_data(2).Ts = 1.5; % s
    CEST_data(2).B1 = 0.5; % uT

    CEST_data(3).raw = loadCESTfromNifti(seqNiiList(8),nx,ny,nz);
    CEST_data(3).TR = 1.6; % s
    CEST_data(3).Ts = 1; % s
    CEST_data(3).B1 = 0.5; % uT

    CEST_data(4).raw = loadCESTfromNifti(seqNiiList(9),nx,ny,nz);
    CEST_data(4).TR = 3; % s
    CEST_data(4).Ts = 2; % s
    CEST_data(4).B1 = 1; % uT

    CEST_data(5).raw = loadCESTfromNifti(seqNiiList(10),nx,ny,nz);
    CEST_data(5).TR = 2.3; % s
    CEST_data(5).Ts = 1.5; % s
    CEST_data(5).B1 = 1; % uT

    CEST_data(6).raw = loadCESTfromNifti(seqNiiList(11),nx,ny,nz);
    CEST_data(6).TR = 1.6; % s
    CEST_data(6).Ts = 1; % s
    CEST_data(6).B1 = 1; % uT
    
    % save
    save(dataDir+"Apro\2D_spm_preproc.mat","T1w","T2w","roi_gm_uint","roi_wm_uint","roi_csf_uint",...
                "WASABIimg", "GRE_FA7", "GRE_FA20","CEST_data");
    fprintf('%.4f s\n', toc);

    
    %% display
    sliceidx = ceil(nz/2);
    Fig = figure(1);
    mosaic(cat(3,T1w(:,:,sliceidx)./max(T1w(:,:,sliceidx),[],'all'),...
        GRE_FA20(:,:,sliceidx)./max(GRE_FA20(:,:,sliceidx),[],'all'), ...
        WASABIimg(:,:,sliceidx,1)./max(WASABIimg(:,:,sliceidx),[],'all'),...
        double(roi_gm_uint(:,:,sliceidx)),...
        double(roi_wm_uint(:,:,sliceidx)),...
        double(roi_csf_uint(:,:,sliceidx))  ),...
        2,3,1,'Top: T1w,GRE FA20, CEST@M0. Buttom: ROI of GM/WM/CSF')
    exportgraphics(Fig, dataDir+"Apro\fig_2D.png", 'BackgroundColor', 'white', 'Resolution', 600);
end

%%
function CESTimg = loadCESTfromNifti(CESTfolder,nx,ny,nz)
    CESTList = dir(CESTfolder+"\*.nii");
    CESTList = CESTList(~ismember({CESTList.name}, {'.', '..'}));
    nf = length(CESTList);
    CESTimg = zeros(nx,ny,nz,nf);
    for i = 1:nf
        CESTimg(:,:,:,i) = read_nii(fullfile(CESTList(i).folder,CESTList(i).name));
    end
end


function img_scaled = read_nii(nii_file)
% READ_NII_SCALED Read a NIfTI (.nii) file and apply rot90.
%
% Usage:
%   [img_scaled, info] = read_nii_scaled('filename.nii');
%
% Input:
%   nii_file - Path to the NIfTI file (.nii or .nii.gz)
%
% Output:
%   img_double - The image array in double

    % Get NIfTI header
    info = niftiinfo(nii_file);
    
    % Read raw image data
    img_raw = niftiread(info);

    % Default slope and intercept
    slope = 1.0;
    inter = 0.0;

    % Try to get scaling factors from raw field (SPM sometimes stores them here)
    if isfield(info, 'raw')
        raw = info.raw;
        if isfield(raw, 'scl_slope') && raw.scl_slope ~= 0
            slope = raw.scl_slope;
        end
        if isfield(raw, 'scl_inter')
            inter = raw.scl_inter;
        end
    elseif isfield(info, 'Slope') && info.Slope ~= 0  % fallback
        slope = info.Slope;
        inter = info.Intercept;
    end

    % Apply scaling
    img_scaled = double(img_raw) * slope + inter;
    img_scaled = rot90(img_scaled);
end