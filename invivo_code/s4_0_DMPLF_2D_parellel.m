% perform fitting on invivo data
%       compare amide, guan, NOE
%       compare MTR, LDA, MPLF (MultiStart1), MPLF (MS-5), MPLF2-ss (MS-5)
% 2025-0203 for invivo data
% 2025-0716 motion correction
clear
% close all
warning off

addpath ..\toolbox\
addpath ..\toolbox\fitting_main_par\
multistart_N = 5; % times for MultiStart


delete(gcp('nocreate'));
parpool(16);

OutputDir = ".\Outprocess_multi_batch";
DirTemp = OutputDir;
cnt = 1;
while(1)
    if exist(DirTemp,'dir')
        DirTemp = OutputDir + "_" + num2str(cnt);
        cnt = cnt + 1;
    else
        break;
    end
end
OutputDir = DirTemp + "\";
mkdir(OutputDir);

%% in-vivo data 
dateDir = dir(".\preprocessedData_multi_batch\QUASS_preproc_*.mat");
% dateDir = dir(".\preprocessedData_multi_B1c\QUASS_preproc_ms_20250522.mat");
dateDir = dateDir(~ismember({dateDir.name}, {'.', '..'}));

for dataidx = 1:length(dateDir)

    folder = dateDir(dataidx).folder;
    dataID = strsplit(dateDir(dataidx).name,'.');
    dataID = strsplit(dataID{1},'.');
    dataID = strrep(dataID, 'QUASSpreproc_','');split(dataID{1},'_');
    fprintf("Processing "+dataID+":\n")

    %% load data
    % load(fileName,"offs","zspec"); % m0 not included
    %     % offs: [nf,1]
    %     %   offs(1): m0
    %     % zspec: [nf, nguanf, nguank, namidef]
    % 
    % [nf,nguanf,nguank,namidef,nMTf,nNOEf] = size(zspec);
    % zspec_vec = reshape(zspec,nf,[]);
    % Npixel = size(zspec_vec,2);
    OutputDirTemp = OutputDir + dataID + "\";
    mkdir(OutputDirTemp);
    
    load(fullfile(dateDir(dataidx).folder,dateDir(dataidx).name), ...
                "T1w","T1map","roi_csf","roi_gm","roi_wm","roi", ...
                "offs_WASABI","offs_CEST", ...
                "WASABI_data","CEST_data");

    offs = offs_CEST;
    [nx,ny,nf] = size(CEST_data(1).zspec);

    for satparaidx = 1:6
        CESTfit_data(satparaidx).TR = CEST_data(satparaidx).TR; % s
        CESTfit_data(satparaidx).Ts = CEST_data(satparaidx).Ts; % s
        CESTfit_data(satparaidx).B1 = CEST_data(satparaidx).B1; % s
        zmap_mea = CEST_data(satparaidx).zspec; % [nx,ny,nf]
        zmap_mea_QUASS = CEST_data(satparaidx).QUASSzspec;
     
        indnz = find(reshape(roi,[],1)~=0); 
    
        zspec_vectemp = reshape(permute(zmap_mea,[3,1,2]),nf,[]); % [nf,Npixel]
        zspec_nz_vec = zspec_vectemp(:,indnz);
        zspec_vectemp = reshape(permute(zmap_mea_QUASS,[3,1,2]),nf,[]); % [nf,Npixel]
        zspec_QUASS_nz_vec = zspec_vectemp(:,indnz);
    
        save(OutputDirTemp+"CESTinput.mat",...
            "T1w","T1map","roi_csf","roi_gm","roi_wm","roi", ...
                    "offs_WASABI","offs_CEST", ...
                    "WASABI_data","CEST_data")
        
        %% DMPLF fitting
        processing_DZ_methods = {@fitting_DMPLF_par};
        output_filenames = {'DMPLF'};
        
        % to be output
        amide_map = zeros(nx,ny);
        guan_map = zeros(nx,ny);
        noe_map = zeros(nx,ny);
    
        for i = 1:length(processing_DZ_methods)
            method = processing_DZ_methods{i};
    
            fit_paratemp = method(offs, zspec_nz_vec, multistart_N);
            
            % amide_map(indnz) = fit_paratemp(5,:);
            % guan_map(indnz) = fit_paratemp(14,:);
            % noe_map(indnz) = fit_paratemp(8,:);
            
            nterm = size(fit_paratemp,1);
            fit_para = zeros(nterm,nx*ny);
            fit_para(:,indnz) = fit_paratemp;
            fit_para = reshape(fit_para,[nterm,nx,ny]);
    
            CESTfit_data(satparaidx).raw_fitpara = permute(fit_para,[2,3,1]);
            % save(OutputDirTemp+"output_raw_"+output_filenames{i}+".mat",'amide_map','guan_map','noe_map','offs','fit_para');
            % fprintf('\t\tsaved to %s\n', "output_raw_"+output_filenames{i}+".mat");
    
            % ============== QUASS ==============
            fit_paratemp = method(offs, zspec_QUASS_nz_vec, multistart_N);
            
            % amide_map(indnz) = fit_paratemp(5,:);
            % guan_map(indnz) = fit_paratemp(14,:);
            % noe_map(indnz) = fit_paratemp(8,:);
            
            nterm = size(fit_paratemp,1);
            fit_para = zeros(nterm,nx*ny);
            fit_para(:,indnz) = fit_paratemp;
            fit_para = reshape(fit_para,[nterm,nx,ny]);
    
            CESTfit_data(satparaidx).QUASS_fitpara = permute(fit_para,[2,3,1]);
            % save(OutputDirTemp+"output_QUASS_"+output_filenames{i}+".mat",'amide_map','guan_map','noe_map','offs','fit_para');
            % fprintf('\t\tsaved to %s\n', "output_QUASS_"+output_filenames{i}+".mat");
    
        end
    end

    save(OutputDirTemp+"output_"+output_filenames{i}+".mat",'CESTfit_data','offs','roi',"roi_csf","roi_gm","roi_wm");
    fprintf('\t\tsaved to %s\n', "output_"+output_filenames{i}+".mat");
    
  
end

delete(gcp('nocreate'));

%%
function CEST_data = motionCorCEST(CEST_data)
% INPUT:
%   CEST_data: [1,3] structure
%       CEST_data(i).zspec: [256,256,nf]
% OUTPUT:
%       CEST_data(i).zspec: [256,256,nf]

    nof = size(CEST_data(1).zspec,3);
    nfra = 3*nof;
    imgSize = size(CEST_data(1).zspec,1);
    imgSeq_ori = cell(nof,3);
    for i = 1:3
        imgArr = squeeze(CEST_data(i).zspec); % [256,256,73]
        for j = 1:nof
            imgSeq_ori{j,i} = imgArr(:,:,j);
        end
    end
    imgSeq_ori = reshape(imgSeq_ori,[],1);
    
    
    %% Motion correction
    calc_dx = zeros(nfra,1);
    calc_dy = zeros(nfra,1);
    calc_rot = zeros(nfra,1);
    baseImg = imgSeq_ori{1};
    
    imgSeq_corrected = cell(nfra,1);
    imgSeq_corrected{1} = baseImg; % Use the first frame as reference
    
    % Use mutual information as registration metric for multimodal images
    [optimizer, metric] = imregconfig('monomodal');optimizer.MaximumIterations = 500;
    % [optimizer, metric] = imregconfig('multimodal');optimizer.MaximumIterations = 500;
    
    % First-pass motion correction
    backNum = 0;
    fprintf('\t1st motion correction: ')
    for i = 2:nfra
        fprintf(repmat('\b',1,backNum));
        backNum = fprintf('%d/%d',i,nfra);
        tform = imregtform(imgSeq_ori{i}, baseImg, 'rigid', optimizer, metric);
            % Compute transforms between frames; find optimal rigid transform tform such that tform @ imgmov = imgref
            % optimizer: optimization algorithm, default regm (regular step gradient descent)
            % metric: similarity metric; use mean square for single-modality; change for different contrasts
        imgSeq_corrected{i} = imwarp(imgSeq_ori{i}, tform, 'OutputView', imref2d(size(baseImg)));
        calc_dx(i) = tform.T(3,1);
        calc_dy(i) = tform.T(3,2);
        calc_rot(i) = atan2d(tform.T(2,1), tform.T(1,1));
    end
    fprintf('\n')
    
    % Second-pass motion correction - residual displacement
    residual_dx = zeros(nfra,1);
    residual_dy = zeros(nfra,1);
    backNum = 0;
    fprintf('\t2nd motion correction: ')
    for i = 2:nfra
        fprintf(repmat('\b',1,backNum));
        backNum = fprintf('%d/%d',i,nfra);
        tform = imregtform(imgSeq_corrected{i}, imgSeq_corrected{1}, 'rigid', optimizer, metric);
        residual_dx(i) = tform.T(3,1);
        residual_dy(i) = tform.T(3,2);
    end
    fprintf('\n')
    
    %% Visualization of results
    figure();set(gcf,'Position', [100 100 600 300]);
    tiledlayout(1,2)
    % 
    % ========================= Before motion correction =========================
    % Motion parameters
    nexttile
    plot(1:nfra, calc_dx, 'r-o', 1:nfra, calc_dy, 'b-s');
    legend('X shift', 'Y shift', 'Location', 'northwest');
    title('Displacement computed for motion correction'); 
    xlabel('Frame index'); ylabel('Displacement (pixels)');
    grid on;
    ylim([-1,2])
    
    % ========================= After motion correction =========================
    imgSeqArr_cor = reshape(cell2mat(imgSeq_corrected'),[imgSize,imgSize,nfra]); % [nx,ny,nf]    

    % Motion parameters
    nexttile
    plot(1:nfra, residual_dx, 'r-o', 1:nfra, residual_dy, 'b-s');
    legend('Residual X shift', 'Residual Y shift', 'Location', 'northwest');
    title('Residual motion after correction'); 
    xlabel('Frame index'); ylabel('Displacement (pixels)');
    grid on;
    ylim([-1,2])
    
    fprintf('\tMean error before registration: %.2f pixels\n', mean(sqrt((calc_dx).^2 + (calc_dy).^2)));
    fprintf('\tMean residual error: %.2f pixels\n', mean(sqrt((residual_dx).^2 + (residual_dy).^2)));

    %% Restore to zspec
    for i = 1:3
        CEST_data(i).zspec = imgSeqArr_cor(:,:,((i-1)*nof+1):(i*nof)); % [256,256,nf]
    end
end