% preprocess .mat data
% Huabin Zhang - huabinz1202@connect.hku.hk, 2024-0627

% Input (2D_spm_preproc.mat):
%   CESTimg: [256,256,18,nf]
%   GRE_FA20/7: [256,256,18]
%   T1w: [256,256,18]
%   roi_csf/gm/wm_uint: [256,256,18], from 0-255

% Ouput (eg. preproc_hc_2025xxxx.mat)
%   anatomy: T1w, T1map [s], roi masks
%   CEST: offs_ppm, CESTM0, db0, zspec


clear
% close all

addpath ..\toolbox\

%% predefine paramter
dateDirList = ["\ms\20250306\"];
dateDirList = "D:\InvivoData\GE\MS\MotionCorr(8MS8HC)" + dateDirList;
parallel_flag = 1;
motioncorList = ["none","multi","mono"];


outputFolder = ".\preprocessedData_motionComp\";
if ~exist(outputFolder,'dir')
    mkdir(outputFolder);
end

if parallel_flag == 1
    delete(gcp('nocreate'));
    parpool(16);
end

for idxMotion = 1:3
    motioncor_flag = motioncorList(idxMotion);

    for fileidx = 1:length(dateDirList)
        dataDir = dateDirList(fileidx);
        
        dataID = strsplit(dataDir,'\');
        dataID = dataID(end-2)+"_"+dataID(end-1);
        fprintf("preprocessing "+dataID+"_"+motioncor_flag+"\n");
        
        %% (1) load data
        fprintf('\tload data from SPM\n')
        load(dataDir+"Apro\2D_spm_preproc.mat","T1w","roi_gm_uint","roi_wm_uint","roi_csf_uint",...
                    "WASABIimg", "GRE_FA7", "GRE_FA20","CEST_data");
     
        %% (2) mask of brain tissue
        % mask from SPM ranges from 0 - 1, probability maps
        maskth = 0.3;
        roi_csf = double(roi_csf_uint >= maskth);
        roi_gm = double(roi_gm_uint >= maskth);
        roi_wm = double(roi_wm_uint >= maskth);
        roi = double(logical(roi_csf + roi_gm + roi_wm));  % merge
    
        %% (3) T1mapping using Dual-Flip Angle (DFA)
        fprintf('\tT1mapping based on DFA\n')
        
        FA1 = 20;
        FA2 = 7;
        img1 = GRE_FA20;
        img2 = GRE_FA7;
        TR_s = 0.03; % TR = 30 ms
        
        R = img1./img2; % img1/img2
        FA1_rad = FA1/180*pi;
        FA2_rad = FA2/180*pi;
        R(roi==0) = 1; 
        Temp = (R*sin(FA2_rad)*cos(FA1_rad)-sin(FA1_rad)*cos(FA2_rad))./(R*sin(FA2_rad)-sin(FA1_rad));
        
        % refine constraint T1 map
        roi(Temp<=1) = 0;
        roi_csf(Temp<=1) = 0;
        roi_gm(Temp<=1) = 0;
        roi_wm(Temp<=1) = 0;
        Temp(Temp<=1) = 1 + eps;
        
        T1map = TR_s./log( Temp ); % in unit [s]
        T1map(roi==0) = 0;
        T1map = min(T1map,10);
    
        %% (4) rearrange Z spectrum of WASABI data
        fprintf('\tWASABI fitting: ')
        tic
        M0idx = 2;
        B0 = 3; % in unit [T]
        gamma_Hz = 42.58; % in unit [Hz/uT]
        load(".\offs_24p_Hz_user18_WASABI.mat","offs");
        
        offs_WASABI = offs/(B0*gamma_Hz); % in unit [ppm], in the order of scanning protocol
        CESTidx = find(abs(offs_WASABI) <= 100); 
        offs_WASABI = offs_WASABI(CESTidx);
        [offs_WASABI, idx_sort_order] = sort(offs_WASABI, 'descend');
        CESTidx = CESTidx(idx_sort_order);
        
        WASABI_data.m0 = WASABIimg(:,:,:,M0idx);
        WASABI_data.raw = WASABIimg(:,:,:,CESTidx);
    
        indnz = find(reshape(roi,[],1)~=0);
        zmap_mea = WASABIimg(:,:,:,CESTidx)./WASABIimg(:,:,:,M0idx).*roi;
        zspec_vec = reshape(permute(zmap_mea,[4,1,2,3]),length(offs_WASABI),[]); % [nf,Npixel]
        zspec_nz_vec = zspec_vec(:,indnz);
        if parallel_flag == 1
            [cnz_vec, afnz_vec, B1nz_vec, db0nz_vec] = calWASABI_par(zspec_nz_vec,offs_WASABI);
        else
            [cnz_vec, afnz_vec, B1nz_vec, db0nz_vec] = calWASABI(zspec_nz_vec,offs_WASABI,1);
        end
        
        map3DTemp = zeros(size(roi));
        map3DTemp(indnz) = cnz_vec; WASABI_data.cmap = map3DTemp;
        map3DTemp(indnz) = afnz_vec; WASABI_data.afmap = map3DTemp;
        map3DTemp(indnz) = B1nz_vec; WASABI_data.B1map = map3DTemp;
        map3DTemp(indnz) = db0nz_vec; WASABI_data.db0map = map3DTemp;
    
        fprintf('%.4f s\n', toc);
    
        %% (5) rearrange Z spectrum of CEST data
        M0idx = 4; % select the 3rd scan as M0 scan
        B0 = 3; % in unit [T]
        gamma_Hz = 42.58; % in unit [Hz/uT]
        load(".\offs_77p_Hz_user17.mat","offs");
        
        offs_CEST = offs/(B0*42.58); % in unit [ppm], in the order of scanning protocol
        CESTidx = find(abs(offs_CEST) <= 100); 
        offs_CEST = offs_CEST(CESTidx);
        [offs_CEST, idx_sort_order] = sort(offs_CEST, 'descend');
        CESTidx = CESTidx(idx_sort_order);
        
        for idx = 1:length(CEST_data)
            CESTimg = CEST_data(idx).raw;
            CEST_data(idx).raw_re = squeeze(CESTimg(:,:,:,[M0idx;CESTidx]));
        end
        CEST_data = rmfield(CEST_data, 'raw');
    
        % motion correction
        switch motioncor_flag
            case "none"
                motionCorData = [];
                fprintf('\tSkip motion correction\n')
            case "multi"
                [CEST_data(1:3),motionCorData] = motionCorCEST(CEST_data(1:3),...
                    outputFolder+"fig_motioncor_multi_"+dataID+".png","multi"); % revise field 'raw_re'
            case "mono"
                [CEST_data(1:3),motionCorData] = motionCorCEST(CEST_data(1:3),...
                    outputFolder+"fig_motioncor_mono_"+dataID+".png","mono");
        end
    
        for idx = 1:length(CEST_data)
            CESTimg = CEST_data(idx).raw_re;
            CEST_data(idx).m0 = CESTimg(:,:,1);
            CEST_data(idx).zspec = CESTimg(:,:,2:end)./(CEST_data(idx).m0+eps).*roi;
        end
        CEST_data = rmfield(CEST_data, 'raw_re');
    
        %% QUASS
        % parallel_flag = 1;
        Tacq = 0; %s
        T1cor_flag = 0;
        B1cor_flag = 0;
        B0cor_flag = 1;
        
        if parallel_flag == 1
            % delete(gcp('nocreate'));
            % parpool(16);
            QUASSpro_methods = @QUASSprocess_par;
        else
            QUASSpro_methods = @QUASSprocess;
        end
        
        fprintf('\tQUASS: \n')
        
        rB1map_norm = WASABI_data.B1map./3.7; % from 0 to 1
        if B0cor_flag == 1
            db0map_vec = reshape(WASABI_data.db0map,1,[]); % ppm
        else
            db0map_vec = reshape(roi * 0,1,[]);
        end
        if T1cor_flag == 1
            R1w_vec = reshape( 1./min(T1map+eps,10),1,[]);
        else
            R1w_vec = reshape(roi * 1,1,[]);
        end
        
        % ============== QUASS ==============
        for i = 1:3
            fprintf("\t\tSat#"+num2str(i)+": ")
            zmap_mea = CEST_data(i).zspec; % [nx,ny,nf]
            [nx,ny,nf] = size(zmap_mea);
            if B1cor_flag == 1
                B1 = CEST_data(i).B1 * rB1map_norm; % nominal
            else
                B1 = CEST_data(i).B1 * roi;
            end
            Ts = CEST_data(i).Ts * roi;
            Trec = (CEST_data(i).TR - Ts - Tacq) * roi;
            
            indnz = find(reshape(roi,[],1)~=0);
        
            % gamma_hz = 267.5153/2/pi; % for protons, in 1e6 Hz/T
            % offs_hz = offs_CEST*3*gamma_hz;
        
            SATpara = [B1(indnz)'; Ts(indnz)'; Trec(indnz)']; % [B1(uT), tsat(s), tdelay(s)]
            
            zspecList = reshape( permute(zmap_mea,[3,1,2]), nf,[]);
            [zspecQUASSList, R1rhoList] = QUASSpro_methods(zspecList(:,indnz), offs_CEST, SATpara, R1w_vec(:,indnz),db0map_vec(:,indnz));
            zmap_mea_QUASS = zeros(nf,nx*ny);
            zmap_mea_QUASS(:,indnz) = zspecQUASSList;
            zmap_mea_QUASS = permute(reshape(zmap_mea_QUASS,nf,nx,ny),[2,3,1]); % [nx,ny,nf]
            
            CEST_data(i).QUASSzspec = zmap_mea_QUASS;
        end
    
        %% save the preprocessed data
        fprintf('\tsave to .mat file\n')
        save(outputFolder+"QUASS_preproc_"+dataID+"_"+motioncor_flag+".mat", ...
                    "T1w","T1map","roi_csf","roi_gm","roi_wm","roi", ...
                    "offs_WASABI","offs_CEST", ...
                    "WASABI_data","CEST_data","motionCorData");
    end
end

if parallel_flag == 1
    delete(gcp('nocreate'));
end
fprintf('done!\n')

%% function
function [cnz_vec, afnz_vec, B1nz_vec, db0nz_vec] = calWASABI(zspecnz_vec,offs,displayflag)
% INPUT:
%   zspec_vec: [nf, Npixel], not zero
%   offs

    if ~exist('displayflag','var') || isempty(displayflag)
        displayflag = 0;
    end
    multistart_N = 5;
    B0 = 3; % main field of scanner
    tp = 5; % ms
    tp_s = tp*10^-3; % s
    GAMMA = 42.576375; %(Hz/uT) = (MHz/T)
    FREQ = B0 * GAMMA; % scanner frequency in MHz
    
    % WASABImodel = @(c,af,B1,db0,xx) abs(c - af * reshape( sin(atan((B1/((FREQ/GAMMA)))./(xx-db0))).^2,[],1) .* ...
    %     reshape(sin(sqrt((B1/((FREQ/GAMMA))).^2+(xx-db0).^2) *FREQ*(2*pi)*tp_s/2).^2,[],1)  );
    WASABImodel = @(par,xx) abs(par(1) - par(2) * reshape( sin(atan((par(3)/((FREQ/GAMMA)))./(xx-par(4)))).^2,[],1) .* ...
        reshape(sin(sqrt((par(3)/((FREQ/GAMMA))).^2+(xx-par(4)).^2) *FREQ*(2*pi)*tp_s/2).^2,[],1)  );
        
    %      c      af     B1    db0       
    iv = [ 1.0    0.9   3.7    0  ];
    lb = [ 0.2    0.5   2.6    -1 ];
    ub = [ 1.0    2.0   4.8    +1 ];
    
    [nf,Npixel] = size(zspecnz_vec);
    cnz_vec = zeros(Npixel,1);
    afnz_vec = zeros(Npixel,1);
    B1nz_vec = zeros(Npixel,1);
    db0nz_vec = zeros(Npixel,1);
    
    runMS = @(prob, n) run(MultiStart('Display','off'), prob, n);
    backNum = 0;
    for idx = 1:Npixel
        if displayflag == 1
            fprintf(repmat('\b',1,backNum));
            backNum = fprintf(' %d/%d',idx,Npixel);
        end

        zspecTemp = zspecnz_vec(:,idx);
        problem = createOptimProblem('lsqcurvefit','x0',iv,'objective',WASABImodel,...
            'lb',lb,'ub',ub,'xdata',offs,'ydata',zspecTemp);
        [par, ~] = runMS(problem, multistart_N);
        
        cnz_vec(idx) = par(1);
        afnz_vec(idx) = par(2);
        B1nz_vec(idx) = par(3);
        db0nz_vec(idx) = par(4);
    end
    fprintf(repmat('\b',1,backNum));

end

function [cnz_vec, afnz_vec, B1nz_vec, db0nz_vec] = calWASABI_par(zspecnz_vec,offs)
% INPUT:
%   zspec_vec: [nf, Npixel], not zero
%   offs

    multistart_N = 5;
    B0 = 3; % main field of scanner
    tp = 5; % ms
    tp_s = tp*10^-3; % s
    GAMMA = 42.576375; %(Hz/uT) = (MHz/T)
    FREQ = B0 * GAMMA; % scanner frequency in MHz
    
    % WASABImodel = @(c,af,B1,db0,xx) abs(c - af * reshape( sin(atan((B1/((FREQ/GAMMA)))./(xx-db0))).^2,[],1) .* ...
    %     reshape(sin(sqrt((B1/((FREQ/GAMMA))).^2+(xx-db0).^2) *FREQ*(2*pi)*tp_s/2).^2,[],1)  );
    WASABImodel = @(par,xx) abs(par(1) - par(2) * reshape( sin(atan((par(3)/((FREQ/GAMMA)))./(xx-par(4)))).^2,[],1) .* ...
        reshape(sin(sqrt((par(3)/((FREQ/GAMMA))).^2+(xx-par(4)).^2) *FREQ*(2*pi)*tp_s/2).^2,[],1)  );
        
    %      c      af     B1    db0       
    iv = [ 1.0    0.9   3.7    0  ];
    lb = [ 0.2    0.5   2.6    -1 ];
    ub = [ 1.0    2.0   4.8    +1 ];
    
    [nf,Npixel] = size(zspecnz_vec);
    cnz_vec = zeros(Npixel,1);
    afnz_vec = zeros(Npixel,1);
    B1nz_vec = zeros(Npixel,1);
    db0nz_vec = zeros(Npixel,1);
    
    runMS = @(prob, n) run(MultiStart('Display','off'), prob, n);
    parfor idx = 1:Npixel
        zspecTemp = zspecnz_vec(:,idx);
        problem = createOptimProblem('lsqcurvefit','x0',iv,'objective',WASABImodel,...
            'lb',lb,'ub',ub,'xdata',offs,'ydata',zspecTemp);
        [par, ~] = runMS(problem, multistart_N);
        
        cnz_vec(idx) = par(1);
        afnz_vec(idx) = par(2);
        B1nz_vec(idx) = par(3);
        db0nz_vec(idx) = par(4);
    end

end

function [CEST_data, motionCorData] = motionCorCEST(CEST_data,figdir,cortype)
% INPUT:
%   CEST_data: [1,3] structure
%       CEST_data(i).raw_re: [256,256,nf]
% OUTPUT:
%       CEST_data(i).raw_re: [256,256,nf]

    nof = size(CEST_data(1).raw_re,3);
    nfra = 3*nof;
    imgSize = size(CEST_data(1).raw_re,1);
    imgSeq_ori = cell(nof,3);
    for i = 1:3
        imgArr = squeeze(CEST_data(i).raw_re); % [256,256,73]
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
    switch cortype
        case "multi"
            [optimizer, metric] = imregconfig('multimodal');optimizer.MaximumIterations = 500;
        case "mono"
            [optimizer, metric] = imregconfig('monomodal');optimizer.MaximumIterations = 500;
    end
    
    % First-pass motion correction
    fprintf('\t1st motion correction: ')
    parfor i = 2:nfra
        % fprintf(repmat('\b',1,backNum));
        % backNum = fprintf('%d/%d',i,nfra);
        tform = imregtform(imgSeq_ori{i}, baseImg, 'rigid', optimizer, metric);
            % Compute transforms between frames; find optimal rigid transform tform such that tform @ imgmov = imgref
            % optimizer: optimization algorithm, default regm (regular step gradient descent)
            % metric: similarity metric; use mean square for single-modality; change for different contrasts
        imgSeq_corrected{i} = imwarp(imgSeq_ori{i}, tform, 'OutputView', imref2d(size(baseImg)));
        calc_dx(i) = tform.T(3,1);
        calc_dy(i) = tform.T(3,2);
        % calc_rot(i) = atan2d(tform.T(2,1), tform.T(1,1));
    end
    fprintf('\n')
    
    % Second-pass motion correction - compute residual displacement
    baseImg_cor = imgSeq_corrected{1};
    residual_dx = zeros(nfra,1);
    residual_dy = zeros(nfra,1);
    fprintf('\t2nd motion correction: ')
    parfor i = 2:nfra
        % fprintf(repmat('\b',1,backNum));
        % backNum = fprintf('%d/%d',i,nfra);
        tform = imregtform(imgSeq_corrected{i}, baseImg_cor, 'rigid', optimizer, metric);
        residual_dx(i) = tform.T(3,1);
        residual_dy(i) = tform.T(3,2);
    end
    fprintf('\n')

    motionCorData.nfra = nfra;
    motionCorData.firstCor_clacdx = calc_dx;
    motionCorData.firstCor_clacdy = calc_dy;
    motionCorData.secondCor_clacdx = residual_dx;
    motionCorData.secondCor_clacdy = residual_dy;
    
    %% Visualization of results
    Fig = figure();set(gcf,'Position', [100 100 600 300]);
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

    exportgraphics(Fig, figdir, 'BackgroundColor', 'white', 'Resolution', 300);
    %% Restore to zspec
    for i = 1:3
        CEST_data(i).raw_re = imgSeqArr_cor(:,:,((i-1)*nof+1):(i*nof)); % [256,256,nf]
    end
end