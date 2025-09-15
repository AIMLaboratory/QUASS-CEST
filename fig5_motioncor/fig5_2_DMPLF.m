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

OutputDir = ".\Outprocess_motionComp";
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
dateDir = dir(".\preprocessedData_motionComp\QUASS_preproc_*.mat");
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

    for satparaidx = 1:3
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
