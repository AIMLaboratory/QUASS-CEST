% 20250617
% use multi-linear regression to assess dependence on pool fraction
%       compare amide, guan, NOE
%       compare MTR, LDA, DMPLF

clear
% close all
warning off

addpath ..\toolbox\
addpath ..\toolbox\fitting_main\
addpath ..\toolbox\fitting_main_par\

parallel_flag = 1;
multistart_N = 5; % times for MultiStart
fileName = ".\simData_random5000_B1_0.5_Ts_2.mat";
outputfolder = ".\ran5000-B1_0.5_Ts_2\";

if ~exist("outputfolder",'dir')
    mkdir(outputfolder)
end

%% parallel calculation option
if parallel_flag == 1
    delete(gcp('nocreate'));
    parpool(16);
    QUASSpro_methods = @QUASSprocess_par;
    processing_DZ_methods = {@fitting_MTR_par,@fitting_LDA_par, @fitting_DMPLF_par};
else
    QUASSpro_methods = @QUASSprocess;
    processing_DZ_methods = {@fitting_MTR,@fitting_LDA, @fitting_DMPLF};
end
output_filenames = {'MTR','LDA','DMPLF';'MTR_QUASS','LDA_QUASS','DMPLF_QUASS'};

%% load data
load(fileName,"offs","zspecList","SATpara"); % m0 not included
SATpara = permute(SATpara,[2,1]);
    % offs: [nf,1]
    % zspec: [nf, Npixel]
    % SATpara: [3,1], [B1(uT), tsat(s), tdelay(s)]
load(fileName,"paraList",'paraInfo');
    % paraList: [4,Npixel]
    %   MTf_relative, Amidef_mM, Guanf_mM, NOEf_mM

MTf_mM = paraList(1,:);
amidef_mM = paraList(2,:); % mM
guanf_mM = paraList(3,:);
noef_mM = paraList(4,:);

% construct predictor variable matrix for linear regression, [nPixel,nvar]
predictVar = [guanf_mM(:),amidef_mM(:),MTf_mM(:),noef_mM(:)];

%% zspec w/wo QUASS
[nf,Npixel] = size(zspecList);
zspec_vec = zspecList; % [nf, nzspec]

% gamma_hz = 267.5153/2/pi; % for protons, in 1e6 Hz/T
% offs_hz = offs*3*gamma_hz;
offs_ppm = offs;
R1w = 1; % Hz
[zspecQUASS_vec, R1rhoList] = QUASSpro_methods(zspec_vec, offs_ppm, SATpara, R1w);
    
%% (1) Delta Z fitting/analysis
for idxalgo = 1:length(processing_DZ_methods)
    method = processing_DZ_methods{idxalgo};
    for idxprepro = 1:2 % QUASS
        if idxprepro == 1
            zspecTemp_vec = zspec_vec;
        else
            zspecTemp_vec = zspecQUASS_vec;
        end
        
        [Zamide_DZ_vec, Zguan_DZ_vec, Znoe_DZ_vec] = method(offs, zspecTemp_vec, multistart_N);
    
        % construct response variable vector, [nPixel,1]
        responVar_guan = Zguan_DZ_vec(:);
        responVar_amide = Zamide_DZ_vec(:);
        responVar_noe = Znoe_DZ_vec(:);
    
        % linear regression
        model_4var_guan = fitlm(predictVar, responVar_guan,'ResponseVar','metric_guan','PredictorVars',{'guanf','amidef','MTf','NOEf'});
        model_4var_amide = fitlm(predictVar, responVar_amide,'ResponseVar','metric_amide','PredictorVars',{'guanf','amidef','MTf','NOEf'});
        model_4var_noe = fitlm(predictVar, responVar_noe,'ResponseVar','metric_noe','PredictorVars',{'guanf','amidef','MTf','NOEf'});
    
        save(outputfolder+"output_"+output_filenames{idxprepro,idxalgo}+".mat",'model_4var_amide','model_4var_guan','model_4var_noe');
    
        fprintf('\t\tsaved to %s\n', outputfolder+"output_"+output_filenames{idxprepro,idxalgo}+".mat");
    end
end

if parallel_flag == 1
    delete(gcp('nocreate'));
end
fprintf('done!\n')
