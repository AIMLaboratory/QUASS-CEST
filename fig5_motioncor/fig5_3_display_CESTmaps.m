
addpath ..\toolbox\

datafolder = ".\Outprocess_motionComp_ms20250316";
preprocfolder = ".\preprocessedData_motionComp_ms20250316";

% datafolder = ".\Outprocess_motionComp_hc20250214_2";
% preprocfolder = ".\preprocessedData_motionComp_hc20250214_2";
dateDir = dir(datafolder);
dateDir = dateDir([dateDir.isdir] & ~ismember({dateDir.name}, {'.', '..','figs_MPLFmaps_hist'}));
if ~exist(datafolder+"\figs_MPLFmaps_hist",'dir')
    mkdir(datafolder+"\figs_MPLFmaps_hist")
end
for idx = 1:length(dateDir)
    dataID = dateDir(idx).name;
    load(preprocfolder+"\"+dataID+".mat","motionCorData")
    load(datafolder+"\"+dataID+"\CESTinput.mat","CEST_data","WASABI_data")
    load(datafolder+"\"+dataID+"\output_DMPLF.mat","CESTfit_data","offs","roi")
    legList = ["2.0s/3.0s", "1.5s/2.3s", "1.0s/1.6s"];

    % B1 = 0.5uT
    histrange = [0.005,0.15];
    colranArr = [[0,0.06];[0,0.05];[0,0.10];[0,0.15]]; % amide, guan, rNOE, MT
    MPLFresultComp_2D(CESTfit_data(1:3),colranArr,datafolder+"\figs_MPLFmaps_hist\figs_B1_0.5uT_maps_"+dataID+".png");
    motionCorDisplay(motionCorData,datafolder+"\figs_MPLFmaps_hist\figs_B1_0.5uT_motion_"+dataID+".png");
end

function motionCorDisplay(motionCorData, pngdir)
% INPUT:
%   motionCorData: structure variable
%       .nfra
%       .firstCor_clac[dx/dy]
%       .secondCor_clac[dx/dy]
    if ~isempty(motionCorData)
        nfra = motionCorData.nfra;
        calc_dx = motionCorData.firstCor_clacdx;
        calc_dy = motionCorData.firstCor_clacdy;
        residual_dx = motionCorData.secondCor_clacdx;
        residual_dy = motionCorData.secondCor_clacdy;
        Fig = figure();set(gcf,'Position', [100 100 600 300]);
        tiledlayout(1,2,'Padding','compact','TileSpacing','compact')
        % 
        % ======================== before motion correction ===================
        nexttile
        plot(1:nfra, calc_dx, 'r-o', 1:nfra, calc_dy, 'b-s');
        title('Before corretion','FontSize',12); 
        xlabel('image index','FontSize',12); ylabel('voxel displacement','FontSize',12);
        grid on;
        ylim([-3,3])

        % ======================== after motion correction ===================
        nexttile
        plot(1:nfra, residual_dx, 'r-o', 1:nfra, residual_dy, 'b-s');
        legend('X-shift', 'Y-shift', 'Location', 'northwest');
        title('After corretion','FontSize',12); 
        xlabel('image index','FontSize',12);
        grid on;
        ylim([-3,3])
    
        exportgraphics(Fig, pngdir, 'BackgroundColor', 'white', 'Resolution', 600);
    end

end

function MPLFresultComp_2D(CESTfit_data,colranArr,pngdir)
    % CESTfit_data(i).QUASS_fitpara:   [nx,ny,npara]
    % CESTfit_data(i).QUASS_fitpara: [nx,ny,npara]

    %% main window
    % Concatenate 3D to 2D
    [nrow,ncol,npara] = size(CESTfit_data(1).QUASS_fitpara);
    
    % main window
    fig = figure('Position', [100, 50, 600, 800]);
    tiledlayout(4,3,'TileSpacing','none','Padding','loose')
    

    %% display
    aximgHandles = gobjects(3,8);

    for preidx = 1:1
        for satidx = 1:3

            aximgHandles(satidx,1) = nexttile(satidx+0*3);imagesc(CESTfit_data(satidx).QUASS_fitpara(:,:,5));axis off;
            colromapTemp = inferno(nrow); colromapTemp(1,:) = 0;
            colormap(aximgHandles(satidx,1),colromapTemp);clim(colranArr(1,:))
            % title('APT');
        
            aximgHandles(satidx,2) = nexttile(satidx+1*3);imagesc(CESTfit_data(satidx).QUASS_fitpara(:,:,14));axis off;
            colromapTemp = jet(nrow); colromapTemp(1,:) = 0;
            colormap(aximgHandles(satidx,2),colromapTemp);clim(colranArr(2,:))
            % title('CEST@2ppm');
        
            aximgHandles(satidx,3) = nexttile(satidx+2*3);imagesc(CESTfit_data(satidx).QUASS_fitpara(:,:,8));axis off;
            colromapTemp = parula(nrow); colromapTemp(1,:) = 0;
            colormap(aximgHandles(satidx,3),colromapTemp);clim(colranArr(3,:))
            % title('NOE');
        
            aximgHandles(satidx,4) = nexttile(satidx+3*3);imagesc(CESTfit_data(satidx).QUASS_fitpara(:,:,11));axis off;
            colromapTemp = hot(nrow); colromapTemp(1,:) = 0;
            colormap(aximgHandles(satidx,4),colromapTemp);clim(colranArr(4,:))
            % title('MT');
        end

    end

    ypos = [0.76,0.55,0.335,0.12];
    for idxpool = 1:4
        ax = nexttile(idxpool*3);
        
        cbh = colorbar(ax);
        cbh.Layout.Tile = 'east';
        
        pos = cbh.Position;
        pos(1) = 0.9152;
        pos(2) = ypos(idxpool); % y position
        pos(3) = 0.02; % width
        pos(4) = 0.15; % height
        set(cbh,'Location','manual','Position', pos); % width can be changed only in 'manual' mode
        cbh.FontSize = 10; 
        cbh.FontWeight = 'bold'; 
    end

    exportgraphics(fig, pngdir, 'BackgroundColor', 'white', 'Resolution', 600);

end
