
addpath ..\toolbox\

datafolder = ".\Outprocess_multi_batch_2";
dateDir = dir(datafolder);
dateDir = dateDir([dateDir.isdir] & ~ismember({dateDir.name}, {'.', '..','figs_MPLFmaps_hist'}));
if ~exist(datafolder+"\figs_MPLFmaps_hist",'dir')
    mkdir(datafolder+"\figs_MPLFmaps_hist")
end
for idx = 1:length(dateDir)
    dataID = dateDir(idx).name;
    load(datafolder+"\"+dataID+"\CESTinput.mat","CEST_data","WASABI_data")
    load(datafolder+"\"+dataID+"\output_DMPLF.mat","CESTfit_data","offs","roi")
    legList = ["2.0s/3.0s", "1.5s/2.3s", "1.0s/1.6s"];

    % B1 = 0.5uT
    histrange = [0.005,0.15];
    colranArr = [[0,0.05];[0,0.05];[0,0.10];[0,0.15]]; % amide, guan, rNOE, MT
    MPLFresultComp_2D(CESTfit_data(1:3),colranArr,datafolder+"\figs_MPLFmaps_hist\figs_B1_0.5uT_maps_"+dataID+".png");
    MPLFHistComp_2D(CESTfit_data(1:3), legList,roi,histrange,datafolder+"\figs_MPLFmaps_hist\figs_B1_0.5uT_hist_"+dataID+".png");

    % B1 = 1.0uT
    histrange = [0.01,0.33];
    colranArr = [[0,0.08];[0,0.07];[0,0.12];[0,0.35]]; % amide, guan, rNOE, MT
    MPLFresultComp_2D(CESTfit_data(4:6),colranArr,datafolder+"\figs_MPLFmaps_hist\figs_B1_1.0uT_maps_"+dataID+".png");
    MPLFHistComp_2D(CESTfit_data(4:6), legList,roi,histrange,datafolder+"\figs_MPLFmaps_hist\figs_B1_1.0uT_hist_"+dataID+".png");

end

function MPLFresultComp_2D(CESTfit_data,colranArr,pngdir)
    % CESTfit_data(i).raw_fitpara:   [nx,ny,npara]
    % CESTfit_data(i).QUASS_fitpara: [nx,ny,npara]


    %% main window
    % Concatenate 3D to 2D
    [nrow,ncol,npara] = size(CESTfit_data(1).raw_fitpara);
    
    % main window
    fig = figure('Position', [100, 50, 1200, 800]);
    tiledlayout(4,6,'TileSpacing','none','Padding','loose')
    

    %% display
    aximgHandles = gobjects(3,8);

    for preidx = 1:2
        for satidx = 1:3
            if preidx == 1
                % (1) show CEST maps
                aximgHandles(satidx,1) = nexttile(satidx+0*6);imagesc(CESTfit_data(satidx).raw_fitpara(:,:,5));axis off;
                colromapTemp = inferno(256); colromapTemp(1,:) = 0;
                colormap(aximgHandles(satidx,1),colromapTemp);clim(colranArr(1,:))
                % title('APT');
            
                aximgHandles(satidx,2) = nexttile(satidx+1*6);imagesc(CESTfit_data(satidx).raw_fitpara(:,:,14));axis off;
                colromapTemp = jet(256); colromapTemp(1,:) = 0;
                colormap(aximgHandles(satidx,2),colromapTemp);clim(colranArr(2,:))
                % title('CEST@2ppm');
            
                aximgHandles(satidx,3) = nexttile(satidx+2*6);imagesc(CESTfit_data(satidx).raw_fitpara(:,:,8));axis off;
                colromapTemp = parula(256); colromapTemp(1,:) = 0;
                colormap(aximgHandles(satidx,3),colromapTemp);clim(colranArr(3,:))
                % title('NOE');
            
                aximgHandles(satidx,4) = nexttile(satidx+3*6);imagesc(CESTfit_data(satidx).raw_fitpara(:,:,11));axis off;
                colromapTemp = hot(256); colromapTemp(1,:) = 0;
                colormap(aximgHandles(satidx,4),colromapTemp);clim(colranArr(4,:))
                % title('MT');
            else
                % (2) show QUASS CEST maps
                aximgHandles(satidx,5) = nexttile(satidx+3+0*6);imagesc(CESTfit_data(satidx).QUASS_fitpara(:,:,5));axis off;
                colromapTemp = inferno(256); colromapTemp(1,:) = 0;
                colormap(aximgHandles(satidx,5),colromapTemp);clim(colranArr(1,:))
                % title('APT');
            
                aximgHandles(satidx,6) = nexttile(satidx+3+1*6);imagesc(CESTfit_data(satidx).QUASS_fitpara(:,:,14));axis off;
                colromapTemp = jet(256); colromapTemp(1,:) = 0;
                colormap(aximgHandles(satidx,6),colromapTemp);clim(colranArr(2,:))
                % title('CEST@2ppm');
            
                aximgHandles(satidx,7) = nexttile(satidx+3+2*6);imagesc(CESTfit_data(satidx).QUASS_fitpara(:,:,8));axis off;
                colromapTemp = parula(256); colromapTemp(1,:) = 0;
                colormap(aximgHandles(satidx,7),colromapTemp);clim(colranArr(3,:))
                % title('NOE');
            
                aximgHandles(satidx,8) = nexttile(satidx+3+3*6);imagesc(CESTfit_data(satidx).QUASS_fitpara(:,:,11));axis off;
                colromapTemp = hot(256); colromapTemp(1,:) = 0;
                colormap(aximgHandles(satidx,8),colromapTemp);clim(colranArr(4,:))
                % title('MT');
            end
        end

    end

    ypos = [0.76,0.55,0.335,0.12];
    for idxpool = 1:4
        ax = nexttile(idxpool*6);
        
        cbh = colorbar(ax);
        cbh.Layout.Tile = 'east';
        
        pos = cbh.Position;
        pos(2) = ypos(idxpool); % y position
        pos(3) = 0.02; % width
        pos(4) = 0.15; % height
        set(cbh,'Location','manual','Position', pos); % width can be changed only in 'manual' mode
        cbh.FontSize = 12; 
        cbh.FontWeight = 'bold'; 
    end

    exportgraphics(fig, pngdir, 'BackgroundColor', 'white', 'Resolution', 600);

end

function MPLFHistComp_2D(CESTfit_data, legList, roi, xrange, pngdir)
    % CESTfit_data(i).raw_fitpara:   [nx,ny,npara]
    % CESTfit_data(i).QUASS_fitpara: [nx,ny,npara]


    %% main window
    
    % main window
    fig = figure('Position', [100, 50, 800, 800]);
    tiledlayout(4,2,'TileSpacing','compact','Padding','compact')
    
    %% display
    poolIdxList = [5,14,8,11]; % amide, guan, NOE, MT
    % xrange = [0.005,0.15];
    leg_on = 1;
    for poolidx = 1:4
        ax = nexttile;
        ylabel_on = 1;
        plot3ImgHisto(ax,legList,xrange,leg_on,ylabel_on,roi,...
                    CESTfit_data(1).raw_fitpara(:,:,poolIdxList(poolidx)),...
                    CESTfit_data(2).raw_fitpara(:,:,poolIdxList(poolidx)),...
                    CESTfit_data(3).raw_fitpara(:,:,poolIdxList(poolidx)));
        leg_on = 0;

        ax = nexttile;
        ylabel_on = 0;
        plot3ImgHisto(ax,legList,xrange,leg_on,ylabel_on,roi,...
                    CESTfit_data(1).QUASS_fitpara(:,:,poolIdxList(poolidx)),...
                    CESTfit_data(2).QUASS_fitpara(:,:,poolIdxList(poolidx)),...
                    CESTfit_data(3).QUASS_fitpara(:,:,poolIdxList(poolidx)));

    end
    exportgraphics(fig, pngdir, 'BackgroundColor', 'white', 'Resolution', 600);

end

%%
function plot3ImgHisto(ax,legList,xrange,leg_on,ylabel_on,roi,image1,image2,image3)
    % ax = gca;
    % image1 = randn(128,128)*50 + 100;
    % image2 = rand(128,128)*255;       
    % image3 = randn(128,128)*30 + 150; 
    idxnz = find(roi(:)~=0);
    pixels1 = image1(idxnz);
    pixels2 = image2(idxnz);
    pixels3 = image3(idxnz);
    
    hold on;
    
    h1 = histogram(ax, pixels1, 'FaceColor', 'r', 'FaceAlpha', 0.5, 'EdgeColor', 'k');
    h2 = histogram(ax, pixels2, 'FaceColor', 'g', 'FaceAlpha', 0.5, 'EdgeColor', 'k');
    h3 = histogram(ax, pixels3, 'FaceColor', 'b', 'FaceAlpha', 0.5, 'EdgeColor', 'k');
    
    hold off;
    xlim(xrange);
    if leg_on == 1
        legend(legList);
    end
    % xlabel('Pixel Value');
    if ylabel_on == 1
        ylabel('Voxel');
    end
    % title('Pixel Value Distribution Comparison');
    
    grid on;
    box on;
end