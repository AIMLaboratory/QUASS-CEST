% load roi of MS lesion from manually drawn result
%   8 MS subjects, but no MS lesion for ms_20250522

MSlist = [...
    "20250206",...
    "20250213",...
    "20250220",...
    "20250227",...
    "20250306",...
    "20250502",...
    "20250520",...
    "20250522"];

tileWidth = 100; 
tileHeight = 100; 

% Compute final figure size (8 rows x 5 columns)
figWidth = 5 * tileWidth;  
figHeight = 8 * tileHeight;

% Create figure and set size
fig = figure;
fig.Units = 'pixels'; 
fig.InnerPosition = [100, 100, figWidth+50, figHeight+50];  % [x, y, width, height]

tiledlayout(8,5,"TileSpacing","none","Padding","loose")

for idx = 1:length(MSlist)
    MSID = MSlist(idx);
    
    load("..\invivo_code\preprocessedData_multi_batch\QUASS_preproc_ms_"+MSID+".mat",'roi_csf','roi_gm','roi_wm','CEST_data','offs_CEST');
    load("..\invivo_code\MSlesionROI_8ms_20250702\roi_lesion_"+MSID+".mat",'roi_lesion','T1w','T2w','M0');

    % refine mask
    roi_csf(roi_lesion==1) = 0;
    roi_gm(roi_lesion==1) = 0;
    roi_wm(roi_lesion==1) = 0;
    roi_gm(roi_wm==1) = 0;
    roi_gm(roi_csf==1) = 0;
    
    % plot
    nexttile
    imshow(T1w,[]);axis equal;axis off;
    
    nexttile
    imshow(T2w,[]);axis equal;axis off;
    
    nexttile
    imshow(M0,[]);axis equal;axis off;
    
    nexttile
    c = jet(4);
    transp = 0.9;
    masks = {roi_csf, roi_gm, roi_wm, roi_lesion};  
    imshow(zeros(size(M0))); % black background
    axis equal;axis off;
    hold on;
    for i = 1:4
        h = imshow(cat(3, c(i,1)*masks{i}, c(i,2)*masks{i}, c(i,3)*masks{i}));
        set(h, 'AlphaData', transp*masks{i});
        axis equal;axis off;
    end
    hold off

    ax = nexttile;
    [labeled_lesion, numLesions] = bwlabeln(roi_lesion);
    imshow(labeled_lesion,[]); 
    colromapTemp = jet(256); colromapTemp(1,:) = 0;
    colormap(ax,colromapTemp);clim([0,inf]);axis equal

end

ax = nexttile(1);title(ax, 'T1w','FontName', 'Times New Roman');
ax = nexttile(2);title(ax, 'T2w','FontName', 'Times New Roman');
ax = nexttile(3);title(ax, 'CEST@M0','FontName', 'Times New Roman');
ax = nexttile(4);title(ax, 'ROI','FontName', 'Times New Roman');
ax = nexttile(5);title(ax, 'Lesion','FontName', 'Times New Roman');

exportgraphics(fig,".\figS1_MS_lesion_raw.png",'BackgroundColor','white',Resolution=600)