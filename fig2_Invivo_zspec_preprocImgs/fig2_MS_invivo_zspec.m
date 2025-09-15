

MSlist = [...
    "20250206",...
    "20250213",...
    "20250220",...
    "20250227",...
    "20250306",...
    "20250502",...
    "20250520",...
    "20250522"];

Fig = figure();set(gcf,"Position",[100,100,900,700])
tiledlayout(3,3,'TileSpacing','compact','Padding','compact')

MSID = "20250220";

load("..\invivo_code\preprocessedData_multi_batch\QUASS_preproc_ms_"+MSID+".mat",...
                'roi_csf','roi_gm','roi_wm','CEST_data','offs_CEST','T1map','WASABI_data');
load("..\invivo_code\MSlesionROI_8ms_20250702\roi_lesion_"+MSID+".mat",'roi_lesion','T1w','T2w','M0');

% refine mask
roi_csf(roi_lesion==1) = 0;
roi_gm(roi_lesion==1) = 0;
roi_wm(roi_lesion==1) = 0;
roi_gm(roi_wm==1) = 0;
roi_gm(roi_csf==1) = 0;

%% figures display
pixelLoc_WM = [173,144];
pixelLoc_GM = [118,82];
pixelLoc_lesion = [98,155];

nexttile(1)
imshow(T1w,[]);axis equal;axis off;axis tight;
title('T1w','FontName', 'Times New Roman','FontSize',13);

nexttile(3)
imshow(T2w,[]);axis equal;axis off;axis tight;
title('T2w','FontName', 'Times New Roman','FontSize',13);

nexttile(2)
c = jet(4);
transp = 0.6;
masks = {roi_csf, roi_gm, roi_wm, roi_lesion};  
imshow(zeros(size(M0))); % black background
axis equal;axis off;axis tight;
hold on;
for i = 1:4
    h = imshow(cat(3, c(i,1)*masks{i}, c(i,2)*masks{i}, c(i,3)*masks{i}));
    set(h, 'AlphaData', transp*masks{i});
    axis equal;axis off;axis tight;
end
r = 2; x = pixelLoc_WM(1); y = pixelLoc_WM(2);
rectangle('Position', [x-r, y-r, 2*r, 2*r], 'Curvature', [1, 1], 'EdgeColor', [1 1 1], 'FaceColor', [1 1 1], 'LineWidth', 2);
r = 2; x = pixelLoc_GM(1); y = pixelLoc_GM(2);
rectangle('Position', [x-r, y-r, 2*r, 2*r], 'Curvature', [1, 1], 'EdgeColor', [1 1 1], 'FaceColor', [1 1 1], 'LineWidth', 2);
r = 2; x = pixelLoc_lesion(1); y = pixelLoc_lesion(2);
rectangle('Position', [x-r, y-r, 2*r, 2*r], 'Curvature', [1, 1], 'EdgeColor', [1 1 1], 'FaceColor', [1 1 1], 'LineWidth', 2);

hold off
title('ROI','FontName', 'Times New Roman','FontSize',13);

nexttile(4)
imagesc(T1map);axis equal;axis off;axis tight;
clim([0,4]);colorbar;%title('T1 map [s]','FontName', 'Times New Roman','FontSize',13)

nexttile(5);
imagesc(WASABI_data.B1map./3.7*100);axis equal;axis off;axis tight;
colorbar;title('rB1 map [%]','FontName', 'Times New Roman','FontSize',13)

nexttile(6);
imagesc(WASABI_data.db0map);axis equal;axis off;axis tight;
colorbar;title('\deltaB0 map [ppm]','FontName', 'Times New Roman','FontSize',13)

%% zspec plot
nexttile(7);
hold on
x = pixelLoc_GM(1);
y = pixelLoc_GM(2);
r = 4; % 8*8 neighborhood
plot(offs_CEST,squeeze(mean(CEST_data(1).zspec(y-r:y+r, x-r:x+r, :), [1, 2])),'r-','LineWidth',1);
plot(offs_CEST,squeeze(mean(CEST_data(2).zspec(y-r:y+r, x-r:x+r, :), [1, 2])),'g-','LineWidth',1);
plot(offs_CEST,squeeze(mean(CEST_data(3).zspec(y-r:y+r, x-r:x+r, :), [1, 2])),'b-','LineWidth',1);
plot(offs_CEST,squeeze(mean(CEST_data(1).QUASSzspec(y-r:y+r, x-r:x+r, :), [1, 2])),'r:','LineWidth',1);
plot(offs_CEST,squeeze(mean(CEST_data(2).QUASSzspec(y-r:y+r, x-r:x+r, :), [1, 2])),'g:','LineWidth',1);
plot(offs_CEST,squeeze(mean(CEST_data(3).QUASSzspec(y-r:y+r, x-r:x+r, :), [1, 2])),'b:','LineWidth',1);
hold off;
set(gca,'XDir','reverse', 'FontWeight', 'bold', 'FontSize', 7)
xlabel('offs [ppm]','FontSize',10);ylabel('Z-value','FontSize',10)
xlim([-6,6]);ylim([0.5,1]);
title('Zspec in GM','FontName', 'Times New Roman','FontSize',13)
lgd = legend('raw\_#1','raw\_#2','raw\_#3','QUASS\_#1','QUASS\_#2','QUASS\_#3',...
    'Location','southeast','fontsize',7);
lgd.Position(1) = lgd.Position(1)- 0.000;
lgd.Position(2) = lgd.Position(2) - 0.00;

nexttile(8);
hold on
x = pixelLoc_WM(1);
y = pixelLoc_WM(2);
% r = 2;
plot(offs_CEST,squeeze(mean(CEST_data(1).zspec(y-r:y+r, x-r:x+r, :), [1, 2])),'r-','LineWidth',1);
plot(offs_CEST,squeeze(mean(CEST_data(2).zspec(y-r:y+r, x-r:x+r, :), [1, 2])),'g-','LineWidth',1);
plot(offs_CEST,squeeze(mean(CEST_data(3).zspec(y-r:y+r, x-r:x+r, :), [1, 2])),'b-','LineWidth',1);
plot(offs_CEST,squeeze(mean(CEST_data(1).QUASSzspec(y-r:y+r, x-r:x+r, :), [1, 2])),'r:','LineWidth',1);
plot(offs_CEST,squeeze(mean(CEST_data(2).QUASSzspec(y-r:y+r, x-r:x+r, :), [1, 2])),'g:','LineWidth',1);
plot(offs_CEST,squeeze(mean(CEST_data(3).QUASSzspec(y-r:y+r, x-r:x+r, :), [1, 2])),'b:','LineWidth',1);
hold off;
set(gca,'XDir','reverse', 'FontWeight', 'bold', 'FontSize', 7)
xlabel('offs [ppm]','FontSize',10);ylabel('Z-value','FontSize',10)
xlim([-6,6]);ylim([0.5,1]);
title('Zspec in WM','FontName', 'Times New Roman','FontSize',13)

nexttile(9);
hold on
x = pixelLoc_lesion(1);
y = pixelLoc_lesion(2);
% r = 2;
plot(offs_CEST,squeeze(mean(CEST_data(1).zspec(y-r:y+r, x-r:x+r, :), [1, 2])),'r-','LineWidth',1);
plot(offs_CEST,squeeze(mean(CEST_data(2).zspec(y-r:y+r, x-r:x+r, :), [1, 2])),'g-','LineWidth',1);
plot(offs_CEST,squeeze(mean(CEST_data(3).zspec(y-r:y+r, x-r:x+r, :), [1, 2])),'b-','LineWidth',1);
plot(offs_CEST,squeeze(mean(CEST_data(1).QUASSzspec(y-r:y+r, x-r:x+r, :), [1, 2])),'r:','LineWidth',1);
plot(offs_CEST,squeeze(mean(CEST_data(2).QUASSzspec(y-r:y+r, x-r:x+r, :), [1, 2])),'g:','LineWidth',1);
plot(offs_CEST,squeeze(mean(CEST_data(3).QUASSzspec(y-r:y+r, x-r:x+r, :), [1, 2])),'b:','LineWidth',1);
hold off;
set(gca,'XDir','reverse', 'FontWeight', 'bold', 'FontSize', 7)
xlabel('offs [ppm]','FontSize',10);ylabel('Z-value','FontSize',10)
xlim([-6,6]);ylim([0.5,1]);
title('Zspec in lesion','FontName', 'Times New Roman','FontSize',13)

exportgraphics(Fig, "fig2_invivo_Zspec_MS_"+MSID+".png", 'BackgroundColor', 'white', 'Resolution', 600);
