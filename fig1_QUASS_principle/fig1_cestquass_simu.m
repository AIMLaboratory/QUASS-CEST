clc; clear
% close all;
% Jianpan Huang - jianpanhuang@outlook.com, 20241118
% modified by Huabin Zhang, 20250225

addpath ..\toolbox\
% demonstrate the capabilities of QUASS in zspec correction

%% simulation
% (0) saturation settings
% [B1(uT), tsat(s), tdelay(s)]
SATpara = [
    [0.5, 0.5, 0.5];...
    [0.5, 1,   1];...
    [0.5, 2,   2];...
    [0.5, 8,   8];...
    ];
SATtags = [
    "0.5s/0.5s";...
    "1s/1s";...
    "2s/2s";...
    "8s/8s";...
    ];  
SATpara = permute(SATpara,[2,1]); % [3,nzspec]
nzspec = size(SATpara,2);

offs = [-300,-300,-300,-300,-5.5:0.1:5.5]';
M0idx = 4;
nf = length(offs);
% offs = [-300;offs(1+randperm(nf-1))];
zspecList = zeros(nf-M0idx,nzspec);
for idxZ = 1:nzspec
    % (1) scanner settings
    B0 = 3;
    gamma = 267.5154109126009;
    gamma_hz = gamma/2/pi;
    
    % (2) saturation rf pulse 
    pulse1_pwr = SATpara(1,idxZ); % in uT
    pulse1_dur = SATpara(2,idxZ); % pulse duration in s, ts
    Td = SATpara(3,idxZ); % tdelay or trecovery in s, td

    pulse1_phase = zeros(size(pulse1_dur)); % pulse duration in s
    pulse1 = [pulse1_pwr*gamma_hz, pulse1_phase, pulse1_dur];
    pulse2_tdelay = 0.00; % pulse delay in s
    pulse2 = [0, 0, pulse2_tdelay];
    pulse_block = {pulse1, pulse2};
    pulse_block_len = length(pulse_block);
    pulse_repeat = 1; % repeat number of pulse cell
    pulse_cell=cell(1,pulse_block_len*pulse_repeat);
    for idxoffs=1:pulse_block_len:size(pulse_cell,2)
        pulse_cell(1,idxoffs:idxoffs+pulse_block_len-1) = pulse_block;
    end
    last_delay = 0; % 1 ends with delay, 0 ends with pulse (no last delay)
    if last_delay == 0
        pulse_cell(end) = [];
    end
    pulse_tpost = 0;
    
    % (3) exchange settings
    %       {name,            t1 [s],   t2 [s],    exch rate [Hz],  dw [ppm],    fraction (0~1)}
    water  = {'water',        1,      0.04       1,               0,           1};
    mt     = {'mt',           1.0,      4.0e-05,   30,              -2.5,        0.1351};
    amide  = {'amide',        1.0,      0.1        50,            3.5,       0.0009009*2};
    guanid = {'guanidine',    1.0,      0.1,       200,            2.0,         0.0009009};
    noe    = {'noe',          1.3,      0.005,     20,              -3.5,        0.0045};
    pools = {water; mt; amide; guanid; noe};
    
    R1w = 1/water{2};
    R2w = 1/water{3};
    
    % (4) BME solver for Z-spectrum
    zspec = zeros(nf,1);
    % =============== use R1w decay for tdelay ==============
    % for n = 1:nf
    %     offs_temp = offs(n);
    %     if n==1
    %         magn0corr=1;
    %     else
    %         magn0corr=(1-zspec(n-1,1))*(1-exp(-R1w*tdelay))+zspec(n-1,1);
    %     end
    %     magn = bmesolver_magn0corr(b0, gamma_hz, pools, pulse_cell, pulse_tpost, offs_temp, magn0corr);
    %     zspec(n,1) = magn(length(pools)*2+1, end, end);
    % end

    for n = 1:nf
        offs_temp = offs(n);

        % (1) CEST saturation, R1rho decay
        if n == 1
            magn = bmesolver(B0, gamma_hz, pools, pulse_cell, pulse_tpost, offs_temp, 0);
        else
            magn = bmesolver(B0, gamma_hz, pools, pulse_cell, pulse_tpost, offs_temp, magn_rec);
        end
        % (2) FA = 90, Mz --> Mx
        magn_spoil = magn*0; % Mx=My=0 due to spoiling after CEST and before FA90 excitation
        magn_spoil(1:length(pools)) = magn(end-length(pools)+1:end); % Mx after FA90 excitation

        % (3) readout, signal propotional to Mxy
        zspec(n,1) = magn_spoil(1); % Mx_water

        % (4) R1 relaxation
        magn_rec = bmesolver(B0, gamma_hz, pools, [], Td, 0, magn_spoil);
        % magn_rec(end:end-3) = 0;
        % magn_rec = 1 - exp(-1*Td);
        
    end

    % [offs,ind] = unique(offs);
    % zspec = zspec(ind);

    zspecList(:,idxZ) = zspec((M0idx+1):end)/zspec(M0idx);
end
offs(1:M0idx) = []; % eliminate M0
clearvars -except zspecList offs B0 R1w R2w SATpara SATtags

%% QUASS process
offs_ppm = offs;
[zspecQUASSList, R1rhoList] = QUASSprocess(zspecList, offs_ppm, SATpara, R1w);

Fig1 = figure();set(gcf,'Position',[150 350 1200 400]);
tiledlayout(1,3,"TileSpacing","compact","Padding","compact")

ax1 = nexttile;
plot(offs_ppm,zspecList,'LineWidth',1);
xlabel('offs [ppm]');set(gca,'xdir','reverse')
ylabel('Z [a.u.]')
title('Raw Zspectra','FontName','Times New Roman')
xlim([-5,5]);ylim([0,1])
set(ax1,'FontSize',14);
legend(SATtags(1),SATtags(2),SATtags(3),SATtags(4),'Location','southeast','fontsize',10)

ax2 = nexttile;
plot(offs_ppm,zspecQUASSList,'LineWidth',1);
xlabel('offs [ppm]');set(gca,'xdir','reverse')
ylabel('Z [a.u.]')
title('QUASS Zspectra','FontName','Times New Roman')
xlim([-5,5]);ylim([0,1])
set(ax2,'FontSize',14);
% legend(SATtags(1),SATtags(2),SATtags(3),SATtags(4),'Location','southeast','fontsize',10)

ax3 = nexttile;
plot(offs_ppm,R1rhoList,'LineWidth',1);
xlim([-5,5])
ylim([1,2])
xlabel('offs [ppm]');set(gca,'xdir','reverse')
title('Fitted R1rho','FontName','Times New Roman')
ylabel('R1rho [s^{-1}]')
set(ax3,'FontSize',14);
% legend(SATtags(1),SATtags(2),SATtags(3),SATtags(4),'Location','southeast','fontsize',10)

exportgraphics(Fig1, ".\fig2_1DZspec_QUASS.png", 'BackgroundColor', 'white', 'Resolution', 600);