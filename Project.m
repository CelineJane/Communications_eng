clc; close all; clear all;

density = 500;
SNR = linspace(0,60,19);

% [ND_MRC_DL ND_ZF_DL ND_MMSE_DL] = NDnointercell(density,'dl');
% [ND_MRC_UL ND_ZF_UL ND_MMSE_UL] = NDnointercell(density,'ul');

[ND_MRC_DL_IC ND_ZF_DL_IC ND_MMSE_DL_IC] = NDintercell(density,'dl');
% [ND_MRC_UL_IC ND_ZF_UL_IC ND_MMSE_UL_IC] = NDintercell(density,'ul');

% [D_MRC_DL D_ZF_DL D_MMSE_DL] = Dnointercell(density,'dl');
% [D_MRC_UL D_ZF_UL D_MMSE_UL] = Dnointercell(density,'ul');
% 
[D_MRC_DL_IC D_ZF_DL_IC D_MMSE_DL_IC] = Dintercell(density,'dl');
% [D_MRC_UL_IC D_ZF_UL_IC D_MMSE_UL_IC] = Dintercell(density,'ul');

% figure
% hold on
% plot(SNR,ND_MRC_UL_IC,'.-.','Color',[0.6350    0.0780    0.1840],'Linewidth',2)
% plot(SNR,ND_ZF_UL_IC,'.-.','Color',[0.3010    0.7450    0.9330],'Linewidth',2)
% plot(SNR,ND_MMSE_UL_IC,'.-.','Color',[ 0.4660    0.6740    0.1880],'Linewidth',2)
% plot(SNR,D_MRC_UL_IC,':','Color',[0.6350    0.0780    0.1840],'Linewidth',2)
% plot(SNR,D_ZF_UL_IC,':','Color',[0.3010    0.7450    0.9330],'Linewidth',2)
% plot(SNR,D_MMSE_UL_IC,':','Color',[ 0.4660    0.6740    0.1880],'Linewidth',2)
% title('MU-MIMO UL Performance with Intercell Interference')
% xlabel('SNR (dB)')
% ylabel('Sum Rate (bits/s/Hz)')
% legend('ND MRC','ND ZF','ND MMSE','D MRC','D ZF','D MMSE','Location','northwest')
% hold off

figure
hold on 
plot(SNR,ND_MRC_DL_IC,'.-.','Color',[0.6350    0.0780    0.1840],'Linewidth',2)
plot(SNR,ND_ZF_DL_IC,'.-.','Color',[0.3010    0.7450    0.9330],'Linewidth',2)
plot(SNR,ND_MMSE_DL_IC,'.-.','Color',[ 0.4660    0.6740    0.1880],'Linewidth',2)
plot(SNR,D_MRC_DL_IC,':','Color',[0.6350    0.0780    0.1840],'Linewidth',2)
plot(SNR,D_ZF_DL_IC,':','Color',[0.3010    0.7450    0.9330],'Linewidth',2)
plot(SNR,D_MMSE_DL_IC,':','Color',[ 0.4660    0.6740    0.1880],'Linewidth',2)
title('MU-MIMO DL Performance with Intercell Interference')
xlabel('SNR (dB)')
ylabel('Sum Rate (bits/s/Hz)')
legend('ND MRC','ND ZF','ND MMSE','D MRC','D ZF','D MMSE','Location','northwest')
hold off

% figure
% hold on
% plot(SNR,ND_MRC_UL,'.-.','Color',[0.6350    0.0780    0.1840],'Linewidth',2)
% plot(SNR,ND_ZF_UL,'.-.','Color',[0.3010    0.7450    0.9330],'Linewidth',2)
% plot(SNR,ND_MMSE_UL,'.-.','Color',[ 0.4660    0.6740    0.1880],'Linewidth',2)
% plot(SNR,D_MRC_UL,':','Color',[0.6350    0.0780    0.1840],'Linewidth',2)
% plot(SNR,D_ZF_UL,':','Color',[0.3010    0.7450    0.9330],'Linewidth',2)
% plot(SNR,D_MMSE_UL,':','Color',[ 0.4660    0.6740    0.1880],'Linewidth',2)
% title('MU-MIMO UL Performance without Intercell Interference')
% xlabel('SNR (dB)')
% ylabel('Sum Rate (bits/s/Hz)')
% legend('ND MRC','ND ZF','ND MMSE','D MRC','D ZF','D MMSE','Location','northwest')
% hold off
% 
% figure
% hold on 
% plot(SNR,ND_MRC_DL,'.-.','Color',[0.6350    0.0780    0.1840],'Linewidth',2)
% plot(SNR,ND_ZF_DL,'.-.','Color',[0.3010    0.7450    0.9330],'Linewidth',2)
% plot(SNR,ND_MMSE_DL,'.-.','Color',[ 0.4660    0.6740    0.1880],'Linewidth',2)
% plot(SNR,D_MRC_DL,':','Color',[0.6350    0.0780    0.1840],'Linewidth',2)
% plot(SNR,D_ZF_DL,':','Color',[0.3010    0.7450    0.9330],'Linewidth',2)
% plot(SNR,D_MMSE_DL,':','Color',[ 0.4660    0.6740    0.1880],'Linewidth',2)
% title('MU-MIMO DL Performance without Intercell Interference')
% xlabel('SNR (dB)')
% ylabel('Sum Rate (bits/s/Hz)')
% legend('ND MRC','ND ZF','ND MMSE','D MRC','D ZF','D MMSE','Location','northwest')
% hold off