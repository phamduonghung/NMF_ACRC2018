%% This script draws figure 6 of paper [1] (see paper for more details)
% [1] 'Cardiac Signal Denoising via an Adaptive Combination of Non-negative Matrix Factorization
% and Contour Representation Computation'. Duong-Hung Pham, Sylvain Meignen,
% Nafissa Dia, Julie Fontecave-Jallon, and Bertrand Rivet.
% Duong Hung PHAM 8/2/2018

main_Evaluations_PHAM;

close all; % clc; 
%mkdir('~/Dropbox/Working/LettreIEEE/Finalresults/MainEva_new1');
%addpath('/users/these/phamdu/Dropbox/Working/LettreIEEE/Finalresults/MainEva_new1');
chemin0 = '~/Dropbox/Papers_PHAM_MEIGNEN/EUSIPCO2018_PCG/figures';
set(0,'DefaultAxesFontSize',18);

SDRgain = zeros(length(participantList),15);
SIRgain = zeros(length(participantList),15);
SARgain = zeros(length(participantList),15);

for j=1:length(participantList)
    for k=[1,3:16]
        if k==1
            SDRgain(j,k) = SDR{j,k}(1)-SDR_Init{j,k}(1); 
            SIRgain(j,k) = SIR{j,k}(1)-SIR_Init{j,k}(1); 
            SARgain(j,k) = SAR{j,k}(1); 
        else
            SDRgain(j,k-1) = SDR{j,k}(1)-SDR_Init{j,k}(1); 
            SIRgain(j,k-1) = SIR{j,k}(1)-SIR_Init{j,k}(1); 
            SARgain(j,k-1) = SAR{j,k}(1); 
        end
    end
end
SDRgain = SDRgain';
SIRgain = SIRgain';
SARgain = SARgain';

FigHandle(1) = figure();  
boxplot(SDRgain,participantList); ylabel('SDR Gain [dB]');
m=1;
if m==1
    for i = 1
     export_fig(FigHandle(i), ... % figure handle
         sprintf('%s/SDR_Gain', chemin0),... % name of output file without extension
         '-painters', ...      % renderer
         '-transparent', ...   % renderer
         '-pdf', ...           % file format
         '-r500' );             % resolution in dpi
    end 
end
FigHandle(2) = figure(); 
boxplot(SIRgain,participantList); ylabel('SIR Gain [dB]');

if m==1
for i = 2
 export_fig(FigHandle(i), ... % figure handle
     sprintf('%s/SIR_Gain', chemin0),... % name of output file without extension
     '-painters', ...      % renderer
     '-transparent', ...   % renderer
     '-pdf', ...           % file format
     '-r500' );             % resolution in dpi
end 
end
FigHandle(3) = figure();
boxplot(SARgain,participantList); ylabel('SAR [dB]');

if m==1
for i = 3
 export_fig(FigHandle(i), ... % figure handle
     sprintf('%s/SAR', chemin0),... % name of output file without extension
     '-painters', ...      % renderer
     '-transparent', ...   % renderer
     '-pdf', ...           % file format
     '-r500' );             % resolution in dpi
end 
end
