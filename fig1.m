%% This script draws figure 1 of paper [1] (see paper for more details)
% [1] 'Cardiac Signal Denoising via an Adaptive Combination of Non-negative Matrix Factorization
% and Contour Representation Computation'. Duong-Hung Pham, Sylvain Meignen,
% Nafissa Dia, Julie Fontecave-Jallon, and Bertrand Rivet.
% Duong Hung PHAM 8/2/2018


clear all;
close all;
%clc;
% addpath(genpath('/users/these/phamdu/Dropbox/Working/LettreIEEE'));
% addpath(genpath('/users/these/phamdu/Dropbox/Working/PCG_Denoisng_Nafissa_test'))
% addpath(genpath('/media/phamdu/USB/PCG_data/SISECFull'))

chemin0 = '~/Dropbox/Papers_PHAM_MEIGNEN/EUSIPCO2018_PCG/figures';
set(0,'DefaultAxesFontSize',12);

for k= 3

threshold = 0.45;
zz = 0; 
X = sprintf('%dth PCG:', k);
disp(X)

%% PCG signal 

myfilenameclean = sprintf('S%d_Clean.mat', k);
myfilenamenoise = sprintf('S%d.mat', k);
myfilenameECG = sprintf('S%d_ECG.mat', k);

sclean = importdata(myfilenameclean);
snoise = importdata(myfilenamenoise);
sECG = importdata(myfilenameECG);

s_clean = sclean.PCG;
s_noise = snoise.x;
s_ECG = sECG;
fs = snoise.fs;

N1=length(s_noise);
index = 200:N1-200;
s_clean = s_clean(1:N1); s_clean1 = s_clean;
s_noise =s_noise(1:N1);
s_ECG = s_ECG(1:N1);

SNRinput(k)= snr(s_clean1(index),s_clean1(index)-s_noise(index))

%% Signal and parameters
nbTest = 1;
%t = (0:N1-1)/fs; 
%% BPF ECG
bECG=fir2(200,[0 5.9 6 33 33.1 fs/2 ]/(fs/2),[0 0 1 1 0 0]);
s_ECG=filtfilt(bECG,1,s_ECG);

% HPF PCG
bPCG=fir2(100,[0 20-eps 20 fs/2 ]/(fs/2),[0 0 1 1]);
s_clean=filtfilt(bPCG,1,s_clean);

%% Number of components for the nnmf
number_nmf_np = 12;
algo = 'als';
number_nmf_p = 2;
number_nmf_e = 1;


%% Window
a=0.05;
nfft=512; %Number of FFT
prec = 10^(-6);
L  = sqrt(nfft/a);
l  = floor(sqrt(-nfft*log(prec)/(a*pi)))+1;
w  = amgauss(2*l+1,l+1,L);
H  = w/norm(w);

t = (0:N1-1)/fs; 
f = (0:nfft/2-1);
%%
make_it_tight = true;
subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.01], [0.2 0.1], [0.15 0.05]);
if ~make_it_tight,  clear subplot;  end
FigHandle(1) = figure(); subplot(3,1,1)
plot(t(index),s_clean1(index),'k');  ylabel('s(t)');
ylim([-1.01*max(abs(s_clean1)) 1.01*max(abs(s_clean1))]);
xlim([0 max(t)]); set(gca,'xtick',[]);

FigHandle(2) = figure(); 
subplot(3,1,1)
plot(t(index),s_noise(index),'b');  ylabel('x(t)');
ylim([-1.01*max(abs(s_noise)) 1.01*max(abs(s_noise))]);
xlim([0 max(t)]); set(gca,'xtick',[]);

FigHandle(3) = figure(); 
subplot(3,1,1)
plot(t(index),s_ECG(index),'r'); 
ylabel('ecg(t)'); 
xlabel('time [sec]');   xlim([0 max(t)]); 

FigHandle(4) = figure();
ha1 = subplot(3,1,1);
plot(t(index),s_clean1(index),'k'); 
legend('(a)','Location','southeast');
xlim([0 max(t)]); set(ha1,'xtick',[]);

ha2 = subplot(3,1,2);
plot(t(index),s_noise(index),'b'); 
legend('(b)','Location','southeast'); 
ylim([-1.01*max(abs(s_noise)) 1.01*max(abs(s_noise))]);
 xlim([0 max(t)]); set(ha2,'xtick',[]);

ha3 = subplot(3,1,3);
plot(t(index),s_ECG(index),'r'); 
legend('(c)','Location','southeast'); 
xlabel('time [sec]');  xlim([0 max(t)]);  
%%%%%%%%%%%%%%%%%%%%% print Figures
% for i = 1:4
%  export_fig(FigHandle(i), ... % figure handle
%      sprintf('%s/fig1%d', chemin0,i),... % name of output file without extension
%      '-painters', ...      % renderer
%      '-transparent', ...   % renderer
%      '-pdf', ...           % file format
%      '-r500' );             % resolution in dpi
% end
end
 