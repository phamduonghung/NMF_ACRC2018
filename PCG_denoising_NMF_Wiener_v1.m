% Combining wiennering fuction  

clear all;
close all;
%clc;
addpath(genpath('/users/these/phamdu/Dropbox/Working/LettreIEEE'));
addpath(genpath('/users/these/phamdu/Dropbox/Working/PCG_Denoisng_Nafissa_test'))
addpath(genpath('/media/phamdu/USB/PCG_data/SISECFull'))

mkdir('~/Dropbox/Working/LettreIEEE/Participants/3_NEW');
chemin0 = '~/Dropbox/Working/LettreIEEE/Participants/3_NEW';
mkdir('~/Dropbox/Working/LettreIEEE/Participants/1_NMF');
chemin1 = '~/Dropbox/Working/LettreIEEE/Participants/1_NMF';

set(0,'DefaultAxesFontSize',18);

for k= 4

threshold = 0.45
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

%Normalize length
% s_clean = s_clean(1:length(s_noise)); 
% s_ECG = s_ECG(1:length(s_noise));

%% Signal and parameters
nbTest = 1;
%t = (0:N1-1)/fs; 
%% BPF ECG
bECG=fir2(200,[0 5.9 6 33 33.1 fs/2 ]/(fs/2),[0 0 1 1 0 0]);
s_ECG=filtfilt(bECG,1,s_ECG);

%% HPF PCG

bPCG=fir2(100,[0 20-eps 20 fs/2 ]/(fs/2),[0 0 1 1]);
s_clean=filtfilt(bPCG,1,s_clean);

%% Number of components for the nnmf
number_nmf_np = 12;
algo = 'als';
number_nmf_p = 2;
number_nmf_e = 1;

temps=(0:length(s_noise)-1)/fs;

%% Window
a=0.05;
nfft=512; %Number of FFT
prec = 10^(-6);
%[H L] = roundgauss(Nx, prec); 
L  = sqrt(nfft/a);
l  = floor(sqrt(-nfft*log(prec)/(a*pi)))+1;
w  = amgauss(2*l+1,l+1,L);
H  = w/norm(w);

%% spectrogram 
%noising PCG
%[snp,fnp,tnp,pnp] = spectrogram(s_noise, hamming(w), w-1, nfft,fs);
%Clean PCG
%[sp,fp,tp,pp] = spectrogram(s_clean, hamming(w), w-1, nfft,fs); 
%Synchronous ECG
%[se,fe,te,pe] = spectrogram(s_ECG, hamming(w), w-1, nfft,fs); 

%New computation 
%noising PCG
%[snp,tnp,fnp] = tfrstft(s_noise,1:N1,nfft,H,0);
[snp,reassign_np,dx_np] = get_tfrs(s_noise,nfft,H,L);

%Clean PCG
%[sp,tp,fp] = tfrstft(s_clean,1:N1,nfft,H,0); sp=sp(1:nfft/2,:);
[sp,reassign_p,dx_p] = get_tfrs(s_clean,nfft,H,L);

%Synchronous ECG
%[se,te,fe] = tfrstft(s_ECG,1:N1,nfft,H,0); se=se(1:nfft/2,:);
[se,reassign_e,dx_e] = get_tfrs(s_ECG,nfft,H,L);

%max(max(s_e-se(1:nfft/2,:)))
%pause

%% nnmf
%noising PCP
[Wnp, Hnp] = nnmf(abs(snp).^2, number_nmf_np,  'algorithm', algo);
%Clean PCG
[Wp, Hp] = nnmf(abs(sp).^2, number_nmf_p,  'algorithm', algo);
%Synchronous ECG
[We, He] = nnmf(abs(se).^2, number_nmf_e, 'algorithm', algo);

%for indTest = 1:nbTest
   % fprintf(['indextest= ' num2str(indTest) '\n'])
%% Method Dia 
  
[sEst1, ~, ~,~] = methodm2(0.75, number_nmf_np, Wnp, Hnp, He, fs, snp, H,nfft);
[~, index_sig_m2, flag_default_current,Hf] = methodm2(threshold, number_nmf_np, Wnp, Hnp, He, fs, snp, H,nfft);


%% Improvement
q=snp.*Hf; 
dx=dx_np.*Hf;
reassign=reassign_np;
%reassign(1,:,:) = squeeze(reassign_np(1,:,:)).*Hf; reassign(2,:,:) = squeeze(reassign_np(2,:,:)).*Hf;  

q1=q(1:80,:);
dx1=dx(1:80,:);
reassign1=reassign(:,1:80,:); 

t = (0:N1-1)/fs; 
f = (0:nfft/2-1);
f1 = (0:79);
if zz==1
    % Display spectrogram
    FigHandle(1) =figure; set(FigHandle(1),'units','normalized','outerposition',[0 0 1 1]);        
    imagesc(t,f,abs(sp)); title('clean STFT');
    set(gca,'YDir','normal');xlabel('time');ylabel('frequency'); 

    FigHandle(2) =figure; set(FigHandle(2),'units','normalized','outerposition',[0 0 1 1]);
    imagesc(t,f,abs(snp)); title('initial noisy STFT');
    set(gca,'YDir','normal');xlabel('time');ylabel('frequency'); 

    FigHandle(22) =figure; set(FigHandle(22),'units','normalized','outerposition',[0 0 1 1]);
    imagesc(t,f,abs(q)); title('noisy STFT');
    set(gca,'YDir','normal');xlabel('time');ylabel('frequency'); 

    FigHandle(3) =figure; set(FigHandle(3),'units','normalized','outerposition',[0 0 1 1]);
    imagesc(t,f1,abs(q1)); title('Low-passed STFT');
    set(gca,'YDir','normal');xlabel('time');ylabel('frequency');
end

%% Compute the basins of attraction associated to the different ridges/modes
[contours, rzeros, ~, ~,basins, ~] = get_contour_basins(q1,dx1,reassign1,1);

%% Display ridges, zeros and basins
NumC = round(3.5*N1/1024); % Number of modes
if zz==1
    FigHandle(4) = figure();  set(FigHandle(4),'units','normalized','outerposition',[0 0 1 1]);
    imagesc(t,f1,basins); set(gca,'clim',[0 NumC],'ydir','normal');
    hold on
    for i=1:NumC
        plot(t(contours{i}.x),f1(contours{i}.y),'k-','linewidth',2)
    end
    set(gca,'YDir','normal');
    tmp = find(~rzeros); [tmp1 tmp2] = ind2sub(size(rzeros),tmp); plot(t(tmp2),f(tmp1),'r.','MarkerSize',15);
    xlabel('time');ylabel('frequency'); %ylim([0 200]); set(gca,'XTick',0:0.2*max(t):max(t));xlim([0 max(t)]); xtickformat('%.2f')

    FigHandle(5) = figure();  set(FigHandle(5),'units','normalized','outerposition',[0 0 1 1]);
    imagesc(t,f1,basins); set(gca,'clim',[0 NumC],'ydir','normal');
    hold on
    for i=1:NumC
        plot(t(contours{i}.x),f1(contours{i}.y),'r-','linewidth',2)
    end
    set(gca,'YDir','normal');
end
%% Denoised STFT
basins1 = NaN(size(q));
basins1(1:80,:) = basins;

q_denoised = zeros(size(q));
for i=1:NumC
    q_denoised = q_denoised+ q.*(basins1==i);    
end
if zz==1
    FigHandle(6) =figure; set(FigHandle(6),'units','normalized','outerposition',[0 0 1 1]);
    imagesc(t,f,abs(q_denoised)); title('denoised STFT');
    set(gca,'YDir','normal');xlabel('time');ylabel('frequency');
end
%% Reconstruction
srec = get_resynth(q_denoised,nfft,H);
sEst = srec;

%% plot the original signal
% figure()
% plot(t(1:length(s_noise)), s_noise, 'b')
% xlim([0 max(t)])
% set(gca, 'FontName', 'Times New Roman', 'FontSize', 14)
% xlabel('Time, s')
% ylabel('Signal amplitude')
% title('Original and reconstructed signal')
% % plot the resynthesized signal 
% hold on
% plot(t, sEst, '-.r')
% legend('Original signal', 'Reconstructed signal')
% figure()
% plot(s_noise(index)-sEst(index)); legend('error')
% pause

%%
if zz==1
    FigHandle(7) = figure(); set(FigHandle(7),'units','normalized','outerposition',[0 0 1 1]);
    subplot(2,2,1)
    plot(t(index),s_noise(index),'b--'); 
    legend('noise','Location','southeast'); xlabel('time');
    subplot(2,2,2)
    plot(t(index),sEst(index),'--',t(index),s_clean1(index)); 
    legend('total reconstruction','clean','Location','southeast'); 
    subplot(2,2,3)
    plot(t(index),sEst(index)); 
    legend('total reconstruction','Location','southeast');        
    subplot(2,2,4)
    plot(t(index),s_clean1(index)-sEst(index),'r--'); 
    legend('error','Location','southeast'); xlabel('time'); 
end
%close all
 
% Computes output SNRs 
SNRoutputold(k)= snr(s_clean1(index),s_clean1(index)-sEst1(index))
SNR_PHAMold(k) = SNRoutputold(k)-SNRinput(k)

SNRoutputnew(k)= snr(s_clean1(index),s_clean1(index)-sEst(index))
SNR_PHAMnew(k) = SNRoutputnew(k)-SNRinput(k)

Gain(k) = SNR_PHAMnew(k)-SNR_PHAMold(k)
pause 

%% Save data  
nameFile = sprintf('S%d_Result.mat',k);
matfile = fullfile(chemin0, nameFile);
save(matfile, 'sEst');

% nameFile = sprintf('S%d_PHAMnew.mat',k);
% matfile = fullfile(chemin0, nameFile);
% save(matfile, 'sEst');

sEst=sEst1;
nameFile = sprintf('S%d_Result.mat',k);
matfile = fullfile(chemin1, nameFile);
save(matfile, 'sEst');
clearvars reassign
%end
end
save(chemin0,'SNR_PHAMnew')
 