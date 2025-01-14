clc; close all; clear all

k= [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16];
for i = 1:length(k)
    X = sprintf('%dth PCG:', k(i));
    disp(X);

    %% PCG signal 

    myfilenameclean = sprintf('S%d_Clean.mat', k(i));
    myfilenamenoise = sprintf('S%d.mat', k(i));
    myfilenameECG = sprintf('S%d_ECG.mat', k(i));

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

    SNRinput(k(i))= snr(s_clean1(index),s_clean1(index)-s_noise(index));
    
    %CRC
    myfilenamecleapoutput0 = sprintf('S%d_Result.mat', k(i));
   % FolderName = 'C:\Users\Justin\Dropbox\Working\LettreIEEE\Participants\CRC';
     FolderName = '/users/these/phamdu/Dropbox/Working/LettreIEEE/Participants/CRC';
    File       = fullfile(FolderName, myfilenamecleapoutput0);
    load(File);
    sEst0 =  sEst; 
    
    %NMF
    myfilenamecleapoutput = sprintf('S%d_Result.mat', k(i));
    FolderName = '/users/these/phamdu/Dropbox/Working/LettreIEEE/Participants/NMF';
    %FolderName = 'C:\Users\Justin\Dropbox\Working\LettreIEEE\Participants\NMF';
    File       = fullfile(FolderName, myfilenamecleapoutput);
    load(File);
    sEst1 =  sEst; 

    % NEW
    myfilenamecleapoutput1 = sprintf('S%d_Result.mat', k(i));
    FolderName = '/users/these/phamdu/Dropbox/Working/LettreIEEE/Participants/NCCOM';
    %FolderName = 'C:\Users\Justin\Dropbox\Working\LettreIEEE\Participants\NMF';
    File       = fullfile(FolderName, myfilenamecleapoutput1);
    load(File);
   
    SNRoutputold0(k(i))= snr(s_clean1(index),s_clean1(index)-sEst0(index));
    SNR_CRC(k(i)) = SNRoutputold0(k(i))-SNRinput(k(i))
    
    SNRoutputold(k(i))= snr(s_clean1(index),s_clean1(index)-sEst1(index));
    SNR_PHAMold(k(i)) = SNRoutputold(k(i))-SNRinput(k(i))

    SNRoutputnew(k(i))= snr(s_clean1(index),s_clean1(index)-sEst(index));
    SNR_PHAMnew(k(i)) = SNRoutputnew(k(i))-SNRinput(k(i))

    Gain(k(i)) = SNR_PHAMnew(k(i))-SNR_PHAMold(k(i))
       
end
SNR_CRC( :,~any(SNR_CRC,1) ) = [];    
SNR_PHAMold( :,~any(SNR_PHAMold,1) ) = [];
SNR_PHAMnew( :,~any(SNR_PHAMnew,1) ) = [];

FigHandle(1) = figure(); 
     boxplot([SNR_PHAMold',SNR_CRC',SNR_PHAMnew'],'Labels',{'NMF','CRC','NCCOM'}) %5 IMF kept
    ylabel('SNR Gain [dB]')
