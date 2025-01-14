function [ye, index_sig, flag_default,Hf] = methodm2 (threshold, number_nmf_np, Wnp, Hnp, He, fs, xf, H,nfft)

%% Intercorr for the choice of index_sig

[index_sig, index_noise, flag_default] = select_threshold(threshold, He, Hnp, fs, number_nmf_np, 1);

%% Wiener filtering after nnmf
%Compute the power spectral density of signal and noise
psds = Wnp(:, index_sig)*Hnp(index_sig, :);
psdn = Wnp(:, index_noise)*Hnp(index_noise, :);

psds = max(psds, 1e-6*max(psds(:)));
psdn = max(psdn, 1e-6*max(psdn(:)));

%Applying Wiener filtering based on psd
Hf = psds./(psds + psdn); %filter transfsr function
N1 = size(xf,2);

%[xf,t_stft,f_stft] = tfrstft(x,1:N1,nfft,H,0);

%[ye,t_istft] = tfristft(xf.*Hf,1:N1,H,0); ye =real(ye); % it corresponds
%to function stft
[ye] = get_resynth(xf.*Hf,nfft,H); 
% figure()
% plot(ye); hold on; plot(ye1,'.-')
%[ye,Hwiene,xf,Hf] = wiener_filtering_nonstatsig(x,psds,psdn, H, L, nfft);
%ye = ye(1:length(x)); 
