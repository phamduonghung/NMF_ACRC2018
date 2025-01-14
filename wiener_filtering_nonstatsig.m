function [y,Hwien,xf,Hf] = wiener_filtering_nonstatsig(x,psds,psdn, H, L, nfft)

% wiener_filtering_nonstatsig - perform Wiener filtering
%   x is the noisy signal
%   psds is the power spectral density (psd) of the original signal
%   psdn is the psd of the noise
%   w : window length
%   h : hop size
%   nfft : number of points for discrete fourier transform
%   fs : sample frequency

% FFT-based wiener filtering
Hf = psds./(psds + psdn); %filter transfsr function

% compute convolution
%old
% [xf, f, t_stft] = stft(x, w, h, nfft, fs); %stft : short-time fourier transform
% [y, t_istft] = istft(xf.*Hf, h, nfft, fs); %istft : inverse short-time fourier transform
%new
N1 = length(x);

%        t=200+(-128:127); sig=[fmconst(200,0.2);fmconst(200,0.4)]; 
%        h=hamming(57); tfr=tfrstft(sig,t,256,h,1);
%        sigsyn=tfristft(tfr,t,h,1);
%        plot(t,abs(sigsyn-sig(t)))

[xf,t_stft,f] = tfrstft(x,1:N1,nfft,H,0);
[y,t_istft] = tfristft(xf.*Hf,1:N1,H,0); y =real(y);

%[xf,~,~] = get_tfrs(x,nfft,H,L); % get TF representations
%[y] = get_resynth(xf.*Hf,nfft,H);

%[xf, f, t_stft] = stft_new(x, w, h, nfft, fs);
%[y, t_istft] = istft_new(xf.*Hf, w, h, nfft, fs);

% % plot the original signal
% figure()
% plot(t_stft(1:length(x)), x, 'b')
% xlim([0 max(t_stft)])
% set(gca, 'FontName', 'Times New Roman', 'FontSize', 14)
% xlabel('Time, s')
% ylabel('Signal amplitude')
% title('Original and reconstructed signal')
% % plot the resynthesized signal 
% hold on
% plot(t_istft, y, '-.r')
% legend('Original signal', 'Reconstructed signal')
% figure
% plot((x(200:length(x)-200))-y(200:length(x)-200)); legend('error')
% %pause
if nargout>1
    Hwien = fftshift( ifft(Hf) ); %wiener filter transfsr function
end