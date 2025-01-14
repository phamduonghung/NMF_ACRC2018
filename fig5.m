%% This script draws figure 5 of paper [1] (see paper for more details)
% [1] 'Cardiac Signal Denoising via an Adaptive Combination of Non-negative Matrix Factorization
% and Contour Representation Computation'. Duong-Hung Pham, Sylvain Meignen,
% Nafissa Dia, Julie Fontecave-Jallon, and Bertrand Rivet.
% Duong Hung PHAM 8/2/2018

clear all;
close all;
%clc; 

addpath(genpath('/users/these/phamdu/Dropbox/Working/LettreIEEE'));
addpath(genpath('/users/these/phamdu/Dropbox/Working/PCG_Denoisng_Nafissa_test'))
addpath(genpath('/media/phamdu/USB/PCG_data/SISECFull'))

chemin0 = '~/Dropbox/Papers_PHAM_MEIGNEN/EUSIPCO2018_PCG/figures';


set(0,'DefaultAxesFontSize',18);
threshold = 0:0.1:1;
numTest = 5;

SNR_NMF = zeros(numTest,length(threshold));
SNR_NCCOM = zeros(numTest,length(threshold));

for k= 3
    for indTest = 1:numTest      
        i
        for p=1:length(threshold)
            p
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

            N1=8192;%length(s_noise);
            index = 200:N1-200;
            s_clean = s_clean(1:N1); s_clean1 = s_clean;
            s_noise =s_noise(1:N1);
            s_ECG = s_ECG(1:N1);

            SNRinput= snr(s_clean1(index),s_clean1(index)-s_noise(index))

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
            L  = sqrt(nfft/a);
            l  = floor(sqrt(-nfft*log(prec)/(a*pi)))+1;
            w  = amgauss(2*l+1,l+1,L);
            H  = w/norm(w);

            %% spectrogram 
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

            [sEst1, ~, ~,~] = methodm2(threshold(p), number_nmf_np, Wnp, Hnp, He, fs, snp, H,nfft);
            [~, index_sig_m2, flag_default_current,Hf] = methodm2(threshold(p), number_nmf_np, Wnp, Hnp, He, fs, snp, H,nfft);


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
            SNRoutputold= snr(s_clean1(index),s_clean1(index)-sEst1(index));
            SNR_PHAMold = SNRoutputold-SNRinput;

            SNRoutputnew= snr(s_clean1(index),s_clean1(index)-sEst(index));
            SNR_PHAMnew = SNRoutputnew-SNRinput;

            Gain = SNR_PHAMnew-SNR_PHAMold;

            SNR_NMF(indTest,p) = SNR_PHAMold;
            SNR_NCCOM(indTest, p) = SNR_PHAMnew;


        %     %% Save data  
        %     nameFile = sprintf('S%d_Result.mat',k);
        %     matfile = fullfile(chemin0, nameFile);
        %     save(matfile, 'sEst');
        % 
        %     % nameFile = sprintf('S%d_PHAMnew.mat',k);
        %     % matfile = fullfile(chemin0, nameFile);
        %     % save(matfile, 'sEst');
        % 
        %     sEst=sEst1;
        %     nameFile = sprintf('S%d_Result.mat',k);
        %     matfile = fullfile(chemin1, nameFile);
        %     save(matfile, 'sEst');
        %     clearvars reassign
        end
    end
        SNR_NMF1 = mean(SNR_NMF);
        SNR_NCCOM1 = mean(SNR_NCCOM);
        
        FigHandle(1) = figure(); 
        plot(threshold,SNR_NMF1,'b-+',threshold,SNR_NCCOM1,'r-*','Linewidth',1.5);  
        legend('NMF','NCCOM','location','southeast'); 
        xlabel('\lambda'); 
        ylabel('SDR Gain [dB]'); 
        ylim([0 12])
        set(gca,'XTick',0:0.2:1)

        export_fig(FigHandle(1), ... % figure handle
            sprintf('%s/Sensibility_lamda', chemin0),... % name of output file without extension
         '-painters', ...      % renderer
         '-transparent', ...   % renderer
         '-pdf', ...           % file format
         '-r500' );             % resolution in dpi
end
 