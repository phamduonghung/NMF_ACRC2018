function [index_sig, index_noise, flag_default] = select_threshold(threshold, He, Hnp, fs, number_nmf_np, index_figure)

j=0;
flag_default=0;
index_sig=[];
figure(index_figure); 
for i=1:number_nmf_np
    [intercorr, tim] = intercorr_nd(He, Hnp(i, :), fs); %[intercorr(i, :), tim(i, :)] pour afficher les valeurs d'intercorr
    %intercorr = intercorr(i, :);
    [m(i),I] = max(abs(intercorr((tim<1)&(tim>-1))));      
    if (m(i)>=threshold)
        j = j+1;
        intTau0(j) = m(i);
        index_sig(j)=i; %index for signal
    end
    %subplot(2,number_nmf_np,2*i-1);
%     figure(i)
%     subplot(1,2,1)
%     plot(tim, intercorr); %xlim([-2 2]); ylim([0 1]); legend('xcorr of He and Hnp'); grid on 
%     %subplot(2,number_nmf_np,2*i);
%     subplot(1,2,2)
%     plot(Hnp(i, :)); hold on;
%     plot(He); legend('Hn','He');
    subplot(2,number_nmf_np,2*i-1);plot(tim, intercorr); xlim([-2 2]); ylim([0 1]); grid on %plot(tim(i, :), intercorr(i, :))
    subplot(2,number_nmf_np,2*i);plot(Hnp(i, :)); hold on;
    plot(He); close 
end
if isempty(index_sig)
    flag_default=1;
    index_sig = find(max(m)==m);
end
index_noise = setdiff(1:number_nmf_np, index_sig);