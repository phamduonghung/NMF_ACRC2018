clearvars
clc

%participantList = {'1_BOUGRINE' ; '2_IHSSENE';'3_DIA';'4_PHAM_max';'5_EMD';'6_HT';'7_IHT'}; 
participantList = {'EMD';'NMF'; 'CRC'; 'NCCOM'}; 

% 1_BOUGRINE
% 2_IHSSENE


for indParticipant = 1:length(participantList)
    participant = participantList{indParticipant};
    participant
    for indSubject = [1 3:16]%16]
        indSubject
        if (indSubject == 1) || (indSubject == 2) || (indSubject == 3) || (indSubject == 4)
            % Loading reference signals
            load(['SISECFull/S' num2str(indSubject) '_Parameters'], 'noiseAmbiant', 'PCG', 'fs', 'x');
            PCGRef          = PCG(1:fix(length(PCG)/fs)*fs);
            noiseAmbiantRef = noiseAmbiant(1:fix(length(noiseAmbiant)/fs)*fs);
            clear PCG noiseAmbiant

            % Loading participant data
            load(['Participants/' participant '/S' num2str(indSubject) '_Result']);

            % Performance evaluation
            [SDR{indParticipant,indSubject},SIR{indParticipant,indSubject},SAR{indParticipant,indSubject},perm{indParticipant,indSubject}]=bss_eval_sources([sEst(:) noiseAmbiantRef(:)].', [PCGRef(:) noiseAmbiantRef(:)].');
            [SDR_Init{indParticipant,indSubject},SIR_Init{indParticipant,indSubject},SAR_Init{indParticipant,indSubject},perm_Init{indParticipant,indSubject}]=bss_eval_sources([x(:) noiseAmbiantRef(:)].', [PCGRef(:) noiseAmbiantRef(:)].');

        elseif (indSubject == 5) || (indSubject == 6) || (indSubject == 7) || (indSubject == 8)
            % Loading reference signals
            load(['SISECFull/S' num2str(indSubject) '_Parameters'], 'noiseAmbiant', 'PCG', 'fs', 'x', 'periodicNoise');
            PCGRef           = PCG(1:fix(length(PCG)/fs)*fs);
            noiseAmbiantRef  = noiseAmbiant(1:fix(length(noiseAmbiant)/fs)*fs);
            periodicNoiseRef = periodicNoise(1:fix(length(periodicNoise)/fs)*fs);
            clear PCG noiseAmbiant periodicNoise

            % Loading participant data
            load(['Participants/' participant '/S' num2str(indSubject) '_Result']);

            % Performance evaluation
            [SDR{indParticipant,indSubject},SIR{indParticipant,indSubject},SAR{indParticipant,indSubject},perm{indParticipant,indSubject}]=bss_eval_sources([sEst(:) noiseAmbiantRef(:) periodicNoiseRef(:)].', ...
                [PCGRef(:) noiseAmbiantRef(:) periodicNoiseRef(:)].');
            [SDR_Init{indParticipant,indSubject},SIR_Init{indParticipant,indSubject},SAR_Init{indParticipant,indSubject},perm_Init{indParticipant,indSubject}]=bss_eval_sources([x(:) noiseAmbiantRef(:) periodicNoiseRef(:)].', ...
                [PCGRef(:) noiseAmbiantRef(:) periodicNoiseRef(:)].');

        elseif (indSubject == 9) || (indSubject == 10) || (indSubject == 11) || (indSubject == 12)
            % Loading reference signals
            load(['SISECFull/S' num2str(indSubject) '_Parameters'], 'noiseAmbiant', 'PCG', 'fs', 'x', 'periodicNoise');
            PCGRef           = PCG(1:fix(length(PCG)/fs)*fs);
            noiseAmbiantRef  = noiseAmbiant(1:fix(length(noiseAmbiant)/fs)*fs);
            periodicNoiseRef = periodicNoise(1:fix(length(periodicNoise)/fs)*fs);
            clear PCG noiseAmbiant periodicNoise

            % Loading participant data
            load(['Participants/' participant '/S' num2str(indSubject) '_Result']);

            % Performance evaluation
            [SDR{indParticipant,indSubject},SIR{indParticipant,indSubject},SAR{indParticipant,indSubject},perm{indParticipant,indSubject}]=bss_eval_sources([sEst(:) noiseAmbiantRef(:) periodicNoiseRef(:)].', ...
                [PCGRef(:) noiseAmbiantRef(:) periodicNoiseRef(:)].');
            [SDR_Init{indParticipant,indSubject},SIR_Init{indParticipant,indSubject},SAR_Init{indParticipant,indSubject},perm_Init{indParticipant,indSubject}]=bss_eval_sources([x(:) noiseAmbiantRef(:) periodicNoiseRef(:)].', ...
                [PCGRef(:) noiseAmbiantRef(:) periodicNoiseRef(:)].');

        elseif (indSubject == 13) || (indSubject == 14) || (indSubject == 15) || (indSubject == 16)
             % Loading reference signals
            load(['SISECFull/S' num2str(indSubject) '_Parameters'], 'noiseAmbiant', 'PCG', 'fs', 'x', 'periodicNoise', 'coughNoise');
            PCGRef           = PCG(1:fix(length(PCG)/fs)*fs);
            noiseAmbiantRef  = noiseAmbiant(1:fix(length(noiseAmbiant)/fs)*fs);
            periodicNoiseRef = periodicNoise(1:fix(length(periodicNoise)/fs)*fs);
            coughNoiseRef    = coughNoise(1:fix(length(coughNoise)/fs)*fs);
            clear PCG noiseAmbiant periodicNoise

            % Loading participant data
            load(['Participants/' participant '/S' num2str(indSubject) '_Result']);

            % Performance evaluation
            [SDR{indParticipant,indSubject},SIR{indParticipant,indSubject},SAR{indParticipant,indSubject},perm{indParticipant,indSubject}]=bss_eval_sources([sEst(:) noiseAmbiantRef(:) periodicNoiseRef(:) coughNoiseRef(:)].', ...
                [PCGRef(:) noiseAmbiantRef(:) periodicNoiseRef(:) coughNoiseRef(:)].');
            [SDR_Init{indParticipant,indSubject},SIR_Init{indParticipant,indSubject},SAR_Init{indParticipant,indSubject},perm_Init{indParticipant,indSubject}]=bss_eval_sources([x(:) noiseAmbiantRef(:) periodicNoiseRef(:) coughNoiseRef(:)].', ...
                [PCGRef(:) noiseAmbiantRef(:) periodicNoiseRef(:) coughNoiseRef(:)].');
        end

%         time = (0:length(PCGRef)-1) / fs;
%         figure(1); clf
%         axs(1) = subplot(3,1,1);
%             plot(time, x)
%             grid on
% 
%         axs(2) = subplot(3,1,2);
%             plot(time, PCGRef)
%             grid on
% 
%         axs(3) = subplot(3,1,3);
%             plot(time, sEst)
% 
%             linkaxes(axs, 'x')
%             drawnow

    end

end


    