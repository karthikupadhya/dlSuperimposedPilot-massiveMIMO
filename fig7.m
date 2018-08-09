clc
clear all
close all

numBSAntennaRange = 150:150:600; % Range of M to test
pilotReuseFactorRange = [1,3,7]; % Range of pilot reuse factors to test
precoderTypeRange = {'zf','mr'}; % Compute throughputs for MR and ZF
pilotTypeRange    = {'regular','staggered','superimposed'}; % Test for regular, SP, and staggered pilot


for gg = 1:numel(precoderTypeRange)
    precoderType = precoderTypeRange{gg};    
    for pp = 1:numel(pilotTypeRange)
        pilotType = pilotTypeRange{pp};
        for rr = 1:numel(pilotReuseFactorRange)
            pilotReuseFactor = pilotReuseFactorRange(rr);
            
            switch pilotType
                case {'staggered','superimposed'}
                    if pilotReuseFactor ~= 7
                        continue; % Superimposed and staggered pilots are tested only for r^{SP} = 7
                    end
            end
            for mm = 1:length(numBSAntennaRange)
                numBSAntenna = numBSAntennaRange(mm);
                achRate(mm) = code7(numBSAntenna,pilotType,precoderType,pilotReuseFactor); % Compute achievable rate
            end
            
            %%% - Plotting - %%%
            switch strcat(pilotType,num2str(pilotReuseFactor));
                case 'regular1'
                    colorStr = 'r';
                case 'regular3'
                    colorStr = 'c';
                case 'regular7'
                    colorStr = 'm';
                case 'staggered7'
                    colorStr = 'b';
                case 'superimposed7'
                    colorStr = 'g';
            end
            
            switch precoderType
                case 'mr'
                    lineStyleStr = '--';
                case 'zf'
                    lineStyleStr = '-';
            end
            
            figure(1); plot(numBSAntennaRange,achRate,'Color',colorStr,'LineStyle',lineStyleStr);
            hold all
            drawnow
            %%%%%%%%%%%%%%%%%%%%%
        end
    end
end

figure(1); legend('RP : $r^{RP} = 1$','RP : $r^{RP} = 3$','RP : $r^{RP} = 7$','Staggered : $r^{SP} = 7$','SP : $r^{SP} = 7$');