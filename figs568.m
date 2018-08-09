% Simulation code for "Downlink Performance of Superimposed Pilots in Massive MIMO Systems"
% by Karthik Upadhya, Sergiy Vorobyov, and Mikko Vehkaperä
% Accepted for publication in IEEE Transactions on Wireless Communications.
%
% This code generates Figures 5-8.

clc
clear all
close all

% - Simulation Parameters - %
numTrial        = 1e4; % Number of user locations in the monte carlo simulation
ulSnrDb         = 10;  % Uplink SNR in dB
dlSnrDb         = 10;  % Downlink SNR in dB
numSymbolUlTimeSlot   = 35;  % Number of symbols in UL time slot

% - System Parameters - %
numCellTotal        = 91; % Number of cells in network
numUser             = 5;  % Number of users per cell
enPowerControl      = 1;  % Enable statistical channel inversion power control

% - Initialize placeholder variables - %
spInterferencePowerTheo=zeros(numTrial,numUser);
regInterferencePowerTheo=zeros(numTrial,numUser);
stagInterferencePowerTheo = zeros(numTrial,numUser);
regSigPowerTheo = regInterferencePowerTheo;
stagSigPowerTheo = stagInterferencePowerTheo;
spSigPowerTheo = spInterferencePowerTheo;

pilotTypeRange = {'staggered','superimposed','regular'}; % Range of pilot types to be tested
pilotReuseFactorRange = [1,3,7]; % Range of pilot reuse factors to be tested

% - Channel setup - %
channelParameter.pathLossCoeff          = 3; % Large scale path-loss coefficient
channelParameter.cellRadius             = 1; % Cell radius in km
channelParameter.d0                     = 0.1;% Users are located at least d0 metres away from the BS
%%%%%%%%%%%%%%%%%%%%%%

% - Simulation Begins - %
for pp = 1:numel(pilotTypeRange)
    pilotType = pilotTypeRange{pp}; % Set pilot type
    for rr = 1:numel(pilotReuseFactorRange)
        pilotReuseFactor = pilotReuseFactorRange(rr); % Set pilot reuse factor
        
        % - For SP and staggered pilot, calculate rate only for r^{sp} = 7 - %
        switch lower(pilotType)
            case {'staggered','superimposed'}
                if pilotReuseFactor ~= 7
                    continue
                end
        end
        
        channelParameter.sharedPilotCellIdx= genSharedPilotCellIndex(numCellTotal,pilotReuseFactor); % Compute cell indices of cells that share the same pilots
        % Here, channelParameter.sharedPilotCellIdx{jj}(ll) == 1 implies
        % that cell jj and ll share the same pilots.
        
        numBSAntennaRange = [50,100,300,600,1000,2000,10000,1e5,1e6,1e7,1e8,1e9,1e10,1e11,1e12]; % Range of number of BS antennas
        for mm = 1:length(numBSAntennaRange)
            numBSAntenna = numBSAntennaRange(mm); % Set number of BS antennas
            for ii = 1:numTrial
                rng(ii); % Uncomment to get the same results as that in the paper
                
                % - Generate user locations - %
                center  = genCellCenters(numCellTotal) .* channelParameter.cellRadius; % Generate the location of BSs in network
                msLocations         = cell(numCellTotal,1); 
                [msLocations{:}]    = deal(zeros(1,numUser));
                userTxGain          = zeros(numCellTotal,numUser);
                for ll = 1:numCellTotal
                    msLocations{ll}  = genHexSample(numUser,channelParameter.d0) * channelParameter.cellRadius + center(ll); % Generate user locations in cell ll.
                    userTxGain(ll,:) = -(channelParameter.pathLossCoeff * 10 * log10(channelParameter.d0./abs(msLocations{ll}-center(ll)))); % Compensate for large-scale path-loss
                end
                
                for ll = 1:numCellTotal % Index of Reference BS
                    for jj = 1:numCellTotal
                        rxPowerInDb = channelParameter.pathLossCoeff * 10 * log10(channelParameter.d0./abs(msLocations{jj}-center(ll))) + userTxGain(jj,:); % Received power of users in dB
                        betaVal     = 10.^(rxPowerInDb/10);
                        channelParameter.betaVal{ll}(jj,:) = betaVal; % The user transmit is subsumed into the large-scale path-loss coefficients (unlike in the paper). In the paper, both are denoted with different symbols
                    end
                end                
                
                channelParameter.ulNoiseVar         = 10^(-ulSnrDb/10); 
                channelParameter.dlNoiseVar         = 10^(-dlSnrDb/10);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                switch pilotType
                    case 'superimposed'
                        % Superimposed pilots %
                        
                        [rhoVal,lambdaVal]                       = rhoLambdaCalc(channelParameter,numSymbolUlTimeSlot,numCellTotal,numUser,numBSAntenna); % Generate values of rho and lambda according to Eqn (75)                        
                        [spSigPowerTheo(ii,:),spInterferencePowerTheo(ii,:)] = spTheoDlInterferenceCalcSpeedup(channelParameter,numBSAntenna,numSymbolUlTimeSlot,numCellTotal,numUser,rhoVal,lambdaVal); % Compute SINR as in Eqn (16)
                        
                        % [spSigPowerTheo(ii,:),spInterferencePowerTheo(ii,:)] = spTheoDlInterferenceCalc(channelParameter,numBSAntenna,numSymbolUlTimeSlot,numCellTotal,numUser,rhoVal,lambdaVal);
                        % This functions also computes the SINR, but is readable. Function titled spTheoDlInterferenceCalcSpeedup has been vectorized for speed
                        
                        spRateTheo{mm}(ii,:)              = log2( 1 + spSigPowerTheo(ii,:)./spInterferencePowerTheo(ii,:)); % Compute spectral efficiency
                        spMse{mm}(ii,:)                   = spMseCalc(channelParameter,numSymbolUlTimeSlot,numCellTotal,numUser,rhoVal,lambdaVal); % Compute MSE
                    case 'staggered'
                        % Staggered Pilots %
                        
                        [rhoVal,lambdaVal]                       = rhoLambdaCalc(channelParameter,numSymbolUlTimeSlot,numCellTotal,numUser,numBSAntenna); % Generate values of rho and lambda according to Eqn (75)
                        
                        % Compute pP, pD, and tau according to Proposition 3 %
                        tau                               = numUser;
                        pP                                = lambdaVal^2 * numSymbolUlTimeSlot/tau; % Compute pP, pD, and tau according to Proposition 3
                        pD                                = rhoVal^2;
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        
                        [stagSigPowerTheo(ii,:),stagInterferencePowerTheo(ii,:)]   = stagTheoDlInterferenceCalcSpeedup(channelParameter,numBSAntenna,tau,numCellTotal,numUser,pD,pP); % Compute SINR as in Eqn (27)
                        
                        % [stagSigPowerTheo(ii,:),stagInterferencePowerTheo(ii,:)]   = stagTheoDlInterferenceCalc(channelParameter,numBSAntenna,tau,numCellTotal,numUser,pD,pP);
                        % This functions also computes the SINR, but is readable. Function titled stagTheoDlInterferenceCalcSpeedup has been vectorized for speed
                        
                        stagRateTheo{mm}(ii,:)            = log2( 1 + stagSigPowerTheo(ii,:)./stagInterferencePowerTheo(ii,:)); % Compute spectral efficiency
                        stagMse{mm}(ii,:)                 = stagMseCalc(channelParameter,tau,numCellTotal,numUser,pD,pP); % Compute MSE
                    case 'regular'
                        % Regular Pilots %
                        
                        [regSigPowerTheo(ii,:),regInterferencePowerTheo(ii,:)]     = regTheoDlInterferenceCalcSpeedup(channelParameter,numBSAntenna,numCellTotal,numUser,numUser * pilotReuseFactor); % Compute SINR as in Eqn (6)
                        
                        % [regSigPowerTheo(ii,:),regInterferencePowerTheo(ii,:)]     = regTheoDlInterferenceCalc(channelParameter,numBSAntenna,numCellTotal,numUser,numUser * pilotReuseFactor);
                        % This functions also computes the SINR, but is readable. Function titled regTheoDlInterferenceCalcSpeedup has been vectorized for speed
                        
                        regRateTheo{mm}(ii,:)              = log2( 1 + regSigPowerTheo(ii,:)./regInterferencePowerTheo(ii,:)); % Compute spectral efficiency
                        regMse{mm}(ii,:)                   = regMseCalc(channelParameter,numCellTotal,numUser,numUser * pilotReuseFactor); % Compute MSE
                end                
            end
            crlb(mm) = mean(1./(numSymbolUlTimeSlot / channelParameter.ulNoiseVar + channelParameter.betaVal{1}(1,:) ));
            switch pilotType
                case 'superimposed'
                    fValRate(mm) = sum(mean(spRateTheo{mm},1)); % Average throughput across user locations
                    fValMse(mm)  = mean(mean(spMse{mm},1)); % Average MSE across user locations
                    
                    if numBSAntenna == 100 
                        [y,n] = hist(spRateTheo{2}(:),1e3); % Generate CDF of throughput for M = 100
                        figure(3); plot(n,cumsum(y) / sum(y));
                        hold all;
                    end
                case 'staggered'
                    % c.f. comments for 'superimposed'. The code is the same.
                    fValRate(mm) = sum(mean(stagRateTheo{mm},1));
                    fValMse(mm)  = mean(mean(stagMse{mm},1));
                    if numBSAntenna == 100
                        [y,n] = hist(stagRateTheo{2}(:),1e3);
                        figure(3); plot(n,cumsum(y) / sum(y));
                        hold all
                    end
                case 'regular'
                    % c.f. comments for 'superimposed'. The code is the same.
                    fValRate(mm) = sum(mean(regRateTheo{mm},1));
                    fValMse(mm)  = mean(mean(regMse{mm},1));
                    if numBSAntenna == 100
                        [y,n] = hist(regRateTheo{2}(:),1e3);
                        figure(3); plot(n,cumsum(y) / sum(y));
                        hold all
                    end
            end
        end
        figure(1);semilogx(numBSAntennaRange,fValRate);
        hold all;
        figure(2);semilogx(numBSAntennaRange,10*log10(fValMse));
        hold all;
        drawnow;
    end
end
figure(1); legend('Staggered Pilot : $r^{SP} = 7$', 'SP : $r^{SP} = 7$', 'RP : $r^{RP} = 1', 'RP : $r^{RP} = 3', 'RP : $r^{RP} = 7');
figure(2); plot(numBSAntennaRange,10*log10(crlb));
figure(2); legend('Staggered Pilot : $r^{SP} = 7$', 'SP : $r^{SP} = 7$', 'RP : $r^{RP} = 1', 'RP : $r^{RP} = 3', 'RP : $r^{RP} = 7','CRLB');
figure(3); legend('Staggered Pilot : $r^{SP} = 7$', 'SP : $r^{SP} = 7$', 'RP : $r^{RP} = 1', 'RP : $r^{RP} = 3', 'RP : $r^{RP} = 7');