function avgAchRate = code7( numBsAntenna , pilotType , precoderType , pilotReuseFactor )
% - Simulation Parameters - %
numUserLocation = 1e3; % Number of realizations of user locations
numTrial        = 10; % Number of channel realizations to average over for each realization of user location
ulSnrDb         = 10; % Uplink SNR
dlSnrDb         = 10; % Downlink SNR
numSymbolUlTimeSlot   = 35;
% - System Parameters - %
numCell             = 91; % Number of cells in network
numUser             = 5; % - Number of users per cell

% - Channel setup - %
channelParameter.pathLossCoeff          = 3; % Large scale path-loss coefficient
channelParameter.cellRadius             = 1; % Cell radius in km
channelParameter.d0                     = 0.1;% Users are located at least d0 metres away from the BS
%%%%%%%%%%%%%%%%%%%%%%

channelParameter.sharedPilotCellIdx= genSharedPilotCellIndex(numCell,pilotReuseFactor); % Compute cell indices of cells that share the same pilots
% Here, channelParameter.sharedPilotCellIdx{jj}(ll) == 1 implies
% that cell jj and ll share the same pilots.

% - Initialize placeholder variables - %
interfTerm1Debug        = cell(numUserLocation,numCell,numUser);
[interfTerm1Debug{:}]   = deal(zeros(numTrial,numUser));
interfTerm2Debug        = cell(numUserLocation,1);
[interfTerm2Debug{:}]   = deal(zeros(numTrial,numUser));
txPrecoderPowerDebug    = cell(numUserLocation,numCell);
[txPrecoderPowerDebug{:}] = deal(zeros(numTrial,numUser));


% - Simulation Begins - %
for nn = 1:numUserLocation
    % rng(nn); % Uncomment to reproduce the results in the paper
    
    % - Generate user locations - %
    center  = genCellCenters(numCell) .* channelParameter.cellRadius; % Generate the location of BSs in network
    msLocations         = cell(numCell,1);
    [msLocations{:}]    = deal(zeros(1,numUser));
    userTxGain          = zeros(numCell,numUser);
    for ll = 1:numCell
        msLocations{ll}  = genHexSample(numUser,channelParameter.d0) * channelParameter.cellRadius + center(ll); % Generate user locations in cell ll.
        userTxGain(ll,:) = -(channelParameter.pathLossCoeff * 10 * log10(channelParameter.d0./abs(msLocations{ll}-center(ll)))); % Compensate for large-scale path-loss
    end
    
    for ll = 1:numCell % Index of Reference BS
        for jj = 1:numCell
            rxPowerInDb = channelParameter.pathLossCoeff * 10 * log10(channelParameter.d0./abs(msLocations{jj}-center(ll))) + userTxGain(jj,:); % Received power of users in dB
            betaVal     = 10.^(rxPowerInDb/10);
            channelParameter.betaVal{ll}(jj,:) = betaVal; % The user transmit is subsumed into the large-scale path-loss coefficients (unlike in the paper). In the paper, both are denoted with different symbols
        end
    end
    
    channelParameter.ulNoiseVar         = 10^(-ulSnrDb/10);
    channelParameter.dlNoiseVar         = 10^(-dlSnrDb/10);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [rhoVal,lambdaVal]             = rhoLambdaCalc(channelParameter,numSymbolUlTimeSlot,numCell,numUser,numBsAntenna);
    for ii = 1:numTrial
        for jj = 1:numCell
            G                   = cNormrnd(0,1,numBsAntenna,numCell*numUser);
            for ll = 1:numCell
                D                                       = (channelParameter.betaVal{jj}(ll,:));
                channelParameter.channelMatrix{jj,ll}   = G(:,(ll-1)*numUser+1:ll*numUser) * sparse(diag(sqrt(D)));
            end
        end
        
        %%%%% - Uplink Stage - %%%%%
        for refCell = 1:numCell
            channelMatrix   = channelParameter.channelMatrix(refCell,:);
            switch lower(pilotType)
                case 'regular'
                    channelEstimate{refCell} = zeros(numBsAntenna,numUser);
                    for ll = 1:numCell
                        if channelParameter.sharedPilotCellIdx{refCell}(ll)
                            channelEstimate{refCell} = channelEstimate{refCell} + channelMatrix{ll};
                        end
                    end
                    channelEstimate{refCell} = channelEstimate{refCell} + cNormrnd(0,sqrt(channelParameter.ulNoiseVar/numUser),numBsAntenna,numUser);
                case 'superimposed'
                    txData          = cNormrnd(0,1,numCell*numUser,numUser);
                    channelEstimate{refCell} = zeros(numBsAntenna,numUser);
                    for ll = 1:numCell
                        channelEstimate{refCell} = channelEstimate{refCell} + rhoVal / (sqrt(numSymbolUlTimeSlot) * lambdaVal) * channelMatrix{ll} * txData((ll-1)*numUser+1:ll*numUser,:);
                        if channelParameter.sharedPilotCellIdx{refCell}(ll);
                            channelEstimate{refCell} = channelEstimate{refCell} + channelMatrix{ll};
                        end
                    end
                    channelEstimate{refCell} = channelEstimate{refCell} + cNormrnd(0,sqrt(channelParameter.ulNoiseVar/(lambdaVal^2*numSymbolUlTimeSlot)),numBsAntenna,numUser);
                    
                case 'staggered'
                    tau             = numUser;
                    pP              = lambdaVal^2 * numSymbolUlTimeSlot/tau;
                    pD              = rhoVal^2;
                    txData          = cNormrnd(0,1,numCell*numUser,numUser);
                    channelEstimate{refCell} = zeros(numBsAntenna,numUser);
                    for ll = 1:numCell
                        if channelParameter.sharedPilotCellIdx{refCell}(ll);
                            channelEstimate{refCell} = channelEstimate{refCell} + channelMatrix{ll};
                        else
                            channelEstimate{refCell} = channelEstimate{refCell} + sqrt( pD / (tau * pP) ) * channelMatrix{ll} * txData((ll-1)*numUser+1:ll*numUser,:);
                        end
                    end
                    channelEstimate{refCell} = channelEstimate{refCell} + cNormrnd(0,sqrt(channelParameter.ulNoiseVar/(pP * tau)),numBsAntenna,numUser);
            end
        end        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%% - Downlink Stage - %%%%%
        for ll = 1:numCell
            for kk = 1:numUser
                switch lower(precoderType);
                    case 'mr'
                        txPrecodingVec= channelEstimate{ll}(:,kk); % MR Precoding Vector
                    case 'zf'
                        txPrecodingMat= channelEstimate{ll} / (channelEstimate{ll}' * channelEstimate{ll}); % ZF Precoding vector
                        txPrecodingVec= txPrecodingMat(:,kk);
                end
                interfTerm1Debug{nn,ll,kk}(ii,:) = abs( channelParameter.channelMatrix{ll,1}' * txPrecodingVec ).^2; % Interference term 1 (before expectation) in the denominator of (15)
                txPrecoderPowerDebug{nn,ll}(ii,kk)   = norm(txPrecodingVec,2)^2;
                
                if ll == 1
                    tVec              = channelParameter.channelMatrix{1,1}' * txPrecodingVec;
                    interfTerm2Debug{nn}(ii,kk) = tVec(kk); % Interference term 2 (before expectation) in the denominator of (15) ( or signal power in the numerator of (15) )
                end
            end
        end
    end
end

achRate = zeros(numUserLocation,numUser);
% Compute signal powers, interference powers, and achievable rate %
for nn = 1:numUserLocation
    interfTerm1 = zeros(1,numUser);
    txPrecoderNormalizingFactor = zeros(numCell,numUser);
    for ll = 1:numCell
        for kk = 1:numUser
            txPrecoderNormalizingFactor(ll,kk) = 1./mean(txPrecoderPowerDebug{nn,ll}(:,kk)); % This is the parameter \nu in the paper
            interfTerm1 = interfTerm1 + txPrecoderNormalizingFactor(ll,kk) * mean(interfTerm1Debug{nn,ll,kk},1); % Term 1 in the denominator of (15) after expectation
        end
    end
    
    interfTerm2 = - txPrecoderNormalizingFactor(1,:) .* abs(mean(interfTerm2Debug{nn},1)).^2; % Term 2 in the denominator of (15) after expectation
    sigTerm     = - interfTerm2; % Average signal power power
    interfTerm3 = channelParameter.dlNoiseVar;
    
    sinrSim = sigTerm ./ ( interfTerm1 + interfTerm2 + interfTerm3 );
    achRate(nn,:) = log2(1 + sinrSim);
end

avgAchRate = sum(mean(achRate,1)); % Compute sum rate in the DL
