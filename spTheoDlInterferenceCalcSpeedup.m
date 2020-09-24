function [ sigPowerTheo , interferencePowerTheo ] = spTheoDlInterferenceCalcSpeedup(channelParameter,M,C,L,K,rhoD,rhoP)

betaVal = channelParameter.betaVal;
dlNoiseVar= channelParameter.dlNoiseVar;
ulNoiseVar= channelParameter.ulNoiseVar(1);
sharedPilotCellIdx     = channelParameter.sharedPilotCellIdx;

invTxPowerNormalization = zeros(L,K);
for ll = 1:L
    tVal1       = M * (sharedPilotCellIdx{ll}(:) == 1).' *  betaVal{ll};
    tVal2(1:K)  = M * (rhoD^2 / (C*rhoP^2)) * sum(betaVal{ll}(:));
    tVal3(1:K)  = M * ulNoiseVar / (C * rhoP^2);
    invTxPowerNormalization(ll,:) = tVal1 + tVal2 + tVal3;
end
txPowerNormalization = 1./invTxPowerNormalization;


jj = 1;
for mm = 1:K    
    interfTerm1(mm) = 0;
    for ll = 1:L
        tVal1       = txPowerNormalization(ll,mm) * M^2 * betaVal{ll}(jj,mm)^2 * (sharedPilotCellIdx{ll}(jj) == 1);
        tVal2       = sum( M * betaVal{ll}(jj,mm) * txPowerNormalization(ll,:) .* ((sharedPilotCellIdx{ll}(:) == 1).' * betaVal{ll}));
        tVal3       = sum(txPowerNormalization(ll,:)  * M^2 * rhoD^2 / (C * rhoP^2) * betaVal{ll}(jj,mm)^2);
        tVal4       = sum(txPowerNormalization(ll,:)  * M * rhoD^2 / ( C * rhoP^2 ) * betaVal{ll}(jj,mm) * sum(betaVal{ll}(:)));
        tVal5       = sum(txPowerNormalization(ll,:)  * M * ulNoiseVar / (C * rhoP^2) * betaVal{ll}(jj,mm) );
        interfTerm1(mm) = interfTerm1(mm) + ( tVal1 + tVal2 + tVal3 + tVal4 + tVal5 );
    end
    
    interfTerm2(mm) = - M^2 * txPowerNormalization(jj,mm) * betaVal{jj}(jj,mm)^2;
    
    interfTerm3(mm) = dlNoiseVar;
end

interferencePowerTheo = interfTerm1 + interfTerm2 + interfTerm3;
sigPowerTheo          = M^2 * txPowerNormalization(jj,:) .* betaVal{jj}(jj,:).^2 ;