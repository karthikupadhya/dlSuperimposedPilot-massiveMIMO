function [ sigPowerTheo , interferencePowerTheo ] = stagTheoDlInterferenceCalcSpeedup(channelParameter,M,tau,L,K,pD,pP)

betaVal = channelParameter.betaVal;
dlNoiseVar= channelParameter.dlNoiseVar;
ulNoiseVar= channelParameter.ulNoiseVar(1);
sharedPilotCellIdx     = channelParameter.sharedPilotCellIdx;

invTxPowerNormalization = zeros(L,K);
for ll = 1:L
    tVal1 = M * (sharedPilotCellIdx{ll}(:) == 1)' * betaVal{ll};
    tVal2(1:K) = M * (sharedPilotCellIdx{ll}(:) == 0).' * sum(betaVal{ll},2) * pD / (tau * pP);
    invTxPowerNormalization(ll,:) = tVal1 + tVal2 + ulNoiseVar * M / (tau * pP);        
end
txPowerNormalization = 1./invTxPowerNormalization;

jj = 1;
for mm = 1:K
    interfTerm1(mm) = 0;
    for ll = 1:L
        tVal1 = txPowerNormalization(ll,mm) * M^2 * betaVal{ll}(jj,mm)^2 * (sharedPilotCellIdx{ll}(jj) == 1);
        tVal2 = sum(M * betaVal{ll}(jj,mm)  * txPowerNormalization(ll,:) .* ((sharedPilotCellIdx{ll}(:) == 1).' * betaVal{ll}));
        tVal3 = sum( txPowerNormalization(ll,:) * M * ulNoiseVar * betaVal{ll}(jj,mm) / (pP * tau));
        tVal4 = sum( txPowerNormalization(ll,:) * pD / (tau * pP) * ( M^2 * betaVal{ll}(jj,mm)^2 * (sharedPilotCellIdx{ll}(jj) == 0)  ) );
        tVal5 = sum( txPowerNormalization(ll,:) * pD / (tau * pP) * M * betaVal{ll}(jj,mm) * ( (sharedPilotCellIdx{ll}(:) == 0).' * sum( betaVal{ll} , 2 ) ) );
        interfTerm1(mm) = interfTerm1(mm) + (tVal1 + tVal2 + tVal3 + tVal4 + tVal5);
    end
    
    interfTerm2(mm) = - M^2 * betaVal{jj}(jj,mm)^2 * txPowerNormalization(jj,mm);
    
    interfTerm3(mm) = dlNoiseVar;
end


interferencePowerTheo = interfTerm1 + interfTerm2 + interfTerm3;
sigPowerTheo          = M^2 * betaVal{jj}(jj,:).^2 .* txPowerNormalization(jj,:);