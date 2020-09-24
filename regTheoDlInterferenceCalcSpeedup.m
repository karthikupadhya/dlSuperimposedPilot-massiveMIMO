function [ sigPowerTheo , interferencePowerTheo ] = regTheoDlInterferenceCalcSpeedup(channelParameter,M,L,K,tau)

betaVal = channelParameter.betaVal;
dlNoiseVar= channelParameter.dlNoiseVar;
ulNoiseVar= channelParameter.ulNoiseVar(1);
sharedPilotCellIdx     = channelParameter.sharedPilotCellIdx;

invTxPowerNormalization = zeros(L,K);
for ll = 1:L
    invTxPowerNormalization(ll,:) = M * (sharedPilotCellIdx{ll}(:) == 1).' * betaVal{ll};
end
invTxPowerNormalization = invTxPowerNormalization + ulNoiseVar * M / tau;
txPowerNormalization    = 1./invTxPowerNormalization;

jj = 1;
interfTerm1 = zeros(1,K);
for mm = 1:K
    interfTerm1(mm) = 0;
    for ll = 1:L
        tVal1           = txPowerNormalization(ll,mm) * M^2 * betaVal{ll}(jj,mm)^2 * (sharedPilotCellIdx{ll}(jj)==1);
        tVal2           = sum( M * betaVal{ll}(jj,mm) * txPowerNormalization(ll,:) .* ((sharedPilotCellIdx{ll}(:) == 1).' * betaVal{ll} ) );
        tVal3           = sum(txPowerNormalization(ll,:) * M * ulNoiseVar * betaVal{ll}(jj,mm) / tau) ;
        interfTerm1(mm) = interfTerm1(mm) + tVal1 + tVal2 + tVal3;
    end
end

interfTerm2           = - M^2 * betaVal{jj}(jj,:).^2 .* txPowerNormalization(jj,:);
interfTerm3(1:K)      = dlNoiseVar;    

sigPowerTheo          = M^2 * betaVal{jj}(jj,:).^2 .* txPowerNormalization(jj,:);
interferencePowerTheo = interfTerm1 + interfTerm2 + interfTerm3;