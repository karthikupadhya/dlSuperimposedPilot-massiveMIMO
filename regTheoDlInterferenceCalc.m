% Aalto University
% Author : Karthik Upadhya
% This function computes Equation (6) in the paper.

function [sigPowerTheo , interferencePowerTheo ] = regTheoDlInterferenceCalc(channelParameter,M,L,K,tau)
betaVal = channelParameter.betaVal;
dlNoiseVar= channelParameter.dlNoiseVar;
ulNoiseVar= channelParameter.ulNoiseVar(1);
sharedPilotCellIdx     = channelParameter.sharedPilotCellIdx;

for ll = 1:L
    for kk = 1:K
        invTxPowerNormalization(ll,kk) = 0;
        for pp = 1:L
            if (sharedPilotCellIdx{ll}(pp) == 1)
                invTxPowerNormalization(ll,kk) = invTxPowerNormalization(ll,kk) + M * betaVal{ll}(pp,kk);
            end
        end
    end
end
invTxPowerNormalization = invTxPowerNormalization + ulNoiseVar * M / tau;
txPowerNormalization    = 1./invTxPowerNormalization;

jj = 1;
for mm = 1:K
    interfTerm1(mm) = 0;
    for ll = 1:L
        for kk = 1:K
            tVal1           = M^2 * betaVal{ll}(jj,mm)^2 * (kk == mm) * (sharedPilotCellIdx{ll}(jj)==1);
            tVal2           = 0;
            for pp = 1:L
                if sharedPilotCellIdx{ll}(pp) == 1
                    tVal2 = tVal2 + M * betaVal{ll}(jj,mm) * betaVal{ll}(pp,kk);
                end
            end
            tVal3           = M * ulNoiseVar * betaVal{ll}(jj,mm) / tau ;
            interfTerm1(mm) = interfTerm1(mm) + txPowerNormalization(ll,kk) * (tVal1 + tVal2 + tVal3);
        end
    end
    
    interfTerm2(mm) = - M^2 * betaVal{jj}(jj,mm)^2 .* txPowerNormalization(jj,mm);
    
    interfTerm3(mm) = dlNoiseVar;    
end


interferencePowerTheo = interfTerm1 + interfTerm2 + interfTerm3;
sigPowerTheo          = M^2 * betaVal{jj}(jj,:).^2 .* txPowerNormalization(jj,:);
