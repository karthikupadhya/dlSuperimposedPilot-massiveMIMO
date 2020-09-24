% Aalto University
% Author : Karthik Upadhya
% This function computes Equation (27) in the paper.

function [sigPowerTheo , interferencePowerTheo] = stagTheoDlInterferenceCalc(channelParameter,M,tau,L,K,pD,pP)
betaVal = channelParameter.betaVal;
dlNoiseVar= channelParameter.dlNoiseVar;
ulNoiseVar= channelParameter.ulNoiseVar(1);
sharedPilotCellIdx     = channelParameter.sharedPilotCellIdx;

for ll = 1:L
    for kk = 1:K
        tVal1 = 0;
        for pp = 1:L
            if sharedPilotCellIdx{ll}(pp) == 1
                tVal1 = tVal1 +  M * betaVal{ll}(pp,kk);
            end
        end
        
        tVal2 = 0;
        for pp = 1:L
            if (sharedPilotCellIdx{ll}(pp) == 0) && (pp ~= ll)
                for nn = 1:K
                    tVal2 = tVal2 + M * betaVal{ll}(pp,nn) * pD / (tau * pP);
                end
            end
        end        
        invTxPowerNormalization(ll,kk) = tVal1 + tVal2 + ulNoiseVar * M / (tau * pP);        
    end
end
txPowerNormalization = 1./invTxPowerNormalization;

jj = 1;
for mm = 1:K
    interfTerm1(mm) = 0;
    for ll = 1:L
        for kk = 1:K
            tVal1 = M^2 * betaVal{ll}(jj,mm)^2 * (sharedPilotCellIdx{ll}(jj) == 1) * (kk == mm);
            
            tVal2 = 0;
            for pp = 1:L
                if sharedPilotCellIdx{ll}(pp) == 1
                    tVal2 = tVal2 + M * betaVal{ll}(jj,mm) * betaVal{ll}(pp,kk);
                end
            end
            
            tVal3 = M * ulNoiseVar * betaVal{ll}(jj,mm) / (pP * tau);
            
            tVal4 = 0;
            for pp = 1:L
                if (sharedPilotCellIdx{ll}(pp) == 0) && (pp~=ll)
                    for qq = 1:K
                        tVal4 = tVal4 + pD / (tau * pP) * ( M^2 * betaVal{ll}(jj,mm)^2 * (pp == jj) * (qq == mm) + M * betaVal{ll}(jj,mm) * betaVal{ll}(pp,qq) );
                    end
                end
            end
            
            interfTerm1(mm) = interfTerm1(mm) + txPowerNormalization(ll,kk) * (tVal1 + tVal2 + tVal3 + tVal4);
        end
    end
    
    interfTerm2(mm) = - M^2 * betaVal{jj}(jj,mm)^2 * txPowerNormalization(jj,mm);
    
    interfTerm3(mm) = dlNoiseVar;
end

interferencePowerTheo = interfTerm1 + interfTerm2 + interfTerm3;
sigPowerTheo          = M^2 * betaVal{jj}(jj,:).^2 .* txPowerNormalization(jj,:);