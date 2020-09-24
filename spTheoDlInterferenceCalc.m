% Aalto University
% Author : Karthik Upadhya
% This function computes Equation (16) in the paper.

function [sigPowerTheo,interferencePowerTheo] = spTheoDlInterferenceCalc(channelParameter,M,C,L,K,rhoD,rhoP)

betaVal = channelParameter.betaVal;
dlNoiseVar= channelParameter.dlNoiseVar;
ulNoiseVar= channelParameter.ulNoiseVar(1);
sharedPilotCellIdx     = channelParameter.sharedPilotCellIdx;

for ll = 1:L
    for kk = 1:K
        
        tVal1 = 0;
        for nn = 1:L
            if sharedPilotCellIdx{ll}(nn) == 1
                tVal1 = tVal1 + M * betaVal{ll}(nn,kk);
            end
        end
        
        tVal2 = 0;
        for nn = 1:L
            for pp = 1:K
                tVal2 = tVal2 + M * (rhoD^2 / (C*rhoP^2)) * betaVal{ll}(nn,pp);
            end
        end
            
        tVal3 = M * ulNoiseVar / (C * rhoP^2);
        
        invTxPowerNormalization(ll,kk) = tVal1 + tVal2 + tVal3;
    end
end
txPowerNormalization = 1./invTxPowerNormalization;

jj = 1;
for mm = 1:K    
    interfTerm1(mm) = 0;
    for ll = 1:L
        for kk = 1:K
            tVal1 = M^2 * betaVal{ll}(jj,mm)^2 * (kk == mm) * (sharedPilotCellIdx{ll}(jj) == 1);
            
            tVal2 = 0;
            for nn = 1:L
                if sharedPilotCellIdx{ll}(nn) == 1
                    tVal2 = tVal2 + M * betaVal{ll}(jj,mm) * betaVal{ll}(nn,kk);
                end
            end
            
            tVal3 = M^2 * rhoD^2 / (C * rhoP^2) * betaVal{ll}(jj,mm)^2;
            
            tVal4 = 0;
            for nn = 1:L
                for pp = 1:K
                    tVal4 = tVal4 + M * rhoD^2 / ( C * rhoP^2 ) * betaVal{ll}(jj,mm) * betaVal{ll}(nn,pp);
                end
            end
            
            tVal5 = M * ulNoiseVar / (C * rhoP^2) * betaVal{ll}(jj,mm);
            
            interfTerm1(mm) = interfTerm1(mm) + txPowerNormalization(ll,kk) * ( tVal1 + tVal2 + tVal3 + tVal4 + tVal5 );
        end
    end
    
    interfTerm2(mm) = - M^2 * txPowerNormalization(jj,mm) * betaVal{jj}(jj,mm)^2;
    
    interfTerm3(mm) = dlNoiseVar;
end


interferencePowerTheo = interfTerm1 + interfTerm2 + interfTerm3;
sigPowerTheo          = M^2 * txPowerNormalization(jj,:) .* betaVal{jj}(jj,:).^2 ;