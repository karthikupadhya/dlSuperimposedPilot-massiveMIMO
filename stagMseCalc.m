function stagMse = stagMseCalc(channelParameter,C,L,K,pD,pP)

betaVal                 = channelParameter.betaVal;
sharedPilotCellIdx      = channelParameter.sharedPilotCellIdx;
ulNoiseVar              = channelParameter.ulNoiseVar(1);

jj = 1;
for mm = 1:K
    
    term1 = 0;
    for ll = 1:L
        if (ll ~= jj) && (sharedPilotCellIdx{jj}(ll) == 1);
            term1 = term1 + betaVal{jj}(ll,mm);
        end
    end
    
    term2 = 0;
    for ll = 1:L
        for kk = 1:K
            if (ll ~= jj) && (sharedPilotCellIdx{jj}(ll) ~= 1);
                term2 = term2 + betaVal{jj}(ll,kk) * pD / ( C * pP );
            end
        end
    end
    
    term3 = ulNoiseVar / ( pP * C );
    
    stagMse(mm) = term1 + term2 + term3;
end