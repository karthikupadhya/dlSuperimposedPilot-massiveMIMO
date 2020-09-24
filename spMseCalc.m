function spMse = spMseCalc(channelParameter,C,L,K,rhoD,rhoP)

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
            term2 = term2 + betaVal{jj}(ll,kk) * rhoD^2 / ( C * rhoP^2 );
        end
    end
    
    term3 = ulNoiseVar / ( rhoP^2 * C );
    
    spMse(mm) = term1 + term2 + term3;
end