function regMse = regMseCalc(channelParameter,L,K,tau)

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
    
    term2 = ulNoiseVar / tau;
    
    regMse(mm) = term1 + term2;
end