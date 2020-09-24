function crlb = crlbCalc(channelParameter,C,L,K,rhoD,rhoP)

rxPower = channelParameter.rxPower;
sharedPilotCellIdx = channelParameter.sharedPilotCellIdx;
ulNoiseVar = channelParameter.ulNoiseVar(1);

jj = 1;
for mm = 1:K
    
    termNr = 0;
    termDr = 0;
    for ll = 1:L
        if sharedPilotCellIdx{jj}(ll) == 1;
            if ll ~= jj
                termNr = termNr + 1/( rhoD^2 + ulNoiseVar / (C * rxPower{jj}(ll,mm)) );
            end
            termDr = termDr + 1/( rhoD^2 + ulNoiseVar / (C * rxPower{jj}(ll,mm)) );
        end                
    end
    crlb(mm) = 1/(C * rhoD^2 / ulNoiseVar + 1 / rxPower{jj}(jj,mm)) * (1 + rhoP^2 * termNr) / ( 1+ rhoP^2 * termDr);
end