function [rhoD,rhoP] = rhoLambdaCalc(channelParameter,numSymbolUlTimeSlot,numCell,numUser,numBsAntenna)

Cu = numSymbolUlTimeSlot;
sharedPilotCellIdx = channelParameter.sharedPilotCellIdx;
noiseVar = channelParameter.ulNoiseVar(1);
M  = numBsAntenna;
for jj = 1:numCell
    for mm = 1:numUser
        betaVal = channelParameter.betaVal{jj};
        betaValVec= vec(betaVal.');
        
        alpha1(jj,mm)  = 1/(Cu-1) * sum( betaValVec.^2 );
        alpha2(jj,mm)  = 1/(Cu-1) * sum( betaValVec ).^2;
        
        alpha3(jj,mm)  = Cu/(Cu-1) * sum(betaValVec) * sum( betaVal(sharedPilotCellIdx{jj} == 1,mm) );
        
        alpha4(jj,mm)  =  noiseVar * sum(betaValVec) / (Cu-1) + noiseVar^2 / Cu;
        
        alpha5(jj,mm)  =  noiseVar * sum( betaVal(sharedPilotCellIdx{jj} == 1,mm) );
        
        alpha6(jj,mm)  =  noiseVar / Cu * sum(betaValVec);        
    end
end

alpha1Val = mean(vec(alpha1));
alpha2Val = mean(vec(alpha2));
alpha3Val = mean(vec(alpha3));
alpha4Val = mean(vec(alpha4));
alpha5Val = mean(vec(alpha5));
alpha6Val = mean(vec(alpha6));

kappaVal  = sqrt( (alpha1Val + alpha2Val / M + alpha4Val/M + alpha6Val/M)/(alpha3Val + alpha4Val + alpha5Val) );

rhoD      = sqrt(1 / ( 1 + sqrt(M) * kappaVal ));
rhoP      = sqrt( 1 - rhoD^2 );
