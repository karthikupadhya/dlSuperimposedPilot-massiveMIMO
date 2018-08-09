function [idx] = genSharedPilotCellIndex(numCell,pilotReuseRatio)

numTier     = ceil( (-3 + sqrt(9 + 12 * (numCell - 1))) / 6 ) ;

cellCenter = genCellCenters(numCell);

switch pilotReuseRatio
    case 1
        ii = 1;
        jj = 0;
    case 3
        ii = 1;
        jj = 1;
    case 7
        ii = 2;
        jj = 1;
end

%%% - Generate pattern - %%%
U = [0,1;1,0;1,-1;0,-1;-1,0;-1,1].';
C = [3/2,0;sqrt(3)/2,sqrt(3)];
R = [cos(pi/3) -sin(pi/3) ; sin(pi/3) cos(pi/3)];
for kk = 1:6
    N(:,kk) = C * U(:,kk) * ii + R * C * U(:,kk) * jj;
end
reusePattern = [0;complex(N(1,:),N(2,:)).'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% - Replicate pattern - %%%
for nn = 1:numTier
    if nn == 1
        replicatePattern{nn} = reusePattern;
    else        
        for kk = 1:numel(replicatePattern{nn-1})
            replicatePattern{nn}(:,kk) = replicatePattern{nn-1}(kk) + reusePattern;
        end
    end
end
replicatePattern{1} = vec(replicatePattern{nn});

for ll = 1:numCell
    reuseCellCenter{ll} = replicatePattern{1} + cellCenter(ll);    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% - Generate index - %%%
if pilotReuseRatio > 7
    error('Code built for pilotReuseRatio <= 7')
end
for ll = 1:7 % Code prepared for reuse ratio uptill 7
    idx{ll} = zeros(length(cellCenter),1);
    for kk = 1:length(reuseCellCenter{ll})
        idx{ll}( abs( reuseCellCenter{ll}(kk) - cellCenter ) < 1e-6 ) = 1;
    end
end

for ll = 1:numCell
    for pp = 1:7
        if idx{pp}(ll) == 1;
            idx{ll} = idx{pp};
            break
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%
end