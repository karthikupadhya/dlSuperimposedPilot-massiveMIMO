% Author : Karthik Upadhya
% Generate centers of tesselated hexagons

function cellCenter = genCellCenters(numCell)

numTier = ceil( (-3 + sqrt(9 + 12 * (numCell - 1))) / 6 ) ;

for nn = 1:numTier
    s = [0,nn,-nn];
    b = [0,1,-1];
    for cc = 1: 6 * nn
        while(1);
            if sum(abs(s + b))/2 == nn
                cellCenterAxial(3*(nn-1)*nn + cc,:) = s + b;
                break;
            else
                b = [-b(2),b(1) + b(2),-b(1)];
            end
        end
        s = s + b;
    end
end

cellCenter = sqrt(3) * sqrt(3)/2 * cellCenterAxial(:,1) + 1i * ( sqrt(3)/2 * cellCenterAxial(:,1) + sqrt(3) * cellCenterAxial(:,2) );
cellCenter = [0;cellCenter];
cellCenter = cellCenter(1:numCell);
end