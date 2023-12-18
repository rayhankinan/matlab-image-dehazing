function atmLight = computeAtmLight(A, numLevels)
% Quad-tree decomposition of dark channel is used for estimation
% of atmospheric light as given in reference paper [2]

[dm, dn, ~] = size(A);
Q = [1, 1, dm, dn];

% Heuristic: Window size for finding spatial minimum value for dark channel
winSize = ceil(min(size(A,1),size(A,2)) / 400 * 15);

% Create dark channel from non-overlapping patches
minFiltStrel = strel('square', winSize);
SE = min(A, [], 3);
darkChannel = imerode(SE, minFiltStrel);

% Iterate over quad-tree levels
for ii = 1:numLevels
    % Quadrants indices matrix
    quadrantIndex = ([
        Q(1), Q(2), (Q(1)+Q(3))/2, (Q(2)+Q(4))/2;
        Q(1),((Q(2)+Q(4))/2)+1, (Q(1)+Q(3))/2, Q(4);
        ((Q(1)+Q(3))/2)+1, Q(2), Q(3), ((Q(2)+Q(4))/2);
        ((Q(3)+Q(1))/2)+1, ((Q(4)+Q(2))/2)+1, Q(3), Q(4)
        ]);
    quadrantIndex = round(quadrantIndex);
    
    % Decomposition of dark channel into four quadrants
    firstQuadrant = darkChannel(quadrantIndex(1,1):quadrantIndex(1,3), quadrantIndex(1,2):quadrantIndex(1,4));
    secondQuadrant = darkChannel(quadrantIndex(2,1):quadrantIndex(2,3), quadrantIndex(2,2):quadrantIndex(2,4));
    thirdQuadrant = darkChannel(quadrantIndex(3,1):quadrantIndex(3,3), quadrantIndex(3,2):quadrantIndex(3,4));
    fourthQuadrant = darkChannel(quadrantIndex(4,1):quadrantIndex(4,3), quadrantIndex(4,2):quadrantIndex(4,4));
    
    % Computation of mean for each quadrant
    mu(1) = mean(firstQuadrant(:));
    mu(2) = mean(secondQuadrant(:));
    mu(3) = mean(thirdQuadrant(:));
    mu(4) = mean(fourthQuadrant(:));
    
    % Selecting maximum average intensity quadrant
    [~, ind] = max(mu);
    Q = quadrantIndex(ind, :);
end

% Selecting bright image pixels based on final decomposed quadrant
img = A(Q(1):Q(3), Q(2):Q(4), :);
[mm, nn, pp] = size(img);
brightIm = ones(mm, nn, pp);

% Minimum Euclidean distance based bright pixel estimation (= atmLight)
equiDist = sqrt((abs(brightIm - img)) .^ 2);
equiDistImage = sum(equiDist, 3);
equiDistVector = equiDistImage(:);
imageVector = reshape(img, mm * nn, pp);
[~, index] = min(equiDistVector);
atmLight = imageVector(index, :);

end