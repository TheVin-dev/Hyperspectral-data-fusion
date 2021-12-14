function [fixedMatchedPoints,movingMatchedPoints] = extractPoints(MOVING,FIXED) 
% Feature-based techniques require license to Computer Vision Toolbox
% Normalize MOVING image

% Get linear indices to finite valued data
finiteIdx = isfinite(MOVING(:));

% Replace NaN values with 0
MOVING(isnan(MOVING)) = 0;

% Replace Inf values with 1
MOVING(MOVING==Inf) = 1;

% Replace -Inf values with 0
MOVING(MOVING==-Inf) = 0;

% Normalize input data to range in [0,1].
MOVINGmin = min(MOVING(:));
MOVINGmax = max(MOVING(:));
if isequal(MOVINGmax,MOVINGmin)
    MOVING = 0*MOVING;
else
    MOVING(finiteIdx) = (MOVING(finiteIdx) - MOVINGmin) ./ (MOVINGmax - MOVINGmin);
end

% Default spatial referencing objects
fixedRefObj = imref2d(size(FIXED));
movingRefObj = imref2d(size(MOVING));

% Detect SURF features
fixedPoints = detectSURFFeatures(FIXED,'MetricThreshold',500.000000,'NumOctaves',4,'NumScaleLevels',6);
movingPoints = detectSURFFeatures(MOVING,'MetricThreshold',500.000000,'NumOctaves',4,'NumScaleLevels',6);

% Extract features
[fixedFeatures,fixedValidPoints] = extractFeatures(FIXED,fixedPoints,'Upright',true);
[movingFeatures,movingValidPoints] = extractFeatures(MOVING,movingPoints,'Upright',true);

% Match features
indexPairs = matchFeatures(fixedFeatures,movingFeatures,'MatchThreshold',18.499985,'MaxRatio',0.185000);
fixedMatchedPoints = fixedValidPoints(indexPairs(:,1));
movingMatchedPoints = movingValidPoints(indexPairs(:,2));
MOVINGREG.FixedMatchedFeatures = fixedMatchedPoints;
MOVINGREG.MovingMatchedFeatures = movingMatchedPoints;



end 