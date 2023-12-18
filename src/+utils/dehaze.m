function [correctedA, thicknessOfHaze, atmLight] = dehaze(A, epsilon, omega, minimumTransmissionValue)
% SimpleDCP
% This function computes dark-channel prior, and refines it using guided
% filter, only across channel elements are considered for dark channel
% estimation

% Check if input image is RGB
isRGB = size(A, 3) == 3;

% 1. Estimate atmospheric light
atmLight = utils.computeAtmLight(A, 5, 30);

if isRGB
    atmLight = reshape(atmLight, [1 1 3]);
end

% 2. Estimate transmission t(x)
normI = min(A, [] , 3);
minI = 1 - normI;

% 3. Use guided filtering to refine the transmission map
rawTransmissionMap = utils.algImguidedFilter(minI, minI, [5 5], epsilon);
transitionMap = min(1, max(0, rawTransmissionMap));

% Thickness of haze in input image is second output of imreducehaze.
% Thickness Map does not depends on amount value.
thicknessOfHaze = 1 - transitionMap;

% Omega value is set to no more than 1, to leave some of haze in restored image for
% natural appearance of dehazed scene
omegaTransitionMap = 1 - omega * (1 - transitionMap);

% This lower bound preserves a small amount of haze in dense haze regions
t0 = cast(minimumTransmissionValue, 'like', A);

% Recover scene radiance
rawRadianceMap = atmLight + (A - atmLight) ./ max(omegaTransitionMap, t0);
radianceMap = min(1, max(0, rawRadianceMap));

% Global stretching of radiance map (Optional)
correctedA = utils.globalStretching(radianceMap, 0.75, [0.001, 0.999], 0.8);

% Reshape atmLight to 1 x 3 vector if input image is RGB
if isRGB
    atmLight = double(reshape(atmLight, [1 3]));
else
    atmLight = double(atmLight);
end

end