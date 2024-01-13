function [MDL sigmas] = calcMDL(MTMs)
s = size(MTMs);
lambdaCount = s(1);
totalModes = s(2);
sigmas = zeros(lambdaCount,totalModes,'single');
%MDL = zeros(1,lambdaPoints,'single');
for lambdaIdx=1:(lambdaCount)
    %SVD
    J0 = squeeze(MTMs(lambdaIdx,:,:));
    [U S V] = svd(J0);
    sigmas(lambdaIdx,:) = diag(S);
end
s2 = (sigmas.^2).';
MDL = 10.*log10(max(s2)./min(s2));