function [BW, T] = niblackThreshold(I, windowSize, k)
%NIBLACKTHRESHOLD  Niblack's local adaptive threshold, implemented from scratch.
%
%   [BW,T] = niblackThreshold(I, windowSize, k)
%
% MATLAB has NO built-in Niblack.  imbinarize(...,'adaptive') uses Bradley's
% method (local MEAN with a sensitivity offset) which is a different formula,
% so we implement Niblack properly here:
%
%       T(x,y) = mu(x,y) + k * sigma(x,y)
%
% where mu and sigma are the mean and standard deviation of the intensities in
% a windowSize x windowSize neighbourhood centred on (x,y).
%
% Reference: W. Niblack, "An Introduction to Digital Image Processing",
%            Prentice-Hall, 1986, pp. 115-116.
%
% INPUTS
%   I           MxN double, expected on [0,1]
%   windowSize  odd integer, side of the local window (default 25)
%   k           Niblack's constant (default 0.2).  Niblack's original text
%               uses k = -0.2 for extracting DARK text from a LIGHT page.  Our
%               feature of interest (the wave area) is BRIGHT on a darker map,
%               so we use k > 0, which pushes the threshold ABOVE the local
%               mean and keeps only pixels that stand out locally.
%
% OUTPUTS
%   BW  MxN logical, true where I > T
%   T   MxN double, the threshold surface itself (useful to plot)
%
% HOW THE LOCAL STATISTICS ARE COMPUTED
%   E[x]  = imfilter(I, ones(w)/w^2)                 -> local mean
%   E[x^2]= imfilter(I.^2, ones(w)/w^2)              -> local mean of squares
%   var   = E[x^2] - E[x]^2                          -> the computational
%                                                       formula for variance
%   A box filter of ones/w^2 is separable and O(1) per pixel in MATLAB's
%   implementation, so this is far cheaper than calling std2 in a loop.

if nargin < 2 || isempty(windowSize), windowSize = 25;  end
if nargin < 3 || isempty(k),          k          = 0.2; end
if mod(windowSize,2) == 0, windowSize = windowSize + 1; end   % force odd

I = double(I);
h = ones(windowSize) / windowSize^2;

mu  = imfilter(I,    h, 'replicate', 'corr');
mu2 = imfilter(I.^2, h, 'replicate', 'corr');

% max(...,0) guards against tiny negative values from floating-point round-off
% in flat regions where E[x^2] and E[x]^2 are numerically identical.
sigma = sqrt(max(mu2 - mu.^2, 0));

T  = mu + k * sigma;
BW = I > T;
end
