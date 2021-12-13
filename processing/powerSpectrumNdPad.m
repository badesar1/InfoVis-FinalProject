function [radialavg, linearavg, rad] = powerspectrum(original, denoised, sigma, nbins)
    % function to compute the power spectrum 
    
    % required inputs - 2D original image, denoised image, noisemap
    % optional inputs - nbins (number of histogram bins over which to
    % average the power spectrum)
    % make sure the crop original, denoised, and sigma to a square grid
    
    % outputs: radialavg - the power spectrum, 
    % linearavg - the autocorrelation
    % rad - normalized radius from k-space center (if you want the radius
    % in voxel units then return the variable "interval" instead
    
    
    sx = size(original,1);
    sy = size(original,2);
    A = numel(original);
    %A = sx*sy;
    
    % compute normalized residual and psd images
    padsize = size(original)*2-1;
    %padsize = 0;
    ep = (denoised - original)./(sigma);
    %epf = fftshift(fft2(ep,padsize(1), padsize(2)));
    epf = fftshift(fft2(ep));
    %epf = fftshift(fft2(ep,sx*2-1,sy*2-1));
    gamma  = (conj(epf).*epf)./A;
    igamma = ifftshift(ifft2(gamma));
    
    % get radii, gamma must be square
    [u, v, w] = size(gamma);
    [x, y, z] = meshgrid(-floor(u/2):floor(u/2)+rem(u,2)-1,...
        -floor(v/2):floor(v/2)+rem(v,2)-1,...
        -floor(w/2):floor(w/2)+rem(w,2)-1);
    R = floor(sqrt(x.^2 + y.^2 + z.^2)); 

    % bin setup
    rad = unique(R(:));
    if nargin < 4
        nbins = length(rad);
    end
    minr  = min(R(:));
    maxr  = max(R(:));
    delta = (maxr-minr)/nbins;
    interval = linspace(minr, maxr, nbins)-delta/2;
    
    ahist = zeros(1,nbins);
    whist = zeros(1,nbins);
    rhist = zeros(1,nbins);
    for i=1:nbins
        % which bin are we in and what are the corresponding radii indices
%         bin = find(interval < rad(i), 1, 'last') + 1;
%         inds = find(R(:)>=interval(bin-1) & R(:)<interval(bin));
%         
        inds = find(R==rad(i));
        % gamma and radius histograms
        ahist(i) = ahist(i) + sum(igamma(inds));
        whist(i) = whist(i) + sum(gamma(inds));
        rhist(i) = rhist(i) + length(inds);
    end
    radialavg = whist./rhist;
    linearavg = ahist./rhist;
    rad = rad./maxr.*sqrt(2)/2;
end