function maxData = computeRunningMaxAdaptive(x, G, Gprime, U, cfg)
%COMPUTERUNNINGMAXADAPTIVE Adaptively resolve streamwise maxima.

error('computeRunningMaxAdaptive:NotImplemented', ...
    ['Adaptive maximum detection has not yet been implemented. ', ...
     'Use cfg.maxDetection.method = ''grid'' for now.']);

maxData = struct();
end
