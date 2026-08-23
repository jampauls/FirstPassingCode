function [x, B, H, eThresh, meta] = runOWNSWorkflow(cfg)
%RUNOWNSWORKFLOW Interface to the user-specific OWNS propagator.

error('runOWNSWorkflow:NotImplemented', ...
    ['The OWNS interface has not yet been implemented. ', ...
     'Use cfg.dataSource = ''synthetic'' or ''matfile'' for now.']);

% Output declarations for MATLAB code analysis:
x = [];
B = {};
H = {};
eThresh = [];
meta = struct();
end