function [x, B, H, eThresh, meta] = makeSyntheticProblem(cfg)
%MAKESYNTHETICPROBLEM Construct a rotating quadratic-form test problem.
%
%   [x, B, H, eThresh, meta] = makeSyntheticProblem(cfg)
%
% The synthetic problem has:
%
%   q(x) = B(x) w,       w ~ N(0,I_2)
%   H(x) = I
%   eThresh(x) = 1
%
% Therefore:
%
%   G(x) = B(x)'*B(x).
%
% Exact values of G'(x) are supplied for verification of continuous
% streamwise maximum detection.
%
% MATLAB version: R2020b

arguments
    cfg struct
end

% Permit the test resolution to be overridden.
Nx = 301;
xMax = 20;

if isfield(cfg, 'synthetic')
    if isfield(cfg.synthetic, 'Nx')
        Nx = cfg.synthetic.Nx;
    end

    if isfield(cfg.synthetic, 'xMax')
        xMax = cfg.synthetic.xMax;
    end
end

x = linspace(0, xMax, Nx).';

nq = 2;
r = 2;

B = cell(Nx, 1);
H = cell(Nx, 1);
eThresh = ones(Nx, 1);

GExact = cell(Nx, 1);
GprimeExact = cell(Nx, 1);

lambdaHistory = zeros(Nx, 2);
lambdaPrimeHistory = zeros(Nx, 2);
thetaHistory = zeros(Nx, 1);
thetaPrimeHistory = zeros(Nx, 1);

% Generator of planar rotations:
%
%   R'(theta) = R(theta)*J*thetaPrime.
J = [0, -1; 1, 0];

for n = 1:Nx
    xn = x(n);

    % ---------------------------------------------------------------------
    % First eigenvalue and derivative
    % ---------------------------------------------------------------------

    z11 = (xn - 6.0) / 1.10;
    z12 = (xn - 14.5) / 1.45;

    term11 = 0.72 * exp(-0.5 * z11^2);
    term12 = 0.30 * exp(-0.5 * z12^2);

    lambda1 = 0.025 + term11 + term12;

    lambda1Prime = ...
        -term11 * (xn - 6.0) / (1.10^2) ...
        -term12 * (xn - 14.5) / (1.45^2);

    % ---------------------------------------------------------------------
    % Second eigenvalue and derivative
    % ---------------------------------------------------------------------

    z21 = (xn - 9.5) / 1.80;
    z22 = (xn - 14.0) / 1.10;

    term21 = 0.18 * exp(-0.5 * z21^2);
    term22 = 0.52 * exp(-0.5 * z22^2);

    lambda2 = 0.015 + term21 + term22;

    lambda2Prime = ...
        -term21 * (xn - 9.5) / (1.80^2) ...
        -term22 * (xn - 14.0) / (1.10^2);

    % ---------------------------------------------------------------------
    % Rotation angle and derivative
    % ---------------------------------------------------------------------

    eta = (xn - 10.0) / 2.0;

    theta = ...
        0.10 * xn ...
        + 0.45 * tanh(eta);

    thetaPrime = ...
        0.10 ...
        + (0.45 / 2.0) * sech(eta)^2;

    R = [cos(theta), -sin(theta); ...
         sin(theta),  cos(theta)];

    Lambda = diag([lambda1, lambda2]);
    LambdaPrime = diag([lambda1Prime, lambda2Prime]);

    % Exact normalized energy matrix.
    Gn = R * Lambda * R.';

    % Since:
    %
    %   R' = thetaPrime * R*J,
    %
    % the matrix derivative is:
    %
    %   G' = R * [Lambda'
    %             + thetaPrime*(J*Lambda - Lambda*J)] * R'.
    Gpn = R * ( ...
        LambdaPrime ...
        + thetaPrime * (J * Lambda - Lambda * J)) * R.';

    % Remove roundoff-scale asymmetry.
    Gn = 0.5 * (Gn + Gn.');
    Gpn = 0.5 * (Gpn + Gpn.');

    % Choose B so that B'*B = G.
    B{n} = sqrt(Lambda) * R.';

    H{n} = eye(nq);

    GExact{n} = Gn;
    GprimeExact{n} = Gpn;

    lambdaHistory(n, :) = [lambda1, lambda2];
    lambdaPrimeHistory(n, :) = ...
        [lambda1Prime, lambda2Prime];

    thetaHistory(n) = theta;
    thetaPrimeHistory(n) = thetaPrime;
end

meta = struct();

meta.name = 'Rotating two-dimensional transient-growth test';
meta.description = [ ...
    'Two-dimensional Gaussian quadratic-form process with ', ...
    'nonmonotone eigenvalues and rotating eigendirections.'];

meta.nq = nq;
meta.r = r;
meta.Nx = Nx;

meta.GExact = GExact;
meta.GprimeExact = GprimeExact;

meta.lambdaHistory = lambdaHistory;
meta.lambdaPrimeHistory = lambdaPrimeHistory;
meta.thetaHistory = thetaHistory;
meta.thetaPrimeHistory = thetaPrimeHistory;

meta.hasGprime = true;
meta.gaussianConvention = 'real';
meta.factorOrientation = 'synthesis';

meta.cfgAtGeneration = cfg;

end
