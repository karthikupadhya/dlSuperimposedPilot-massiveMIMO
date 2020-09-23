function output = cNormrnd(varargin)

mu      = varargin{1};
sigma   = varargin{2};

if nargin >= 3
    M   = varargin{3};
else
    M   = 1;
end

if nargin >= 4
    N   = varargin{4};
else
    N   = 1;
end

output = complex(normrnd(real(mu),sigma/sqrt(2),M,N) , normrnd(imag(mu),sigma/sqrt(2),M,N));