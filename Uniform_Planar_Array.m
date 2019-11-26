function [AF, AF_dB, AV] = Uniform_Planar_Array(phi, theta, f, d, weights)
%Uniform_Planar_Array Calculate array factor of linear antenna array for a given
% theta [rad] based on an array defined by:
%   frequency f [Hz],
%   element spacing d [m],
%   and element complex weight

% Here a zero value for theta corresponds to the normal vector to the plane
% that the array makes
% i.e. a vector with theta zero points directly in perpendicular to the
% array face and a theta of 90 degrees if parallel/in plane to the face of the array

% Note that this planar array has a reference point that is at the very end
% of the array in the corner and not in the center

% Also note that this planar array does NOT need to be square
% it is M elements in the x and N elements in the y direction

% Since k^2 = (2*pi*f)^2*mu*epsi then k = 2*pi*f*sqrt(mu_r*mu_0*epsi_r*epsi_0)
% and c = lambda * f => f = c / Lambda and c = 1/sqrt(mu_r*mu_0*epsi_r*epsi_0)
% k = 2*pi*(c / lambda)*sqrt(mu_r*mu_0*epsi_r*epsi_0) => k = (2*pi)/lambda

% Physical constants
c = 299792458;

% Derived values
M = size(weights, 1);
N = size(weights, 2);
lambda = c / f;
k = (2 * pi) / lambda;

PHI_X = (k * d * sin(theta) * cos(phi));
PHI_Y = (k * d * sin(theta) * sin(phi));

m = 1:1:M; n = 1:1:N;
M_VECTOR = exp(j * PHI_X * (m-1));
N_VECTOR = exp(j * PHI_Y * (n-1));
AV = weights .* ((M_VECTOR.')*(N_VECTOR));
AF = sum(sum(AV));

AF_dB = 10*log10(AF*AF');

end