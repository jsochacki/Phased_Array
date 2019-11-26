function [AF, AF_dB, AV] = Uniform_Linear_Array(theta, f, d, weights)
%Uniform_Linear_Array Calculate array factor of linear antenna array for a given
% theta [rad] based on an array defined by:
%   frequency f [Hz],
%   element spacing d [m],
%   and element complex weight

% Here a zero value for theta corresponds to the line that the array makes
% i.e. a vector with theta zero points directly in line with the array 
% i.e. parallel to the length of the array and
% a theta of 90 degrees if perpendicular to the length of the array

% Note that this linear array has a reference point that is at the very end
% of the array and not in the center

% Since k^2 = (2*pi*f)^2*mu*epsi then k = 2*pi*f*sqrt(mu_r*mu_0*epsi_r*epsi_0)
% and c = lambda * f => f = c / Lambda and c = 1/sqrt(mu_r*mu_0*epsi_r*epsi_0)
% k = 2*pi*(c / lambda)*sqrt(mu_r*mu_0*epsi_r*epsi_0) => k = (2*pi)/lambda

if (size(weights, 2) < size(weights, 1)), weights = weights.';, end;

% Physical constants
c = 299792458;

% Derived values
N = length(weights);
lambda = c / f;
k = (2 * pi) / lambda;
PHI = (k * d * cos(theta));

n = 1:1:N;
AV = weights .* exp(j * PHI * (n-1));
AF = sum(AV);

AF_dB = 10*log10(AF*AF');

%For sanity check if you like
%(sin((N/2)*phi)/sin((1/2)*phi))^2 - AF*AF' < 1e-14
%20*log10(sin((N/2)*phi)/sin((1/2)*phi)) - AF_dB < 1e-14

end