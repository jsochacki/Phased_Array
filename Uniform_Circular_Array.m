function [AF, AF_dB, AV] = Uniform_Circular_Array(phi, theta, f, a, weights)
%Uniform_Planar_Array Calculate array factor of linear antenna array for a given
% theta [rad] based on an array defined by:
%   frequency f [Hz],
%   array radius [m],
%   and element complex weight

% Here a zero value for theta corresponds to the normal vector to the plane
% that the circular array makes
% i.e. a vector with theta zero points directly in perpendicular to the circular
% array face and a theta of 90 degrees if parallel/in plane to the circular face of the array

% Note that this circular array has a reference point that is at the very
% center of the circle and the first element in the circle is element 1
% i.e. phi_n = 0 at element 1 and phi_n = ((2*pi)/N)*(N-1) at element N

% Note that this array does not have an element in the middle

% Since k^2 = (2*pi*f)^2*mu*epsi then k = 2*pi*f*sqrt(mu_r*mu_0*epsi_r*epsi_0)
% and c = lambda * f => f = c / Lambda and c = 1/sqrt(mu_r*mu_0*epsi_r*epsi_0)
% k = 2*pi*(c / lambda)*sqrt(mu_r*mu_0*epsi_r*epsi_0) => k = (2*pi)/lambda

% Physical constants
c = 299792458;

% Derived values
N = length(weights);
lambda = c / f;
k = (2 * pi) / lambda;
delta_phi = (2 * pi) / N;

n = 0:1:N-1;
PHI = k * a * sin(theta);
AV = weights .* exp(j * PHI * cos(phi - (delta_phi * n)));
AF = sum(AV);

AF_dB = 10*log10(AF*AF');

end