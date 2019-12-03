clear all
clc

frequency = 28e9;
d_a = 0.005;
weights = ones(1, 10);
array_type = 'circular';

[phi_theta, phi_theta_dB, phi, theta] = generate_phi_theta(frequency, d_a, weights, array_type);

patternCustom(phi_theta_dB', theta, phi);
patternCustom(phi_theta_dB', theta, phi,'CoordinateSystem','rectangular');
patternCustom(phi_theta_dB', theta, phi,'CoordinateSystem','polar','Slice', ...
              'phi','SliceValue',[45 90 180 359]);
patternCustom(phi_theta_dB', theta, phi,'CoordinateSystem','rectangular', ...
              'Slice','phi','SliceValue',[45 90 180 359]);
 
%converts the pattern from phi - theta coordinates to az-el coordinates.
[pattern_azel,az,el] = phitheta2azelpat(phi_theta_dB, phi, theta);

freqVector  = [24 28].*1e9;        % Frequency range for element pattern 
antenna = phased.CustomAntennaElement('FrequencyVector', freqVector, ...
                                      'AzimuthAngles', az, ...
                                      'ElevationAngles', el, ...
                                      'MagnitudePattern', pattern_azel, ...
                                      'PhasePattern', zeros(size(pattern_azel)));

fmax = freqVector(end);
pattern(antenna,fmax,'Type','powerdb')
view(90,0)

%more here
%https://www.mathworks.com/help/antenna/examples/visualize-custom-radiation-patterns.html
%https://www.mathworks.com/help/phased/examples/antenna-array-analysis-with-custom-radiation-pattern.html