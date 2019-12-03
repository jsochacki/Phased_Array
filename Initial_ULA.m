clear all
clc

relative_filename = './Radiation_Patterns/DF1_3GHz_Gain_3D_pt.csv';

%Import the custom pattern (can add phase too if you want eve though it isn't in there now)
%Input columns need to be phi, theta, pattern_dB pattern_phase
[pattern_phitheta, phi, theta] = HFSS_pattern_import(relative_filename);

%converts the pattern from phi - theta coordinates to az-el coordinates.
[pattern_azel,az,el] = phitheta2azelpat(pattern_phitheta, phi, theta);


freqVector  = [2.75 3.25].*1e9;        % Frequency range for element pattern 
antenna = phased.CustomAntennaElement('FrequencyVector', freqVector, ...
                                      'AzimuthAngles', az, ...
                                      'ElevationAngles', el, ...
                                      'MagnitudePattern', pattern_azel, ...
                                      'PhasePattern', zeros(size(pattern_azel)));

fmax = freqVector(end);
pattern(antenna,fmax,'Type','powerdb')
view(90,0)


c = 299792458;
lambda = c / fmax;
array = phased.UCA('Element',antenna,'NumElements',4,'Radius',lambda)

pattern(array,fmax,'PropagationSpeed',c,'Type','powerdb',...
    'CoordinateSystem','UV');
 
pattern(array,fmax,-1:0.01:1,0,'PropagationSpeed',c,'Type','powerdb', ...
    'CoordinateSystem','UV')
axis([-1 1 -50 0]);



f0 = 3e9;
scanPhi = -30:30;

steeringvec = phased.SteeringVector('SensorArray',array,...
    'PropagationSpeed',c);            

arrayresp = phased.ArrayResponse('WeightsInputPort',true,...
    'SensorArray',array);
                             
scanTheta = zeros(1,numel(scanPhi));
scanAngles  = [scanPhi;scanTheta];
weights = steeringvec(f0,scanAngles);       % Calculate Weights for steering

az = -180:180;
el = zeros(1,numel(az));
ang_pairs = [az;el];

arrayresp2d = ArrayResponseDemo2DPolarScope;             % Initialize scope

for t = 1:length(scanTheta)                        % Calculate response
    w = weights(:,t);
    temp_out = abs(arrayresp(f0,ang_pairs,w)); 
    temp_out = temp_out/max(temp_out);
    arrayresp2d(temp_out);                         % Display on scope
end


% Derive the elevation pattern.
el_ang = -90:90;
arrayresp = phased.ArrayResponse('SensorArray',array, ...
    'PropagationSpeed',c);
el_pat = abs(arrayresp(fmax,el_ang));  % elevation pattern


% Plot the radar vertical diagram.
freespace_rng = 100;  % in km
ant_height = 20;      % in m

radarvcd(fmax,freespace_rng,ant_height,...
    'HeightUnit','m','RangeUnit','km',...
    'AntennaPattern',el_pat/max(el_pat),'PatternAngles',el_ang.');
 
 

subplot(211)
pattern(array,fc,-180:180,0,'CoordinateSystem','rectangular', ...
    'PropagationSpeed',c,'Type','powerdb')
title('Without steering')
subplot(212)
pattern(array,fc,-180:180,0,'CoordinateSystem','rectangular', ...
    'PropagationSpeed',c,'Type','powerdb','Weights',sv)
title('With steering')