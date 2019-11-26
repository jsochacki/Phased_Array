clear all
clc

relative_filename = './Radiation_Patterns/DF1_3GHz_Gain_3D_pt.csv';

%Import the custom pattern (can add phase too if you want eve though it isn't in there now)
%Input columns need to be phi, theta, pattern_dB pattern_phase
%Here Boresight is at phi = 0 and theta = 90 and phi sweeps from boresight
%to back to boresignt around the wider portion of the beam while theta is
%from horizon to horizon over the skinny portion of the beam
[pattern_phitheta, phi, theta] = HFSS_pattern_import(relative_filename);
theta = theta - 90;

%NOTE THAT FOR THIS ELEMENT THE BEAM IS POINTING AT THETA = 90 (THETA = 0
%AFTER SUBTRACTING 90 LIKE DONE SO ABOVE) AND THAT PHI SWEEPS FROM
%BORESIGHT AROUND THE BACK OF THE BEAM AND NOT AROUND THE BEAM IN A CIRCLE
%ABOUT BORESIGHT LIKE ONE MIGHT EXPECT

freq_range = [24e9 28e9];   % [start_frequency stop_frequency] [Hz]
freq_points = 101;          % No. of frequency points
d = 0.005;                  % Element spacing [m]
weights = ones(10, 10);    % Amplitude weights for N elements from 0 to 1  

freq_step = (freq_range(2) - freq_range(1)) / (freq_points - 1); 
frequency_vector = freq_range(1):freq_step:freq_range(2); 
theta_vector_degrees = -90:1:90; 
theta_vector_rads = theta_vector_degrees * (pi / 180); 

phi_val = 0;
for n = 1:1:length(frequency_vector)
   frequency = frequency_vector(n);
   
   % Calculate Array Factor
   for nn = 1:1:length(theta_vector_rads)
      % extract element characteristic
      %This way we are placing antenna theta in the array theta direction
      Pattern = pattern_phitheta(find(theta == theta_vector_degrees(nn)), phi_val + 1);
      
      [AF, AF_dB, AV] = Uniform_Planar_Array(phi_val*(pi / 180), theta_vector_rads(nn), frequency, d, weights);
      weight_vs_angle{n}{nn} = {theta_vector_rads(nn), AV};
      %Combine for full characteristics
      %Note that if you have complex pattern info you need to multiply the
      %pattern response and the array factor in the linear domain and then 
      %convert to dB
      F(n, nn) = AF_dB + Pattern;
   end
end

weights = weight_vs_angle{end}{45}{2}'.';

figure(1);
surf(theta_vector_degrees, frequency_vector, F);
ax = gca;
ax.YAxis.TickLabelFormat = '%,.1g';
set(ax, 'FontSize',12)
rotate3d on;
xlim([-90 90]);
ylabel('\fontsize{14}Frequency / GHz');
xlabel('\fontsize{14}Angle / °');
zlabel('\fontsize{14}Level / dBm');

figure(2); 
plot(theta_vector_degrees, F(end,:)); 