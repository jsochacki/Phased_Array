clear all
clc

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
   
   % create omnidirectional characteristic 
   iPattern = zeros(1, length(theta_vector_degrees)); 
   
   % Calculate Array Factor
   for nn = 1:1:length(theta_vector_rads)
      [AF, AF_dB, AV] = Uniform_Planar_Array(phi_val*(pi / 180), theta_vector_rads(nn), frequency, d, weights);
      weight_vs_angle{n}{nn} = {theta_vector_rads(nn), AV};      
      %Combine for full characteristics
      F(n, nn) = AF_dB + iPattern(nn);
   end
end 

weights = weight_vs_angle{end}{90}{2}'.';

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