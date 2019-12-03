clear all
clc

freq_range = [24e9 28e9];   % [start_frequency stop_frequency] [Hz]
freq_points = 101;          % No. of frequency points
d = 0.005;                  % Element spacing [m]
weights = [1 1 1 1];    % Amplitude weights for N elements from 0 to 1  

freq_step = (freq_range(2) - freq_range(1)) / (freq_points - 1); 
frequency_vector = freq_range(1):freq_step:freq_range(2); 
theta_vector_degrees = -90:1:90; 
theta_vector_rads = theta_vector_degrees * (pi / 180); 

for n = 1:1:length(frequency_vector)
   frequency = frequency_vector(n);
   
   % create omnidirectional characteristic 
   iPattern = zeros(1, length(theta_vector_degrees)); 
   
   % Calculate Array Factor
   for nn = 1:1:length(theta_vector_rads)
      [AF, AF_dB, AV] = Uniform_Linear_Array(theta_vector_rads(nn), frequency, d, weights);
      weight_vs_angle{n}{nn} = {theta_vector_rads(nn), AV};
      %Combine for full characteristics
      F(n, nn) = AF_dB + iPattern(nn);
   end
end 

weights = weight_vs_angle{end}{45}{2}'.';
% %weights = weights .* exp(j*45*(pi/180).*(0:1:3))
% c = 299792458; lambda = c / freq_range(end); k1 = (2 * pi)*(30*(pi/180)) / lambda;
% weights = weights .* exp(j*k1.*(0:1:3))

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