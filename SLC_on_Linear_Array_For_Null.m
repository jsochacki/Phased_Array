clear all
clc

beam_thetas = [0 50 70];
null_thetas = [35 60 80 90];

freq_range = [24e9 28e9];   % [start_frequency stop_frequency] [Hz]
freq_points = 101;          % No. of frequency points
d = 0.005;                  % Element spacing [m]
weights = ones(1, 16);      % Amplitude weights for N elements from 0 to 1  

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

figure(2);
plot(theta_vector_degrees, F(end,:));

%Calculate weights for sidelobe canceller
for n = 1:1:length(beam_thetas)
   weights_desired_direction(:, :, n) = weight_vs_angle{end}{find(theta_vector_degrees == beam_thetas(n))}{2}';
end
%Can just sum desired direction vectors and normalize for total beam response
weights_desired_direction = sum(weights_desired_direction, 3);
weights_desired_direction = weights_desired_direction ./ abs(weights_desired_direction);

%Need to itteratively subtract out null vectors, cannot do cumulatively
%like desired vectors, however, you can only effectively place one deep null
%with this method which will be the last one that you subtract.  The other
%nulls may or may not be there but if they are there they will not be there
%not like that final one you place

%weights = desired_weights - null_weights * ((null_weights' * desired_weights) / (null_weights' * null_weights));
%What is happening here is from a paper about optimization but it can be
%thought of as the following
%(null_weights' * desired_weights) is the difference vector between null_weights and desired_weights
% i.e. is the angle between the two vectors and also difference in
% magnitude
%(null_weights' * null_weights) is the magnitude of the null_weights
% so null_weights * ((null_weights' * desired_weights) / (null_weights' * null_weights));
% is doing (null_weights / magnitude of null_weights) * difference in
% magnitude and angle between null_weights and desired_weights
% again this is actually derives as the LMS value in the paper
% 5.0  SIDELOBE CANCELLATION 
% but this is one way to think of how I have applied it here vs how they
% apply it / derive its application in the paper
weights = weights_desired_direction;
for n = 1:1:length(null_thetas)
   null_weights = weight_vs_angle{end}{find(theta_vector_degrees == null_thetas(n))}{2}';
   weights = weights - null_weights * ((null_weights' * weights) / (null_weights' * null_weights));
   weights_null_direction(:, :, n) = null_weights;
end

%Apply
for n = 1:1:length(frequency_vector)
   frequency = frequency_vector(n);
   
   % create omnidirectional characteristic 
   iPattern = zeros(1, length(theta_vector_degrees)); 
   
   % Calculate Array Factor
   for nn = 1:1:length(theta_vector_rads)
      [AF, AF_dB, AV] = Uniform_Linear_Array(theta_vector_rads(nn), frequency, d, weights);
      %Combine for full characteristics
      F(n, nn) = AF_dB + iPattern(nn);
   end
end 

figure(2);
hold on
plot(theta_vector_degrees, F(end,:));