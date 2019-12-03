clear all
clc

freq_range = [24e9 28e9];   % [start_frequency stop_frequency] [Hz]
freq_points = 101;          % No. of frequency points
d = 0.005353;                  % Element spacing [m]
N_elements_per_side = 10;
weights = ones(N_elements_per_side, N_elements_per_side);    % Amplitude weights for N elements from 0 to 1  

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

figure(2); 
plot(theta_vector_degrees, F(end,:)); 
hold on

pel = findpeaks(F(end, :), theta_vector_degrees, 'NPeaks', 2, 'SortStr', 'descend');
sllopt = pel(1) - pel(2);

% Array Thining using the Natural Selection Genetic Algorithm
% References
% Randy L. Haupt, Thinned Arrays Using Genetic Algorithms,
%     IEEE Transactions on Antennas and Propagation, Vol 42, No 7, 1994
% Randy L. Haupt, An Introduction to Genetic Algorithms for Electromagnetics,
%     IEEE antennas and Propagation Magazine, Vol 37, No 2, 1995
% Harry L. Van Trees, Optimum Array Processing, Wiley-Interscience, 2002


% Set random number generator for repeatibility
rng_state = rng(2013, 'twister');

% Set active / inactive weight matrix for only 1/4 of the array since it is
% symmetric about both axes using uniform distribution
rows = N_elements_per_side / 2;
columns = N_elements_per_side / 2;
candidates_per_generation = 200;
w0 = double(rand(rows, columns, candidates_per_generation, 'double') > 0.5);
% Set the inner N_elements_per_side / 4 elements to active and the outer
% ones to random activity, this is so that we have a decent generation to
% start with
w0(1:idivide(rows, int32(2)), 1:idivide(columns , int32(2)), :) = 1;

%Single Iteration Example
% Pick one random candidate, plot the beam pattern, and compute the sidelobe
% level
wtemp = w0(:, :, randi(candidates_per_generation));
weights = [fliplr(flipud(wtemp)) flipud(wtemp); fliplr(wtemp) wtemp];

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

figure(2); 
plot(theta_vector_degrees, F(end,:));

pel = findpeaks(F(end, :), theta_vector_degrees, 'NPeaks', 2, 'SortStr', 'descend');
sllopt = pel(1) - pel(2);

fillrate = (sum(weights(:)) / power(N_elements_per_side, 2)) * 100; %fill percent

%Actual implementation use
clear sllopt;
generations = 2;
for iteration = 1:1:generations
   %Calculate the sidelobe level for each candidate from the current
   %generation
   for candidate = 1:1:candidates_per_generation
      % Pick one random candidate, plot the beam pattern, and compute the sidelobe
      % level
      wtemp = w0(:, :, candidate);
      weights = [fliplr(flipud(wtemp)) flipud(wtemp); fliplr(wtemp) wtemp];

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

      pel = findpeaks(F(end, :), theta_vector_degrees, 'NPeaks', 2, 'SortStr', 'descend');
      sllopt(candidate) = pel(1) - pel(2);
   end
   
   %Now sort the performance of all candidates
   [~, idx] = sort(sllopt, 'descend');
   
   %Now remove the half of the current generation that has the worse
   %performance.
   w0 = w0(:, :, idx(1:1:(candidates_per_generation / 2)));
   
   %Finally perform random mutation of candidates that performed best and
   %add create new generation.  All this does is take the best performes
   %and perform a random mutation and add those to the new generation.
   %First it takes all weights from the row mutation start idx onward and does a
   %flipud therby mutating the rows and then takes all rows from the column
   %mutation start idx and does a fliplr therby mutating the columns
   row_mutation_start_idx = randi(rows);
   column_mutation_start_idx = randi(columns);
   w0(row_mutation_start_idx:1:end , :, ((candidates_per_generation / 2) + 1):1:candidates_per_generation) ...
      = flipud(w0(row_mutation_start_idx:1:end , :, 1:1:(candidates_per_generation / 2)));
   w0(column_mutation_start_idx:1:end , :, ((candidates_per_generation / 2) + 1):1:candidates_per_generation) ...
      = fliplr(w0(column_mutation_start_idx:1:end , :, 1:1:(candidates_per_generation / 2)));
end

wtemp = w0(:, :, 1);
optimum_weights = [fliplr(flipud(wtemp)) flipud(wtemp); fliplr(wtemp) wtemp];
fillrate = (sum(optimum_weights(:)) / power(N_elements_per_side, 2)) * 100; %fill percent

rng(rng_state);

for n = 1:1:length(frequency_vector)
   frequency = frequency_vector(n);

   % create omnidirectional characteristic 
   iPattern = zeros(1, length(theta_vector_degrees)); 

   % Calculate Array Factor
   for nn = 1:1:length(theta_vector_rads)
      [AF, AF_dB, AV] = Uniform_Planar_Array(phi_val*(pi / 180), theta_vector_rads(nn), frequency, d, optimum_weights);
      weight_vs_angle{n}{nn} = {theta_vector_rads(nn), AV};      
      %Combine for full characteristics
      F(n, nn) = AF_dB + iPattern(nn);
   end
end 

figure(2); 
plot(theta_vector_degrees, F(end,:));

pel = findpeaks(F(end, :), theta_vector_degrees, 'NPeaks', 2, 'SortStr', 'descend');
sllopt = sllopt(idx(1));

fillrate = (sum(optimum_weights(:)) / power(N_elements_per_side, 2)) * 100; %fill percent