clear all
clc

beam_theta = 0;
Goal_Beam_Power = 24;
Goal_SLL_Value = -30;
beam_width_in_Degrees = 20;  %Must be even

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

mask = [(max(F(end,:)) + Goal_SLL_Value)*ones(1, (181 - (beam_width_in_Degrees + 1)) / 2) max(F(end,:))*ones(1, beam_width_in_Degrees + 1) (max(F(end,:)) + Goal_SLL_Value)*ones(1, (181 - (beam_width_in_Degrees + 1)) / 2)];
slmask = ~(mask == max(F(end,:)));

figure(2); 
plot(theta_vector_degrees, F(end,:));
hold on
plot(theta_vector_degrees, mask,'k');

%Now do PSO algorithm
particles_per_value = 64;

% PSO paramters
Max_iteration = 60;
values = size(weights, 1) * size(weights, 2);
magnitude_max = 1;
velocity_max = 0.4;
inertia_max = 0.9;
inertia_min = 0.2;
c1 = 2; c2 = 2;
cg_curve = zeros(1,Max_iteration, values);

particle_grid_steps = sqrt(particles_per_value);
x = linspace(-1/sqrt(2), 1/sqrt(2), particle_grid_steps);

% Initializations
for value = 1:1:values
   for n = 1:1:particle_grid_steps
      for nn = 1:1:particle_grid_steps
         particles_value(n, nn, value) = x(n) + j * x(nn);
      end
   end
end

particles_velocities = velocity_max * (rand(particle_grid_steps, particle_grid_steps, values) ...
                     + (j * rand(particle_grid_steps, particle_grid_steps, values)));
particles_value_personal_best = rand(particle_grid_steps, particle_grid_steps, values) ...
                              + (j * rand(particle_grid_steps, particle_grid_steps, values));
particle_personal_best_objective_function = 1e3*ones(particle_grid_steps, particle_grid_steps, values);

swarm_global_best_value = ones(1, 1, values);
swarm_global_best_objective_function = 1e3*ones(1, 1, values);

for itteration = 1:1:Max_iteration % main loop
   
   %Calculate objective function for each particle
   for value = 1:1:values
      for n = 1:1:particle_grid_steps
         for nn = 1:1:particle_grid_steps

            %update active weights
            weights = permute(swarm_global_best_value, [3 1 2]).';
            weights(1, value) = particles_value(n, nn, value);
            
            %Compute objective function START
            frequency = frequency_vector(end);

            % create omnidirectional characteristic 
            iPattern = zeros(1, length(theta_vector_degrees)); 

            % Calculate Array Factor
            for nnn = 1:1:length(theta_vector_rads)
               [AF, AF_dB, AV] = Uniform_Linear_Array(theta_vector_rads(nnn), frequency, d, weights);
               %Combine for full characteristics
               F(end, nnn) = AF_dB + iPattern(nnn);
            end

            %Remove -Inf values
            temp = F(end, find(F(end, :).*slmask ~= 0)); ind = find(temp == -Inf); temp(ind) = -100;
            SLL = sum(temp - Goal_SLL_Value) / sum(slmask);
            temp = F(end, find(F(end, :).*(~slmask) ~= 0)); ind = find(temp == -Inf); temp(ind) = -100;
            Beam_Power = sum(Goal_Beam_Power - temp) / sum(~slmask);
            %Compute objective function END            
            
            objective_function_value = SLL + Beam_Power;
            if objective_function_value < particle_personal_best_objective_function(n, nn, value)
               particle_personal_best_objective_function(n, nn, value) = objective_function_value;
               particles_value_personal_best(n, nn, value) = particles_value(n, nn, value);
            end
            if(objective_function_value < swarm_global_best_objective_function(1, 1, value))
               swarm_global_best_objective_function(1, 1, value) = objective_function_value;
               swarm_global_best_value(1, 1, value) = particles_value(n, nn, value);
            end
         end
      end
   end

   %Update the inertia weight
   inertial_weight = inertia_max - (itteration * ( (inertia_max - inertia_min) / Max_iteration));

   %Update the Velocity and Position of particles
   particles_velocities = inertial_weight * particles_velocities + ...  % inertia
      c1 * rand(particle_grid_steps, particle_grid_steps, values) .* (particles_value_personal_best - particles_value) +  ...   % congnitive
      c2 * rand(particle_grid_steps, particle_grid_steps, values) .* (swarm_global_best_value - particles_value);  % social

   index = find(particles_velocities > velocity_max);
   particles_velocities(index) = velocity_max * rand;
   
   index = find(particles_velocities < -velocity_max);
   particles_velocities(index) = -velocity_max * rand;

   particles_value = particles_value + particles_velocities;
   
   magnitude = sqrt(particles_value .* conj(particles_value));
   index = find(magnitude > magnitude_max);
   particles_value(index) = particles_value(index) ./ magnitude(index);

   cg_curve(itteration, 1:1:values) = swarm_global_best_objective_function;
   swarm_global_best_objective_function
end
figure(1)
hold on
for value = 1:1:values
   plot(cg_curve(:,value))
end
xlabel('Iteration')


weights = permute(swarm_global_best_value, [3 1 2]).';

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