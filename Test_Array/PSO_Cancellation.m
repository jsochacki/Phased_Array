clear all
clc

%In this case we want all tx to rx values to be 0
%Since the tx ports are 1-5 and the rx ports are 6-10
%we want to make sure that S(6,7,8,9,10)1, S(6,7,8,9,10)2, etc....
%are all 0.  Since the matrix has input ports vs columns and output ports
%vs rows as below

%s11 s12 s13 s14 s15 s16 s17 s18 s19 s110
%s21 s22 s23 s24 s25 s26 s27 s28 s29 s210
%s31 s32 s33 s34 s35 s36 s37 s38 s39 s310
%s41 s42 s43 s44 s45 s46 s47 s48 s49 s410
%s51 s52 s53 s54 s55 s56 s57 s58 s59 s510
%s61 s62 s63 s64 s65 s66 s67 s68 s69 s610
%s71 s72 s73 s74 s75 s76 s77 s78 s79 s710
%s81 s82 s83 s84 s85 s86 s87 s88 s89 s810
%s91 s92 s93 s94 s95 s96 s97 s98 s99 s910
%s101 s102 s103 s104 s105 s106 s107 s108 s109 s1010

%The inputs I have are to the weights on the upper left diagonal only
%So woud defactor set x6-x10 = 0
% x = [x1 x2 x3 x4 x5 0 0 0 0 0].';
%The output vector y results from s * x and since x is tx only that also
%further masks off our s matrix effectively to the following

%s11 s12 s13 s14 s15 0   0   0   0   0
%s21 s22 s23 s24 s25 0   0   0   0   0
%s31 s32 s33 s34 s35 0   0   0   0   0
%s41 s42 s43 s44 s45 0   0   0   0   0
%s51 s52 s53 s54 s55 0   0   0   0   0
%s61 s62 s63 s64 s65 0   0   0   0   0
%s71 s72 s73 s74 s75 0   0   0   0   0
%s81 s82 s83 s84 s85 0   0   0   0   0
%s91 s92 s93 s94 s95 0   0   0   0   0
%s101 s102 s103 s104 s105 0   0   0   0   0

%So effectively all we are reallly solving for is the left half of the
%sparameter matrix so we can truncate the square matrix to a rectangular as
%such

%| s11 s12 s13 s14 s15 |    *   | x1 | = |  y1  |
%| s21 s22 s23 s24 s25 |        | x2 |   |  y2  |
%| s31 s32 s33 s34 s35 |        | x3 |   |  y3  |
%| s41 s42 s43 s44 s45 |        | x4 |   |  y4  |
%| s51 s52 s53 s54 s55 |        | x5 |   |  y5  |
%| s61 s62 s63 s64 s65 |                 |  y6  |
%| s71 s72 s73 s74 s75 |                 |  y7  |
%| s81 s82 s83 s84 s85 |                 |  y8  |
%| s91 s92 s93 s94 s95 |                 |  y9  |
%| s101 s102 s103 s104 s105 |            |  y10 |

%The goal is to minize the y6 - y10 terms and we don't care about the
%y1 - y5 terms so again we can effectively reduce the matrix to the
%following

%| s61 s62 s63 s64 s65 |    *   | x1 | = |  y6  |
%| s71 s72 s73 s74 s75 |        | x2 |   |  y7  |
%| s81 s82 s83 s84 s85 |        | x3 |   |  y8  |
%| s91 s92 s93 s94 s95 |        | x4 |   |  y9  |
%| s101 s102 s103 s104 s105 |   | x5 |   |  y10  |

%And lastly we can say that we want to just minimize all values equally so
%we will say that min( y*y' ) is the fitness function that we want to
%evaluate/minimize which is just to total power from the tramsmit array
%into the receive array

%Also note that the manifold file has x different cell for each frequency
%Once in a manifold you have a single element per column and rows are vs
%phi

[FileName, sp, freq_GHz] = ReadTouchstone();
sparameter_freq_vector = freq_GHz * 1e9;

Near_Field_Struct = load('./Near_Field_Data/5ELx2 1meter NearField Tx Only.mat');
Far_Field_Struct = load('./Far_Field_Data/5ELx2 FarField Tx Only.mat');

beam_theta = 90;
null_phis = [0];
transmit_phis = [-180:5:15 15:5:175];
Max_SINR = 70;

freq_range = [5e8 6e9];     % [start_frequency stop_frequency] [Hz]
freq_step = 50e6;           % Frequency Step size [Hz]
weights = ones(1, 5);       % Amplitude weights for N elements from 0 to 1  

Starting_Element_Number = 6;
Ending_Element_Number = 10;
S = sp(Starting_Element_Number:1:Ending_Element_Number, 1:1:length(weights), :);

freq_points = ((freq_range(2) - freq_range(1)) / freq_step) + 1; 
frequency_vector = freq_range(1):freq_step:freq_range(2); 
phi_vector_degrees = -180:1:179; 
phi_vector_rads = phi_vector_degrees * (pi / 180); 

N = zeros(length(frequency_vector), length(phi_vector_degrees));
F = zeros(length(frequency_vector), length(phi_vector_degrees));
for n = 1:1:length(frequency_vector)
   frequency = frequency_vector(n);
   
   near_field_findex = find(Near_Field_Struct.flist_MHz == (frequency/1e6));
   near_field_vs_phi = Near_Field_Struct.manifold{near_field_findex};
   temp = near_field_vs_phi(ismember(Near_Field_Struct.phi, phi_vector_degrees), :).*weights;
   AF = sum(temp, 2);
   N(n, :) = 10*log10(AF.*conj(AF));
   
   far_field_findex = find(Far_Field_Struct.flist_MHz == (frequency/1e6));
   far_field_vs_phi = Far_Field_Struct.manifold{far_field_findex};
   temp = far_field_vs_phi(ismember(Far_Field_Struct.phi, phi_vector_degrees), :).*weights;
   AF = sum(temp, 2);
   F(n, :) = 10*log10(AF.*conj(AF));
end

analysis_and_plot_frequency_index = find(frequency_vector == 1e9);
analysis_and_plot_frequency = frequency_vector(analysis_and_plot_frequency_index);

figure(2);
plot(phi_vector_degrees, N(analysis_and_plot_frequency_index,:), phi_vector_degrees, F(analysis_and_plot_frequency_index,:));
figure(3)
polarplot(phi_vector_degrees*pi/180, N(analysis_and_plot_frequency_index,:));
figure(4)
polarplot(phi_vector_degrees*pi/180, F(analysis_and_plot_frequency_index,:));

FF_Max = max(F(analysis_and_plot_frequency_index,:));


%Now do PSO algorithm
particles_per_value = 64;

% PSO paramters
Max_iteration = 200;
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
            frequency = analysis_and_plot_frequency;

            near_field_findex = find(Near_Field_Struct.flist_MHz == (frequency/1e6));
            near_field_vs_phi = Near_Field_Struct.manifold{near_field_findex};
            temp = near_field_vs_phi(ismember(Near_Field_Struct.phi, null_phis), :).*weights;
            AF = sum(temp, 2);
            NF = 10*log10(sum(AF.*conj(AF))/length(null_phis));

            far_field_findex = find(Far_Field_Struct.flist_MHz == (frequency/1e6));
            far_field_vs_phi = Far_Field_Struct.manifold{far_field_findex};
            temp = far_field_vs_phi(ismember(Far_Field_Struct.phi, transmit_phis), :).*weights;
            AF = sum(temp, 2);
            FF = 10*log10(sum(AF.*conj(AF))/length(transmit_phis));

            sparameter_frequency_index = find(sparameter_freq_vector == frequency);

            temp = S(:, :, sparameter_frequency_index) * weights.';
            MSP = 10*log10((temp'*temp) / length(temp));

            Power_Loss = 10*log10(length(weights) / (weights*weights'));

            % Use Near Fields (Modified for SINR)
            % Just need to limit the near field null effect so it doesn't
            % dominate the cost function
            MNF = ((FF - Max_SINR) * (NF < (FF - Max_SINR))) + (NF * ~(NF < (FF - Max_SINR)));
            MPL = (Power_Loss * (3 > Power_Loss)) + (Power_Loss*10 * ~(3 > Power_Loss));
            % Dont use Sparameter results, they skew the cost function and
            % have no bearing on nulling performance.  They do reflect good
            % nulling once weights are computed but do not do well to adapt
            % for good weights
            %objective_function_value = MSP + MNF + MPL + (FF_Max - FF);
            objective_function_value = MNF + MPL + (FF_Max - FF);
            % Use S parameters and Near Fields
            %objective_function_value = MSP + NF + Signal_Power + (FF_Max - FF);
            % Use Near Fields only
            %objective_function_value = NF + Signal_Power;
            % Use Far Fields only (Creats null in far filed but it shifts
            % in the near filed and fills in so this will not work)
            %objective_function_value = Signal_Power + FF;
            % Use S parameters only (Does not create null like we would
            % hope so we can ignore these in optimization algorithm)
            %objective_function_value = MSP + Signal_Power;


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

weights.*weights'.' %Magnitudes
angle(weights)*180/pi %Phases
y = S(:, :, sparameter_frequency_index) * weights.';
10*log10(y.*y'.') %S parameters
10*log10(weights*weights' / length(weights)) %power loss

N = zeros(length(frequency_vector), length(phi_vector_degrees));
F = zeros(length(frequency_vector), length(phi_vector_degrees));
for n = 1:1:length(frequency_vector)
   frequency = frequency_vector(n);
   
   near_field_findex = find(Near_Field_Struct.flist_MHz == (frequency/1e6));
   near_field_vs_phi = Near_Field_Struct.manifold{near_field_findex};
   temp = near_field_vs_phi(ismember(Near_Field_Struct.phi, phi_vector_degrees), :).*weights;
   AF = sum(temp, 2);
   N(n, :) = 10*log10(AF.*conj(AF));
   
   far_field_findex = find(Far_Field_Struct.flist_MHz == (frequency/1e6));
   far_field_vs_phi = Far_Field_Struct.manifold{far_field_findex};
   temp = far_field_vs_phi(ismember(Far_Field_Struct.phi, phi_vector_degrees), :).*weights;
   AF = sum(temp, 2);
   F(n, :) = 10*log10(AF.*conj(AF));
end

figure(2);
hold on
plot(phi_vector_degrees, N(analysis_and_plot_frequency_index,:), phi_vector_degrees, F(analysis_and_plot_frequency_index,:));
figure(3)
hold on
polarplot(phi_vector_degrees*pi/180, N(analysis_and_plot_frequency_index,:));
ax = gca;
ax.RLim = [0 FF_Max+2];
figure(4)
hold on
polarplot(phi_vector_degrees*pi/180, F(analysis_and_plot_frequency_index,:));
polarplot(phi_vector_degrees*pi/180, N(analysis_and_plot_frequency_index,:));
ax = gca;
ax.RLim = [0 FF_Max+2];