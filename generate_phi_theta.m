function [phi_theta, phi_theta_dB, phi, theta] = generate_phi_theta(frequency, d_a, weights, array_type, element_phi_theta)

phi = 0:1:359;
theta = 0:1:179;
phi_theta = zeros(length(theta), length(phi));

switch lower(array_type)
   case 'linear'
      for n = 1:1:180
         for nn = 1:1:360
            [phi_theta(n, nn), phi_theta_dB(n, nn), ~] = ...
               Uniform_Linear_Array(n * (pi / 180), frequency, d_a, weights);
         end
      end
   case 'planar'
      for n = 1:1:180
         for nn = 1:1:360
            [phi_theta(n, nn), phi_theta_dB(n, nn), ~] = ...
               Uniform_Planar_Array(nn * (pi / 180), n * (pi / 180), frequency, d_a, weights);
         end
      end
   case 'circular'
      for n = 1:1:180
         for nn = 1:1:360
            [phi_theta(n, nn), phi_theta_dB(n, nn), ~] = ...
               Uniform_Circular_Array(nn * (pi / 180), n * (pi / 180), frequency, d_a, weights);
         end
      end
   otherwise
      phi_theta=[]; phi=[]; theta=[];
      error('The array geometry %s is not currently supported.', array_type);
end
      
switch nargin
   case 4
   case 5
      phi_theta = phi_theta .* element_phi_theta;
      phi_theta_dB = phi_theta_dB + 10*log10(element_phi_theta .* element_phi_theta'.')
   otherwise
      phi_theta=[]; phi=[]; theta=[];
      error('Not enough input arguements provided');
end
      
end