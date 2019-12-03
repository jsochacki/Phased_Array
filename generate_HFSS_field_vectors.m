function [phi_theta_vector, all_theta, all_phi] = generate_HFSS_field_vectors(phi_theta, phi, theta)

all_theta = repelem(theta, length(phi));
all_phi = repelem(phi, length(theta));
index = 1;
for n = 0:1:(length(theta) - 1)
   for nn = 0:1:(length(phi) - 1)
      phi_theta_vector(index) = phi_theta(n + 1, nn + 1);
      index = index + 1;
   end
end

end