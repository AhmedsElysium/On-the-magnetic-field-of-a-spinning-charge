radius = 1;                 
charge_density = 1;         
omega = 1;                  

[X, Y, Z] = meshgrid(-2:0.2:2, -2:0.2:2, -2:0.2:2);

% Magnetic field calculation using Biot-Savart law
mu_0 = 4 * pi * 1e-7;  
B = zeros(size(X, 1), size(X, 2), size(X, 3), 3);
for theta = 0:pi/36:2*pi
    for phi = 0:pi/18:pi
        r_prime = [radius * sin(phi) * cos(theta), radius * sin(phi) * sin(theta), radius * cos(phi)];
        
        dl = radius * sin(phi) * (omega * [-sin(theta), cos(theta), 0]);
        
        r_obs = cat(4, X - r_prime(1), Y - r_prime(2), Z - r_prime(3));
        
        r = sqrt(sum(r_obs.^2, 4));
        
        cross_prod = zeros(size(X, 1), size(X, 2), size(X, 3), 3);
        for i = 1:3
            cross_prod(:,:,:,i) = dl(2) * r_obs(:,:,:,3) - dl(3) * r_obs(:,:,:,2);
        end
        cross_prod(:,:,:,2) = dl(3) * r_obs(:,:,:,1) - dl(1) * r_obs(:,:,:,3);
        cross_prod(:,:,:,3) = dl(1) * r_obs(:,:,:,2) - dl(2) * r_obs(:,:,:,1);
        
        dB = (mu_0 / (4 * pi)) * (cross_prod ./ (r.^3));
        
        B = B + dB;
    end
end

theta_sphere = 0:pi/36:2*pi;
phi_sphere = 0:pi/18:pi;
[theta_grid, phi_grid] = meshgrid(theta_sphere, phi_sphere);

% Add sphere to represent the shell
[x_sphere, y_sphere, z_sphere] = sphere(500);  % Adjust resolution
hold on;
surf(x_sphere, y_sphere, z_sphere, 'FaceAlpha', 1, 'EdgeColor', 'none');


[x_stream, y_stream, z_stream] = meshgrid(-2:0.5:2, -2:0.5:2, -2:0.5:2);

% Plot magnetic field lines
streamline(X, Y, Z, B(:,:,:,1), B(:,:,:,2), B(:,:,:,3), x_stream, y_stream, z_stream);

colormap('white');
set(gca, 'Color', 'black');

title('Magnetic Field Outside Rotating Charged Sphere');
xlabel('X-axis');
ylabel('Y-axis');
zlabel('Z-axis');
axis equal;
view(3);
camlight;
hold off;
