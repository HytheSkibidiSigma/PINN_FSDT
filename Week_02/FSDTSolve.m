clear; clc; close all;

a = 5; b = 5; 
Lx = 2*a; Ly = 2*b;   
h = 1;               
n_layer = 2;          
t = h / n_layer;      
layup = [0, 90];     
q0 = 1e-3;           

E1 = 25;
E2 = E1 / 25;
G12 = 0.5 * E2;
G13 = G12;
G23 = 0.2 * E2;
mu12 = 0.25;
mu21 = mu12 * E2 / E1;
Ks = 5/6;             

Q11 = E1 / (1 - mu12 * mu21);
Q12 = mu12 * E2 / (1 - mu12 * mu21);
Q22 = E2 / (1 - mu12 * mu21);
Q66 = G12;
Q = [Q11, Q12, 0; Q12, Q22, 0; 0, 0, Q66];

A = zeros(3,3); B = zeros(3,3); D = zeros(3,3);
A44 = 0; A55 = 0;

for i = 1:n_layer
    phi = layup(i);
    m = cosd(phi); n = sind(phi);
    
    T = [m^2, n^2, m*n; n^2, m^2, -m*n; -2*m*n, 2*m*n, m^2-n^2];
    Q_bar = T' * Q * T;
    
    z1 = (i - 1 - n_layer/2) * t;
    z2 = (i - n_layer/2) * t;
    
    A = A + Q_bar * (z2 - z1);
    B = B + Q_bar * (z2^2 - z1^2) / 2;
    D = D + Q_bar * (z2^3 - z1^3) / 3;
    
    if phi == 0
        Q44_bar = G23; Q55_bar = G13;
    else
        Q44_bar = G13; Q55_bar = G23;
    end
    A44 = A44 + Q44_bar * (z2 - z1);
    A55 = A55 + Q55_bar * (z2 - z1);
end

m_terms = 100;
n_terms = 100; 

[X_grid, Y_grid] = meshgrid(linspace(0, Lx, 50), linspace(0, Ly, 50));
W_domain = zeros(50, 50);
w_center = 0; 

for m = 1:2:m_terms
    for n = 1:2:n_terms
        alpha_m = m * pi / Lx;
        beta_n = n * pi / Ly;
        K = zeros(5,5);
        K(1,1) = A(1,1)*alpha_m^2 + A(3,3)*beta_n^2; 
        K(1,2) = (A(1,2) + A(3,3))*alpha_m*beta_n;
        K(1,4) = B(1,1)*alpha_m^2 + B(3,3)*beta_n^2; 
        K(1,5) = (B(1,2) + B(3,3))*alpha_m*beta_n;
        K(2,1) = K(1,2); 
        K(2,2) = A(3,3)*alpha_m^2 + A(2,2)*beta_n^2;
        K(2,4) = K(1,5); 
        K(2,5) = B(3,3)*alpha_m^2 + B(2,2)*beta_n^2;
        K(3,3) = Ks*A55*alpha_m^2 + Ks*A44*beta_n^2; 
        K(3,4) = Ks*A55*alpha_m; 
        K(3,5) = Ks*A44*beta_n;
        K(4,1) = K(1,4); K(4,2) = K(2,4); K(4,3) = K(3,4);
        K(4,4) = D(1,1)*alpha_m^2 + D(3,3)*beta_n^2 + Ks*A55; 
        K(4,5) = (D(1,2) + D(3,3))*alpha_m*beta_n;
        K(5,1) = K(1,5); K(5,2) = K(2,5); K(5,3) = K(3,5); K(5,4) = K(4,5);
        K(5,5) = D(3,3)*alpha_m^2 + D(2,2)*beta_n^2 + Ks*A44;
        
        F = [0; 0; (16*q0) / (pi^2*m*n); 0; 0];

        d = K \ F;
        W_mn = d(3);
        
        W_domain = W_domain + W_mn * sin(alpha_m * X_grid) .* sin(beta_n * Y_grid);
        w_center = w_center + W_mn * sin(alpha_m * Lx/2) * sin(beta_n * Ly/2);
    end
end

W_bar_matrix = W_domain * 100 * (h^3) * E2 / (q0 * Lx^4);
w_bar_center = w_center * 100 * (h^3) * E2 / (q0 * Lx^4);

fprintf('Độ võng lớn nhất: %.4f\n', w_bar_center);

figure('Color', 'w', 'Position', [100, 100, 800, 600]);
surf(X_grid - a, Y_grid - b, W_bar_matrix, 'EdgeColor', 'none');
colormap(jet); 
cb = colorbar;
cb.Label.String = 'Độ võng';
cb.Label.FontSize = 12;

title('Mặt cong biến dạng FSDT', 'FontSize', 14);
xlabel('Trục X (mm)', 'FontSize', 12);
ylabel('Trục Y (mm)', 'FontSize', 12);
zlabel('Độ võng W_{bar}', 'FontSize', 12);
axis tight; grid on; view(-35, 45);