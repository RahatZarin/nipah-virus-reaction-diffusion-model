clear; clc; close all;

%% ------------------ PARAMETERS ------------------ %%
params = struct( ...
    'Psi_H', 5285, ...
    'Psi_F', 35, ...
    'Psi_B', 0.411, ...
    'beta1', 1.39897e-10, ...
    'beta2', 3.49641e-5, ...
    'beta3', 4e-6, ...
    'beta4', 2.53342e-4, ...
    'beta5', 1.55344e-5, ...
    'mu_H', 3.91e-5, ...
    'mu_F', 1.5632e-2, ...
    'mu_B', 1.3699e-4, ...
    'delta_H', 4.36999e-2, ...
    'delta_B', 6.22043e-2, ...
    'sigma_H', 2.25266e-2, ...
    'gamma_H', 1.53737e-3 ...
);

%% ------------------ SPATIAL & TEMPORAL SETUP ------------------ %%
L = 20; Tmax = 60; Nx = 100; Nt = 5000;
dx = L / (Nx - 1); dt = Tmax / Nt;
x = linspace(0, L, Nx); t = linspace(0, Tmax, Nt);
[Tgrid, Xgrid] = meshgrid(t, x);

%% ------------------ DIFFUSION COEFFICIENTS ------------------ %%
%dSH = 0.0; dIH = 0.0; dRH = 0.0; dSF = 0.0; dIF = 0.0; dSB = 0.0; dIB = 0.0;
dSH = 0.02; dIH = 0.02; dRH = 0.02; dSF = 0.05; dIF = 0.05; dSB = 0.5;  dIB = 0.5; 


chi = struct( ...
    'SH', dSH * dt / dx^2, ...
    'IH', dIH * dt / dx^2, ...
    'RH', dRH * dt / dx^2, ...
    'SF', dSF * dt / dx^2, ...
    'IF', dIF * dt / dx^2, ...
    'SB', dSB * dt / dx^2, ...
    'IB', dIB * dt / dx^2 ...
);

% ------------------ INITIAL CONDITIONS (Multiple Gaussians) ------------------ %
% ------------------ INITIAL CONDITIONS (Multiple Gaussians - Fixed) ------------------ %
sigma = 2.5;

% Outbreak centers and amplitudes for each infected group
centers = [L/4, 3*L/4];  % Spatial locations of outbreaks

% Amplitudes (total infected at each center)
amplitudes_IH = [25, 22];   % Infected Humans
amplitudes_IF = [20, 19];   % Infected Farmers
amplitudes_IB = [40, 30];   % Infected Bats

% Initialize arrays
S_H(:,1) = 135179953 * ones(Nx,1);
I_H(:,1) = zeros(Nx,1);
R_H(:,1) = 30 * ones(Nx,1);

S_F(:,1) = 2200 * ones(Nx,1);
I_F(:,1) = zeros(Nx,1);

S_B(:,1) = 2930 * ones(Nx,1);
I_B(:,1) = zeros(Nx,1);

% Apply Gaussian sum for Infected compartments
for i = 1:length(centers)
    I_H(:,1) = I_H(:,1) + amplitudes_IH(i) * exp(-((x - centers(i)).^2) / (2 * sigma^2))';
    I_F(:,1) = I_F(:,1) + amplitudes_IF(i) * exp(-((x - centers(i)).^2) / (2 * sigma^2))';
    I_B(:,1) = I_B(:,1) + amplitudes_IB(i) * exp(-((x - centers(i)).^2) / (2 * sigma^2))';
end
figure;
plot(x, I_H(:,1), 'r', x, I_F(:,1), 'm--', x, I_B(:,1), 'k:', 'LineWidth', 2);
legend('I_H', 'I_F', 'I_B');
xlabel('Space'); ylabel('Initial Infection');
title('Multiple Gaussian Initial Conditions'); grid on;


% %% ------------------ INITIAL CONDITIONS (Single Gaussian) ------------------ %%
% sigma = 2.5;
% S_H = zeros(Nx, Nt); I_H = zeros(Nx, Nt); R_H = zeros(Nx, Nt);
% S_F = zeros(Nx, Nt); I_F = zeros(Nx, Nt);
% S_B = zeros(Nx, Nt); I_B = zeros(Nx, Nt);
% 
% % Human compartments
% S_H(:,1) = 135179953 * ones(Nx,1);
% I_H(:,1) = 47 * exp(-((x - L/2).^2) / (2 * sigma^2));
% R_H(:,1) = 30 * ones(size(x));  % Recovered stays uniform
% 
% % Farmer compartments
% S_F(:,1) = 2200 * ones(Nx,1);
% I_F(:,1) = 39 * exp(-((x - L/2).^2) / (2 * sigma^2));
% 
% % Bat compartments
% S_B(:,1) = 2930 * ones(Nx,1);
% I_B(:,1) = 70 * exp(-((x - L/2).^2) / (2 * sigma^2));

% figure('Name','Initial Conditions - Nipah Virus Model','NumberTitle','off', 'Position', [100, 100, 1000, 900]);
% 
% subplot(4,2,1);
% plot(x, S_H(:,1), 'b', 'LineWidth', 2); grid on;
% title('S_H: Susceptible Humans'); xlabel('Space (x)'); ylabel('Population');
% 
% subplot(4,2,2);
% plot(x, I_H(:,1), 'r', 'LineWidth', 2); grid on;
% title('I_H: Infected Humans'); xlabel('Space (x)'); ylabel('Population');
% 
% subplot(4,2,3);
% plot(x, R_H(:,1), 'g', 'LineWidth', 2); grid on;
% title('R_H: Recovered Humans'); xlabel('Space (x)'); ylabel('Population');
% 
% subplot(4,2,4);
% plot(x, S_F(:,1), 'c', 'LineWidth', 2); grid on;
% title('S_F: Susceptible Farmers'); xlabel('Space (x)'); ylabel('Population');
% 
% subplot(4,2,5);
% plot(x, I_F(:,1), 'm', 'LineWidth', 2); grid on;
% title('I_F: Infected Farmers'); xlabel('Space (x)'); ylabel('Population');
% 
% subplot(4,2,6);
% plot(x, S_B(:,1), 'k', 'LineWidth', 2); grid on;
% title('S_B: Susceptible Bats'); xlabel('Space (x)'); ylabel('Population');
% 
% subplot(4,2,7);
% plot(x, I_B(:,1), 'y', 'LineWidth', 2); grid on;
% title('I_B: Infected Bats'); xlabel('Space (x)'); ylabel('Population');
% 
% sgtitle('Initial Spatial Distributions of All Compartments (t = 0)', 'FontSize', 14, 'FontWeight', 'bold');

% figure;
% plot(x, I_H(:,1), 'r', x, I_F(:,1), 'm--', x, I_B(:,1), 'k:', 'LineWidth', 2.5);
% legend('I_H', 'I_F', 'I_B');
% xlabel('Space'); ylabel('Initial Infection');
% title('Single Gaussian Initial Conditions'); grid on;

%% ------------------ MAIN TIME LOOP ------------------ %%
for n = 1:Nt-1
    SH = S_H(:,n); IH = I_H(:,n); RH = R_H(:,n);
    SF = S_F(:,n); IF = I_F(:,n);
    SB = S_B(:,n); IB = I_B(:,n);

    SH_half = SH + dt * (params.Psi_H - params.beta1*SH.*IH - params.beta2*SH.*IF - params.mu_H*SH + params.gamma_H*RH);
    IH_half = IH + dt * (params.beta1*SH.*IH + params.beta2*SH.*IF - (params.mu_H + params.delta_H + params.sigma_H)*IH);
    RH_half = RH + dt * (params.sigma_H*IH - (params.mu_H + params.gamma_H)*RH);

    SF_half = SF + dt * (params.Psi_F - params.beta3*SF.*IF - params.beta4*SF.*IB - params.mu_F*SF);
    IF_half = IF + dt * (params.beta3*SF.*IF + params.beta4*SF.*IB - params.mu_F*IF);

    SB_half = SB + dt * (params.Psi_B - params.beta5*SB.*IB - params.mu_B*SB);
    IB_half = IB + dt * (params.beta5*SB.*IB - (params.mu_B + params.delta_B)*IB);

    S_H(:,n+1) = crank_nicolson_neumann(SH_half, chi.SH, Nx);
    I_H(:,n+1) = crank_nicolson_neumann(IH_half, chi.IH, Nx);
    R_H(:,n+1) = crank_nicolson_neumann(RH_half, chi.RH, Nx);
    S_F(:,n+1) = crank_nicolson_neumann(SF_half, chi.SF, Nx);
    I_F(:,n+1) = crank_nicolson_neumann(IF_half, chi.IF, Nx);
    S_B(:,n+1) = crank_nicolson_neumann(SB_half, chi.SB, Nx);
    I_B(:,n+1) = crank_nicolson_neumann(IB_half, chi.IB, Nx);
end

%% ------------------ PLOT: TIME SERIES AT CENTER ------------------ %%
center = round(Nx/2);

figure;
plot(t, S_H(center,:), 'b', 'DisplayName','Susceptible Human','linewidth', 1.5); hold on;
plot(t, I_H(center,:), 'r', 'DisplayName','Infected Human','linewidth', 1.5);
plot(t, R_H(center,:), 'g', 'DisplayName','Recovered Human','linewidth', 1.5);
plot(t, S_F(center,:), 'c', 'DisplayName','Susceptible Farmer','linewidth', 1.5); 
plot(t, I_F(center,:), 'm', 'DisplayName','Infected Farmer','linewidth', 1.5);
plot(t, S_B(center,:), 'k', 'DisplayName','Susceptible Bat','linewidth', 1.5); 
plot(t, I_B(center,:), 'y', 'DisplayName','Infected Bat','linewidth', 1.5);
xlabel('Time'); ylabel('Population');
%title('Bat Dynamics at Spatial Center');
legend; grid on;

% figure;
% plot(t, S_H(center,:), 'b', 'DisplayName','Susceptible Human'); %hold on;
% figure;
% plot(t, I_H(center,:), 'r', 'DisplayName','Infected Human');
% plot(t, R_H(center,:), 'g', 'DisplayName','Recovered Human');
% xlabel('Time'); ylabel('Population');
% title('Human Dynamics at Spatial Center'); legend; grid on;
% 
% figure;
% plot(t, S_F(center,:), 'c', 'DisplayName','Susceptible Farmer'); %hold on;
% plot(t, I_F(center,:), 'm', 'DisplayName','Infected Farmer');
% xlabel('Time'); ylabel('Population');
% title('Farmer Dynamics at Spatial Center'); legend; grid on;
% 
% figure;
% plot(t, S_B(center,:), 'k', 'DisplayName','Susceptible Bat'); %hold on;
% plot(t, I_B(center,:), 'y', 'DisplayName','Infected Bat');
% xlabel('Time'); ylabel('Population');
% title('Bat Dynamics at Spatial Center'); legend; grid on;


% %% ------------------ SURFACE PLOTS ------------------ %%
%% ------------------ SURFACE PLOTS (Enhanced) ------------------ %%
titles = {'S_H (Susceptible Humans)', 'I_H (Infected Humans)', 'R_H (Recovered Humans)', ...
          'S_F (Susceptible Farmers)', 'I_F (Infected Farmers)', ...
          'S_B (Susceptible Bats)', 'I_B (Infected Bats)'};

data_list = {S_H, I_H, R_H, S_F, I_F, S_B, I_B};
zlabels = {'S_H(t,X)', 'I_H(t,X)', 'R_H(t,X)', ...
           'S_F(t,X)', 'I_F(t,X)', 'S_B(t,X)', 'I_B(t,X)'};

colormaps = {'parula', 'hot', 'summer', 'cool', 'spring', 'autumn', 'winter'};

for i = 1:length(data_list)
    figure('Name', titles{i}, 'NumberTitle', 'off');
    surf(Tgrid, Xgrid, data_list{i}, 'EdgeColor', 'none');
    xlabel('Time', 'FontSize', 12); ylabel('Space', 'FontSize', 12); zlabel(zlabels{i}, 'FontSize', 12);
    %title(['Spatial-Temporal Plot of ', titles{i}], 'FontSize', 14, 'FontWeight', 'bold');
    title([titles{i}], 'FontSize', 14, 'FontWeight', 'bold');
    colormap(colormaps{i}); colorbar;
    shading interp;
    view([135 30]);  % Change to [az el] to rotate
    lighting gouraud;
end

% figure; surf(Tgrid, Xgrid, S_H, 'EdgeColor', 'none');
% xlabel('Time'); ylabel('Space'); zlabel('S_H(t,X)'); title('Susceptible Humans'); colorbar;
% 
% figure; surf(Tgrid, Xgrid, I_H, 'EdgeColor', 'none');
% xlabel('Time'); ylabel('Space'); zlabel('I_H(t,X)'); title('Infected Humans'); colorbar;
% 
% figure; surf(Tgrid, Xgrid, R_H, 'EdgeColor', 'none');
% xlabel('Time'); ylabel('Space'); zlabel('R_H(t,X)'); title('Recovered Humans'); colorbar;
% 
% figure; surf(Tgrid, Xgrid, S_F, 'EdgeColor', 'none');
% xlabel('Time'); ylabel('Space'); zlabel('S_F(t,X)'); title('Susceptible Farmers'); colorbar;
% 
% figure; surf(Tgrid, Xgrid, I_F, 'EdgeColor', 'none');
% xlabel('Time'); ylabel('Space'); zlabel('I_F(t,X)'); title('Infected Farmers'); colorbar;
% 
% figure; surf(Tgrid, Xgrid, S_B, 'EdgeColor', 'none');
% xlabel('Time'); ylabel('Space'); zlabel('S_B(t,X)'); title('Susceptible Bats'); colorbar;
% 
% figure; surf(Tgrid, Xgrid, I_B, 'EdgeColor', 'none');
% xlabel('Time'); ylabel('Space'); zlabel('I_B(t,X)'); title('Infected Bats'); colorbar;

%% ------------------ CRANK-NICOLSON (NEUMANN) ------------------ %%
function u_new = crank_nicolson_neumann(u_half, chi, Nx)
    u_new = zeros(Nx,1);
    main = (1 + chi) * ones(Nx,1);
    upper = (-chi/2) * ones(Nx-1,1);
    lower = (-chi/2) * ones(Nx-1,1);
    A = diag(main) + diag(upper,1) + diag(lower,-1);

    % Neumann boundary adjustments
    A(1,2) = -chi; A(1,1) = 1 + chi;
    A(end,end-1) = -chi; A(end,end) = 1 + chi;

    rhs = u_half;
    for i = 2:Nx-1
        rhs(i) = chi/2*u_half(i-1) + (1 - chi)*u_half(i) + chi/2*u_half(i+1);
    end
    rhs(1)   = (1 - chi)*u_half(1) + chi*u_half(2);      % left Neumann
    rhs(end) = (1 - chi)*u_half(end) + chi*u_half(end-1);% right Neumann

    u_new = A \ rhs;
end
