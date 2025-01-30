clear all;
close all;
clc;

% Ask user if they want to use default values or enter custom values
choice = input('Do you want to use default values? (1 for Yes, 0 for No): ');

if choice == 1
    % Default values
    P1 = 1;          % Initial pressure (bar)
    T1 = 300;        % Initial temperature (K)
    r = 16;          % Compression ratio
    alpha = 2;       % Pressure ratio
    beta = 1.5;      % Cut-off ratio
    gamma = 1.4;     % Specific heat ratio
else
    % Get user input
    P1 = input('Enter initial pressure (bar): ');
    T1 = input('Enter initial temperature (K): ');
    r = input('Enter compression ratio: ');
    alpha = input('Enter pressure ratio: ');
    beta = input('Enter cut-off ratio: ');
    gamma = input('Enter specific heat ratio: ');
end

% Constants
R = 0.287;          % Gas constant (kJ/kg-K)

% Process calculations
V1 = R*T1/P1;       % Initial volume
V2 = V1/r;          % Volume after compression
V3 = V2;            % Volume at constant volume heat addition
V4 = V3*beta;       % Volume after constant pressure heat addition
V5 = V1;            % Final volume

% Pressure calculations
P2 = P1*r^gamma;    % Pressure after compression
P3 = P2*alpha;      % Pressure after constant volume heat addition
P4 = P3;            % Pressure during constant pressure heat addition
P5 = P4*(V4/V5)^gamma; % Final pressure

% Temperature calculations
T2 = T1*r^(gamma-1);
T3 = T2*alpha;
T4 = T3*beta;
T5 = T4*(V4/V5)^(gamma-1);

% Generate points for plotting
n_points = 100;

% Process 1-2 (Isentropic compression)
V12 = linspace(V1, V2, n_points);
P12 = P1*(V1./V12).^gamma;

% Process 2-3 (Constant volume heat addition)
V23 = linspace(V2, V3, n_points);
P23 = linspace(P2, P3, n_points);

% Process 3-4 (Constant pressure heat addition)
V34 = linspace(V3, V4, n_points);
P34 = ones(1, n_points)*P3;

% Process 4-5 (Isentropic expansion)
V45 = linspace(V4, V5, n_points);
P45 = P4*(V4./V45).^gamma;

% Plot PV diagram
figure(1)
plot(V12, P12, 'b', 'LineWidth', 2)
hold on
plot(V23, P23, 'r', 'LineWidth', 2)
plot(V34, P34, 'g', 'LineWidth', 2)
plot(V45, P45, 'm', 'LineWidth', 2)
grid on
xlabel('Volume (m^3/kg)')
ylabel('Pressure (bar)')
title('P-V Diagram of Dual Cycle')
legend('1-2 Compression', '2-3 Constant Volume', '3-4 Constant Pressure', '4-5 Expansion')

% Calculate efficiency for different compression ratios
r_range = 4:0.5:24;
efficiency = zeros(size(r_range));

for i = 1:length(r_range)
    r_current = r_range(i);
    
    % Calculate efficiency
    term1 = 1/r_current^(gamma-1);
    term2 = (alpha-1)/(gamma-1);
    term3 = alpha*beta*(1-1/r_current^(gamma-1));
    term4 = gamma*alpha*(beta-1);
    
    efficiency(i) = 1 - (term1*(1 + term2 + term3 - term4))/(alpha*(beta-1) + (alpha-1)/(gamma-1));
end

% Plot efficiency curve
figure(2)
plot(r_range, efficiency*100, 'b', 'LineWidth', 2)
grid on
xlabel('Compression Ratio')
ylabel('Thermal Efficiency (%)')
title('Thermal Efficiency vs Compression Ratio')

% Display efficiency for the chosen compression ratio
if choice == 1
    fprintf('\nFor default compression ratio r = %.1f:\n', r);
else
    fprintf('\nFor chosen compression ratio r = %.1f:\n', r);
end
current_efficiency = interp1(r_range, efficiency, r)*100;
fprintf('Thermal Efficiency = %.2f%%\n', current_efficiency);