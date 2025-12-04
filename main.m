% Synchronization of a dynamical network of Lorenz oscillators
% ------------------------------------------------------------
% This script simulates the synchronization of 6 Lorenz oscillators 
% Synchronization error data is saved in "sync_data.xlsx".

clc; clear; close all;

%% Lorenz oscillator parameters
sigma = 10;
rho   = 25;
beta  = 8/3;

% Coupling strength (analytically derived)
a = -sigma + (beta*(beta+1)*(rho+sigma)^2) / (16*(beta-1));
% a=1;

%% Define connectivity digraphs
% Connectivity for first time interval

tail1 = [1 2 3 4 4 5 4 5 6];
head1 = [2 3 1 1 2 1 5 6 4];
G1 = digraph(tail1, head1);

N         = G1.numnodes;  % Number of oscillators (nodes)
numStates = 3;   % State dimension of each oscillator

%% Simulation settings
data_length = 10;
t_end = 2;
tspan = linspace(0, t_end, data_length);

% Assign coupling strengths
G1 = SyncCouplingAssign(G1,a);
figure
plot(G1,'EdgeLabel',G1.Edges.Weight,'Layout','circle')

%% Initial conditions
x_mean = 0; 
x_std  = 12;
P      = diag([1, 0, 0]);
X0     = x_mean + x_std*rand(1, numStates*N);

%% Simulate coupled Lorenz systems
[X, t] = SimulateCoupledSystems(@LorenzOscillator, tspan, X0, G1, P);
X0 = X(end,:);

%% Visualization settings
state_indices   = 1:numStates;
state_index_all = 1:numStates*N;

colors    = lines(N);
linestyle = {'-','--','-.',':','-','--','-.',':','-','--'};
lw        = [2*ones(1,4) 2*ones(1,4) 2 2];
markers   = {'none','none','none','none','*','*','o','o','.','.'};

%% Data storage setup
E         = zeros(N, length(t));   % Synchronization error matrix
cols      = 'ABCDEFG';         % Excel column labels
filename  = 'synchronization_data.xlsx';
range_end = length(t) + 1;

%% Plot synchronization errors
figure; hold on; grid on;

X_circ = [X(:,4:end) X(:,1:3)];
X_diff = X-X_circ;

for i = 1:N
    
    % Extract state indices for oscillator i
    slice_i   = (i-1)*numStates + state_indices;
    % Synchronization error (L2 distance from others)
    e = vecnorm(X_diff(:,slice_i), 2, 2).';
    E(i,:) = e;

    % Plot error trajectory
    plot(t, e, ...
        'Color', colors(i,:), ...
        'LineWidth', lw(i), ...
        'LineStyle', linestyle{i}, ...
        'Marker', markers{i}, ...
        'MarkerFaceColor', 'none', ...
        'DisplayName', sprintf('Systems %d and %d', i,mod(i+1,6)));

    % Save to Excel
    range_str = strcat(cols(i+1),'2:',cols(i+1),string(range_end));
    writematrix(e.', filename,'Sheet',1,'Range',range_str);
    writematrix("e"+i, filename,'Sheet',1,'Range',cols(i+1)+"1");
end

xlabel('Time');
ylabel('Synchronization error (L2 norm)');
title('Synchronization of 10 Lorenz Oscillators');
legend show;
hold off;

%% Save time vector to Excel
writematrix(t, filename, 'Sheet', 1, 'Range', strcat(cols(1),'2:',cols(1),string(range_end)));
writematrix('t', filename, 'Sheet', 1, 'Range', 'A1');
disp(strcat(cols(1),'2:',cols(1),string(range_end)))
