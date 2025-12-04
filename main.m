% Synchronization of a dynamical network of Lorenz oscillators
% ------------------------------------------------------------
% Purpose:
%   Simulate a network of N identical Lorenz oscillators coupled via a directed graph.
%   Coupling strengths are assigned via SyncCouplingAssign to guarantee synchronization.
%   After simulation, compute pairwise synchronization error (state-space L2 distance)
%   between each oscillator and its "next" in a ring (1–2, 2–3, …, N–1–N, N–1).
%   Plot the error vs time with distinct line/marker styles for each pair.
%
% Note:
%   P (projection matrix) determines which state variables (e.g. x-component) are coupled.
%   Here P = diag([1 0 0]) → only the first coordinate is used in coupling.

clc; clear; close all;

%% === System parameters (Lorenz oscillator) ===
sigma = 10;
rho   = 25;
beta  = 8/3;

% Analytical coupling strength derived (sufficient for synchronization)
a = -sigma + (beta*(beta+1)*(rho+sigma)^2) / (16*(beta-1));
% a = 1;  % (alternate manual coupling strength)

%% === Define connectivity digraph ===
% Define directed edges via tail → head lists
tail1 = [1 2 3 4 4 5 4 5 6];
head1 = [2 3 1 1 2 1 5 6 4];
G = digraph(tail1, head1);

N = G.numnodes;       % Number of oscillators
numStates = 3;        % Dimension of each oscillator's state

%% === Assign coupling strengths ===
G = SyncCouplingAssign(G, a);
figure;
plot(G, 'EdgeLabel', G.Edges.Weight, 'Layout', 'circle');
title('Coupling digraph with assigned weights');

%% === Projection matrix for coupling ===
P = diag([1, 0, 0]);   % couple only first coordinate (e.g. x-variable)

%% === Simulation settings ===
data_length = 10;
t_end = 2;
tspan = linspace(0, t_end, data_length);

%% === Initial conditions ===
x_mean = 0;
x_std  = 12;
X0 = x_mean + x_std * rand(1, numStates * N);

%% === Simulate coupled Lorenz oscillators ===
[X, t] = SimulateCoupledSystems(@LorenzOscillator, tspan, X0, G, P);

%% === Compute synchronization errors (pairwise in a ring) ===
syncError = zeros(N, length(t));

figure; hold on; grid on;
colors = lines(N);

% Define arrays of marker types and line-styles to cycle through
markerList    = {'o', '+', '*', 's', 'd', 'v', '^', '>', '<', 'p', 'h'};
lineStyleList = {'-', '--', ':', '-.'};

for i = 1:N
    % Determine "next" oscillator index in ring (with wrap-around)
    next = i + 1;
    if next > N
        next = 1;
    end

    % Indices of state components in X for oscillator i and for next oscillator:
    idx_i    = (i-1)*numStates + (1:numStates);
    idx_next = (next-1)*numStates + (1:numStates);

    % Extract state trajectories:
    Xi    = X(:, idx_i);
    Xnext = X(:, idx_next);

    % Difference in full state space at each time:
    D = Xi - Xnext;            % size: [ length(t) × numStates ]

    % Compute Euclidean norm (L2) of difference at each time:
    e_i = vecnorm(D, 2, 2);    % returns [length(t) × 1] vector  :contentReference[oaicite:0]{index=0}

    syncError(i, :) = e_i;

    % Choose style for this line
    m  = markerList{ mod(i-1, numel(markerList)) + 1 };
    ls = lineStyleList{ mod(i-1, numel(lineStyleList)) + 1 };
    c  = colors(i, :);

    plot(t, e_i, ...
         'LineStyle', ls, ...
         'Marker',   m, ...
         'Color',    c, ...
         'LineWidth', 2, ...
         'MarkerSize', 6, ...
         'DisplayName', sprintf('Osc %d vs %d', i, next) ...
    );
end

xlabel('Time');
ylabel('Synchronization error (L2 norm in state space)');
title(sprintf('Synchronization error: ring-coupled %d Lorenz oscillators', N));
legend('Location','bestoutside');
hold off;

%% === Save synchronization error data to Excel ===
filename = 'synchronization_data.xlsx';
cols = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ';
range_end = length(t) + 1;

% First column: time vector
writematrix(t.', filename, 'Sheet', 1, 'Range', strcat(cols(1), '2:', cols(1), string(range_end)));
writematrix('t', filename, 'Sheet', 1, 'Range', strcat(cols(1), '1'));

for i = 1:N
    % Write error data for each oscillator-pair
    col_letter = cols(i+1);
    range_str = strcat(col_letter, '2:', col_letter, string(range_end));
    writematrix(syncError(i,:).', filename, 'Sheet', 1, 'Range', range_str);
    writematrix( sprintf('e_%d', i), filename, 'Sheet', 1, 'Range', strcat(col_letter, '1') );
end

disp(['Time and synchronization-error data saved to ', filename]);
