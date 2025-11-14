% Synchronization of a dynamical network of Lorenz oscillators
% ------------------------------------------------------------
% This script simulates the synchronization of 10 Lorenz oscillators 
% over two time intervals with different connectivity digraphs (G1 and G2).
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
% tail1 = [1 2 3 3 4 8 8 7 8 6 7 9 10 5 5 7];
% head1 = [2 3 1 4 1 1 3 3 6 7 8 10 5 9 8 9];

tail1 = [1 2 3 4 4 5 4 5 6];
head1 = [2 3 1 1 2 1 5 6 4];
G1 = digraph(tail1, head1);

% Connectivity for second time interval
% tail2 = [1 2 2 3 3 4 1 4 2 3 4 5 5 6 7 8 8 8 7 9 10 1];
% head2 = [2 1 3 1 4 1 5 5 6 6 7 6 7 8 5 6 7 9 9 10 9 10];

tail2 = [1 2 3 4 5 6 3 4];
head2 = [2 3 4 1 6 5 6 5];
G2    = digraph(tail2, head2);

% Connectivity for second time interval
% tail3 = [1 2 8 8 8 2 9 7 6 4 4 4 3 10 5];
% head3 = [2 8 1 9 7 6 7 6 4 9 3 5 10 5 3];

tail3 = [1 2 3 4 5 6 2 6];
head3 = [2 1 4 3 6 5 3 1];
G3    = digraph(tail3, head3);

% figure(1)
% plot(G1,'Layout','circle')
% figure(2)
% plot(G2,'Layout','circle')
% figure(3)
% plot(G3,'Layout','circle')

N         = G1.numnodes;  % Number of oscillators (nodes)
numStates = 3;   % State dimension of each oscillator

%% Simulation settings
data_length1 = 10;
data_length2 = 10;
data_length3 = 10;
t_end1 = 1;                  % End of first interval
t_end2 = 1;                  % End of second interval
t_end3 = 1;                  % End of third interval
tspan1 = linspace(0, t_end1, data_length1);
tspan2 = linspace(0, t_end2, data_length2);
tspan3 = linspace(0, t_end3, data_length3);

% Assign coupling strengths
G1 = SyncCouplingAssign(G1,a);
figure
plot(G1,'EdgeLabel',G1.Edges.Weight,'Layout','circle')
G2 = SyncCouplingAssign(G2,a);
figure
plot(G2,'EdgeLabel',G2.Edges.Weight,'Layout','circle')
G3 = SyncCouplingAssign(G3,a);
figure
plot(G3,'EdgeLabel',G3.Edges.Weight,'Layout','circle')

% figure(1)
% plot(G1,'EdgeLabel',G1.Edges.Weight,'Layout','circle')
% figure(2)
% plot(G2,'EdgeLabel',G2.Edges.Weight,'Layout','circle')
% figure(3)
% plot(G3,'EdgeLabel',G3.Edges.Weight,'Layout','circle')

%% Initial conditions
x_mean = 5; 
x_std  = 20;
P      = diag([1, 0, 0]);                    % Projection matrix
X0     = x_mean + x_std*rand(1, numStates*N);

%% Simulate coupled Lorenz systems
% First interval (graph G1)
[X1, t1] = SimulateCoupledSystems(@LorenzOscillator, tspan1, X0, G1, P);

% Use final state from first interval as initial state for second
X0 = X1(end,:);

% Third interval (graph G2)
[X2, t2] = SimulateCoupledSystems(@LorenzOscillator, tspan2, X0, G2, P);

% Use final state from second interval as initial state for third
X0 = X2(end,:);

% Second interval (graph G3)
[X3, t3] = SimulateCoupledSystems(@LorenzOscillator, tspan3, X0, G3, P);


% Concatenate results from both intervals
X = [X1;X2(2:end,:);X3(2:end,:)];
t = [t1;t2(2:end,:)+t1(end);t3(2:end,:)+t2(end)+t3(end)];

%% Visualization settings
state_indices   = 1:numStates;
state_index_all = 1:numStates*N;

colors    = lines(N);
linestyle = {'-','--','-.',':','-','--','-.',':','-','--'};
lw        = [2*ones(1,4) 1.5*ones(1,4) 1 1];
markers   = {'none','none','none','none','*','*','o','o','.','.'};

%% Data storage setup
E         = zeros(N, length(t));   % Synchronization error matrix
cols      = 'ABCDEFG';         % Excel column labels
filename  = 'synchronization_data.xlsx';
range_end = length(t) + 1;

%% Plot synchronization errors
figure; hold on; grid on;

for i = 1:N
    % Extract state indices for oscillator i
    slice_i   = (i-1)*numStates + state_indices;
    slice_rem = setdiff(state_index_all, slice_i);

    % Synchronization error (L2 distance from others)
    e = vecnorm(repmat(X(:,slice_i),1,N-1) - X(:,slice_rem), 2, 2).';
    E(i,:) = e;

    % Plot error trajectory
    plot(t, e, ...
        'Color', colors(i,:), ...
        'LineWidth', lw(i), ...
        'LineStyle', linestyle{i}, ...
        'Marker', markers{i}, ...
        'MarkerFaceColor', 'none', ...
        'DisplayName', sprintf('System %d', i));

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
