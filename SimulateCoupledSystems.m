% SimulateCoupledSystems integrates a network of coupled oscillators.
%
%   [X,t] = SimulateCoupledSystems(systemDynamics,tSpan,X0,G,P)
%
%   Inputs:
%     systemDynamics – function handle for oscillator dynamics (dx/dt = f(x))
%     tSpan          – time vector for simulation
%     X0             – initial condition vector for all oscillators
%     G              – digraph defining network connectivity
%     P              – projection matrix (coupling applied to specific states)
%
%   Outputs:
%     X – simulated state trajectories (rows: time, cols: state variables)
%     t – time vector returned by ode45

function [X,t] = SimulateCoupledSystems(systemDynamics,tspan,X0,G,P)   

    %% Build weighted Laplacian matrix from graph
    % Weighted adjacency matrix (transpose because digraph adjacency
    % gives row→col edges, but we want incoming edges in rows)
    A = full(adjacency(G, G.Edges.Weight)).';

    % In-degree Laplacian (row sums = in-degree weights)
    L = diag(sum(A,2)) - A;

    %% Simulate network dynamics with ode45
    [t, X] = ode45(@(t,X) CoupledDynamics(t,X,systemDynamics,L,P),tspan,X0);  
end
