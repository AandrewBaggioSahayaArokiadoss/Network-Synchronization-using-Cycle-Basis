function G = SyncCouplingAssign(G,a)
% Robust version: supports numeric, string, or categorical node names in EndNodes.

    if ~isa(G,'digraph')
        error('G must be a digraph object.');
    end
    if ~isscalar(a) || a <= 0
        error('a must be a positive scalar.');
    end
    
    G = NegativeImbalanceVector(G,a);
    
    G = CycleBasisVector(G,a);

    figure
    plot(G,'Layout','circle','EdgeLabel',G.Edges.Weight);
    title('Directed graph with edge weights and vertex imbalance');
end
