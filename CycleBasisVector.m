% CycleBasisVector updates edge weights of a strongly connected digraph
% so that all edges have nonzero weight while preserving vertex imbalances.
%
%   Inputs:
%     G – strongly connected digraph
%
%   Outputs:
%     G – updated digraph with edge_weight property modified
%
function G = CycleBasisVector(G)

%% Initialize temporary weights
null_weight = zeros(G.numedges,1);

% Strongly Connected Component (SCC) decomposition

nodeSCC = conncomp(G,'Type','strong');
K = max(nodeSCC);

% Extract edge end nodes
E = G.Edges.EndNodes;

vert_imb = -incidence(G)*G.Edges.Weight;


for comp = 1: K
    nodesInComp = find(nodeSCC == comp);
    if isempty(nodesInComp)
        continue;
    else
        % Now find the edges whose **both** endpoints lie in this component:
        % Assuming G.Edges.EndNodes is [M×2] table of node‐IDs (numeric indices)
        % Create a logical mask: source in component & target in component
        edgesInComp_mask = ismember(E(:,1),nodesInComp) & ismember(E(:,2),nodesInComp);
        % edgesInComp = find(mask);                      % indices of edges
        % in this component
    end
    % Continue until all edges have been assigned nonzero weights in this
    % SCC

    while any(null_weight<1 & edgesInComp_mask)

        % --- Step 1: Pick a random edge among those still unweighted ---
        zeroIdx = find(null_weight<1 & edgesInComp_mask);
        eidx    = zeroIdx(randi(numel(zeroIdx)));

        % --- Step 2: Get source (tail) and target (head) nodes of edge ---
        uv = E(eidx,:);
        t  = uv(1);   % tail/source
        h  = uv(2);   % head/target

        % --- Step 3: Find a directed path from head back to tail ---
        % shortestpath can return the sequence of edge indices directly
        [~,~,edgePath] = shortestpath(G,h,t);

        % --- Step 4: Build cycle and update weights ---

        if ~isempty(edgePath)
            C = [eidx, edgePath(:)'];            % cycle edges
            null_weight(C) = null_weight(C) + 1; % increment weights
        end
    end
end

    

    P_sum = ;

     %% Update edge weights of G
     if
         G.Edges.Weight = G.Edges.Weight + null_weight;
     else
         G.Edges.Weight = G.Edges.Weight + null_weight;
     end
end
