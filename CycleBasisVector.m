function G = CycleBasisVector(G, a)
    arguments
        G {mustBeA(G,"digraph")} 
        a (1,1) double {mustBePositive}
    end

    %% Ensure weight property exists
    if ~ismember('Weight', G.Edges.Properties.VariableNames)
        G.Edges.Weight = zeros(G.numedges, 1);
    end

    % %% Initialize temporary weights
    null_weight = zeros(G.numedges,1);

    % Strongly Connected Component (SCC) decomposition
    nodeSCC = conncomp(G,'Type','strong');
    K = max(nodeSCC);

    % Extract edge end nodes
    E = G.Edges.EndNodes;

    for comp = 1:K
        nodesInComp = find(nodeSCC == comp);
        if isempty(nodesInComp)
            continue;
        else
            % find edges whose both endpoints are in this component
            edgesInComp_mask = ismember(E(:,1),nodesInComp) & ismember(E(:,2),nodesInComp);
        end

        %% --- NEW CODE BLOCK: compute D_root, P_sum, scaling_factor ---
        % Restrict to edges within this component
        % eidxComp = find(edgesInComp_mask);
        % For each vertex in the component compute (sum of outgoing weights within comp) minus (sum of incoming weights within comp)
        w = G.Edges.Weight;   % existing weight vector (must exist)
        diff_vect = zeros(numel(nodesInComp),1);
        for i = 1:numel(nodesInComp)
            v = nodesInComp(i);
            % outgoing edges from v to nodesInComp
            outMask = (E(:,1)==v) & ismember(E(:,2),nodesInComp);
            inMask  = (E(:,2)==v) & ismember(E(:,1),nodesInComp);
            sum_out = sum( w(outMask));
            sum_in  = sum( w(inMask));
            diff_vect(i) = sum_out - sum_in;
        end

        % pick the maximal positive difference (if any)
        posDiffs = diff_vect(diff_vect>0);
        D_root = max(posDiffs);

        % Check if this SCC has incoming edges from other SCCs
        % i.e., edges whose head is in this SCC and whose tail is *not* in this SCC
        incoming_from_other_mask = ismember(E(:,2),nodesInComp) & ~ismember(E(:,1),nodesInComp);
        % m_comp = sum(edgesInComp_mask>0);
        % null_weight=(zeros(numel(m_comp),1));
        if ~any(incoming_from_other_mask)
            % no incoming from other SCCs
            P_sum = D_root / (a);
            numVerts = numel(nodesInComp);
            scaling_factor = 2 * a * P_sum * (1 + P_sum) / numVerts;
        else
            % has incoming from other SCCs
            scaling_factor = 1;
        end

        % Continue until all edges in this SCC have been assigned non‐zero weights
        while any(null_weight<1 & edgesInComp_mask)
            zeroIdx = find(null_weight<1 & edgesInComp_mask);
            eidx    = zeroIdx(randi(numel(zeroIdx)));
            uv = E(eidx,:);
            t  = uv(1);
            h  = uv(2);
            [~,~,edgePath] = shortestpath(G,h,t);
            if ~isempty(edgePath)
                C = [eidx, edgePath(:)'];
                null_weight(C) = null_weight(C) + scaling_factor;
            end
        end
    end
    %% Update edge weights of G
    G.Edges.Weight = G.Edges.Weight + null_weight;
    figure
    plot(G,'EdgeLabel',null_weight,'Layout','circle')
end
