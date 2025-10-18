function Gout = SyncCouplingAssign(G,a)
% Robust version: supports numeric, string, or categorical node names in EndNodes.

    if ~isa(G,'digraph')
        error('G must be a digraph object.');
    end
    if ~isscalar(a) || a <= 0
        error('a must be a positive scalar.');
    end

    M = numedges(G);
    N = numnodes(G);

    % Initialize or clear existing weights
    if ~ismember('Weight', G.Edges.Properties.VariableNames)
        G.Edges.Weight = zeros(M,1);
    else
        G.Edges.Weight(:) = 0;
    end

    % Strongly connected components
    nodeSCC = conncomp(G,'Type','strong');
    K = max(nodeSCC);

    % Extract edge end nodes
    endNodes = G.Edges.EndNodes;

    % Convert to numeric node indices robustly
    sList = zeros(M,1);
    tList = zeros(M,1);
    for e = 1:M
        if isnumeric(endNodes)
            src = endNodes(e,1);
            tgt = endNodes(e,2);
        elseif isstring(endNodes) || ischar(endNodes) || iscategorical(endNodes)
            src = findnode(G, string(endNodes(e,1)));
            tgt = findnode(G, string(endNodes(e,2)));
        elseif iscell(endNodes)
            src = findnode(G, endNodes{e,1});
            tgt = findnode(G, endNodes{e,2});
        else
            error('Unsupported EndNodes type: %s', class(endNodes));
        end
        sList(e) = src;
        tList(e) = tgt;
    end

    % Build incoming/outgoing edge lists
    incomingEdges = cell(N,1);
    outgoingEdges = cell(N,1);
    for e = 1:M
        incomingEdges{tList(e)} = [incomingEdges{tList(e)}; e];
        outgoingEdges{sList(e)} = [outgoingEdges{sList(e)}; e];
    end

    % Process each SCC
    for comp = 1:K
        nodesInComp = find(nodeSCC == comp);
        if isempty(nodesInComp)
            continue;
        end

        % --- Choose representative (chosen) vertex
        candidateNodes = [];
        for v = nodesInComp(:).'
            incE = incomingEdges{v};
            if isempty(incE), continue; end
            if any(nodeSCC(sList(incE)) ~= comp)
                candidateNodes(end+1) = v; %#ok<AGROW>
            end
        end

        if ~isempty(candidateNodes)
            chosenVertex = candidateNodes(randi(numel(candidateNodes)));
        else
            chosenVertex = nodesInComp(randi(numel(nodesInComp)));
        end

        % --- Aggregate incidence vectors
        sumVec = zeros(M,1);
        for v = nodesInComp(:).'
            if v == chosenVertex, continue; end
            pth = shortestpath(G, chosenVertex, v);
            if isempty(pth)
                error('No path found inside SCC %d from %d to %d.', comp, chosenVertex, v);
            end
            L = numel(pth)-1;
            for s = 1:L
                u = pth(s); w = pth(s+1);
                eIdx = find(sList==u & tList==w, 1);
                if isempty(eIdx)
                    error('Edge inconsistency: %d -> %d not found.', u, w);
                end
                sumVec(eIdx) = sumVec(eIdx) + (L-s+1);
            end
        end

        % Scale aggregated vector
        G.Edges.Weight = G.Edges.Weight + 2*a*sumVec;

        % --- Assign weights for incoming edges from other SCCs
        incE_chosen = incomingEdges{chosenVertex};
        incFromOther = incE_chosen(nodeSCC(sList(incE_chosen)) ~= comp);
        if ~isempty(incFromOther)
            outgoingSum = sum(G.Edges.Weight(outgoingEdges{chosenVertex}));
            incomingSum = sum(G.Edges.Weight(incE_chosen));
            imbalanceChosen = outgoingSum - incomingSum;
            totalWeight = (imbalanceChosen/2 + a);
            G.Edges.Weight(incFromOther) = G.Edges.Weight(incFromOther) + totalWeight/numel(incFromOther);
        end

        otherVerts = setdiff(nodesInComp, chosenVertex);
        for v = otherVerts(:).'
            incE_v = incomingEdges{v};
            if isempty(incE_v), continue; end
            incFromOther_v = incE_v(nodeSCC(sList(incE_v)) ~= comp);
            if ~isempty(incFromOther_v)
                G.Edges.Weight(incFromOther_v) = G.Edges.Weight(incFromOther_v) + a/numel(incFromOther_v);
            end
        end

        % --- Verify imbalance condition
        imbalances = zeros(numel(nodesInComp),1);
        for k = 1:numel(nodesInComp)
            v = nodesInComp(k);
            outgoingSum = sum(G.Edges.Weight(outgoingEdges{v}));
            incomingSum = sum(G.Edges.Weight(incomingEdges{v}));
            imbalances(k) = outgoingSum - incomingSum;
        end
        chosenIdx = find(nodesInComp == chosenVertex);
        if imbalances(chosenIdx) > 0 && all(imbalances(setdiff(1:numel(nodesInComp),chosenIdx)) <= 0)
            fprintf('SCC %d: chosen vertex %d OK (imbalance %.4g)\n', comp, chosenVertex, imbalances(chosenIdx));
        else
            warning('SCC %d imbalance check failed.', comp);
        end
    end

    % --- Final output + visualization
    Gout = G;
    figure;
    p = plot(Gout,'Layout','layered','EdgeLabel',Gout.Edges.Weight,'NodeLabel',imbalances);
    allImb = zeros(N,1);
    for v = 1:N
        outE = outgoingEdges{v};
        inE  = incomingEdges{v};
        allImb(v) = sum(Gout.Edges.Weight(outE)) - sum(Gout.Edges.Weight(inE));
    end
    p.NodeCData = allImb;
    colorbar;
    title('Directed graph with edge weights and vertex imbalance');
end
