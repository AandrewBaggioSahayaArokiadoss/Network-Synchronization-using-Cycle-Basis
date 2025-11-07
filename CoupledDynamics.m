% This function returns dX/dt for the entire network by evaluating the
% state of each system

function dXdt = CoupledDynamics(t,X,systemDynamics,L,P)

N = length(L);
stateDim = size(X,1)/N;
X = reshape(X,stateDim,N); % Reshape into state matrix [stateDim x N]
couplingTerm = L * X.'; % Diffusive coupling term

dXdt = zeros(stateDim,N);
for i = 1:N
    dXdt(:, i) = systemDynamics(t, X(:, i)) - P*couplingTerm(i,:).' - P* X(:, i);
end
dXdt = dXdt(:); % Flatten back to column vector

end
