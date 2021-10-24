function X = binlineq(T)
% X = binlineq(T)
% solve a homogenous linear equation over the field F_2
% i.e. returns a nontrivial solution of equation: T * X = 0
% if it's impossible, returns a trivial solution
    i = 1;
    N = size(T, 1);
    inds = zeros(1, N);
    % transform matrix into upper triangular form
    while i <= N
        non0 = find(T(i, :), 1);
        if isempty(non0)
            T(i, :) = [];
            N = N - 1;
            continue;
        end
        inds(i) = non0;
        toUpd = [false(i, 1); (T(i+1 : end, non0) == 1)];
        T(toUpd, :) = xor(T(toUpd, :), T(i, :));
        i = i + 1;
    end
    inds = inds(1:N);
    T1 = T(:, inds);
    X = zeros(size(T, 2), 1);
    % backward propagation
    nonbase = setdiff(1:size(T, 2), inds);
    if isempty(nonbase), return; end
    [b, initInds] = chooseInit(T, nonbase); 
    X(initInds) = 1;
    X(inds) = backprop(T1, b);
end

function x = backprop(T1, b)
    N = size(T1, 1);
    x = zeros(N, 1);
    for i = N : -1 : 1
        x(i) = b(i);
        b = xor(b, T1(:, i) * x(i));
    end
end

function [vec, inds] = chooseInit(T, nonbaseInds)
    T0 = T(:, nonbaseInds);
    % try to find square-free numbers
    hasOdd = any(mod(T0(2:end, :), 2), 1);
    oddInds = find(hasOdd);
    if isempty(oddInds)
        % try to find square-free or negative numbers
        hasOddFull = any(mod(T0, 2), 1);
        oddInds = find(hasOddFull);
    end
    if isempty(oddInds)
        inds = nonbaseInds;
    else
        inds = nonbaseInds(oddInds(1));
    end
    vec = mod(sum(T(:, inds), 2), 2);
end
