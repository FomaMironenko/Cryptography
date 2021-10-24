function [s, t] = factorization(n)
% [s, t] = factorization(n)
    if mod(n, 2) == 0, s = 2; t = n/2; return; end
    if mod(n, 3) == 0, s = 3; t = n/3; return; end
    Zp = ModuloRing(n);
    computer = struct(...
        'B', uint64([0, 0]),... % last two B's from iterations
        'T', [],            ... % a dividors matrix
        'FB', [],           ... % a dynamically-updated factor base
        'x', [],            ... % a remainder of a continuous fraction
        'Bres', uint64([]), ... % a sequence of b_i
        'Cres', int64([])   ... % a sequence of squares of b_i in Zp
    );
    for k = 1 : ceil(n^(1/4))
        % continue iterations
        computer = continuedFracFB(computer, Zp, 10);
        T = computer.T;
        B = computer.Bres;
        C = computer.Cres;
        if numel(C) ~= numel(unique(C))
            % all elements of B are unique due to continuedFracFB function
            [C, inds] = sort(C);
            B = B(inds);
            d = diff(C);
            i = find(d == 0, 1);
            b1 = B(i);
            b2 = B(i+1);
        else
            FB = computer.FB;
            assert(all(T(:) >= 0));
            assert(all(round(T(:)) == T(:)));
            % solve a linear system
            X = binlineq(mod(T, 2));
            if all(X == 0), continue; end
            P = sum(T(:, logical(X)), 2);
            assert(all(mod(P, 2) == 0));
            % try to extract a divider
            b1 = 1;
            b2 = 1;
            for i = find(X)'
                b1 = Zp.mul(b1, B(i));
            end
            for j = 1 : numel(P)
                b2 = Zp.mul( b2, Zp.pow(FB(j), P(j)/2) );
            end
            assert(Zp.pow(b1, 2) == Zp.pow(b2, 2));
            if b1 == b2 || b1 == n - b2, continue; end
        end
        % luck luck, divider found
        s = gcd(b1 + b2, n);
        assert(s ~= 1);
        t = n / s;
        return;
    end
    s = [];
    t = [];
end

function computer = continuedFracFB(computer, Zp, Niters)
    n = Zp.base;
    % init data
    if isempty(computer.x)
        x = sqrt(double(n));
        a0 = uint64(floor(x));
        x = x - double(a0);
        B = uint64([1, a0]);
        Bres = uint64([]);
        Cres = int64([]);
        FB = -1;
        T = zeros(numel(FB), 0);
    else
        x = computer.x;
        B = computer.B;
        Bres = computer.Bres;
        Cres = computer.Cres;
        FB = computer.FB;
        T = computer.T;
    end
    % new iterations
    for i = 1 : Niters
        assert(x ~= 0);
        % compute new values
        A = floor(1 / x);
        B = [
            Zp.add( Zp.mul( A, B(1) ), B(2) ),...
            B(1) ...
        ];
        c = Zp.min(Zp.pow(B(1), 2));
        if ismember(B(1), Bres) || ismember(n - B(1), Bres), continue; end
        % update factor base
        fact = factor(abs(c));
        newFact = setdiff(fact, [1, FB]);
        FB = [FB, newFact]; %#ok <ARGOV>
        T = [T; zeros(numel(newFact), size(T, 2))]; %#ok <ARGOV>
        % add new element resulting sequence
        Bres = [Bres, B(1)];   %#ok <ARGOV>
        Cres = [Cres, c];      %#ok <ARGOV>
        T(1, end+1) = (c < 0); %#ok <ARGOV>
        for j = 2 : numel(FB)
            T(j, end) = sum(FB(j) == fact);
        end
        if numel(Cres) >= numel(FB) + 1, break; end
        x = 1 / x - double(A);
    end
    % return data to computer
    computer.x = x;
    computer.B = B;
    computer.Bres = Bres;
    computer.Cres = Cres;
    computer.FB = FB;
    computer.T = T;
end
