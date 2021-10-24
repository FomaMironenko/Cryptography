function [s, t] = factorization(n)
% [s, t] = factorization(n)
    if mod(n, 2) == 0, s = 2; t = n/2; return; end
    if mod(n, 3) == 0, s = 3; t = n/3; return; end
    Zp = ModuloRing(n);
    start = max(2, floor(n ^ (1/3)));
    for y = start : 1 : 10*sqrt(n)
        [B, ~, FB, T] = continuedFracFB(y, Zp);
        assert(all(T(:) >= 0));
        assert(all(round(T(:)) == T(:)));
        X = binlineq(mod(T, 2));
        if all(X == 0), continue; end
        P = sum(T(:, logical(X)), 2);
        assert(all(mod(P, 2) == 0));
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
        s = gcd(b1 + b2, n);
        assert(s ~= 1);
        t = n / s;
        return;
    end
    s = [];
    t = [];
end

function n = max_pow(a, p) %#ok <UNUSD>
    n = 0;
    while mod(a, p) == 0
        a = a / p;
        n = n + 1;
    end
end

function [B, C, FB] = randGenFB(y, Zp) %#ok <UNUSD>
    FB = [-1, primes(y)];
    B = [];
    C = [];
    for i = 1 : y^2
        b = round(2 + rand() * (Zp.base - 3));
        c = Zp.min(Zp.pow(b, 2));
        fact = factor(abs(c));
        if all(ismember(fact, FB)) && ~ismember(b, B)
            B = [B, b]; %#ok <ARGOV>
            C = [C, c]; %#ok <ARGOV>
            if numel(C) >= numel(FB) + 1, return; end
        end
    end
end

function [Bres, Cres, FB, T] = continuedFracFB(k, Zp)
    n = Zp.base;
    x = sqrt(double(n));
    a0 = uint64(floor(x));
    A = uint64(zeros(k + 2, 1));
    B = uint64(zeros(k + 2, 1));
    A(1:2) = [1, a0];
    B(1:2) = [1, a0];
    Bres = [];
    Cres = [];
    x = x - double(a0);
    FB = [-1, primes(sqrt(double(n)))];
    T = zeros(numel(FB), 0);
    for i = 3 : k
        assert(x ~= 0);
        A(i) = floor(1 / x);
        B(i) = Zp.add(...
            Zp.mul(A(i), B(i-1)), ...
            B(i-2) ...
        );
        c = Zp.min(Zp.pow(B(i), 2));
        fact = factor(abs(c));
        if (all(ismember(fact, FB)) || all(fact == 1)) && ...
                ~ismember(c, Cres)
            Bres = [Bres, B(i)];   %#ok <ARGOV>
            Cres = [Cres, c];      %#ok <ARGOV>
            T(1, end+1) = (c < 0); %#ok <ARGOV>
            for j = 2 : numel(FB)
                T(j, end) = sum(FB(j) == fact);
            end
        end
        if numel(Cres) > numel(FB) + 1, return; end
        x = 1 / x - double(A(i));
    end
end

% function [s, t] = factorization(n)
% % [s, t] = factorization(n)
%     rng(1);
%     if mod(n, 2) == 0, s = 2; t = n/2; return; end
%     if mod(n, 3) == 0, s = 3; t = n/3; return; end
%     Zp = ModuloRing(n);
%     start = max(2, floor(n ^ (1/3)));
%     for y = start : 1 : 10*sqrt(n)
%         [B, C, FB] = randGenFB(y, Zp);
%         FB = sort(unique(FB));
%         C = Zp.min(C);
%         % T is a multipliers matrix, i.e. j-th column represents
%         % a decomposition of C(j) into multipliers from factor base FB
%         T = zeros(numel(FB), numel(C));
%         T(1, :) = (C < 0);
%         for i = 2 : numel(FB)
%             for j = 1 : numel(C)
%                 T(i, j) = max_pow(C(j), FB(i));
%             end
%         end
%         assert(all(T(:) >= 0));
%         assert(all(round(T(:)) == T(:)));
%         X = binlineq(mod(T, 2));
%         if all(X == 0), continue; end
%         P = sum(T(:, logical(X)), 2);
%         assert(all(mod(P, 2) == 0));
%         b1 = 1;
%         b2 = 1;
%         for i = find(X)'
%             b1 = Zp.mul(b1, B(i));
%         end
%         for j = 1 : numel(P)
%             b2 = Zp.mul( b2, Zp.pow(FB(j), P(j)/2) );
%         end
%         assert(Zp.pow(b1, 2) == Zp.pow(b2, 2));
%         if b1 == b2 || b1 == n - b2, continue; end
%         s = gcd(b1 + b2, n);
%         assert(s ~= 1);
%         t = n / s;
%         return;
%     end
%     s = [];
%     t = [];
% end
