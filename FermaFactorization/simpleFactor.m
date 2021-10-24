function [a, b] = simpleFactor(n)
    if mod(n, 2) == 0
        a = 2;
        b = n / 2;
        return;
    end
    t = ceil(sqrt(n));
    while t <= n
        s = t^2 - n;
        if isSquare(s)
            s = sqrt(s);
            a = t + s;
            b = t - s;
            return;
        end
        t = t + 1;
    end
    assert(false);
end

