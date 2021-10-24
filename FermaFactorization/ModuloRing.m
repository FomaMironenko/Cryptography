classdef ModuloRing
% ModuloRing
% each object introduces an interface for modulo ring with 
% base: \prop base
    
    properties (SetAccess = immutable)
        base uint64;
    end
    
    methods
        function this = ModuloRing(n)
            assert(n > 0 && round(n) == n);
            this.base = uint64(n);
        end
        
        function c = rep(this, a)
            c = uint64(mod(double(a), double(this.base)));
        end
        
        function c = min(this, a)
            b = this.rep(a);
            c = int64(zeros(numel(b), 1));
            c(b > this.base/2) = -int64(this.base - b(b > this.base/2));
            c(b <= this.base/2) = b(b <= this.base/2);
        end
        
        function c = add(this, a, b)
            a = this.rep(a);
            b = this.rep(b);
            c = mod(a + b, this.base);
        end
        
        function c = mul(this, a, b)
            a = this.rep(a);
            b = this.rep(b);
            c = mod(a .* b, this.base);
        end
        
        function c = inv(this, a)
            [g, x, ~] = gcd(double(a), -double(this.base));
            if g ~= 1, c = 0; return; end
            c = this.rep(x);
        end
        
        function c = div(this, a, b)
            c = this.mul(a, this.inv(b));
        end
        
        function c = pow(this, a, n)
            if n == 0
                c = 1;
            elseif mod(n, 2) == 0
                tmp = this.pow(a, n / 2);
                c = this.mul(tmp, tmp);
            else
                c = this.mul(a, this.pow(a, n - 1));
            end
        end
    end
end
