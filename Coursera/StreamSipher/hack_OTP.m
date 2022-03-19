c1 = 'attack at dawn';
s1 = '6c73d5240a948c86981bc294814d';
s1 = reshape(s1, 2, [])';

m = uint8(c1);
mk = uint8(hex2dec(s1));
k = bitxor(m, mk);

c2 = 'attack at dusk';
m2 = uint8(double(c2));

mk2 = bitxor(m2, k);
s2 = lower(dec2hex(mk2))
reshape(s2', 1, [])
    