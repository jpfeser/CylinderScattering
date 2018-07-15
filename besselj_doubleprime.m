function f = besselj_doubleprime(n,z)

f = -1./z.*besselj_prime(n,z)+(((n./z).^2)-1).*besselj(n,z);

end
