function f = besselh_doubleprime(n,z)

f = -1./z.*besselh_prime(n,z)+((n./z).^2-1).*besselh(n,z);

end


