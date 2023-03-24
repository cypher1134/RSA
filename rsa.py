import random

def gcd(a,b):
  while b != 0:
      a,b = b, a % b
  return a

def xgcd(a,b):
    """Extended GCD:
    Returns (gcd, x, y) where gcd = ax+by."""
    prevx, x = 1, 0;  prevy, y = 0, 1
    while b:
        q, r = divmod(a,b)
        x, prevx = prevx - q*x, x  
        y, prevy = prevy - q*y, y
        a, b = b, r
    return a, prevx, prevy

def invmod(a, m):
    g, x, _ = xgcd(a, m)
    if g!= 1:
        raise ValueError('a n\'est pas inversible modulo m')
    return x % m

def rsa_chiffrement (x,N,e):
    return pow(x, e, N)

def rsa_dechiffrement (y,p,q,d):
    n = p*q
    return pow(y, d, n)

# Retourne s tel que s % n1 == a1 et s % n2 == a2
def crt2(a1, a2, n1, n2):
    _, x, y = xgcd(n1, n2)
    m = n1 * n2
    s1 = (a1 * n2 * y + a2 * n1 * x) % m
    s2 = (s1 - a1) * invmod(n1, n2) % n2
    return s1, s2

def rsa_dechiffrement_crt(y, p, q, up, uq, dp, dq, N):
    mp = pow(y, dp, p)
    mq = pow(y, dq, q)
    s = crt2(mp, mq, p, q)
    return (s * up * q + uq * p) % N

#### Wiener
def cfrac(a, b):
    coeffs = []
    while b != 0:
        q = a // b
        r = a % b
        coeffs.append(q)
        a = b
        b = r
    return coeffs



def reduite(L):
    n = [] # Nominators
    d = [] # Denominators

    for i in range(len(L)):
        if i == 0:
            ni = L[i]
            di = 1
        elif i == 1:
            ni = L[i]*L[i-1] + 1
            di = L[i]
        else: # i > 1
            ni = L[i]*n[i-1] + n[i-2]
            di = L[i]*d[i-1] + d[i-2]

        n.append(ni)
        d.append(di)
        yield (ni, di)


def Wiener(m, c, N, e):

    cf_expansion = cfrac(e, N)
    convergents = reduite(cf_expansion) 
    
    for pk, pd in convergents: # pk - possible k, pd - possible d
        if pk == 0:
            continue;
        if pow(c,pd,N)==m:
            return pd

    print('[-] Wiener\'s Attack failed; Could not factor N')
    return None


### Generation de premiers

def is_probable_prime(N, nbases=20):
    """
    True if N is a strong pseudoprime for nbases random bases b < N.
    Uses the Miller--Rabin primality test.
    """

    def miller(a, n):
        """
        Returns True if a proves that n is composite, False if n is probably prime in base n
        """

        def decompose(i, k=0):
            """
            decompose(n) returns (s,d) st. n = 2**s * d, d odd
            """
            if i % 2 == 0:
                return decompose(i // 2, k + 1)
            else:
                return (k, i)

        (s, d) = decompose(n - 1)
        x = pow(a, d, n)
        if (x == 1) or (x == n - 1):
            return False
        while s > 1:
            x = pow(x, 2, n)
            if x == n - 1:
                return False
            s -= 1
        return True

    if N == 2:
        return True
    for i in range(nbases):
        
        a = random.randint(2, N - 1)
        if miller(a, N):
            return False
    return True


def random_probable_prime(bits):
    """
    Returns a probable prime number with the given number of bits.
    Remarque : on est sur qu'un premier existe par le postulat de Bertrand
    """
    n = 1 << bits
    import random
    p = random.randint(n, 2 * n - 1)
    while (not (is_probable_prime(p))):
        p = random.randint(n, 2 * n - 1)
    return p
