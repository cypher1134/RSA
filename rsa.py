def rsa_chiffrement (x,N,e):
    return

def rsa_dechiffrement (y,p,q,d):
    return

# Retourne s tel que s % n1 == a1 et s % n2 == a2
def crt2 (a1,a2,n1,n2):
    return

def rsa_dechiffrement_crt (y,p,q,up,uq,dp,dq,N):
    return

#### Wiener
def cfrac(a,b):
    return

def reduite(L):
    return

def Wiener(m,c,N,e):
    return


### Generation de premiers
import random
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
        import random
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
