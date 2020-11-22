import numpy as np
from sympy.ntheory.factor_ import totient as eulerphi

def circmat(v, m):
    '''This function produces a circulant matrix of the type
    that is used in the function lfsrlength'''
    if len(v) < (2 * m - 1):
        raise ValueError("length of vector must be < 2 * m - 1")
    return np.array([v[i:i+m] for i in range(m)])

def coinc(txt, n):
    """This function returns the number of matches between txt and txt 
    shifted by n"""
    if isinstance(txt, str):
        txt = np.array(list(txt))
    return np.sum(np.roll(txt, n) == np.array(txt))
    
def convm(x, p):
    x = np.array(x)
    if isinstance(x, np.ndarray):
        x = x.reshape((len(x), 1))
    X = np.zeros((len(x) + p - 1, p))
    for i in range(p):
        X[:,i][i:i + len(x)] = x.T
    return X

def corr(v):
    fvec=(.082, .015, .028, .043, .127, .022, .020, .061, .070, .002, .008, .040, .024, .067, .075, .019, .001, .060, .063, .091, .028, .010, .023, .001, .020, .001)
    return (circulant(fvec).T * v)[:,2]

def crt(a, m):
    """This function solves the Chinese Remainder Theorem problem:
    x= a(1) mod m(1)
    x= a(2) mod m(2)
    ...
    x= a(r) mod m(r)
    The values for a and m should be a vector of the same dimension"""
    r = len(a)
    M = np.prod(m)
    x = 0
    for i in range(r):
        x += a[i] * M/m[i] * invmodn(M/m[i], m[i])
        x %= M
    return int(x)

def invmodn(b, n):
    """This function calculates the inverse of an element b mod n
    uses extended Euclidian"""
    n0 = n
    b0 = b
    t0 = 0
    t = 1

    q = n0 // b0
    r = n0 - q * b0
    while r > 0:
       temp = t0 - q * t
       if (temp >= 0):
          temp = temp % n
       if (temp < 0):
          temp = n - (-temp % n)
       t0 = t
       t = temp
       n0 = b0
       b0 = r
       q = n0 // b0
       r = n0 - q * b0

    if b0 !=1:
       return None
    else:
       return t % n

# eulerphi goes here, already defined

def frequency(seq):
    """This function tabulates the amount of times each character 'a'-'z' 
    occurs in txt."""
    counts = [0 for _ in range(26)]
    for c in seq:
        counts[c2n(c)] += 1
    return counts

# I'm forgoing correct spacing here for the slickness of the indentation. don't judge me
golay=np.array([[1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 1, 1, 0, 0, 0, 1, 0],\
                [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1, 0, 0, 0, 1],\
                [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 1, 0, 1, 1, 1, 0, 0, 0],\
                [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 1, 0, 1, 1, 1, 0, 0],\
                [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1, 0],\
                [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1],\
                [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 1, 1, 0, 1, 1],\
                [0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 0, 1, 1, 0, 1],\
                [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 1, 0, 1, 1, 0],\
                [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 0, 1, 1],\
                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1, 0, 0, 0, 1, 0, 1],\
                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]])

golayb = golay[:, 12:23]

hammingpc=np.array([[0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0],\
                    [1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 0, 1, 0, 0],\
                    [0, 1, 1, 1, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 0],\
                    [1, 0, 1, 1, 1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 1]])
def int2text(x):
    """This function takes the number in the variable x and converts 
    it to a character y. 
    The convention for this function is
      a  <-- 0
      b  <-- 1
      and so on..."""
    return chr(x + ord('a'))

def int2text1(x):
    """This function takes the number in the variable x and converts 
    it to a character y. 
    The convention for this function is
      a  <-- 0
      b  <-- 1
      and so on..."""
    s = ''
    i = 0
    while x:
        s += int2text((x % 100) - 1)
        x //= 100
    return s[::-1]

# still not super sure what this one does
def interppoly(x, f, m):
    """This function calculates the interpolating polynomial that goes 
    through the points (x_j,f_j) in the vectors x and f
    The modulus is m
    The method used in this program is the Newton form of the interpolant."""
    if len(x) != len(f):
       raise ValueError('Vectors are not the same length')
    
    n = len(x)
    x = np.array(x)
    f = np.array(f)
    y = f
    temp = 0
    for k in range(1, n):
        for j in range(n - 1, k - 1, -1):
            temp = (y[j]-y[j-1]) % m
            temp1 = invmodn(x[j]-x[j-k], m)
            y[j] = (temp*temp1) % m
          
    y[1] = y[1] % m

    for k in range(n - 1, 0, -1):
        for j in range(k, n):
            temp = (y[j]*x[k-1]) % m
            y[j-1] = (y[j-1]-temp) % m
    return y
inv_field = invmodn

def khide(p):
    """This function generates a random integer k that is between
    1 and p-2 and is relatively prime to p-1."""
    while True:
       k = 1 + int(np.random.rand() * (p - 2))
       if np.gcd(k, p - 1) == 1:
          break
    return k

def lfsr(c, k, n):
    """This function generates a sequence of n bits produced by the linear feedback
    recurrence relation that is governed by the coefficient vector c.
    The initial values of the bits are given by the vector k"""

    y = [0 for _ in range(n)]

    kln = len(k)
    for j in range(n):
        if j < kln:
            y[j] = k[j]
        else:
            reg = y[j-kln:j]
            y[j] = np.mod(np.dot(reg, c), 2)
    return y

def lfsrlength(v,n):
    """This function tests the vector v of bits to see if it is generated
    by a linear feedback recurrence of length at most n"""

    res = dict()
    for j in range(1, n+1):
       M = circmat(v,j)
       Mdet = int(np.round(np.linalg.det(np.array(M)))) % 2
       res[j] = Mdet
       if j == 38:
           det = np.linalg.det(np.array(M))
    return res

def lfsrsolve(v,n):
    """Given a guess n for the length of the recurrence that generates
    the binary vector v, this function returns the coefficients of the
    recurrence."""
    v = v[:]

    vln = len(v)

    if (vln < 2*n):
       raise ValueError('The vector v needs to be atleast length 2n')

    M = np.array(circmat(v,n))
    Mdet = np.linalg.det(M)

    x = np.array(v[n:2*n]).reshape(-1, 1)
    Minv = np.linalg.inv(M)
    Minv = np.mod(np.round(Minv*Mdet), 2)
    # A note: Technically, the round() function should never show up, but
    # since Matlab does double precision arithmetic to calculate the inverse matrix
    # we need to bring the result back to integer values so we can perform a meaningful
    # mod operation. As long as this routine is not used on huge examples, it should
    # be ok

    y = np.mod(np.matmul(Minv, x), 2)
    return list(y[:,0])

##def multell(p,M,a,b,n):
##% This function prints the Mth multiple of p on the elliptic
##%  curve with coefficients a and b mod n.
##
##z1=M;
##y = [inf inf];
##while z1 != 0:
##    while (np.mod(z1,2) == 0):
##        z1 = (z1/2)
##        p = addell(p, p, a, b, n)
##        if (length(p)==0):
##           y=[];
##           disp('Multell found a factor of n and exited');
##           z1
##           return y
##    z1 = z1 - 1
##    y=addell(y, p, a, b, n)
##    if (length(y)==0)
##       return y

# https://gist.github.com/Ayrx/5884790
def primetest(n, k=30):
    """This function tests a number n for primeness.
    It uses the Miller-Rabin primality (compositeness) test
    z=1 means prime
    z=0 means composite"""
    # Implementation uses the Miller-Rabin Primality Test
    # The optimal number of rounds for this test is 40
    # See http://stackoverflow.com/questions/6325576/how-many-iterations-of-rabin-miller-should-i-use-for-cryptographic-safe-primes
    # for justification

    # If number is even, it's a composite number

    if n == 2:
        return True

    if n % 2 == 0:
        return False

    r, s = 0, n - 1
    while s % 2 == 0:
        r += 1
        s //= 2
    for _ in range(k):
        a = np.random.randint(2, n - 1)
        x = pow(a, s, n)
        if x == 1 or x == n - 1:
            continue
        for _ in range(r - 1):
            x = pow(x, 2, n)
            if x == n - 1:
                break
        else:
            return False
    return True

def _powermod_single(a,z,n):

    #take care of negative exponent
    if (z<0):
        z *= -z
        a = invmodn(a, n)

    return pow(a, z, n)

def powermod(a, z, n):
    """Calculates y = a^z mod n
If a is a matrix, it calculates a(j,k)^z mod for every element in a"""
    if isinstance(a, (int, float)):
        return _powermod_single(a, z, n)
    #this one takes in lists for a and b
    ret = []
    for x in a:
        ret.append(_powermod_single(x, z, n))
    return(ret)


if __name__ == "__main__":
    pass
