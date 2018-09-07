class Poly:
    ### Polynomials are represented as a list with
    ### leading coefficient (highest power of x) at position zero.
    def __init__(self, coeffs):
        assert coeffs.__class__ == list or coeffs.__class__ == tuple
        for a in coeffs:
            if a.__class__ != int:
                raise TypeError ('Poly only currently implemented for ints')
        while len(coeffs) > 0:
            if coeffs[0] == 0:
                coeffs.pop(0)
            else:
                break
        self.c = coeffs
        self.deg = len(coeffs)-1
        self.len = self.deg + 1
        self.content = 'int'
    def __repr__(self):
        def nxt_nz(alist, cur_ix):
            for i,x in enumerate(alist[(cur_ix+1):]):
                if x != 0:
                    return cur_ix + i + 1
            return None
        if len(self.c) == 0: # constant (zero) polynomial
            return "0"
        elif len(self.c) == 1:
            return str(self.c[0]) # constant (non-zero) polynomial
        else: # second or higher
            # First write the leading term:
            p = self.deg
            if p == 1:
                if self.c[0] == 1:
                    rep = "x"
                elif self.c[0] == -1:
                    rep = "-x"
                else:
                    rep = "%dx" % self.c[0]
            else:
                if self.c[0] == 1:
                    rep = "x**%d" % p
                elif self.c[0] == -1:
                    rep = "-x**%d" % p
                else:
                    rep = "%dx**%d" % (self.c[0], self.deg)
            # Find the index of next non-zero term:
            nxt = nxt_nz(self.c, 0)
            while not nxt is None:
                if self.c[nxt] < 0:
                    rep += " - "
                else:
                    rep += " + "
                a = abs(self.c[nxt])
                p = self.deg - nxt
                if p == 0:
                    rep += "%d" % a
                elif p == 1:
                    if a == 1:
                        rep += "x"
                    else:
                        rep += "%dx" % a
                else:
                    if a == 1:
                        rep += "x**%d" % p
                    else:
                        rep += "%dx**%d" % (a,p)
                    
                nxt = nxt_nz(self.c, nxt)
            return rep
    def __mul__(self,g):
        return self.mul(g)
    def __add__(self,g):
        return self.add(g)
    def __neg__(self):
        return self.neg()
    def __sub__(self,g):
        return self.sub(g)

    def mul(self, g):
        '''fft would be faster but float rounds off, this method is slower
        but exact!'''
        # Naive approach
        P = [0,] * (self.len+g.len-1)
        for i,x in enumerate(self.c):
            for j,y in enumerate(g.c):
                P[i+j] += x*y
        return Poly(P)

    def add(self, g):
        nf = self.len
        ng = g.len
        N = max(nf, ng)
        Xf = [0,]*(N-nf) + self.c
        Xg = [0,]*(N-ng) + g.c
        S = []
        for x,y in zip(Xf,Xg):
            S.append(x + y)
        return Poly(S)

    def neg(self):
        return Poly([-x for x in self.c])

    def sub(self, g):
        return self.add(g.neg())


    def gcd(self,g):
        if self.deg > g.deg:
            u = self
            v = g
        else:
            u = g
            v = self
        while v != [0,]:
            t = u.div(v)[1]
            if t.len == 0:
                return u
            while t.c[0] == 0:
                t.pop(0)
            u = v
            v = t
        return u

    def div(self,g):
        if len(self.c) < len(g.c):
            return Poly(()), self
        assert g.c[0] == 1 # This ensures all coefficients of solution are
                         # integers
        n = self.len
        nv = g.len
        r = list(self.c)
        q = []
        for k in range(n-nv+1):
            qa = r[k]
            q.append(qa)
            for j in range(1,nv):
                r[j+k] -= qa*g.c[j]
        return Poly(q), Poly(r[(-nv+1):])

    def rem(self,g):
        return div(self,g)[1]

    def eval(self,x):
        ans = self.c[0]
        for c in self.c[1:]:
            ans = ans * x + c
        return ans


        

class Fraction:
    ### Base methods
    def __init__(self, num, denom):
        self.n = num
        self.d = denom
        self.__reduce__()
    def __repr__(self):
        return str(self.n) + '/' + str(self.d)
    def __reduce__(self): # Fix +/- and possibly reduce num. and denom. by gcd
        if self.d < 0:
            self.n = -self.n
            self.d = -self.d
        if self.n == 0:
            self.d = 1
            return
        g = _gcd(abs(self.n), abs(self.d))
        if g > 1:
            print(g)
            self.n //= g
            self.d //= g

    ### Operations
    def __mul__(self, y):
        ans = Fraction(self.n * y.n, self.d * y.d)
        ans.__reduce__()
        return ans
    def __add__(self, y):
        ans = Fraction(self.n*y.d + y.n*self.d, self.d * y.d)
        ans.__reduce__()
        return ans
    def __sub__(self, y):
        my = y.__neg__()
        return self.__add__(my)
    def __truediv__(self, y):
        return self * y.__rtruediv__(1)
    def __pow__(self, m):
        return Fraction(self.n ** m, self.d ** m)

    ### Comparisons
    def __eq__(self, y):
        if self.n == y.n and self.d == y.d:
            return True
        return False
    def __neq__ (self, y):
        return not self.__eq__(y)
    def __lt__(self, y):
        return self.n * y.d < y.n * self.d
    def __gt__(self, y):
        return self.n * y.d > y.n * self.d
    def __le__(self, y):
        return not self.__gt__(y)
    def __ge__(self, y):
        return not self.__lt__(y)
    def __neg__(self):
        return Fraction(-self.n, self.d)

    ### Allow integers as left operand for some operations:
    def __rtruediv__(self, m):
        return Fraction(m*self.d, self.n)
    def __rmul__(self, m):
        return Fraction(m*self.n, self.d)
    def __radd__(self, m):
        return Fraction(m*self.d + self.n, self.d)
    def __rsub__(self, m):
        return Fraction(m, 1).__add__(self.__neg__())
        
def _gcd(x,y):
    return _eeuclid(x,y)[2]
    
def _eeuclid(x,y):
    '''Extended Euclid algorithm for gcd and inverse.
    This implementation means the results need not be backtracked to find
    inverses and solution to linear Diophantine equations.
    Returns (a,b,g) such that ax + by = g = gcd(x,y)'''
    # Crandall and Pomerance, algo. 2.1.4

    assert x > 0
    assert y > 0
    if x > y:
        x,y=y,x
    olda = 1
    oldb = 0
    oldg = x
    u = 0
    v = 1
    w = y
    g = oldg
    while w > 0:
        q = g // w
        a = u
        b = v
        g = w
        u = olda - q*u
        v = oldb - q*v
        w = oldg - q*w

        olda = a
        oldb = b
        oldg = g
    return (a,b,g)

    
