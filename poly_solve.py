def solve2(b,c):
    '''Returns real solutions to x^2 + b.x + c = 0 using a numerically stable
    version of the quadratic formula (see Numerical Recipes in C, pp. 183-4)'''
    import math
    # we all know that when a is one, x = -b/2 +/- sqrt(b^2-4c)/2
    det = b**2 - 4*c
    if det < 0:
        return []
    elif det == 0:
        return [-b/2,]
    else:
        qt = math.copysign(math.sqrt(det), b) # i.e. sgn(b) * sqrt(det)
        q = -(b+qt)/2
        return [q, c/q]
    
def solve3(b,c,d, ZEROTOL=10**-8):
    '''Returns real solutions to the cubic equation x^3 + b.x^2 + c.x + d = 0.'''
    import pdb
    import math
    def feq0(afloat, TOL=ZEROTOL):
        if abs(afloat) < TOL:
            return True
        return False
    def cuberoot(areal):
        tmp = pow(abs(areal), 1/3)
        if areal < 0:
            return -tmp
        return tmp
    def xfromy(ys,b):
        res = []
        for y in ys:
            res.append(y-b/3)
        return res
    # First reduce to the depressed cubic:
    # y^3 + Ay = B
    A = (3*c-b**2)/3 # -3Q
    B = (9*b*c - 2*b**3 - 27*d) / 27 # -2R
    print('y**3 + %2.f y = %.2f' % (A,B))
    
    # And the associated quadratic in u=t^3:
    # u^2 + Bu + C = 0
    C = -(A**3)/27 # Q^3
    print('u^2 + %.2f u + %.2f = 0' % (B,C))

    # determinant of this is B^2 - 4C
    det = B**2 - 4*C # -4R^2-4Q^3

    print(A,B,C,det)
    if feq0(det):# == 0:
        print('det == 0, two real solutions')
        t3 = -B/2
        t = cuberoot(t3)
        s3 = B/2
        s = cuberoot(s3)
        y = s - t
        return xfromy((y,t), b)
    elif det < 0: # R^2 + Q^3 > 0
        print('det < 0, three real solutions')
        #Q = (b**2 - 3*c)/9
        Q = A / -3
        #R = (2*b**3 - 9*b*c + 27*d)/54
        R = B / -2
        theta = math.acos(R/math.sqrt(Q**3)) # should be valid here
        j = -2*math.sqrt(Q)
        return xfromy((j*math.cos(theta/3),
                       j*math.cos((theta + 2*math.pi)/3),
                       j*math.cos((theta - 2*math.pi)/3)), b)
    else:
        print('det > 0, one real solution')
        print(solve2(B,C))
        u = solve2(B,C)[0] # Just take one solution as it's the same
        print(u)
        t = cuberoot(u)
        print(t)
        s = cuberoot(B+u)
        print(s)
        return xfromy((s-t,),b)
