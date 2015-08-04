# -*- coding: utf-8 -*-
"""
Created on Wed Apr 17 16:25:40 2013

@author: jmccormick
"""

def polyeval(a, x):
    """
    p(x) = polyeval(a, x)
         = a[0] + a[1]x + a[2]x^2 + ....
    """
    p = 0
    a.reverse()
    for coef in a:
        p = p * x + coef
    a.reverse()
    return p
    
def polyderiv(a):
    """
    p'(x) = polyderiv(a)
          = b[0] + b[1]x + b[2]x^2 + .......
    """
    b = []
    for i in range(1, len(a)):
        b.append(i * a[i])
    return b
    
def polyreduce(a, root):
    """
    Given x = r is a root of n'th degree polynomial p(x) = (x-r)q(x),
    divide p(x) by linear factor (x-r) using the same algorithm as
    polynomial evaluation.  Then, return the (n-1)'th degree quotient 
    q(x) = polyreduce(a, r)
         = c[0] + c[1]x + c[2]x^2 +...+ c[n-2]x^{n-2} + c[n-1]x^{n-1}
    """
    c, p = [], 0
    a.reverse()
    for coef in a:
    	p = p * root + coef
	c.append(p)
    a.reverse()
    c.reverse()
    return c[1:]
    
def quadratic(a, b, c=None):
    """
    x^2 + ax + b = 0  (or ax^2 + bx + c = 0)
    By substituting x = y-t and t = a/2, the equation reduces to 
        y^2 + (b-t^2) = 0 
    which has easy solution
        y = +/- sqrt(t^2-b)
    """
    import math, cmath
    if c:		# (ax^2 + bx + c = 0)
	a, b = b / float(a), c / float(a)
    t = a / 2.0
    r = t**2 - b
    if r >= 0:		# real roots
	y1 = math.sqrt(r)
    else:		# complex roots
	y1 = cmath.sqrt(r)
    y2 = -y1
    return y1 - t, y2 - t
    
def cubic(a, b, c, d=None):
    """
    x^3 + ax^2 + bx + c = 0  (or ax^3 + bx^2 + cx + d = 0)
    With substitution x = y-t and t = a/3, the cubic equation reduces to    
        y^3 + py + q = 0,
    where p = b-3t^2 and q = c-bt+2t^3.  Then, one real root y1 = u+v can
    be determined by solving 
        w^2 + qw - (p/3)^3 = 0
    where w = u^3, v^3.  From Vieta's theorem,
        y1 + y2 + y3 = 0
        y1 y2 + y1 y3 + y2 y3 = p
        y1 y2 y3 = -q,
    the other two (real or complex) roots can be obtained by solving
        y^2 + (y1)y + (p+y1^2) = 0
    """
    from math import cos
    if d:			# (ax^3 + bx^2 + cx + d = 0)
	a, b, c = b / float(a), c / float(a), d / float(a)
    t = a / 3.0
    p, q = b - 3 * t**2, c - b * t + 2 * t**3
    u, v = quadratic(q, -(p/3.0)**3)
    if type(u) == type(0j):	# complex cubic root
	r, w = polar(u.real, u.imag)
	y1 = 2 * cbrt(r) * cos(w / 3.0)
    else:			# real root
        y1 = cbrt(u) + cbrt(v)
    y2, y3 = quadratic(y1, p + y1**2)
    return y1 - t, y2 - t, y3 - t