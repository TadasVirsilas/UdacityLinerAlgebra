
from decimal import Decimal

def isclose(a, b, rel_tol=0.0, abs_tol=1e-9):
    a = float(a)
    b= float(b)
    return abs(a-b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)
       
    
class MyDecimal(Decimal):
    def is_near_zero(self, eps=1e-10):
        return abs(self) < eps