from vector import Vector
from plane import Plane
from lin_sys import LinearSystem
from decimal import Decimal

p1 = Plane(Vector(['5.862','1.178','-10.366']),'-8.15')
p2 = Plane(Vector(['-2.391','-0.589','5.183']),'-4.075')

ls = LinearSystem([p1,p2]).compute_rref()
print(ls)
print('->')
print(ls.solve())


p1 = Plane(Vector(['8.631','5.112','-1.816']),'-5.113')
p2 = Plane(Vector(['4.315','11.132','-5.27']),'-6.775')
p3 = Plane(Vector(['-2.158','3.01','-1.727']),'-0.831')

ls = LinearSystem([p1,p2,p3]).compute_rref()
print(ls)
print('->')
print(ls.solve())

p1 = Plane(Vector(['5.262','2.739','-9.878']),'-3.441')
p2 = Plane(Vector(['5.111','6.358','7.638']),'-2.152')
p3 = Plane(Vector(['2.016','-9.92','-1.367']),'-9.278')
p4 = Plane(Vector(['2.167','-13.543','-18.883']),'-10.567')


ls = LinearSystem([p1,p2, p3,p4]).compute_rref()
print(ls)
print('->')
print(ls.solve())

