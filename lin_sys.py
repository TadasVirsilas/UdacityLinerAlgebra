from decimal import Decimal, getcontext
from copy import deepcopy

from vector import Vector
from plane import Plane
from nltk.app.nemo_app import initialFind
from numba.cuda import initialize

getcontext().prec = 30


class LinearSystem(object):

    ALL_PLANES_MUST_BE_IN_SAME_DIM_MSG = 'All planes in the system should live in the same dimension'
    NO_SOLUTIONS_MSG = 'No solutions'
    INF_SOLUTIONS_MSG = 'Infinitely many solutions'

    def __init__(self, planes):
        try:
            d = planes[0].dimension
            for p in planes:
                assert p.dimension == d

            self.planes = planes
            self.dimension = d

        except AssertionError:
            raise Exception(self.ALL_PLANES_MUST_BE_IN_SAME_DIM_MSG)


    def swap_rows(self, row1, row2):
        self.planes[row1], self.planes[row2] = self.planes[row2] , self.planes[row1] 


    def multiply_coefficient_and_row(self, coefficient, row):
        normal_vector = self[row].normal_vector
        k = self[row].constant_term          
        normal_vector = normal_vector * coefficient;
        k = k * coefficient      
        self[row] = Plane(normal_vector,k)        
        


    def add_multiple_times_row_to_row(self, coefficient, row_to_add, row_to_be_added_to):
        planeByCoeficient = self[row_to_add] * coefficient
        self.planes[row_to_be_added_to] = self.planes[row_to_be_added_to] + planeByCoeficient


    def indices_of_first_nonzero_terms_in_each_row(self):
        num_equations = len(self)
        num_variables = self.dimension

        indices = [-1] * num_equations

        for i, p in enumerate(self.planes):
            try:
                indices[i] = p.first_nonzero_index(p.normal_vector)
            except Exception as e:
                if str(e) == Plane.NO_NONZERO_ELTS_FOUND_MSG:
                    continue
                else:
                    raise e

        return indices


    def __len__(self):
        return len(self.planes)


    def __getitem__(self, i):
        return self.planes[i]


    def __setitem__(self, i, x):
        try:
            assert x.dimension == self.dimension
            self.planes[i] = x

        except AssertionError:
            raise Exception(self.ALL_PLANES_MUST_BE_IN_SAME_DIM_MSG)


    def __str__(self):
        ret = 'Linear System:\n'
        temp = ['Equation {}: {}'.format(i + 1, p) for i, p in enumerate(self.planes)]
        ret += '\n'.join(temp)
        return ret
    
    def organize_system(self):
        
        number_of_vars = self.dimension         
        
        
        for row,_ in enumerate(self.planes):            
            for col in range(0,number_of_vars):
                colVal= MyDecimal(self.planes[row].normal_vector[col])
                if (colVal.is_near_zero()):
                    self.swap_with_first_nonzero(row,col)        
    
    
    def swap_with_first_nonzero(self,index,col):
        
        for i in range(index,len(self)):
            colVal = MyDecimal(self[i].normal_vector[col])
            if(not colVal.is_near_zero()):
                self.swap_rows(index, i) 


    """ This function will return a different copy of this system in triangular form"""
    def compute_triangular_form(self):
        system = deepcopy(self)
        system.organize_system()
                
        elimination_info = {'system' : system}
        
        elimination_info['initial_idx'] = 0        
        elimination_info['final_idx'] = len(self.planes)  - 1
        
        self.elimination(elimination_info)
        
        elimination_info['initial_idx'] = len(self.planes)  -1      
        elimination_info['final_idx'] = 0 
        
        self.elimination(elimination_info)
        return system
        
        
    """This function will eliminate all possible variables in the direction going from elimination_info['starting_index'] to elimination_info['final_idx']"""
    def elimination(self,elimination_info):
        
        system = elimination_info['system']
        initial_idx = elimination_info['initial_idx']
        final_idx = elimination_info['final_idx']
                
        
        steps =  1 if initial_idx < final_idx else -1  
                
        
        for i in range(initial_idx,final_idx,steps):
            for j in range(i+steps,final_idx,steps):
                coeficient = system.find_coeficient(system.planes[i],system.planes[j]);
                self.add_multiple_times_row_to_row(coeficient,i,j);
                
    
    """This function will find the coeficient that needs to be multplied by the first equation so that the sum of the 2 eq. will eliminate the first non-zero variable in common"""
    def find_coeficient(self,eq1,eq2):
        first_nonzero = eq1.first_nonzero_index(eq1.normal_vector)
        first_nonzero_is_shared = True if eq2.normal_vector[first_nonzero] > 0 else False
        coeficient = 0
        
        if(first_nonzero_is_shared):
            eq1_coeficient = eq1.normal_vector[first_nonzero]
            eq2_coeficient = eq2.normal_vector[first_nonzero]
            coeficient = eq1_coeficient/eq2_coeficient
            if(eq1_coeficient > 0 and eq2_coeficient > 0):
                coeficient *= -1             
            
        else:
            coeficient = 0
            
        return coeficient
        
            
class MyDecimal(Decimal):
    def is_near_zero(self, eps=1e-10):
        return abs(self) < eps


p0 = Plane(normal_vector=Vector(['1', '1', '1']), constant_term='1')
p1 = Plane(normal_vector=Vector(['0', '1', '0']), constant_term='2')
p2 = Plane(normal_vector=Vector(['1', '1', '-1']), constant_term='3')
p3 = Plane(normal_vector=Vector(['1', '0', '-2']), constant_term='2')

s = LinearSystem([p0, p1, p2, p3])

print(s.indices_of_first_nonzero_terms_in_each_row())
print('{},{},{},{}'.format(s[0], s[1], s[2], s[3]))
print(len(s))
print(s)

s[0] = p1
print(s)

print (MyDecimal('1e-9').is_near_zero())
print (MyDecimal('1e-11').is_near_zero())
