from decimal import Decimal, getcontext
from copy import deepcopy

from vector import Vector
from math_util import MyDecimal
from plane import Plane
from nltk.app.nemo_app import initialFind
from numba.cuda import initialize

getcontext().prec = 30

""""LinearSysten class documentation.

This class is a representation a system of first grade equations. Equations
are represented with the class plane.Plane. The equations cannot have more than
3 independent variables.

Examples:
    You create a system of 2 planes like following: 
    
    p1 = new plane.Plane(normal_vector= vector.Vector(['2','2','4'], constant_term='1')
    
    p2 = new plane.Plane(normal_vector= vector.Vector(['0','0','2'], constant_term='5')
    
    lin_sys = lin_sys.LinearSystem([p1,p2])
        
    
Attributes:
    planes (list[plane.Plane]): The list of planes in this system.
        
    dimension(int): The number of independent variables in the plane equations.
                    Always equals to 3.
                    
"""


class LinearSystem(object):

    ALL_PLANES_MUST_BE_IN_SAME_DIM_MSG = 'All planes in the system should live in the same dimension'
    NO_SOLUTIONS_MSG = 'No solutions'
    INF_SOLUTIONS_MSG = 'Infinitely many solutions'


    def __init__(self, planes):
        """Initialize the system. Saves the provided planes.
            
        """
        try:
            d = planes[0].dimension
            for p in planes:
                assert p.dimension == d

            self.planes = planes
            self.dimension = d

        except AssertionError:
            raise Exception(self.ALL_PLANES_MUST_BE_IN_SAME_DIM_MSG)


    #Provided by Udacity.
    def swap_rows(self, row1, row2):
        
        """Swap the position of 2 equations.
            
            Args:
                row1(int): The position of the first row to swap.
                row2(int): The position of the second row to swap.
        """
        self.planes[row1], self.planes[row2] = self.planes[row2] , self.planes[row1] 


    def multiply_coefficient_and_row(self, coefficient, row):
        """Multiplies a coefficient by each term of the equation in the row position.
        
            This method modifies the object on which is called.
            
            Args:
                coefficient(float): Coefficient to multiply.
                
                row(int): Index of the row to which the coefficient will be multiplied.         
            
        """
        normal_vector = self[row].normal_vector
        k = self[row].constant_term          
        normal_vector = normal_vector * coefficient;
        k = k * coefficient      
        self[row] = Plane(normal_vector,k)        
        


    def add_multiple_times_row_to_row(self, coefficient, row_to_add, row_to_be_added_to):
        """Multiplies a row by a coefficient and add the result to another row.
        
            Multiplies "coefficient" by the equation in "row_to_add" index and add 
            the result to equation in the "row_to_be_added_to" index. Only the 
            equation in the index "row_to_be_added_to" will be modified.
            
            Args:
                coefficient(float): The number that will be multiplied by each term of the equation
                                    in the index "row_to_add".
                                    
                row_to_add(int): The index of the equation that will be multiplied
                                 by "coefficient".
                                    
                row_to_be_added(int): The index of the equation that will be 
                                      modified.
            
            
        """
        planeBycoefficient = self[row_to_add] * coefficient
        self.planes[row_to_be_added_to] = self.planes[row_to_be_added_to] + planeBycoefficient        


    #Provided by Udacity
    def indices_of_first_nonzero_terms_in_each_row(self):
        """Returns a list with the the first non-zero index of the equations for each row.
                    
                    
            Returns:
                list[int]: The ith position of this list will have the index of 
                            the first non-zero coefficient in the ith equation. 
            
        """
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

    #Provided by Udacity
    def __len__(self):
        """Returns the number of equations in the system.
                    
                    
            Returns:
                int: Number of equations in the system. 
            
        """
        return len(self.planes)


    def __getitem__(self, i):
        """Returns the ith equation (plane) in this system.
                    
            Args:
                i(int): Index of the equation. 
                
            Returns:
                plane.Plane: Plane in the ith representing the ith equation of 
                             the system.
            Raises:
                 IndexError: If i is bigger than the number of equations in the
                             system.
        """
        return self.planes[i]

    #Provided by Udacity
    def __setitem__(self, i, x):
        
        """Set the ith equation (plane) in this system.
                    
            Args:
                i(int): Index of the plane. 
                
            
            Raises:
                 IndexError: If i is bigger than the number of equations in the
                             system.
                             
                Exception: If the new plane has more dimensions than the ones
                           already in the system.
        """
        try:
            assert x.dimension == self.dimension
            self.planes[i] = x

        except AssertionError:
            raise Exception(self.ALL_PLANES_MUST_BE_IN_SAME_DIM_MSG)


    def __str__(self):
        """Returns a string representation of this equations system.
             
            Returns:
                str: A string representation of this equations system.
             
        """
        
        ret = 'Linear System:\n'
        temp = ['Equation {}: {}'.format(i + 1, p) for i, p in enumerate(self.planes)]
        ret += '\n'.join(temp)
        return ret
    
    def organize_system(self):
        
        number_of_vars = self.dimension        
        non_zero_idx_per_row = self.indices_of_first_nonzero_terms_in_each_row()
        
        for row,_ in enumerate(self.planes):    
            for row2 in range(row+1,len(self.planes)):
                
                if(non_zero_idx_per_row[row] < 0 and non_zero_idx_per_row[row2] >= 0):
                    self.swap_rows(row, row2)
                    non_zero_idx_per_row = self.indices_of_first_nonzero_terms_in_each_row()
                    continue
                
                if(non_zero_idx_per_row[row2] < 0):                    
                    continue
                
                if(non_zero_idx_per_row[row2] < non_zero_idx_per_row[row]):
                    self.swap_rows(row, row2)
                    non_zero_idx_per_row = self.indices_of_first_nonzero_terms_in_each_row()
                    continue
                
                first_col_where_one_coefficient_is_zero_and_other_no = self.get_first_col_where_one_coefficient_is_zero_and_other_no(row,row2)
                if(first_col_where_one_coefficient_is_zero_and_other_no == None):
                    continue
                 
                row_col_coefficient = MyDecimal(self.planes[row].normal_vector[first_col_where_one_coefficient_is_zero_and_other_no])
                row2_col_coefficient = MyDecimal(self.planes[row2].normal_vector[first_col_where_one_coefficient_is_zero_and_other_no])
                
                if (row2_col_coefficient.is_near_zero() and not row_col_coefficient.is_near_zero() and non_zero_idx_per_row[row2] <= non_zero_idx_per_row[row]):
                    self.swap_rows(row,row2)
                    non_zero_idx_per_row = self.indices_of_first_nonzero_terms_in_each_row()
                    
    
    
    def get_first_col_where_one_coefficient_is_zero_and_other_no(self,row,row2):
                        
                
        for col in range(0,self.dimension):
            coefficient_row = 0 if MyDecimal(self.planes[row].normal_vector[col]).is_near_zero() else self.planes[row].normal_vector[col]  
            coefficient_row2 = 0 if MyDecimal(self.planes[row].normal_vector[col]).is_near_zero() else self.planes[row2].normal_vector[col]
            one_is_zero_but_not_both = coefficient_row * coefficient_row2 == 0 and (coefficient_row + coefficient_row2 != 0) 
            if(one_is_zero_but_not_both):
                return col
            
        return None
            


    
    def compute_triangular_form(self):
        """This function will return a different copy of this system in triangular form.
             
            Returns:
                LinearSystem: A new system equals to this one in triangular form.
             
        """
        
        system = deepcopy(self)
        system.organize_system()
                
        elimination_info = {'system' : system}
        
        elimination_info['initial_idx'] = 0        
        elimination_info['final_idx'] = len(self.planes)
        
        self.elimination(elimination_info)
        
        """elimination_info['initial_idx'] = len(self.planes)  -1      
        elimination_info['final_idx'] = -1
        
        self.elimination(elimination_info)"""
        system.organize_system();
        return system
        
        
    """This function will eliminate all possible variables in the direction going from elimination_info['starting_index'] to elimination_info['final_idx']"""
    def elimination(self,elimination_info):
        
        system = elimination_info['system']
        initial_idx = elimination_info['initial_idx']
        final_idx = elimination_info['final_idx']
                
        #first_non_zeros = system.indices_of_first_nonzero_terms_in_each_row()
        steps =  1 if initial_idx < final_idx else -1  
                
        
        for i in range(initial_idx,final_idx,steps):
            for j in range(i+steps,final_idx,steps):                
                coefficient = system.find_coefficient(system.planes[i],system.planes[j]);
                system.add_multiple_times_row_to_row(coefficient,i,j);
            
            system.organize_system()
                
    
    
    def find_coefficient(self,eq1,eq2):
        
        """Find coefficient to eliminate first non-zero coefficient in eq2. 
        
            Find coefficient that if multiplied by eq1 will produce an equation, that added to eq2, will eliminate eq2's first non-zero coefficient.
            
            Returns:
                LinearSystem: A new system equals to this one in triangular form.
             
        """
        try:            
            first_nonzero = eq1.first_nonzero_index(eq1.normal_vector)        
            first_nonzero_is_shared = True if eq2.first_nonzero_index(eq2.normal_vector) == first_nonzero else False
            coefficient = 0
        
            if(first_nonzero_is_shared):
                eq1_coefficient = eq1.normal_vector[first_nonzero]
                eq2_coefficient = eq2.normal_vector[first_nonzero]
                coefficient = eq2_coefficient/eq1_coefficient
                if(eq1_coefficient * coefficient != eq2_coefficient*-1):
                    coefficient *= -1             
                
            else:
                coefficient = 0
            
            return coefficient
        except Exception:
            return 0
        
            