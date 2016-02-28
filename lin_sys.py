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
        n = system.dimension
        j= 0
        for i,_ in enumerate(system.planes):
            while j < n:
                c = MyDecimal(system.planes[i].normal_vector[j])
                if c.is_near_zero():
                    row_with_nonzero = system.find_idx_with_nonzero(j, i)
                    if row_with_nonzero:
                        system.swap_rows(i, row_with_nonzero)
                    else:
                        j += 1
                        continue;
                        
                system.clear_var(i,j)
                
                j += 1
                break;
                
            
        return system
    
    def solve(self):
        """Returns the solution of this system of equation.
        
        Returns:
            vector.Vector: If there is an unique solution, will return a vector
                representing the x,y,z points of the solution.
            
            bool: False if there is no solution and True if there are many 
                solutions.
        """
        rref = self.compute_rref()
        unique_solution = ['0','0','0']
        first_nonzeros = rref.indices_of_first_nonzero_terms_in_each_row()
        response = None
        one_variable_alone = False
        single_solution = True
        
        for i,plane in enumerate(rref.planes):
            pivot_var_idx = first_nonzeros[i]    
            constant_term = MyDecimal(plane.constant_term)
            if(pivot_var_idx < 0 and not constant_term.is_near_zero() ):
                return False
            
            number_of_vars = rref[i].var_count()
            if(number_of_vars == 1):
                one_variable_alone = True
                unique_solution[pivot_var_idx] = plane.constant_term/plane.normal_vector[pivot_var_idx]
            
            elif(number_of_vars > 1):
                single_solution = False
        
        
       
        if(one_variable_alone and not single_solution):
            response = False
        elif(not one_variable_alone and not single_solution):
            response = True
        else:
            response = unique_solution
        
        return response 
    
    def compute_rref(self):
        """Returns a copy of this system in Reduced Row-Echelon Form:
        
            Reduced Row-Echelon means that if possible every variable will
            have a coefficient of 1 and will be the only variable in each 
            equation.
            
            Returns:
                lin_sys.LinearSystem: Returns a copy of this system in Reduced 
                                      Row-Echelon Form.
        """
        rref = self.compute_triangular_form()
        
        start_idx = len(rref.planes) -1
        end_idx = 0 -1
        steps = -1
        
        first_nonzero_by_row = rref.indices_of_first_nonzero_terms_in_each_row()
        
        for i in range(start_idx, end_idx, steps):
            first_nonzero = first_nonzero_by_row[i]
            has_no_nonzero = first_nonzero < 0 
            if(has_no_nonzero):
                continue;
            rref.coef_to_one(i,first_nonzero)
            rref.remove_var_above(first_nonzero,i)       
        return rref
        
    def remove_var_above(self,var_idx,row):
        """Remove variable in "var_idx" index from all rows above "row"
        
        Args:
            var_idx(int): Index of the variable to eliminate.
            
            row(int): The index of the equation that has the variable as a pivot
                variable.
        
        """    
        for i in range(row-1, -1, -1):                
            coefficient = self.find_coefficient(self.planes[row], self.planes[i], var_idx)
            self.add_multiple_times_row_to_row(coefficient, row, i)
            
        
        
    def coef_to_one(self,row_idx,variable_idx):
        """Puts a coeficient of 1 to the variable in "variable_idx"  index  for the equation in "row_idx" index.
        Args:
            row_idx(int): Index of the equation to modify.
            
            variable_idx(int): The index of the variable. 
        
        """
        coefficient_inverse = Decimal('1.0')/self.planes[row_idx].normal_vector[variable_idx]
        self.multiply_coefficient_and_row(coefficient_inverse,row_idx) 
           
    def find_idx_with_nonzero(self, coeficient_idx, start_idx=0):
        """Find the first equation index where coeficient "coeficient_idx" is not zero.
            
            Args:
                coeficient_idx(int): The index of the coefficient to check.
                
                start_idx(int): Start the search in this row of the system.
            
            Returns:
                bool: Returns False if there is not an equation in this system 
                    with coefficient "coeficient_idx" different than zero.  
                
                int: Returns the index of the equation with the first non-zero 
                     coefficient in variable "coeficient_idx". Default value is
                     zero.
            Raises:
                IndexError: If start_idx is out of the bounds of self.planes array.
                        
                
        """
        for i in range(start_idx, len(self.planes)):
            row_var_value = self.planes[i].normal_vector[coeficient_idx]
            row_var_value = MyDecimal(row_var_value)
            if(not row_var_value.is_near_zero()):
                return i
        
        
        return False
        
    
    def clear_var(self,row,col):
        """Eliminate the variable in column "col" from equations bellow "row"
        
            Args:
                row(int): All equations bellow "row" will get their variable in
                          index "col" cleared.
                
                col(int): The index of the variable to ble cleared.
            
            Raises:
                IndexError: If "row" or "col" params areout of the bounds 
                            of self.planes.            
        """
        for i in range(row+1,len(self.planes)):                      
            coefficient = self.find_coefficient(self.planes[row],self.planes[i], col);
            self.add_multiple_times_row_to_row(coefficient,row,i);
                                                   
    
    
    def find_coefficient(self,eq1,eq2,col):
        
        """Find coefficient to eliminate the "col" variable in eq2. 
        
            Find coefficient that if multiplied by eq1 will produce an equation, 
            that added to eq2, will eliminate eq2's variable in index "col".
            
            Args: 
                eq1(plane.Plane): The first plane.
                
                eq2(plane.Plane): Secod plane.
                
                col(int): The column from which to calculate the coefficient.
                
            Returns:
                int: The coefficient that if multiplied by eq1 will produce an 
                     equation, that added to eq2, will eliminate eq2's variable 
                     in index "col". Zero will be returned if one of the two
                     variables in pos "col" is zero.
                     
             
        """
        try:            
                               
            coefficient = 0
        
            eq1_coefficient = eq1.normal_vector[col]
            eq2_coefficient = eq2.normal_vector[col]
            coefficient = -eq2_coefficient/eq1_coefficient
            
            return coefficient
        except Exception:
            return 0
        
            