import types
from decimal import Decimal, getcontext
from vector import Vector
from copy import deepcopy
from math_util import MyDecimal

getcontext().prec = 30



""""Plane class documentation.

This class is a representation of a 3D plane. The plane is represented by a normal
vector of 3 dimensions and basepoint. The normal vector to the line and the 
constant term of the plane's equation are asked and with them the basepoint is calculated.


Examples:
    The equation of a plane is Ax + By + Cz = K. K is the constant term and the 
    normal vector has the coordinates (A,B). 
    
    To instantiate:: 
        line = Line(['1.6','2','3'],5)
        
    
Attributes:
    normal_vector (Vector): The normal vector of the line.
        
    constant_term(Decimal): The constant term of the line's equation.
        
    basepoint(Vector): The basepoint of the line. It is calculated based on the
    normal_vector and the constant_term. 
    
"""



class Plane(object):

    NO_NONZERO_ELTS_FOUND_MSG = 'No nonzero elements found'

    def __init__(self, normal_vector=None, constant_term=None):
        """Plane class constructor.
        
            Args:
                normal_vector(vector.Vector): The normal vector of the plane.
                
                constant_term(Decimal): The constant term of the plane's equation.
            
        """
        self.dimension = 3

        if not normal_vector:
            all_zeros = ['0']*self.dimension
            normal_vector = Vector(all_zeros)
        self.normal_vector = normal_vector

        if not constant_term:
            constant_term = Decimal('0')
        self.constant_term = Decimal(constant_term)

        self.set_basepoint()

    def is_parallel_to(self,p2):
        """Check if self and p2 are parallel planes.
        
            Args:
                p2(plane.Plane): The plane to check against.
            
            Returns: 
                bool: True if self and p2 are parallel, false otherwise.
                    
            Raises:
                ValueError: If self and p2 normal vectors doesn't have the same 
                dimensions.
        """
        return self.normal_vector.is_parallel_to(p2.normal_vector)
    
    def __eq__(self, other):
        return self.is_same_as(other)
    
    def is_same_as(self,p2):
        """Check if self and p2 are the same plane.
        
            Args:
                p2(plane.Plane): The plane to check against.
            
            Returns: 
                bool: True if self and p2 are the same plane, false otherwise.
                    
            Raises:
                ValueError: If self and p2 normal vectors doesn't have the same 
                dimensions.
        """
        if(self.normal_vector.is_zero() and p2.normal_vector.is_zero()):
            diff = self.constant_term - p2.constant_term
            return MyDecimal(diff).is_near_zero()
        
        atLeastOneZero =  self.normal_vector.is_zero() or p2.normal_vector.is_zero()
        
        if(atLeastOneZero):
            return False
                 
        areParallel = self.is_parallel_to(p2)
        if(not areParallel):
            return False
        vectorBetweenLines = self.basepoint - p2.basepoint
        return vectorBetweenLines.is_orthogonal_to(self.normal_vector) and vectorBetweenLines.is_orthogonal_to(p2.normal_vector)

    def var_count(self):
        """Returns a count of the number of variables present in the equation.
        
        Returns:
            int: Count of variables in this plane's equation.
        
        """
        var_count = 0;
        for _, coefficient in enumerate(self.normal_vector):
            var_coefficient = MyDecimal(coefficient)
            if(not var_coefficient.is_near_zero()):
                var_count += 1
                
        return var_count
        
    def __add__(self,operand):
        """Add a number or a plane to this plane equation. The result is a new plane.
        
            Args:
                operand(plane.Plane): Another plane to be added to this one.
                operand(int): An int to add to this plane's equation.
                operand(float): A float to add to this plane's equation.
                operand(Decimal): A Decimal to add to this plane's equation.
                
            Returns:
                plane.Plane: A plane equals to this plane plus the operand.
                
            Raises:
                TypeError:  Thrown if the operand param is not a valid type.
            
        
        """
        response = deepcopy(self)
        
        
        if(isinstance(operand, Plane)):
            response.normal_vector = response.normal_vector + operand.normal_vector
            response.constant_term =    response.constant_term + operand.constant_term
            

        elif(isinstance(operand,types.numeric_types)):
            response.normal_vector = response.normal_vector + operand
            response.constant_term =    response.constant_term + operand
        else:
            raise TypeError("You can only add numbers and planes.")
        
        response.set_basepoint()
        return response
            
    #Provided by Udacity.
    def set_basepoint(self):
        """Calculates and set a basepoint based on the normal_vector and the constant_term.
            
        """
        try:
            n = self.normal_vector
            c = self.constant_term
            basepoint_coords = ['0']*self.dimension

            initial_index = Plane.first_nonzero_index(n)
            initial_coefficient = n[initial_index]

            basepoint_coords[initial_index] = c/initial_coefficient
            self.basepoint = Vector(basepoint_coords)

        except Exception as e:
            if str(e) == Plane.NO_NONZERO_ELTS_FOUND_MSG:
                self.basepoint = None
            else:
                raise e

    #Provided by Udacity
    def __str__(self):
        """Represent this plane as an equation.
        """
        num_decimal_places = 3

        def write_coefficient(coefficient, is_initial_term=False):
            coefficient = round(coefficient, num_decimal_places)
            if coefficient % 1 == 0:
                coefficient = int(coefficient)

            output = ''

            if coefficient < 0:
                output += '-'
            if coefficient > 0 and not is_initial_term:
                output += '+'

            if not is_initial_term:
                output += ' '

            if abs(coefficient) != 1:
                output += '{}'.format(abs(coefficient))

            return output

        n = self.normal_vector

        try:
            initial_index = Plane.first_nonzero_index(n)
            terms = [write_coefficient(n[i], is_initial_term=(i==initial_index)) + 'x_{}'.format(i+1)
                     for i in range(self.dimension) if round(n[i], num_decimal_places) != 0]
            output = ' '.join(terms)

        except Exception as e:
            if str(e) == self.NO_NONZERO_ELTS_FOUND_MSG:
                output = '0'
            else:
                raise e

        constant = round(self.constant_term, num_decimal_places)
        if constant % 1 == 0:
            constant = int(constant)
        output += ' = {}'.format(constant)

        return output


    #Provided by Udacity
    @staticmethod
    def first_nonzero_index(iterable):
        
        """Returns the index of the first non-zero value in iterable. 
        
            Args:
                iterable(iterable): The iterable to search in.
                
            Returns:
                int: The index of the first non-zero value of iterable.
                
            Raises:
                Exception: If all items are zero.
        """
        for k, item in enumerate(iterable):
            if not MyDecimal(item).is_near_zero():
                return k
        raise Exception(Plane.NO_NONZERO_ELTS_FOUND_MSG)

    def __mul__(self, coefficient):
        response = deepcopy(self)        
        response.normal_vector = response.normal_vector  * coefficient
        response.constant_term = response.constant_term * coefficient
        response.set_basepoint()
        return response