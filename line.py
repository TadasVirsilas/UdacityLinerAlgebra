from decimal import Decimal, getcontext

from vector import Vector
from math_util import MyDecimal

getcontext().prec = 30

""""Line class documentation.

This class is a representation of a 2D line. The line is represented by a normal
vector and basepoint. The normal vector to the line and the constant term of the line's equation
are asked and with them the basepoint is calculated.


Examples:
    The line equation is Ax + By = K. K is the constant term and the 
    normal vector has the coordinates (A,B). 
    
    To instantiate:: 
        line = Line(['1.6','2','3'],5)
        
    
Attributes:
    normal_vector (Vector): The normal vector of the line.
        
    constant_term(Decimal): The constant term of the line's equation.
        
    basepoint(Vector): The basepoint of the line. It is calculated based on the
    normal_vector and the constant_term. 
    
"""


class Line(object):

    NO_NONZERO_ELTS_FOUND_MSG = 'No nonzero elements found'

    #Provided by Udacity
    def __init__(self, normal_vector=None, constant_term=None):
        """Constructor for Line class.
        
        Args:
            normal_vector(vector.Vector): The line's normal vector.
            constant_term(Decimal): The constant_term of the line's equation.
                
            
        """
        self.dimension = 2

        if not normal_vector:
            all_zeros = ['0']*self.dimension
            normal_vector = Vector(all_zeros)
        self.normal_vector = normal_vector

        if not constant_term:
            constant_term = Decimal('0')
        self.constant_term = Decimal(constant_term)

        self.set_basepoint()


    def is_parallel_to(self,l2):
        """Check if self and l2 are parallel.
            
        Args:
            l2(line.Line): We will check if self is parallel to this line. 
        
        Returns:
            bool: True if self and l2 are parallel lines. False otherwise.
            
        Raises: 
            ValueError: If self and v normal vectors doesn't have the same dimensions.
        """
        return self.normal_vector.is_parallel_to(l2.normal_vector)
    
    def is_same_as(self,l2):
        """Check if self and l2 are the same.
            
        Args:
            l2(line.Line): We will check if self and this line are the same. 
        
        Returns:
            bool: True if self and l2 are the same line. False otherwise.
            
        Raises: 
            ValueError: If self and l2 normal vectors doesn't have the same dimensions.
        """
        if(self.normal_vector.is_zero() and l2.normal_vector.is_zero()):
            diff = self.constant_term - l2.constant_term
            return MyDecimal(diff).is_near_zero()
        
        atLeastOneZero =  self.normal_vector.is_zero() or l2.normal_vector.is_zero()
        
        if(atLeastOneZero):
            return False
                 
        areParallel = self.is_parallel_to(l2)
        if(not areParallel):
            return False
        vectorBetweenLines = self.basepoint - l2.basepoint
        return vectorBetweenLines.is_orthogonal_to(self.normal_vector) and vectorBetweenLines.is_orthogonal_to(l2.normal_vector)
        
    def get_intersection_with(self,l2):
        """Get intersection point between self and l2.
        
        Args:
            l2(line.Line): The that will be used to find the intersection with. 
        
        Returns:
            list[Decimal]: Containing x and y coordinates.
            None: If no intersection in common.
            
        Raises: 
            ValueError: If self and v normal vectors doesn't have the same dimensions.        
        """
        areParallels = self.is_parallel_to(l2)
        areSame = self.is_same_as(l2)
        intersection = None
        
        if(areParallels and not areSame):
            intersection = None            
        elif(areSame):
            intersection = self        
        else:
            try :
                A = self.normal_vector[0]
                B = self.normal_vector[1]
                C=l2.normal_vector[0]
                D = l2.normal_vector[1]
                k1 = self.constant_term
                k2 = l2.constant_term
                x = (D*k1 - B*k2) / (A*D-B*C)
                y = (-C*k1 + A*k2) / (A*D-B*C)
                intersection = [x,y]
            except ZeroDivisionError:
                intersection = None
             
        return intersection
    
    
    #Provided by Udacity
    def set_basepoint(self):
        """Calculate and set the basepoint for this line.
        
        """
        try:
            n = self.normal_vector
            c = self.constant_term
            basepoint_coords = ['0']*self.dimension

            initial_index = Line.first_nonzero_index(n)
            initial_coefficient = n[initial_index]

            basepoint_coords[initial_index] = c/initial_coefficient
            self.basepoint = Vector(basepoint_coords)

        except Exception as e:
            if str(e) == Line.NO_NONZERO_ELTS_FOUND_MSG:
                self.basepoint = None
            else:
                raise e

    #Provided by Udacity
    def __str__(self):
        """Returns a string representing the equation of this class.
        
        Returns:
            string: The equation of this line.
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
            initial_index = Line.first_nonzero_index(n)
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
        """Returns the index of the first non-zero value for the iterable.
        
        Returns:
            int: The index first nonzero value of the iterable.
            
        Raises:
            Exception: If all values in the iterable are zero. 
        
        """
        for k, item in enumerate(iterable):
            if not MyDecimal(item).is_near_zero():
                return k
        raise Exception(Line.NO_NONZERO_ELTS_FOUND_MSG)
