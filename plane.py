import types
from decimal import Decimal, getcontext
from vector import Vector
from copy import deepcopy
from math_util import MyDecimal

getcontext().prec = 30


class Plane(object):

    NO_NONZERO_ELTS_FOUND_MSG = 'No nonzero elements found'

    def __init__(self, normal_vector=None, constant_term=None):
        self.dimension = 3

        if not normal_vector:
            all_zeros = ['0']*self.dimension
            normal_vector = Vector(all_zeros)
        self.normal_vector = normal_vector

        if not constant_term:
            constant_term = Decimal('0')
        self.constant_term = Decimal(constant_term)

        self.set_basepoint()

    def is_parallel_to(self,l2):
        return self.normal_vector.is_parallel_to(l2.normal_vector)
    
    def __eq__(self, other):
        return self.is_same_as(other)
    
    def is_same_as(self,l2):
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

    def __add__(self,operand):
        
        response = deepcopy(self)
        
        
        if(isinstance(operand, Plane)):
            response.normal_vector = response.normal_vector + operand.normal_vector
            response.constant_term =    response.constant_term + operand.constant_term
            

        elif(isinstance(operand,types.numeric_types)):
            response.normal_vector = response.normal_vector + operand
            response.constant_term =    response.constant_term + operand
        else:
            raise Exception("You can only add numbers and planes.")
        
        response.set_basepoint()
        return response
            
        
    def set_basepoint(self):
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


    def __str__(self):

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


    
    @staticmethod
    def first_nonzero_index(iterable):
        for k, item in enumerate(iterable):
            if not MyDecimal(item).is_near_zero():
                return k
        raise Exception(Plane.NO_NONZERO_ELTS_FOUND_MSG)

    def __mul__(self, coeficient):
        response = deepcopy(self)        
        response.normal_vector = response.normal_vector  * coeficient
        response.constant_term = response.constant_term * coeficient
        response.set_basepoint()
        return response