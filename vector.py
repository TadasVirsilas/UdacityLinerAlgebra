
"""Vector class documentation.

This class is a representation of a vector as described by lineal algebra.

Example on how to instantiate:
    vector = Vector(['1.6','2','3'])
    vector = Vector([1.6,2,3])
    
Attributes:
    coodinates (tuple): Contains the x,y,z... coordinates of the vector in the
        form of Decimals.
        
    dimension(int): Especifies in how many dimensions (2D,3D,4D) is this vector 
        represented (how many coordinates. If the vector has x, y and z
        coordinates then dimension would be 3). 
"""

import math
import math_util
from math import acos
from decimal import Decimal


class Vector(object):
    
    
    def __iter__(self):
        """Return the iterator for this vector's coordinates.
    
            Returns:
                iteratr: The iterator for the coordinates of this vector. 
        """
        return iter(self.coordinates)
    
    
    def __getitem__(self, key):
        """Return coordinates in key position.
        
            Arguments: 
                key(int): index of the coordinate to return.
                
            Returns:
                Decimal: The coordinate value in the index represented by the parameter
                key.
            
            Raises:
                IndexError: When key is greater than the coordinates length.
                
        
        """
        if(key > len(self.coordinates) - 1):
            raise IndexError
        return self.coordinates[key]
        
    #Provided by Udacity       
    def __init__(self, coordinates):
        """
            Args:
                coordinates (list): A list of numbers representing the coordinates (x,y,z,...).
                
            Example:
                #Instantiating a vector.
                vector = Vector(['1.6','2','3'])
                vector = Vector([1.6,2,3]).

        """
        self._current = -1
        try:
            if not coordinates:
                raise ValueError
            for i, _ in enumerate(coordinates):
                coordinates[i] = Decimal(coordinates[i])
                
            self.coordinates = tuple(coordinates)
            self.dimension = len(coordinates)
        except ValueError:
            raise ValueError('The coordinates must be nonempty')

        except TypeError:
            raise TypeError('The coordinates must be an iterable')


    def __str__(self):
        return 'Vector: {}'.format(self.coordinates)


    def __eq__(self, v):
        
        return self.coordinates == v.coordinates

    def __add__(self, v):        
        response = []
        for i, value in enumerate(self.coordinates):
            response.append(v.coordinates[i] + value)
        return Vector(response)

    def __sub__(self, v):
        neg = v * -1
        return self +neg

    def __mul__(self, v):
        response = list(self.coordinates)
        for i, _ in enumerate(response):
            response[i] = self.coordinates[i] * v;
        return Vector(response)

    def dot(self, v):
        """Returns the dot product between this instance and another vector.
        
            Parameters:
                v(Vector): The other vector to perform the dot product with.
                
            Returns:
                Decimal: The dot vector
            
            Raises:
                ValueError: If the two vectors doesn't have the same length.
            
        """
        if(len(self.coordinates) != len(v.coordinates)):
            raise ValueError("Vectors should have same length")
        result = 0
        for i, _ in enumerate(self.coordinates):
            result += self.coordinates[i] * v.coordinates[i]
        return result

   

    def get_unit_vector(self):
        """Returns a new instance with the the value of this instance's unit vector.
            
            Returns:
                vector.Vector: The module vector of this instance.
                
            Raises:
                ZeroDivisionError: Thrown if the module of self is zero.
        """
        module = self.module()
        return self * (1 / module)
    
    def module(self):
        """Returns the module of this instance.
            
            Returns:
                float
        """
        response = 0 
        for _, val in enumerate(self.coordinates):
            response += val ** 2
        return math.sqrt(response)

    def get_projection_on(self, v):
        """Get the project of this vector on to another vector.
        
            Args:
                v(vector.Vector): The vector on to this instance will be 
                projected.
                
            Returns: 
                vector.Vector: A vector instance that represents the projection
                    of this instance on the vector v.
                    
            Raises:
                ValueError: If self and v doesn't have the same dimensions.
        """
        unitV = v.get_unit_vector()        
        return unitV * self.dot(unitV)

    def is_parallel_to(self, v):
        """        
            Args:
                v(vector.Vector): The vector that we will test if it's
                    parallel to self.
            
            Returns:
                bool: True if this instance and v are parallel vectors.
                
            Raises: 
                ValueError: If self and v doesn't have the same dimensions.
        """
        dotProduct = abs(self.dot(v))
        modulesMultiplacation = self.module() * v.module()
        return (self.is_zero() or v.is_zero() or math_util.isclose(dotProduct, modulesMultiplacation))
    
    def get_projection_parallel_to(self, v):
        """Gets the horizontal component of the projection of self into v.
            
            Args:
                v(vector.Vector): The vector where this vector will be projected on.
            
            Returns:
                vector.Vector : A new vector that is the projection of self on v.
                
            Raises: 
                ValueError: If self and v doesn't have the same dimensions.
            
        """
        unitV = v.get_unit_vector()
        return unitV * self.dot(unitV)
    
    def get_projection_orthogonal_to(self, v):
        """Gets the orthogonal projection of self onto v.
      
            Args:
                v(vector.Vector):  The vector orthogonal the this vector projection.
                
            Returns:
                vector.Vector: : A new vector that is orthogonal to v and if its added to the projection will result in the same self vector.
                
            Raises: 
                ValueError: If self and v doesn't have the same dimensions.
          
        """
        horizontalComponent = self.get_projection_parallel_to(v)
        verticalComponent = self -horizontalComponent
        return verticalComponent 
    
    def angle_with(self, v, inDegrees=False):
        """Returns the angle betwee self and v.
        
            Args:
                v(vector.Vector): The vector with which we are going to find the angle difference
                                
            Returns:
                float: The angle between self and v. The angle is in degrees
                if the param inDegrees is True, otherwise its in radians.
                
            Raises:
                ValueError: If the two vectors doesn't have the same length.
            
                
        """
        dotProduct = self.get_unit_vector().dot(v.get_unit_vector())
        angleInRadians = acos(dotProduct)
        if inDegrees:
            degrees_per_radian = 180 / math.pi
            return angleInRadians * degrees_per_radian
        return angleInRadians;
        
        
    def is_orthogonal_to(self, v, tolerance=1e-10):
        """Returns True if self and v are orthogonal.
            Args:
                v(vector.Vector): Will check if this vector and self are orthogonal
                    to each other.
                    
                tolerance(float): The minimum value that is considered zero.
            Returns:
                bool: True if self and v are orthogonal. False otherwise.
                
            Raises:
                ValueError: If the two vectors doesn't have the same length.

        """
        return self.dot(v) <= tolerance

    def is_zero(self, tolerance=1e-10):
        """Return True if self is the zero vector.
            
            Args:
                tolerance(float): Minimum value considered as zero.
                
            Returns:
                bool: True if module is zero.
        """
        return self.module() <= tolerance
     
    
    def cross_product(self, v):     
        """Returns a new vector representing the cross-product between self and v.
        
            Args: 
                v(vector.Vector): The other vector with which the cross product
                    will be performed.
                
            Returns: 
                vector.Vector: Returns a vector that is orthogonal vector to self.
                
            Raises:
                Exception: When the vectors are not specifically in 3 dimensions.
        """
        try:
            x1, y1, z1 = self.coordinates
            x2, y2, z2 = v.coordinates
            
            response = [y1 * z2 - y2 * z1,
                      - (x1 * z2 - x2 * z1),
                        x1 * y2 - x2 * y1
                        ]
            return Vector(response)
        except ValueError as e:
            msg = str(e)
            if msg == "need more than 2 values to unpack":
                v1 = Vector(self.coordinates + (0,))
                v2 = Vector(v.coordinates + (0,))
                return v1.cross_product(v2)
            elif (msg == 'too many values to unpack' or
                  msg == 'need more than 1 value to unpack'):
                raise Exception(self.ONLY_DEFINED_IN_MORE_THAN_TREE_DIM)
            else:
                raise e
        
        
    def area_of_parallelogram(self, v):
        """Return the area of the parallelogram formed by self and v.
        
            Args: 
                v(vector.Vector): Vector that will produce the parallelogram with self.
                
            Returns:
                float: The area of the parallelogram.
                            
        """
        angeleWith = self.angle_with(v)
        return self.module() * v.module() * math.sin(angeleWith)
    
    def area_of_triangle(self, v):
        """Return the area of one of the 2 triangles of the parallelogram. 
        
            Args:
                v(vector.Vector): Vector that will produce the parallelogram with self.
                
            Returns: The area of the parallelogram.
            
        """
        return self.area_of_parallelogram(v) / 2
       
