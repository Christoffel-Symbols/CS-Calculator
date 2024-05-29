"""
Py-Christoffel Symbols Calculator

"""

import sympy as sym


class PyCSC:
    """
    PyCSC is a class that calculates the following from a metric tensor:
        - Christoffel Symbols (First and Second kinds)
        - Riemann Tensor (First and Seconds kinds)
        - Ricci Tensor (covariant version)
        - Ricci Scalar 
        - Einstein Tensor (covariant version)
    
    """

    def __init__(self, coordinates):
        """
        Create a PyCSC instance.

        Parameters
        ----------
            coordinates  :: list
              List of coordinates specified by the user.

        Attributes
        ----------
            num_coordinates :: int
              Number of dimensions in the space-time model. 

            coordinate_list :: list
              List of coordinates according to `coordinates`.

        Returns
        -------
            `PyCSC` instance
        
        """

        # Check input
        if isinstance(coordinates,list):
            # Check for duplicates
            if len(coordinates) != len(set(coordinates)):
                raise ValueError(
                    "Repeated coordinates not allowed."
                )

            for coord in coordinates:
                if not isinstance(coord, str):
                    raise ValueError(
                        "coordinates parameter should be a list of coordinates with each coordinate represented as a string. For example, coordaintes = ['theta', 'phi']."
                    )
                
            if (len(coordinates) <= 4 and len(coordinates) >= 2):
                self.num_coordinates = len(coordinates)
            else:
                raise ValueError(
                    "Number of coordinates can be maximum 4 (1 time coordinate and 3 space coordinates, i.e., t,x,y,z) or minimum 2 (both space coordinates, i.e., x,y)."
                )
        else:
            raise TypeError(
                f"Invalid type {coordinates}" + " coordinates parameter"
            )

        symbols_string = ''
        for coord in coordinates:
            symbols_string += coord + ' '
        
        self.coordinate_list = sym.symbols(symbols_string)
