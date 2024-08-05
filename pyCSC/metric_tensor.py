"""
Metric tensor class.
"""


import sympy as sym
from sympy import Matrix
from IPython.display import display
from sympy.parsing.sympy_parser import parse_expr

from .pycsc import PyCSC

class MetricTensor:
    """
    MetricTensor class.

    """

    def __init__(self, matrix, PyCSCObj, variable_values={}, config='ll'):
        """
        Create a MetricTensor instance.
        
        Parameters
        ----------
            matrix :: str
                String representation of the metric tensor.

            PyCSCObj :: `pycsc.pycsc.PyCSC` object
                The `pycsc.pycsc.PyCSC` object containing the information regarding the coordinates.

            variable_values :: dict
                Variable values that will be substitued into the metric tensor.

            config :: str
                Configuration of the Metric Tensor: contravariant or covariant.

        Attributes
        ----------

            PyCSCObj :: `pycsc.pycsc.PyCSC` object
                The `pycsc.pycsc.PyCSC` object containing the information regarding the coordinates.

            tensor :: `sympy.matrices.dense.MutableDenseMatrix`
                metric tensor.

        Returns
        -------
            None
        
        """

        if not isinstance(PyCSCObj, PyCSC):
            raise TypeError(
                "PyCSCObj must be a `pycsc.pycsc.PyCSC` object"
            )

        self.PyCSCObj = PyCSCObj

        if isinstance(config, str):
            if (config != 'll') and (config != 'ul'):
                raise ValueError(
                    "Not a valid config"
                )
            self.config = config
        else:
            raise TypeError(
                f"Invalid type {type(config)}" + " for config parameter"
            )

        self.tensor = None

        if isinstance(matrix, str):
            nested_list = parse_expr(matrix, local_dict=locals())
            if Matrix(nested_list).shape == (self.PyCSCObj.num_coordinates,self.PyCSCObj.num_coordinates):
                self.tensor = parse_expr(matrix, local_dict=locals())
            else:
                raise Exception(
                    f"Dimensions of matrix should be according to the number of coordinates specified"
                )
        else:
            raise TypeError(
                f"Invalud type {type(matrix)}" " for matrix parameter. Must be a string representation of a nested list" + "(i.e., '[[r,0],[0,sin(theta)]]' "
            )
        
        allowed_variables = ['alpha', 'delta', 'epsilon']
        if not isinstance(variable_values, dict):
            raise TypeError(
                "variable_values must be type `dict`"
            )
        else:
            if len(variable_values) > 3:
                raise ValueError(
                    'Maximum 3 variables can de defined'
                )

            for key in variable_values:
                if key not in allowed_variables:
                    raise Exception(
                        "Allowed variable are 'alpha', 'delta', and 'epsilon'"
                    )
                if not isinstance(variable_values[key],str):
                    raise TypeError(
                        f"Invalid type {type(variable_values[key])}" + f" for {key}"
                    )
                if not variable_values[key]:
                    raise Exception(
                        "Value cannot be empty strings"
                    )
                
        variables_dict = dict()
        if 'alpha' in variable_values.keys():
            alpha = sym.Symbol('alpha')
            variables_dict[alpha] = parse_expr(variable_values['alpha'], local_dict=locals())
        if 'delta' in variable_values.keys():
            delta = sym.Symbol('delta')
            variables_dict[delta] = parse_expr(variable_values['delta'], local_dict=locals())
        if 'epsilon' in variable_values.keys():
            epsilon = sym.Symbol('epsilon')
            variables_dict[epsilon] = parse_expr(variable_values['epsilon'], local_dict=locals())

        for dummy_index1 in range(self.PyCSCObj.num_coordinates):
            for dummy_index2 in range(self.PyCSCObj.num_coordinates):
                if self.tensor[dummy_index1][dummy_index2] != 0:
                    try:
                        # Need to parse it again when doing non-notebook local development
                        self.tensor[dummy_index1][dummy_index2] = parse_expr(self.tensor[dummy_index1][dummy_index2], local_dict=locals()).subs(variables_dict)
                    except:
                        raise Exception(
                            "Could not substitute variable values" + "into the metric"
                        )
                    

        self.tensor = Matrix(self.tensor)

        display(self.tensor)

    def change_config(self):
        """
        Changes the configuration of the metric tensor.

        Parameters
        ----------
            None

        Attributes
        ----------
            None

        Returns
        -------
            tensor.inv() :: sympy.matrices.dense.MutableDenseMatrix
                Inverse of the metric.
        """

        return self.tensor.inv()
