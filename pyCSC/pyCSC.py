"""
Py-Christoffel Symbols Calculator

"""

import sympy as sym
from sympy import init_printing
from sympy import Matrix
from IPython.display import display, Math
init_printing()


class PyCSC:
    """
    PyCSC is a class that calculates the following from a metric tensor:
        - Christoffel Symbols Second kind
        - Non-zero components of Christoffel Symbols First kind
        - Non-zero components of Riemann Tensor
        - Ricci Tensor
        - Ricci Scalar
        - Einstein Tensor
    
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
            
            variables_dict :: dict
              Mapping variable parameters (i.e., alpha, delta, epsilon) to the values
              specified in the `variables_values` parameter in `PyCSC.metric_tensor()`.
            
            metric :: `sympy.matrices.dense.MutableDenseMatrix`
              Metric tensor used for computing Christoffel Symbols. Components can 
              consist of variable and constant parameters.
              
            contra_metric :: `sympy.matrices.dense.MutableDenseMatrix`
              Inverse of the Metric tensor.

            riemann_dict :: dict
              Non-zero components of the Riemannian tensor.

            christoffel_fk :: dict
              Non-zero components of Christoffel Symbols First kind.
            
            christoffel_sk :: list of `sympy.matrices.dense.MutableDenseMatrix`
              Christoffel Symbols of second kind.
            
            ricci_tensor :: `sympy.matrices.dense.MutableDenseMatrix`
              Ricci tensor.
            
            ricci_scalar :: `sympy.core.add.Add`
              Ricci scalar.

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
                "`coordinates` parameter must be a list of coordinates. For example, coordinates = ['r', 'theta', 'phi']."
            )

        symbols_string = ''
        for coord in coordinates:
            symbols_string += coord + ' '
        
        self.coordinate_list = sym.symbols(symbols_string)

        self.variables_dict = dict() 

        self.metric = None
        self.contra_metric = None

        # Intializing the output
        self.christoffel_sk = None
        self.riemann_dict = None
        self.christoffel_fk = None
        self.ricci_tensor = None
        self.ricci_scalar = None

    def metric_tensor(self, matrix, variable_values={}):
        """
        Metric tensor.

        Parameters
        ----------
            matrix :: str
              A matrix with dimensions [`PyCSC.num_coordiantes`,`PyCSC.num_coordiantes`] 
              that represents the components of the metric tensor.
    
            variable_values :: dict()
              Variable parameters that will be substituted into the ,metric tensor before calculating christoffel symbols.

        Attributes
        ----------
            metric :: `sympy.matrices.dense.MutableDenseMatrix`
              Metric tensor used for computing Christoffel Symbols. Components can 
              consist of variable and constant parameters.

            variables_dict :: dict
              Mapping variable parameters (i.e., alpha, delta, epsilon) to the values
              specified in the `variables_values` parameter in `PyCSC.metric_tensor()`.

            contra_metric :: `sympy.matrices.dense.MutableDenseMatrix`
              Inverse of the metric tensor.

        Returns
        -------
            None
        
        """

        # Check input
        if isinstance(matrix, str):
           if Matrix(sym.parsing.sympy_parser.parse_expr(matrix)).shape == (self.num_coordinates,self.num_coordinates):
               self.metric = sym.parsing.sympy_parser.parse_expr(matrix, locals())
           else:
               raise ValueError(
                   f" `matrix` shape must be ({self.num_coordinates},{self.num_coordinates})"
               )

        else:
            raise TypeError(
                'matrix must be an n-dimensional list inside a str'
            )

        allowed_keys = ['alpha', 'delta', 'epsilon']
        if not isinstance(variable_values, dict):
            raise TypeError(
                "'variable_values' must be Type `dict` with three distinct keys: 'alpha', 'delta' and 'epsilon'"
            )
        else:
            if len(variable_values) > 3:
                raise ValueError(
                    'Maximum number of 3 variable parameters are allowed.'
                )
            
            if len(variable_values) !=0:

                for key in variable_values:
                    if key not in allowed_keys:
                        raise ValueError(
                            "Allowed keys in the parameter 'variable_values' are 'alpha', 'delta', 'epsilon'"
                        )
                    
                    if not isinstance(variable_values[key], str):
                        raise TypeError(
                            "Value for each key in the parameter 'variable_values' should be a mathematical expression inside a string (i.e., variable_values = {'alpha': 'x**2 + 2'})"
                        )
                    if not variable_values[key]:
                        raise ValueError(
                            "Value for each key in the parameter 'variable_values' cannot be an empty string."
                        )        
                 
        # substitute variables' definition into the metric tensor
        if 'alpha' in variable_values.keys():
            alpha = sym.symbols('alpha')
            self.variables_dict[alpha] = sym.parsing.sympy_parser.parse_expr(variable_values['alpha'])
            
        if 'delta' in variable_values.keys():
            delta = sym.symbols('delta')
            self.variables_dict[delta] = sym.parsing.sympy_parser.parse_expr(variable_values['delta'])
            
        if 'epsilon' in variable_values.keys():
            epsilon = sym.symbols('epsilon')
            self.variables_dict[epsilon] = sym.parsing.sympy_parser.parse_expr(variable_values['epsilon'])

        
        for i in range(self.num_coordinates):
            for j in range(self.num_coordinates):
                    if self.metric[i][j] == 0:
                        pass
                    else:
                        try:
                            self.metric[i][j] = sym.parsing.sympy_parser.parse_expr(self.metric[i][j]).subs(self.variables_dict)
                        except:
                            print('could not substitute')
            
        self.metric = sym.Matrix(self.metric)
        print('Metric Tensor')
        display(self.metric)
        self.contra_metric = self.metric.inv()

    def calculate_christoffel_symbol_fk(self, show_symbols=True):
        """
        Calculate Christoffel Symbols of first kind.

        Parameters
        ----------
            show_symbols :: Bool
              If True, display the non-zero components of the Christoffel symbols
              of first kind.

        Attributes
        ----------
            christoffel_fk :: dict
              Non-zero components of the Christoffel Symbols First kind.

        Returns
        -------
            None
            
        """

        if not isinstance(show_symbols, bool):
            raise TypeError(
                'show_symbols can either be True or False.'
            )

        # Check if the metric is there or not
        if self.metric is None:
            raise Exception(
                'Please specify the metric tensor using `PyCSC.metric_tensor()`' +
                ' before calculating Christoffel Symbols of first kind.'
            )
        
        self.christoffel_fk = dict()

        for a in range(self.num_coordinates):
            for b in range(self.num_coordinates):
                for c in range(b+1):
                    ans = (1/2) * (self.metric[a,c].diff(self.coordinate_list[b]) - self.metric[b,c].diff(self.coordinate_list[a]) + self.metric[b,a].diff(self.coordinate_list[c]))
                    try:
                        self.christoffel_fk[str(a) + str(b) + str(c)] = sym.simplify(ans)
                    except:
                        self.christoffel_fk[str(a) + str(b) + str(c)] = ans

        if show_symbols:    
            for key in self.christoffel_fk:
                Gamma = sym.symbols(f'Gamma_{key}')
                ans = self.christoffel_fk[key]
                if ans != 0:
                    display(Math(sym.latex(Gamma) + ' = ' + sym.latex(ans)))

    def calculate_christoffel_symbol(self, show_symbols=True):
        """
        Calculate Christoffel Symbols of second kind

        Parameters
        ----------
            show_symbols :: Bool
              If True, display the Christoffel symbols of second kind.

        Attributes
        ----------
            christoffel_sk :: list of `sympy.matrices.dense.MutableDenseMatrix`  
              A list that stores matrices of Christoffel symbols of second kind.

        Returns
        -------
            None
        
        """

        # Check pre-requisites
        if self.metric is None:
            raise Exception(
                'Please specify the metric tensor using `PyCSC.metric_tensor()`' +
                ' before calculating Christoffel Symbols of first kind.'
            )

        # Check inputs
        if not isinstance(show_symbols, bool):
            raise TypeError(
                'show_symbols can either be True or False.'
            )

        self.christoffel_sk = []
        # Create a list of Christoffel symbols matrices full of zeros
        it = 0 
        while it < self.num_coordinates:
            self.christoffel_sk.append(sym.zeros(self.num_coordinates,self.num_coordinates))
            it += 1
        
        for a in range(self.num_coordinates):
            for b in range(self.num_coordinates):
                # Using the torsion free condition to reduce number of iterations.
                for c in range(b+1):
                    for i in range(self.num_coordinates):
                        self.christoffel_sk[a][b,c] += (1/2) * self.contra_metric[a,i] * (self.metric[i,c].diff(self.coordinate_list[b]) + self.metric[b,i].diff(self.coordinate_list[c]) - self.metric[b,c].diff(self.coordinate_list[i]))
                    if (b == c):
                        pass
                    else:
                        self.christoffel_sk[a][c,b] = self.christoffel_sk[a][b,c]

        if show_symbols:
            for i in range(self.num_coordinates):
                Gamma = sym.symbols(f'Gamma^{i}_mu_nu')
                display(Math(sym.latex(Gamma) + ' = ' + sym.latex(self.christoffel_sk[i])))

    def calculate_riemann_tensor(self, show_tensor = True):
        """
        Calculate the Riemannian Tensor.

        Parameters
        ----------
            show_tensor :: bool
              If True, show the non-zero components of the Riemann Tensor.

        Attributes
        ----------
            riemann_dict :: dict
              Non-zero components of the Riemannian tensor.

        Returns
        -------
            None
    
        """

        # check input
        if not isinstance(show_tensor,bool):
            raise TypeError(
                'show_tensor can either be True or False.'
            )

        # check prerequisites
        if self.christoffel_sk is None:
            raise Exception(
                'Please calculate Christoffel symbols of second kind before calculating'
                + ' the Riemann Tensor.'
            )
        
        self.riemann_dict = dict()
        
        for a in range(self.num_coordinates):
            for b in range(self.num_coordinates):
                for c in range(self.num_coordinates):
                    for d in range(c+1):
                        key = str(a)+str(b)+str(c)+str(d)
                        key2 = str(a)+str(b)+str(d)+str(c)

                        if key in self.riemann_dict:
                            pass
                        else:                        
                            summation1 = 0
                            summation2 = 0
                            
                            for k in range(self.num_coordinates):
                                summation1 += self.christoffel_sk[a][c,k]*self.christoffel_sk[k][d,b]
                                summation2 += self.christoffel_sk[a][d,k]*self.christoffel_sk[k][c,b]
                                
                                ans = self.christoffel_sk[a][d,b].diff(self.coordinate_list[c]) - self.christoffel_sk[a][c,b].diff(self.coordinate_list[d]) + summation1 - summation2
                                
                                try:
                                    final_ans = sym.simplify(ans)
                                except:
                                    final_ans = ans

                                self.riemann_dict[key] = final_ans
                                self.riemann_dict[key2] = -final_ans

        if show_tensor:    
            for key in self.riemann_dict:
                R = sym.symbols(f'R^{key[0]}_{key[1:]}')
                ans = self.riemann_dict[key]
                if ans != 0:
                    display(Math(sym.latex(R) + ' = ' + sym.latex(ans)))
    
    def calculate_ricci_tensor(self, show_tensor=True):
        """
        Calculate Ricci Tensor. It is a symmetric tensor with respect to its indices R_{ab} = R_{ba}.

        Parameters
        ----------
            show_tensor :: bool
              If True, display the Ricci Tensor.

        Attributes
        ----------
            ricci_tensor :: `sympy.matrices.dense.MutableDenseMatrix`
                Ricci tensor.

        Returns
        -------
            None
        
        """

        
        # Check pre-requisites
        if self.riemann_dict is None:
            raise Exception(
                'Please calculate the Riemann Tensor before calculating' +
                ' the Ricci Tensor.'
            )
    
        # check input
        if not isinstance(show_tensor,bool):
            raise TypeError(
                'show_tensor can either be True or False.'
            )

        self.ricci_tensor = sym.zeros(self.num_coordinates, self.num_coordinates)

        for a in range(self.num_coordinates):
            for b in range(a+1):

                summation = 0
                for k in range(self.num_coordinates):
                    key = str(k)+str(a)+str(k)+str(b)
                    summation += self.riemann_dict[key]

                try:
                    ans = sym.simplify(summation)
                except:
                    ans = summation
                self.ricci_tensor[a,b] = ans
                if a==b:
                    pass
                else:
                    self.ricci_tensor[b,a] = ans

        if show_tensor:
            R = sym.symbols(f'R_mu_nu')
            display(Math(sym.latex(R) + ' = ' + sym.latex(self.ricci_tensor)))
    

    def calculate_ricci_scalar(self, show_scalar=True):
        """
        Calculate the Ricci

        Parameters
        ----------
            show_scalar :: bool
              If True, display the value of the Ricci scalar.

        Attributes
        ----------
            ricci_scalar :: `sympy.core.add.Add`
              Ricci scalar

        Returns
        -------
            None
        
        """

        # check input
        if not isinstance(show_scalar,bool):
            raise TypeError(
                'show_scalar can either be True or False.'
            )

        # prerequisites
        if self.ricci_tensor is None:
            raise Exception(
                'Please calculate Ricci Tensor before calculating the Ricci scalar.'
            )

        self.ricci_scalar = 0

        for k in range(self.num_coordinates):
            for m in range(self.num_coordinates):
                self.ricci_scalar += self.contra_metric[k,m] * self.ricci_tensor[k,m]

        if show_scalar:
            try:
                display(sym.simplify(self.ricci_scalar))
            except:
                display(self.ricci_scalar)
    
    def calculate_einstein_tensor(self, show_tensor=False):
        """
        Calculate Einstein Tensor. It is a symmetric tensor with respect to its indices. G_{ab} = G_{ba}

        Parameters
        ----------
            show_tensor :: bool
              If True, display the Einstein Tensor.

        Attributes
        ----------
            None

        Returns
        -------
            einstein_tensor :: `sympy.matrices.dense.MutableDenseMatrix`
              Einstein Tensor
        
        """

        if not isinstance(show_tensor, bool):
            raise TypeError(
                'show_tensor can either be True or False'
            )

        if self.ricci_scalar is None:
            raise Exception(
                'Please calculate Ricci scalar before calculating Einstein Tensor.'
            )
        
        einstein_tensor = sym.zeros(self.num_coordinates, self.num_coordinates)

        for a in range(self.num_coordinates):
            for b in range(a+1):
                try:
                    ans = sym.simplify(self.ricci_tensor[a,b] - (1/2)*self.metric[a,b]*self.ricci_scalar)
                except:
                    ans = self.ricci_tensor[a,b] - (1/2)*self.metric[a,b]*self.ricci_scalar
                einstein_tensor[a,b] = ans
                if a==b:
                    pass
                else:
                    einstein_tensor[b,a] = ans

        if show_tensor:
            G = sym.symbols(f'G_mu_nu')
            display(Math(sym.latex(G) + ' = ' + sym.latex(einstein_tensor)))
        
        return einstein_tensor
