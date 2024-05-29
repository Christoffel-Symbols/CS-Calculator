"""
Christoffel symbols of the given space-time model.
"""


import sympy as sym
from IPython.display import display, Math
from copy import copy


from .pycsc import PyCSC
from .metric_tensor import MetricTensor


class ChristoffelSymbols:
    """
    Create a ChristoffelSymbols instance.
    """

    def __init__(self, PyCSCObj, MetricTensorObj):
        """
        
        Parameters
        ----------
            PyCSCObj :: `pycsc.pycsc.PyCSC` object
                The `pycsc.pycsc.PyCSC` object containing the information regarding the coordinates.

            MetricTensorObj :: `pycsc.metric_tensor.MetricTensorObj` object
                The `pycsc.metric_tensor.MetricTensorObj` object containing the information regarding the metric tensor.

        Attributes
        ----------
            PyCSCObj :: `pycsc.pycsc.PyCSC` object
                The `pycsc.pycsc.PyCSC` object containing the information regarding the coordinates.

            MetricTensorObj :: `pycsc.metric_tensor.MetricTensorObj` object
                The `pycsc.metric_tensor.MetricTensorObj` object containing the information regarding the metric tensor.

            christoffel_fk :: dict
                Dictionary consisting values of Christoffel Symbols of first kind.

            cleaner_christoffel_fk :: dict
                Cleaner version of christoffel_fk with one key amongst its symmetric counterparts and non-zero values.

            christoffel_sk :: list
                List of `sympy.matrices.dense.MutableDenseMatrix` representing christoffel Symbols of second kind.

        Returns
        -------
            None

        """

        if not isinstance(PyCSCObj, PyCSC):
            raise TypeError(
                "PyCSCObj must be a `pycsc.pycsc.PyCSC` object"
            )
        
        if not isinstance(MetricTensorObj, MetricTensor):
            raise TypeError(
                "MetricTensorObj must be a `pyCSC.metric_tensor.MetricTensor` object"
            )

        self.MetricTensorObj = MetricTensorObj
        self.PyCSCObj = PyCSCObj
        self.christoffel_fk = None
        self.cleaner_christoffel_fk = None
        self.christoffel_sk = None


    def calculate(self, config='ull', simplify=True, show=True):
        """
        Calculates christoffel symbols of both first and second kinds.

        Parameters
        ----------
            config :: str
                Configuration of christoffel symbols asked to calculate.

            simplify :: bool
                If True, then simplify sympy expressions.
                if False, skip simplification.

            show :: Bool
                If True, then show the calculated values.
                If False, do not show the calculate values.

        Attributes
        ----------
            christoffel_fk :: dict
                Dictionary consisting of values of Christoffel Symbols of first kind.

            christoffel_sk :: list
                List of `sympy.matrices.dense.MutableDenseMatrix` representing christoffel Symbols of second kind.

            cleaner_christoffel_fk :: dict
                Cleaner version of christoffel_fk with one key amongst its symmetric counterparts and non-zero values.

            config :: str
                Configuration of christoffel symbols asked to calculate.

        Returns
            if config is 'lll':
                christoffel_fk :: dict
                    Dictionary consisting of values of Christoffel Symbols of first kind.

            if config is 'ull':
                christoffel_sk :: list
                    List of `sympy.matrices.dense.MutableDenseMatrix` representing christoffel Symbols of second kind.
        """


        if isinstance(config, str):
            if (config != 'ull' and config != 'lll'):
                raise ValueError(
                     "Not a valid config"
                )
            else:
                self.config = config
        else:
            raise TypeError(
                f"Invalid type {type(config)}" + " for config parameter"
            )
        
        if not isinstance(show, bool):
            raise TypeError(
                f"invalid type {type(show)}" + " for show parameter"
            )
        

        if self.MetricTensorObj.tensor is None:
            raise Exception(
                "`pyCSC.christoffel_symbols.ChristoffelSymbols.calculate` called before" + " initializing `pyCSC.metric_tensor.MetricTensor` class"
            )

        if self.MetricTensorObj.config == 'll':
            covariant_metric = self.MetricTensorObj.tensor # ll config
            contra_metric = self.MetricTensorObj.change_config()
        else:
            contra_metric = self.MetricTensorObj.tensor # uu config
            covariant_metric = self.MetricTensorObj.change_config()

        iter = self.PyCSCObj.num_coordinates
        coordinate_list = self.PyCSCObj.coordinate_list

        if self.config == 'ull':
            self.christoffel_sk = []
            for _ in range(iter):
                # Need to initialize christoffel_sk as a list of zero matrices everytime a user decides to call calculate method.
                self.christoffel_sk.append(sym.zeros(iter,iter))

            for a in range(iter):
                for b in range(iter):
                    #Using torsion free condition to reduce the 3rd iteration's loop count
                    for c in range(b+1):
                        ans = 0
                        for i in range(iter):
                            ans += (1/2) * contra_metric[a,i] * (sym.diff(covariant_metric[i,c],coordinate_list[b]) + sym.diff(covariant_metric[b,i],coordinate_list[c]) - sym.diff(covariant_metric[b,c],coordinate_list[i]))
                        
                        if simplify:
                            try:
                                ans = sym.simplify(sym.trigsimp(sym.expand_trig(sym.factor(ans))), rational=True)
                            except:
                                print("Simplification not possible")
                        self.christoffel_sk[a][b,c] = ans

                        if (b!=c):
                            self.christoffel_sk[a][c,b] = self.christoffel_sk[a][b,c]

            if show:
                for i in range(self.PyCSCObj.num_coordinates):
                    Gamma = sym.symbols(f'Gamma^{i}_mu_nu')
                    display(Math(sym.latex(Gamma) + ' = ' + sym.latex(self.christoffel_sk[i])))
        
            return self.christoffel_sk # list of matrices

        else:
            # covariant version of christoffel symbols
            self.christoffel_fk = dict()
            self.cleaner_christoffel_fk = dict()

            for a in range(iter):
                for b in range(iter):
                    for c in range(b+1):
                        key = str(a) + str(b) + str(c)
                        symm_key = str(a) + str(c) + str(b)
                        ans = (1/2) * (sym.diff(covariant_metric[a,c],coordinate_list[b]) - sym.diff(covariant_metric[b,c],coordinate_list[a]) + sym.diff(covariant_metric[b,a],coordinate_list[c]))

                        if simplify:                    
                            try:
                                ans = sym.simplify(sym.trigsimp(sym.expand_trig(sym.factor(ans))), rational=True)
                            except:
                                print('Simplification not possible')

                        self.christoffel_fk[key] = ans
                        if b != c:
                            self.christoffel_fk[symm_key] = ans

                        if ans !=0:
                            self.cleaner_christoffel_fk[key] = ans


            if show:
                # why shallow copy work
                # <https://stackoverflow.com/questions/3975376/why-updating-shallow-copy-dictionary-doesnt-update-original-dictionary>
                for key,value in self.cleaner_christoffel_fk.items():
                    Gamma = sym.symbols(f'Gamma_{key}')
                    display(Math(sym.latex(Gamma) + ' = ' + sym.latex(value)))

            return self.christoffel_fk
        
    # Avoid using this method! It takes 3 times more time than the calculate() method.
        
    # def change_config(self, show=True):
    #     """
        
    #     """

    #     if self.MetricTensorObj.config == 'll':
    #         covariant_metric = self.MetricTensorObj.tensor # ll config
    #         contra_metric = self.MetricTensorObj.change_config()
    #     else:
    #         contra_metric = self.MetricTensorObj.tensor # uu config
    #         covariant_metric = self.MetricTensorObj.change_config()

    #     if self.config == 'ull':
        
    #         # Error handling
    #         if self.christoffel_sk == []:
    #             raise Exception(
    #                 "`pyCSC.ChristoffelSymbols.change_config` called before"
    #                  + "`pyCSC.ChristoffelSymbols.calculate`"
    #             )
            
    #         self.christoffel_fk = dict()
            
    #         for a in range(self.PyCSCObj.num_coordinates):
    #             for b in range(self.PyCSCObj.num_coordinates):
    #                 for c in range(b+1):
    #                     key = str(a) + str(b) + str(c)
    #                     ans = 0
    #                     for i in range(self.PyCSCObj.num_coordinates):
    #                         ans += covariant_metric[a,i]*self.christoffel_sk[i][b,c]
                        
    #                     if ans != 0:
    #                         try:
    #                             self.christoffel_fk[key] = sym.simplify(ans)
    #                         except:
    #                             self.christoffel_fk[key] = ans

    #         if show:
    #                 for key in self.christoffel_fk:
    #                     Gamma = sym.symbols(f'Gamma_{key}')
    #                     value = self.christoffel_fk[key]
    #                     # earlier we would check whether value was non zero, now we won't because keys with value  = 0 won't be added to the dict
    #                     display(Math(sym.latex(Gamma) + ' = ' + sym.latex(value)))

    #         return self.christoffel_fk
    #     else:

    #         if self.christoffel_fk is None:
    #             raise Exception(
    #                 "`pyCSC.ChristoffelSymbols.change_config` called before"
    #                  + "`pyCSC.ChristoffelSymbols.calculate`"
    #             )
            
    #         for _ in range(self.PyCSCObj.num_coordinates):
    #             self.christoffel_sk.append(sym.zeros(self.PyCSCObj.num_coordinates,self.PyCSCObj.num_coordinates))

    #         for a in range(self.PyCSCObj.num_coordinates):
    #             for b in range(self.PyCSCObj.num_coordinates):

    #                 # last two indices are symmetric with respect to torsion free condition
    #                 for c in range(b+1):
    #                     ans = 0
    #                     for i in range(self.PyCSCObj.num_coordinates):
    #                         key = str(i) + str(b) + str(c)
    #                         symmetry_key = str(i) + str(c) + str(b) 
    #                         try:
    #                             # torsion free condition
    #                             ans += contra_metric[a,i]*self.christoffel_fk[key or symmetry_key]
    #                         except:
    #                             # If there's a keyError then that means value for that key is 0.
    #                             ans += 0
    #                     try:
    #                         self.christoffel_sk[a][b,c] = sym.simplify(ans)
    #                     except:
    #                         self.christoffel_sk[a][b,c] = ans

    #                     if (b!=c):
    #                         self.christoffel_sk[a][c,b] = self.christoffel_sk[a][b,c]

    #         if show:
    #             for i in range(self.PyCSCObj.num_coordinates):
    #                 Gamma = sym.symbols(f'Gamma^{i}_mu_nu')
    #                 display(Math(sym.latex(Gamma) + ' = ' + sym.latex(self.christoffel_sk[i])))

    #         return self.christoffel_sk
    


        



