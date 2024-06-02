"""
Riemann Tensor
"""

import sympy as sym
from IPython.display import display, Math
from .pycsc import PyCSC
from .metric_tensor import MetricTensor
from .christoffel_symbols import ChristoffelSymbols

class RiemannTensor:
    """
    `pyCSC.riemann_tensor.RiemannTensor` class.
    """

    def __init__(self, PyCSCObj, MetricTensorObj, ChristoffelSymbolsObj):
        """
        Create a RiemannTensor instance.

        Parameters
        ----------
            PyCSCObj :: `pycsc.pycsc.PyCSC` object
                The `pycsc.pycsc.PyCSC` object containing the information regarding the coordinates.

            MetricTensorObj :: `pycsc.metric_tensor.MetricTensorObj` object
                The `pycsc.metric_tensor.MetricTensorObj` object containing the information regarding the metric tensor.
            
            ChristoffelSymbolsObjs :: `pycsc.christoffel_symbols.ChristoffelSymbolsObj` object
                The `pycsc.christoffel_symbols.ChristoffelSymbolsObj` object containing the information regarding the christoffel symbols of first and second kinds.

        Attributes
        ----------
            riemann_fk :: dict
                Dictionary consisting values of Riemann Tensor of first kind.

            riemann_sk :: dict
                Dictionary consisting values of Riemann Tensor of second kind.

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
                "MetricTensorObj must be a `pycsc.metric_tensor.MetricTensor` object"
            )

        if not isinstance(ChristoffelSymbolsObj, ChristoffelSymbols):
            raise TypeError(
                "ChristoffelSymbolsObj must be a `pyCSC.christoffel_symbols.ChristoffelSymbols` object"
            )
        
        self.MetricTensorObj = MetricTensorObj
        self.PyCSCObj = PyCSCObj
        self.ChristoffelSymbolsObj = ChristoffelSymbolsObj
        self.riemann_fk = None
        self.cleaner_riemann_fk = None
        self.riemann_sk = None
        self.cleaner_riemann_sk = None


    def calculate(self, config='llll', simplify=True, show=True):
        """
        Calculates Riemann Tensor of both first and second kinds. This method always calculates the covariant version of Riemann Tensor indifferent to the configuration specified by the user because it's quicker to calculate contravariant Riemann Tensor by contracting covariant Riemann Tensor.

        Parameters
        ----------
            config :: str
                Configuration of Riemann Tensor specified.

            simplify :: bool
                If True, then simplify sympy expressions.
                if False, skip simplification.

            show :: bool
                If True, then shows non-zero values and skips symmetric counterparts of indices.
                If False, do not show the calculated values.
        """

        if isinstance(config, str):
            if (config != 'ulll' and config != 'llll'):
                raise ValueError(
                    "Not a valid config."
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
        
        if not isinstance(simplify, bool):
            raise TypeError(
                f"invalid type {type(simplify)}" + " for simplify parameter"
            )
        
        if self.ChristoffelSymbolsObj.christoffel_fk is None and self.ChristoffelSymbolsObj.christoffel_sk is None:
            raise Exception(
                "`pyCSC.riemann_tensor.RiemannTensor.calculate` called before" + " `pyCSC.christoffel_symbols.ChristoffelSymbolsObj.calculate`"
            )
    
        christoffel_fk = self.ChristoffelSymbolsObj.christoffel_fk
        christoffel_sk = self.ChristoffelSymbolsObj.christoffel_sk

        coordinate_list = self.PyCSCObj.coordinate_list

        if self.config == 'llll':
            #From Rindler
            #R_{abcd} = T{abd,c} - T{abc,d} + T{i,bc}*T{iad} - T{i,bd}*T{iac}

            self.riemann_fk = dict()
            self.cleaner_riemann_fk = dict()
            
            for a in range(self.PyCSCObj.num_coordinates):
                for b in range(self.PyCSCObj.num_coordinates):
                    for c in range(self.PyCSCObj.num_coordinates):
                        for d in range(self.PyCSCObj.num_coordinates):

                            key = str(a) + str(b) + str(c) + str(d)

                            if (a == b) or (c==d):
                                self.riemann_fk[key] = 0

                            # For a given key, 3 more keys are filled
                            elif key not in self.riemann_fk:

                                key2 = str(a) + str(b) + str(d) + str(c)
                                key3 = str(b) + str(a) + str(c) + str(d)
                                key4 = str(c) + str(d) + str(a) + str(b)

                                # Bianchi identity
                                key5  = str(a) + str(c) + str(d) + str(b)
                                key6 = str(a) + str(d) + str(b) + str(c)

                                # simultaneously switching the pair of indices
                                key7 = str(b) + str(a) + str(d) + str(c)

                                # first switch the pair of indices and then switching order in them
                                key8 = str(c) + str(d) + str(b) + str(a)
                                key9 = str(d) + str(c) + str(a) + str(b)

                                # switch individual pairs and then flip the pairs
                                key10 = str(d) + str(c) + str(b) + str(a)

                                minor_key1 = str(a) + str(b) + str(d)
                                term1 = christoffel_fk[minor_key1]

                                minor_key2 = str(a) + str(b) + str(c)
                                term2 = christoffel_fk[minor_key2]

                                summation1 = 0
                                summation2 = 0


                                for i in range(self.PyCSCObj.num_coordinates):
                                    cs_key1 = str(i) + str(a) + str(d) #christoffel symbols key 1
                                    cs_key2 = str(i) + str(a) + str(c) #christoffel symbols key 2
                                    summation1 += christoffel_sk[i][b,c]*christoffel_fk[cs_key1]
                                    summation2 += christoffel_sk[i][b,d]*christoffel_fk[cs_key2]

                                ans = sym.diff(term1,coordinate_list[c]) - sym.diff(term2,coordinate_list[d]) + summation1 - summation2
                                
                                if simplify:
                                    try:
                                        ans = sym.simplify(sym.trigsimp(sym.expand_trig(sym.factor(ans))), rational=True)
                                    except:
                                        # set simplified_failed = True
                                        print('Simplification not possible')

                                    
                                # Come back to this to eliminiate symmetric keys
                                self.riemann_fk[key] = ans
                                self.riemann_fk[key2] = -ans
                                self.riemann_fk[key3] = -ans
                                self.riemann_fk[key4] = ans
                                self.riemann_fk[key7] = ans
                                self.riemann_fk[key8] = -ans
                                self.riemann_fk[key9] = -ans
                                self.riemann_fk[key10] = ans

                                if key5 in self.riemann_fk:
                                    self.riemann_fk[key6] = -(ans + self.riemann_fk[key5])

                                elif key6 in self.riemann_fk:
                                    self.riemann_fk[key5] = -(ans + self.riemann_fk[key6])

                                if ans !=0:
                                    self.cleaner_riemann_fk[key] = ans


            if show:
                for key, value in self.cleaner_riemann_fk.items():
                    R = sym.symbols(f'R_{key}')
                    display(Math(sym.latex(R) + ' = ' + sym.latex(value)))

            return self.riemann_fk
        
        
        else:

            self.riemann_sk = dict()
            self.cleaner_riemann_sk = dict()

            # Calculating using Christoffel symbols.
            for a in range(self.PyCSCObj.num_coordinates):
                for b in range(self.PyCSCObj.num_coordinates):
                    for c in range(self.PyCSCObj.num_coordinates):
                        for d in range(c+1):
                            
                            key = str(a)+str(b)+str(c)+str(d)
                            key2 = str(a)+str(b)+str(d)+str(c)
          
                            summation1 = 0
                            summation2 = 0
                            
                            for k in range(self.PyCSCObj.num_coordinates):
                                summation1 += christoffel_sk[a][c,k]*christoffel_sk[k][d,b]
                                summation2 += christoffel_sk[a][d,k]*christoffel_sk[k][c,b]
                                
                                ans = sym.diff(christoffel_sk[a][d,b],coordinate_list[c]) - sym.diff(christoffel_sk[a][c,b],coordinate_list[d]) + summation1 - summation2
                                
                            if simplify:
                                try:
                                    ans = sym.simplify(sym.trigsimp(sym.expand_trig(sym.factor(ans))), rational=True)
                                except:
                                    print('Simplification not possible')

                            self.riemann_sk[key] = ans
                            if c!=d:
                                self.riemann_sk[key2] = - ans

                            if ans!=0:
                                self.cleaner_riemann_sk[key] = ans

            if show:
                for key,value in self.cleaner_riemann_sk.items():
                    R = sym.symbols(f'R^{key[0]}_{key[1:]}')
                    display(Math(sym.latex(R) + ' = ' + sym.latex(value)))

            return self.riemann_sk