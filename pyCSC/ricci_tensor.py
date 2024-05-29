"""
Ricci Tensor
"""

import sympy as sym
from IPython.display import display, Math

from .metric_tensor import MetricTensor
from .pycsc import PyCSC
from .riemann_tensor import RiemannTensor

class RicciTensor:
    """
    `pycsc.ricci_tensor.RicciTensor` class
    """

    def __init__(self, PyCSCObj,MetricTensorObj, RiemannTensorObj):
        """
        Create a RicciTensor instance.

        Parameters
        ----------
            PyCSCObj :: `pycsc.pycsc.PyCSC` object
                The `pycsc.pycsc.PyCSC` object containing the information regarding the coordinates.

            MetricTensorObj :: `pycsc.metric_tensor.MetricTensorObj` object
                The `pycsc.metric_tensor.MetricTensorObj` object containing the information regarding the metric tensor.

            RiemannTensorObj :: `pycsc.riemann_tensor.RiemannTensorObj` object
                The `pycsc.riemann_tensor.RiemannTensorObj` object containing the information regarding the Riemann Tensor of first and second kinds.

        Attributes
        ----------
            PyCSCObj :: `pycsc.pycsc.PyCSC` object
                The `pycsc.pycsc.PyCSC` object containing the information regarding the coordinates.

            MetricTensorObj :: `pycsc.metric_tensor.MetricTensorObj` object
                The `pycsc.metric_tensor.MetricTensorObj` object containing the information regarding the metric tensor.

            RiemannTensorObj :: `pycsc.riemann_tensor.RiemannTensorObj` object
                The `pycsc.riemann_tensor.RiemannTensorObj` object containing the information regarding the Riemann Tensor of first and second kinds.

            ricci_fk :: `sympy.matrices.dense.MutableDenseMatrix`
                Ricci tensor of first kind.

            ricci_sk :: `sympy.matrices.dense.MutableDenseMatrix`
                Ricci tensor of first kind.

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
        
        if not isinstance(RiemannTensorObj, RiemannTensor):
            raise TypeError(
                "RiemannTensorObj must be a `pyCSC.riemann_tensor.RiemannTensor` object"
            )

        self.PyCSCObj = PyCSCObj
        self.MetricTensorObj = MetricTensorObj
        self.RiemannTensorObj = RiemannTensorObj
        

        self.ricci_fk = None
        self.ricci_sk = None

    def calculate(self, show=True, config='ll', simplify=True):
        """
        R_{ab} = g^{cd}*R_{dacb}
        Calculate Ricci Tensor from Riemann Tensor.

        Parameters
        ----------
            show :: bool
                If True, then print Ricci Tensor
                If False, skip printing

            config :: str
                Configuration of Ricci Tensor

            simplify :: bool
                If True, simplify sympy expressions
                If False, skip simplification

        Attributes
        ----------
            ricci_fk :: `sympy.matrices.dense.MutableDenseMatrix`
                Ricci tensor of first kind.

        Returns
        -------
            ricci_fk :: `sympy.matrices.dense.MutableDenseMatrix`
                Ricci tensor of first kind.
        """

        if isinstance(config, str):
            if (config != 'll' and config != 'uu'):
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
                f"invalid type {type(simplify)}" + " for show parameter"
            )
        
        if self.RiemannTensorObj.riemann_fk is None:
            raise Exception(
                "`pycsc.riemann_tensor.RiemannTensor.calculate_from_riemann`" + " method called before `pyscs.riemann_tensor.RiemannTensor.calculate_from_symbols`"
            )
        
        self.ricci_fk = sym.zeros(self.PyCSCObj.num_coordinates, self.PyCSCObj.num_coordinates)

        if self.config == 'll':

            if self.MetricTensorObj.config == 'll':
                contra_metric = self.MetricTensorObj.change_config()
            else:
                contra_metric = self.MetricTensorObj.ricci_fk

            riemann_fk = self.RiemannTensorObj.riemann_fk

            for a in range(self.PyCSCObj.num_coordinates):
                for b in range(a+1):
                    
                    ans = 0
                    for c in range(self.PyCSCObj.num_coordinates):
                        for d in range(self.PyCSCObj.num_coordinates):
                            ans += contra_metric[c,d]*riemann_fk[str(d) + str(a) + str(c) + str(b)]

                    if simplify:
                        try:
                            ans = sym.simplify(sym.trigsimp(sym.expand_trig(sym.factor(ans))), rational=True)
                        except:
                            print('Simplification skipped')
                    
                    self.ricci_fk[a,b] = ans
                    if a!=b:
                        self.ricci_fk[b,a] = ans

            if show:
                R = sym.symbols(f'R_mu_nu')
                display(Math(sym.latex(R) + ' = ' + sym.latex(self.ricci_fk)))

            return self.ricci_fk

        else:
            raise Exception(
                "Contravariant version not supported yet"
            )
