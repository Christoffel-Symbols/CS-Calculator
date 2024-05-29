"""
Ricci Scalar

"""

import sympy as sym
from IPython.display import display, Math

from .metric_tensor import MetricTensor
from .pycsc import PyCSC
from .riemann_tensor import RiemannTensor
from .christoffel_symbols import ChristoffelSymbols
from .ricci_tensor import RicciTensor

class RicciScalar:
    """
    `pycsc.ricci_scalar.RicciScalar`class
    
    """

    def __init__(self, PyCSCObj, MetricTensorObj, RicciTensorObj):
        """
        Create a `pycsc.ricci_scalar.RicciScalar` instance.

        Parameters
        ----------
            PyCSCObj :: `pycsc.pycsc.PyCSC` object
                The `pycsc.pycsc.PyCSC` object containing the information regarding the coordinates.

            MetricTensorObj :: `pycsc.metric_tensor.MetricTensorObj` object
                The `pycsc.metric_tensor.MetricTensorObj` object containing the information regarding the metric tensor.

            RicciTensorObj :: `pycsc.ricci_tensor.RicciTensorObj` object
                The `pycsc.ricci_tensor.RicciTensorObj` object containing the information regarding the Ricci Tensor of first and second kinds.

        Attributes
        ----------
            PyCSCObj :: `pycsc.pycsc.PyCSC` object
                The `pycsc.pycsc.PyCSC` object containing the information regarding the coordinates.

            MetricTensorObj :: `pycsc.metric_tensor.MetricTensorObj` object
                The `pycsc.metric_tensor.MetricTensorObj` object containing the information regarding the metric tensor.

            RicciTensorObj :: `pycsc.ricci_tensor.RicciTensorObj` object
                The `pycsc.ricci_tensor.RicciTensorObj` object containing the information regarding the Ricci Tensor of first and second kinds.

            scalar :: `sympy.core.xxx.xxx`
                Ricci scalar

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
    
        if not isinstance(RicciTensorObj, RicciTensor):
            raise TypeError(
            "RicciTensorObj must be a `pycsc.ricci_tensor.RicciTensor` object"
            )
        
        self.PyCSCObj = PyCSCObj
        self.MetricTensorObj = MetricTensorObj
        self.RicciTensorObj = RicciTensorObj

        self.scalar = None

    def calculate(self,show=True,simplify=True):
        """
        Calculate Ricci scalar.

        Parameters
        ----------
            show :: bool
                If True, then print Ricci Tensor
                If False, skip printing

            simplify :: bool
                If True, simplify sympy expressions
                If False, skip simplification

        Attributes
        ----------
            scalar :: `sympy.core.xxx.xxx`
                Ricci scalar

        Returns
        -------
            scalar :: `sympy.core.xxx.xxx`
                Ricci scalar

        """

        if not isinstance(show, bool):
            raise TypeError(
                f"invalid type {type(show)}" + " for show parameter"
            )
        
        if not isinstance(simplify, bool):
            raise TypeError(
                f"invalid type {type(simplify)}" + " for simplify parameter"
            )    
        
        if self.RicciTensorObj.ricci_fk is None:
            raise Exception(
                "`pycsc.ricci_tensor.RicciTensor.ricci_fk is type None`"
            )
        
        ricci_fk = self.RicciTensorObj.ricci_fk

        if self.MetricTensorObj.config == 'll':
            contra_metric = self.MetricTensorObj.change_config()
        else:
            contra_metric = self.MetricTensorObj.tensor

        ans = 0
        for k in range(self.PyCSCObj.num_coordinates):
            for m in range(self.PyCSCObj.num_coordinates):
                ans += contra_metric[k,m] * ricci_fk[k,m]

        if simplify:
            try:
                ans = sym.simplify(sym.trigsimp(sym.expand_trig(sym.factor(ans))), rational=True)
            except:
                print('Simplification not possible.')

        self.scalar = ans

        if show:
            display(self.scalar)

        return self.scalar

        