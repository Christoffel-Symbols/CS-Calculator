"""
Einstein Tensor

"""

import sympy as sym
from IPython.display import display, Math

from .metric_tensor import MetricTensor
from .pycsc import PyCSC
from .ricci_tensor import RicciTensor
from .ricci_scalar import RicciScalar
from .metric_tensor import MetricTensor

class EinsteinTensor:
    """
    `pycsc.einstein_tensor.EinsteinTensor` class
    """

    def __init__(self, PyCSCObj, MetricTensorObj,RicciTensorObj, RicciScalarObj):
        """
        Create `pycsc.einstein_tensor.EinsteinTensor` instance

        Parameters
        ----------
            PyCSCObj :: `pycsc.pycsc.PyCSC` object
                The `pycsc.pycsc.PyCSC` object containing the information regarding the coordinates.

            RicciTensorObj :: `pycsc.ricci_tensor.RicciTensorObj` object
                The `pycsc.ricci_tensor.RicciTensorObj` object containing the information regarding the Ricci Tensor of first and second kinds.

            RicciScalarObj :: `pycsc.ricci_scalar.RicciScalarObj` object
                The `pycsc.ricci_scalar.RicciScalarObj` object containing the information regarding the Ricci scalar.

        Attributes
        ----------
            PyCSCObj :: `pycsc.pycsc.PyCSC` object
                The `pycsc.pycsc.PyCSC` object containing the information regarding the coordinates.

            RicciTensorObj :: `pycsc.ricci_tensor.RicciTensorObj` object
                The `pycsc.ricci_tensor.RicciTensorObj` object containing the information regarding the Ricci Tensor of first and second kinds.

            RicciScalarObj :: `pycsc.ricci_scalar.RicciScalarObj` object
                The `pycsc.ricci_scalar.RicciScalarObj` object containing the information regarding the Ricci scalar.

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
            "MetricTensroObj must be a `pycsc.metric_tensor.MetricTensor` object"
            )
        
        if not isinstance(RicciTensorObj, RicciTensor):
            raise TypeError(
            "RicciTensorObj must be a `pycsc.ricci_tensor.RicciTensor` object"
            )
        
        if not isinstance(RicciScalarObj, RicciScalar):
            raise TypeError(
            "RicciScalarObj must be a `pycsc.ricci_scalar.RicciScalar` object"
            )
        
        
        self.PyCSCObj = PyCSCObj
        self.RicciTensorObj = RicciTensorObj
        self.RicciScalarObj = RicciScalarObj
        self.MetricTensorObj = MetricTensorObj
        self.tensor = None

    def calculate(self,show=True,simplify=True):
        """
        Calculate Einstein Tensor.

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
            tensor :: `sympy.matrices.dense.MutableDenseMatrix`
                Einstein tensor.

        Returns
        -------
            tensor :: `sympy.matrices.dense.MutableDenseMatrix`
                Einstein tensor.
        
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
        
        if self.RicciScalarObj.scalar is None:
            raise Exception(
                "`pycsc.ricci_scalar.RicciScalar.scalar is type None`"
            )

        self.tensor = sym.zeros(self.PyCSCObj.num_coordinates, self.PyCSCObj.num_coordinates)

        ricci_scalar = self.RicciScalarObj.scalar

        ricci_fk = self.RicciTensorObj.ricci_fk

        if self.MetricTensorObj.config == 'll':
            covariant_metric = self.MetricTensorObj.tensor
        else:
            covariant_metric = self.MetricTensorObj.change_config()

        ans = 0
        for a in range(self.PyCSCObj.num_coordinates):
            for b in range(a+1):
                ans = ricci_fk[a,b] - (1/2)*covariant_metric[a,b]*ricci_scalar

                if simplify:
                    try:
                        ans = sym.simplify(sym.trigsimp(sym.expand_trig(sym.factor(ans))), rational=True)
                    except:
                        print('simplification not possible')

                self.tensor[a,b] = ans
                if a!=b:
                    self.tensor[b,a] = ans

        if show:
            G = sym.symbols(f'G_mu_nu')
            display(Math(sym.latex(G) + ' = ' + sym.latex(self.tensor)))

        return self.tensor