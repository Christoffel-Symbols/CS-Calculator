# Python Christoffel Symbols Calculator (PyCSC)

All rights reserved 2024-25.

**This package is actively under development**

This Python package is for calculating Christoffel Symbols and tensors (i.e., Riemann Tensor, Ricci Tensor, etcetera) that are helpful in doing gravitational and relativistic astrophysics. With time, the package is expected to mature, both in terms of efficiency and functionality.

## Table of Contents

1. [Installation](#installtion)
2. This package comes with a Graphical User Interface (GUI) [Christoffel Symbols Calculator frontend](https://github.com/Christoffel-Symbols/Christoffel-Symbols-Calculator-frontend) repository.
3. [Known issues](#known-issues)
4. [Roadmap](#roadmap)
5. [Questions, Issues, Suggestions, and Other Feedback](#questions-issues)

## Installation

- Option 1 (Pip install)

```bash
pip install git+https://github.com/Christoffel-Symbols/Christoffel-Symbols-Calculator.git
```

- Option 2 (Clone the repository)

```bash
git clone https://github.com/Christoffel-Symbols/Christoffel-Symbols-Calculator.git
```

After cloning the repository, run the following command in the Christoffel-Symbols-Calculator folder

```bash
pip install .
```

or

```bash
python setup.py install
```

## Known Issues

- Need to develop a parser that can render pythonnic mathematical expression to latex in real-time. This, in part, can be achieved using shunting yard algorithm.
- Need to validate matrix elements before submitting a request.
- Put a restriction on sending a zero n-dimensional matrix.
- Catch the error when contravariant metric tensor cannot be computed.
- Incoporate latex style text when filling out the metric tensor.
- Remove `real` parameter when generating symbols using `sympy.core.symbols.Symbol` and run tests (scale factor 'at(t)' is causing an issue).
- `alpha`, `delta`, and `epsilon` reserved parameters don't work when a request is sent from the Graphical UI.

## Roadmap

Here are some future plans for the Christoffel Symbols Calculator:

- Include electromagnetism tensor to describe electromagnetic field in spacetime.

## Questions

If you have any feedback, suggestions, or new ideas on how to improve the software, plaese feel free to reach out to me through email ([dhananjhay03@gmail.com](mailto:dhananjhay03@gmail.com)) or the [discussion page](https://github.com/Christoffel-Symbols/Christoffel-Symbols-Calculator/discussions).
