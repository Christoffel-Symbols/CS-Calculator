"""
Python Christoffel Symbols Calculator (PyCSC)
-------------------------------------------------------

A sophisticated and fast Python package to calculate a variety of tensors helpful
in doing gravitational and relativistic astrophysics. 
See the [`Christoffel-Symbols-Calculator-frontend`](https://github.com/Christoffel-Symbols) 
GitHub repository for a Graphical User Interface (GUI) to complement this package.

Author: Dhananjhay Bansal
Contact: dhananjhay03@gmail.com

---


    GNU General Public License v3 (GNU GPLv3)

    (c) 2024. All rights reserved                  

    Dhananjhay Bansal disclaims any warranties, expressed, implied, or statutory, of any kind with rrespect to the software, 
    including without limitation any warranty of merchantability or fitness for a particular purpose. Dhananjhay Bansal shall 
    not be liable in any event for any damages, whether direct or indirect, special or general, consequential or incidental, 
    arising from the use of the software. The name of the Dhananjhay Bansal may be used to endorse or promote products derived 
    from this software without specific prior written permission.                  
                                        

    This file is part of the Christoffel-Symbols project.              

    Christoffel-Symbols is free software:  you can redistribute it and/or modify it under the terms of the GNU General Public 
    License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.                   

    Christoffel-Symbols is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied 
    warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.                        

    You should have received a copy of the GNU General Public License along with Christoffel-Symbols project. 
    If not, see<http://www.gnu.org/licenses/>.

"""

from setuptools import setup

long_description = __doc__.strip()

setup(
    name="PyCSC",
    version="2.1.4",
    description=long_description,
    url="https://github.com/Christoffel-Symbols/Christoffel-Symbols-Calculator",
    author="Dhananjhay Bansal",
    maintainer="Dhananjhay Bansal",
    maintainer_email="dhananjhay03@gmail.com",
    packages=[
        "pyCSC"
    ],
    install_requires=["sympy","IPython"],
    license="GPLv3",
    python_requies=">=3.9",
    platforms=["Ubuntu"] # Tested on Ubuntu and MacOS. Windows likely okay.
)
