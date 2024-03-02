# matrix_coeff_polynomial_arithmetic
## Why the code?
In my work, I get some polynomials of the form (e.g., the [paper](https://doi.org/10.22331/q-2022-11-03-848))
$$F(t, t_1, ..., t_M) = \sum_k C_k t_0^{n_{0,k}} t_1^{n_{1,k}} ... t_M^{n_{M,k}} \exp \left( \sum_{l=0}^{M} \alpha_{k,l} t_l \right).$$
Here, ${C_k}$ are matrices, $n$ are integer numbers, and $\alpha$ are (potentially complex) numbers. I need to perform basic arithmetic operations like addition and (matrix) multiplication of such polynomials. These expressions, for example, can be emerged from calculating the covariance matrix derived from solving two different systems of linear differential equations.

I want the calculation to be reasonably accurate. In principle, this can be easily done in Sympy.
However, in my experience, it can get pretty slow when the number of terms is large.
Therefore, I write my onw implementation in C++. The internal calculation is done in quadruple precision (can be easily changed) with the Boost Multiprecision library.

## Depencencies
* python
* pybind11
* Numpy
* libboost (multiprecision)
* Eigen3

## Compilation and usage
To compile the library, navigate to the `lib/` directory and execute the following command, replacing `<eigen3>` and `<libboost>` with the appropriate paths:

```shell
g++ -O3 -Wall -shared -fPIC $(python3 -m pybind11 --includes) CtexptArith.cpp -o PolyCtexpt$(python3-config --extension-suffix) -std=c++17 -I<eigen3> -I<libboost> -I -march=native -DNDEBUG
```

The library can be called from Python. An example is provided (`example.py`).