# matlab-diamond-norm

This is a MATLAB&reg; function that computes the diamond norm distance between two completely positive trace-preserving (CPTP) superoperators. The calculation is based on a simplified semidefinite program specified in [arXiv:0901.4709](http://arxiv.org/abs/0901.4709) [J. Watrous, [Theory of Computing 5, 11, pp. 217-238 (2009)](http://theoryofcomputing.org/articles/v005a011/)]. More explicitly, it is specified in page 11 of the arXiv article, or page 231 of the published article.

**Note:** the function will not give the correct answer on arbitrary superoperators -- it is only valid on superoperators corresponding to the difference between two CPTP superoperators.

The input is the Choi matrix corresponding to the superoperator in question. The Choi matrix follows the convention described in [J. Watrous, [Theory of Computing 5, 11, pp. 217-238 (2009)](http://theoryofcomputing.org/articles/v005a011/)], i.e., the partial trace over the first half of the Choi matrix yields the zero matrix for a superoperator that is the difference of two CPTP maps, or the identity matrix for a CPTP map (the Jamiolkoski convention, on the other hand, would lead to a matrix proportional to the identity matrix, but with trace 1).

## Installation

Simply add the `src` directory to your MATLAB path.

## Dependencies

This library depends on the [CVX](http://cvxr.com/cvx/) library, and will not function without it. You may find instructions on installing the CVX library [here](http://web.cvxr.com/cvx/doc/install.html).

The test function `dnormtest` also depends on [QIP.m](http://github.com/marcusps/QIP.m), but this library is not necessary to run the `dnorm` function.

## Authors

Marcus P. da Silva (@marcusps)

## License

Copyright 2014 Raytheon BBN Technologies

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
