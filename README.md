# Lenstra–Lenstra–Lovász Reduction Implementation in C++

> ## Contains:
- [l3.cpp](l3.cpp "Code File") (L3 Reduction for Integral Basis)
- [l3_r.cpp](l3_r.cpp "Code File") (L3 Reduction for Rational/Decimal Basis)
- [rational.cpp](rational.cpp "Code File") (Rational numbers implementation)

> ## Instructions to run either of the file:
### Compile:
``` g++ -std=c++14 fileName.cpp -o fileName```
### Run:
``` fileName ```

> ## User Inputs:
Modifiy the following variables in int main() to change user input as per need.

```
std::vector<std::vector<int>> basis;
basis = {{201, 37}, {1648, 297}};

double delta = 0.99;
```

> ## Paper Refered:
- [Lattice Basis Reduction Improved Practical Algorithm and Solving Subset Sum Problem](https://www.researchgate.net/publication/226499611_Lattice_Basis_Reduction_Improved_Practical_Algorithms_and_Solving_Subset_Sum_Problems "Paper Reference")

> ## Purpose:
Part of Independent Project titled "**An Experimental Study of BKZ Algorithm**""
at [IIIT Delhi](https://www.iiitd.ac.in/)
in Monsoon 2022 Semester
under the guidance of Professor [Subhabrata Samajder](https://www.iiitd.ac.in/subhabrata "Profile")

> ## Author:
- [Abhimanyu Gupta](https://github.com/0deadLock0 "GitHub Profile")

> ## Useful Links:
- [Lattice?](https://en.wikipedia.org/wiki/Lattice_(group) "Wikpedia Article")
- [Why short Orthogonal Vectors?](https://www.esat.kuleuven.be/cosic/blog/lattice-reduction/ "Blog Post")
- [Gram-Schmidt Process](https://en.wikipedia.org/wiki/Gram%E2%80%93Schmidt_process "Wikpedia Article")
- [Short Example](https://www.youtube.com/watch?v=XEMEiBcwSKc "Youtube Video")

## Note:
> - The l3_r.cpp is not running properly.
> - The root cause for this is identified to be the imprecisions in convering `double` values to `Rational` values.
> - The current algorithm (based on idea of Continued Fractions) fails to represent a subset of `double` values using small integer `numerator` and `denominator`.
> - During computations, these numerator and denominator values grows rapidly leading to `integer overflow`, eventually causing code termination.
