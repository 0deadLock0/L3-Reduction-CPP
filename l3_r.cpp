#include <iostream>
#include <cassert>
#include <vector>
#include <utility>
#include <tuple>
#include <cmath>
#include <limits>

namespace Rationals
{
	int64_t calculate_lcm(int64_t, int64_t);
	int64_t calculate_gcd(int64_t, int64_t);

	class Rational
	{
		private:
			int64_t numerator;
			int64_t denominator;
			void simplify()
			{
				int64_t gcd = calculate_gcd(this->numerator, this->denominator);
				this->numerator /= gcd;
				this->denominator /= gcd;

				if(this->denominator < 0)
				{
					this->numerator *= -1;
					this->denominator *= -1;
				}
			}
		public:
			Rational() : Rational(0)
			{
			}
			Rational(int64_t n, int64_t d)
			{
				this->numerator = n;
				this->denominator = d;
				this->simplify();
			}
			Rational(int n, int d) : Rational((int64_t)n, (int64_t)d)
			{
			}
			Rational(int64_t r) : Rational(r, (int64_t)1)
			{
			}
			Rational(int r) : Rational((int64_t)r)
			{
			}
			Rational(double r, long double precision = 1e-5, int64_t cycles = 10)
			{
				//Reference: https://gist.github.com/mikeando/7073d62385a34a61a6f7
				//Reference: https://stackoverflow.com/a/64828741/12512406

				double eps = 1e-15;

				int sign  = r > 0 ? 1 : -1;
				r = r * sign;

				int64_t Aprev[2] = {1, 0};
				int64_t Bprev[2] = {0, 1};

				double x = r;
				double dif = precision + 1;
				while(cycles > 0 && dif - precision > std::numeric_limits<double>::epsilon())
				{
					int64_t n = std::floor(x);
					x-=n;
					x=1/x;

					int64_t A = Aprev[0] + n * Aprev[1];
					Aprev[0] = Aprev[1];
					Aprev[1] = A;

					int64_t B = Bprev[0] + n * Bprev[1];
					Bprev[0] = Bprev[1];
					Bprev[1] = B;

					if(A < 0  ||  B < 0)
						break;

					this->numerator = Bprev[1];
					this->denominator = Aprev[1];

					double approx = (double)B/(double)A;
					double dif = std::abs(approx - r);

					--cycles;
				}

				this->numerator *= sign;

				this->simplify();
			}
			Rational(const Rational& r) : Rational(r.numerator, r.denominator)
			{
			}

			int64_t get_numerator() const
			{
				return this->numerator;
			}
			int64_t get_denominator() const
			{
				return this->denominator;
			}

			Rational& operator=(const Rational& r2)
			{
				this->numerator = r2.numerator;
				this->denominator = r2.denominator;

				return *this;
			}
			Rational& operator+=(const Rational& r2)
			{
				int64_t lcm = calculate_lcm(this->denominator, r2.denominator);
				int64_t term1 = lcm / this->denominator * this->numerator;
				int64_t term2 = lcm / r2.denominator * r2.numerator;

				this->numerator = term1 + term2;
				this->denominator = lcm;
				this->simplify();

				return *this;
			}
			Rational& operator-=(const Rational& r2)
			{
				Rational r(-r2.numerator, r2.denominator);
				*this += r;

				return *this;
			}
			Rational& operator*=(const Rational& r2)
			{
				int64_t gcd1 = calculate_gcd(this->numerator, r2.denominator);
				int64_t gcd2 = calculate_gcd(this->denominator, r2.numerator);

				this->numerator /= gcd1;
				this->numerator *= (r2.numerator / gcd2);

				this->denominator /= gcd2;
				this->denominator *= (r2.denominator / gcd1);

				this->simplify();

				return *this;
			}
			Rational& operator/=(const Rational& r2)
			{
				Rational r(r2.denominator, r2.numerator);
				*this *= r;

				return *this;
			}

			explicit operator double() const
			{
				return (double)(this->numerator) / this->denominator;
			}


			friend std::istream& operator>>(std::istream& , Rational&);
	};

	std::istream& operator>>(std::istream& in, Rational& r)
	{
		int64_t n,d;
		char c;

		in >> n >> c >> d;

		if( c != '/' )
			throw std::invalid_argument( "Ill format encountered while reading" );

		r.numerator = n;
		r.denominator = d;
		r.simplify();

		return in;
	}
	std::ostream& operator<<(std::ostream& out, const Rational& r)
	{
		out << r.get_numerator() << "/" << r.get_denominator();
		return out;
	}

	Rational operator+(const Rational& r1, const Rational& r2)
	{
		Rational r = r1;
		r += r2;

		return r;
	}
	Rational operator-(const Rational& r1, const Rational& r2)
	{
		Rational r = r1;
		r -= r2;

		return r;
	}
	Rational operator*(const Rational& r1, const Rational& r2)
	{
		Rational r = r1;
		r *= r2;

		return r;
	}
	Rational operator/(const Rational& r1, const Rational& r2)
	{
		Rational r = r1;
		r /= r2;

		return r;
	}

	bool operator<(const Rational& r1, const Rational& r2)
	{
		int64_t lcm = calculate_lcm(r1.get_denominator(), r2.get_denominator());
		int64_t term1 = lcm / r1.get_denominator() * r1.get_numerator();
		int64_t term2 = lcm / r2.get_denominator() * r2.get_numerator();

		return term1 < term2;
	}
	bool operator>(const Rational& r1, const Rational& r2)
	{
		return r2 < r1;
	}
	bool operator<=(const Rational& r1, const Rational& r2)
	{
		return !(r2 < r1);
	}
	bool operator>=(const Rational& r1, const Rational& r2)
	{
		return !(r1 < r2);
	}
	bool operator==(const Rational& r1, const Rational& r2)
	{
		return !(r1 < r2) && !(r2 < r1);
	}
	bool operator!=(const Rational& r1, const Rational& r2)
	{
		return !(r1 == r2);
	}


	int64_t calculate_lcm(int64_t a, int64_t b)
	{
		return a * (b / calculate_gcd(a, b));
	}
	int64_t calculate_gcd(int64_t a, int64_t b)
	{
		while(b != 0)
		{
			int64_t r = a % b;
			a = b;
			b = r;
		}
		return a;
	}

	Rational floor(const Rational& r)
	{
		int64_t rem = r.get_numerator() % r.get_denominator();
		if( rem < 0 )
			rem = r.get_denominator() + rem;
		Rational r2(r.get_numerator() - rem, r.get_denominator());
		return r2;
	}
	Rational ceil(const Rational& r)
	{
		Rational r2(r.get_numerator() - 1, r.get_denominator());
		return floor(r2) + 1;
	}
	Rational round(const Rational& r)
	{
		Rational r2 = r - Rational(1, 2);
		return ceil(r2);
	}

	Rational abs(const Rational& r)
	{
		Rational abs_r(std::abs(r.get_numerator()), r.get_denominator());
		return abs_r;
	}

	Rational pow(Rational num, int exponent)
	{
		if( exponent < 0 )
		{
			exponent *= -1;
			num = 1 / num;
		}

		Rational result(1);
		while(exponent > 0)
		{
			if(exponent & 1)
				result *= num;
			num *= num;
			exponent >>= 1;
		}
		return result;
	}
}

using namespace Rationals;

template<typename T>
std::ostream& operator<<(std::ostream &out, const std::vector<T> &vec)
{
	size_t size = vec.size();
	out << "(";
	for(size_t i = 0 ; i < size ; ++i)
	{
		out << vec[i];
		if(i != size - 1)
			out << ",";
	}
	out << ")";
	return out;
}

std::vector<Rational> subtract_vectors(const std::vector<Rational>& vector1, const std::vector<Rational>& vector2)
{
	const int m1 = vector1.size();
	const int m2 = vector2.size();
	assert(m1 == m2);

	std::vector<Rational> result(m1);
	for(int i = 0 ; i < m1 ; ++i)
		result[i] = vector1[i] - vector2[i];

	return result;
}

std::vector<Rational> add_vectors(const std::vector<Rational>& vector1, const std::vector<Rational>& vector2)
{
	const int m1 = vector1.size();
	const int m2 = vector2.size();
	assert(m1 == m2);

	std::vector<Rational> result(m1);
	for(int i = 0 ; i < m1 ; ++i)
		result[i] = vector1[i] + vector2[i];

	return result;
}

std::vector<Rational> multipy_scaler_to_vector(Rational scaler_value, const std::vector<Rational>& vector)
{
	const int m = vector.size();
	std::vector<Rational> scaled_vector(m);
	for(int i = 0 ; i < m ; ++i)
		scaled_vector[i] = scaler_value * vector[i];

	return scaled_vector;
}

//Assuming Matirx is Square
template<typename T1, typename T2>
T1 calculate_determinant(const std::vector<std::vector<T2>>& matrix)
{
	const int n = matrix.size();
	if(n == 1)
		return matrix[0][0];
	else if(n == 2)
		return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];
	else
	{
		T1 determinant = 0;
		for(int col_i = 0, sign = 1 ; col_i < n ; ++col_i, sign = -sign)
		{
			std::vector<std::vector<T2>> sub_matrix(n - 1, std::vector<T2>(n - 1));
			for(int row_i = 1, subrow = 0 ; row_i < n ; ++row_i, ++subrow)
			{
				for(int col_j = 0, subcol = 0 ; col_j < n ; ++col_j)
				{
					if(col_j == col_i)
						continue;
					sub_matrix[subrow][subcol] = matrix[row_i][col_j];
					++subcol;
				}
			}

			T1 sub_determinant = matrix[0][col_i];
			if(sub_determinant != 0)
				sub_determinant *= sign * calculate_determinant<T1, T2>(sub_matrix);
			determinant += sub_determinant;
		}

		return determinant;
	}
}

template<typename T1, typename T2>
Rational calculate_inner_product(const std::vector<T1>& vector1, const std::vector<T2>& vector2)
{
	const int m1 = vector1.size();
	const int m2 = vector2.size();
	assert(m1 == m2);

	Rational inner_product = 0;
	for(int i = 0 ; i < m1 ; ++i)
	{
		Rational terms_product = vector1[i] * vector2[i];
		inner_product += terms_product;
	}

	return inner_product;
}

template<typename T>
Rational calculate_norm(const std::vector<T>& vector)
{
	Rational norm = std::sqrt(double(calculate_inner_product(vector, vector)));
	return norm;
}

template<typename T1, typename T2>
Rational calculate_gram_schmidt_coefficient(const std::vector<T1>& vector, const std::vector<T2>& normalised_vector)
{
	Rational numerator = calculate_inner_product(vector, normalised_vector);
	Rational denominator = calculate_inner_product(normalised_vector, normalised_vector);
	Rational coefficient = numerator / denominator;

	return coefficient;
}

template<typename T>
std::tuple<std::vector<std::vector<Rational>>, std::vector<std::vector<Rational>>, std::vector<Rational>> calculate_gram_schmidt_coefficients_and_normalised_vector_norms(const std::vector<std::vector<T>>& vectors)
{
	const int m = vectors.size();

	std::vector<std::vector<Rational>> normalised_vectors(m, std::vector<Rational>(m));

	std::vector<std::vector<Rational>> gram_schmidt_coefficients(m, std::vector<Rational>(m));
	for(int i = 0 ; i < m ; ++i)
		gram_schmidt_coefficients[i][i] = 1;

	for(int i = 0 ; i < m ; ++i)
	{
		for(int j = 0 ; j < m ; ++j)
			normalised_vectors[i][j] = vectors[i][j];
		for(int j = 0 ; j < i ; ++j)
		{
			gram_schmidt_coefficients[i][j] = calculate_gram_schmidt_coefficient(vectors[i], normalised_vectors[j]);
			std::vector<Rational> projection_vector = multipy_scaler_to_vector(gram_schmidt_coefficients[i][j], normalised_vectors[j]);

			normalised_vectors[i] = subtract_vectors(normalised_vectors[i], projection_vector);
		}

	}

	std::vector<Rational> normalised_vector_norms_squared(m);
	for(int i = 0 ; i < m ; ++i)
		normalised_vector_norms_squared[i] = calculate_inner_product(normalised_vectors[i], normalised_vectors[i]);

	return std::make_tuple(gram_schmidt_coefficients, normalised_vectors, normalised_vector_norms_squared);
}

void size_reduce_basis_k_and_update_gram_schmidt_coefficients(int k, std::vector<std::vector<Rational>>* reduced_basis, std::vector<std::vector<Rational>>* gram_schmidt_coefficients)
{
	// std::cout << *gram_schmidt_coefficients << std::endl;

	const int m = reduced_basis->size();

	for(int j = k - 1 ; j >= 0 ; --j)
	{
		if(abs((*gram_schmidt_coefficients)[k][j]) > Rational(1, 2))
		{
			std::vector<Rational> scaled_vector = multipy_scaler_to_vector(round((*gram_schmidt_coefficients)[k][j]), (*reduced_basis)[j]);
			(*reduced_basis)[k] = subtract_vectors((*reduced_basis)[k], scaled_vector);

			for(int i = 0 ; i < m ; ++i)
			{
				(*gram_schmidt_coefficients)[k][i] -= round((*gram_schmidt_coefficients)[k][j]) * (*gram_schmidt_coefficients)[j][i];
			}
		}
	}

	// std::cout << *gram_schmidt_coefficients << std::endl;
}

std::vector<std::vector<Rational>> l3_reduction(const std::vector<std::vector<Rational>>& basis, Rational delta)
{
	const int m = basis.size();
	std::vector<std::vector<Rational>> reduced_basis(m, std::vector<Rational>(m));
	for(int i = 0 ; i < m ; ++i)
		for(int j = 0 ; j < m ; ++j)
			reduced_basis[i][j] = basis[i][j];

	std::vector<std::vector<Rational>> gram_schmidt_coefficients;
	std::vector<Rational> normalised_vector_norms_squared;
	std::tie(gram_schmidt_coefficients, std::ignore, normalised_vector_norms_squared) = calculate_gram_schmidt_coefficients_and_normalised_vector_norms(reduced_basis);

	std::cout << gram_schmidt_coefficients << std::endl ;
	std::cout << normalised_vector_norms_squared << std::endl;

	int k = 1; //1 as 0 base indexing
	while( k < m )
	{
		size_reduce_basis_k_and_update_gram_schmidt_coefficients(k, &reduced_basis, &gram_schmidt_coefficients);
		if( delta * normalised_vector_norms_squared[k - 1]  > normalised_vector_norms_squared[k] + pow(gram_schmidt_coefficients[k][k-1], 2) * normalised_vector_norms_squared[k - 1] )
		{
			std::swap(reduced_basis[k], reduced_basis[k - 1]);
			k = std::max(k - 1, 1);

			std::tie(gram_schmidt_coefficients, std::ignore, normalised_vector_norms_squared) = calculate_gram_schmidt_coefficients_and_normalised_vector_norms(reduced_basis);
		}
		else
			++k;

		// std::cout << "k=" << k << std::endl;
		// std::cout << reduced_basis << std::endl;
		// std::cout << gram_schmidt_coefficients << std::endl;
	}

	return reduced_basis;
}


//|u_(i,j)| <= 1/2 for 1 <= j < i <= n
bool check_condition1(const std::vector<std::vector<Rational>>& gram_schmidt_coefficients)
{
	const int n = gram_schmidt_coefficients.size();
	for(int i = 0 ; i < n ; ++i)
	{
		for(int j = 0 ; j < i ; ++j)
		{
			if( gram_schmidt_coefficients[i][j] > 1.0/2 )
				return false;
		}
	}

	return true;
}

// delta * || ~(b_i) ||^2 <= || u_(i+1, i) * ~(b_i) + ~(b_i+1) ||^2
bool check_condition2(const std::vector<std::vector<Rational>>& gram_schmidt_coefficients, const std::vector<std::vector<Rational>>& normalised_vectors, const Rational& delta)
{
	const int n = normalised_vectors.size();
	for(int i = 0 ; i < n - 1; ++i)
	{
		Rational left_side_term = delta * calculate_inner_product(normalised_vectors[i], normalised_vectors[i]);
		std::vector<Rational> right_side_term1 = multipy_scaler_to_vector(gram_schmidt_coefficients[i + 1][i], normalised_vectors[i]);
		std::vector<Rational> right_side_term2 = add_vectors(right_side_term1, normalised_vectors[i + 1]);
		Rational right_side_term = calculate_inner_product(right_side_term2, right_side_term2);

		if( left_side_term > right_side_term )
			return false;
	}

	return true;
}

//det(lattice) <= PI_i_1_to_n || b_i(reduced_basis) || <= det(lattice) * (PI_i_1_to_n (1 + (aplha/4) * ((alpha^(i - 1) - 1 ) / (alpha -1)))) ^ 0.5
bool check_condition3(const std::vector<std::vector<Rational>>& lattice, const std::vector<std::vector<Rational>>& reduced_basis, const Rational& alpha)
{//Code needs to be looked for correctness

	const int n = lattice.size();
	Rational determinant = calculate_determinant<Rational, Rational>(lattice);
	Rational abs_determinant = abs(determinant);
	Rational norm_product = 1;
	for(int i = 0 ; i < n ; ++i)
		norm_product *= calculate_norm(reduced_basis[i]);

	if( abs_determinant > norm_product )
		return false;

	Rational product = 1;
	for(int i = 0 ; i < n ; ++i) //0-based indexing
	{
		Rational term = 1 + (alpha / 4.0) * ((pow(alpha, i) - 1) / (alpha - 1));
		// std::cout << i << " -> " << alpha << " | " << alpha - 1 << " " << pow(alpha, i) << " : " << term << "\n";
		product *= term;
	}

	Rational right_side_term = abs_determinant * std::sqrt(double(product));

	// std::cout << alpha << "\n";
	// std::cout << abs_determinant << " " << norm_product << " " << right_side_term << "\n";

	if( norm_product > right_side_term )
		return false;

	return true;
}

//|| b_j || <= alpha ^ ((i - j)/2) * (1 + (aplha/4) * ((alpha^(i - 1) - 1 ) / (alpha -1))) ^ 0.5 * || ~(b_i) || for all i > j
bool check_condition4(const std::vector<std::vector<Rational>>& normalised_vectors, const std::vector<std::vector<Rational>>& reduced_basis, const Rational& alpha)
{
	const int n = normalised_vectors.size();

	std::vector<Rational> normalised_vectors_norm(n);
	for(int i = 0 ; i < n ; ++i)
		normalised_vectors_norm[i] = calculate_norm(normalised_vectors[i]);

	std::vector<Rational> reduced_basis_norm(n);
		for(int i = 0 ; i < n ; ++i)
			reduced_basis_norm[i] = calculate_norm(reduced_basis[i]);

	for(int j = 0 ; j < n ; ++j)
	{
		for(int i = j + 1 ; i < n ; ++i)
		{
			Rational term1 = std::pow(double(alpha), (i - j) / 2.0);
			Rational term2 = 1 + (alpha / 4.0) * ((pow(alpha, i) - 1) / (alpha - 1));
			Rational term3 = normalised_vectors_norm[i];

			Rational term = term1 * std::sqrt(double(term2)) * term3;

			if(reduced_basis_norm[j] > term)
				return false;
		}
	}

	return true;
}

// || b_1 || <= alpha^((n - 1) / 4) * det(lattice) ^ 1 / n
bool check_condition5(const std::vector<std::vector<Rational>>& lattice, const std::vector<std::vector<Rational>>& reduced_basis, const Rational& alpha)
{
	const int n = lattice.size();
	Rational determinant = calculate_determinant<Rational, Rational>(lattice);

	Rational reduced_basis_0_norm = calculate_norm(reduced_basis[0]);
	Rational right_side = std::pow(double(alpha), (n - 1) / 4.0) * pow(determinant, n);

	if(reduced_basis_0_norm > right_side)
		return false;

	return true;
}

// || b_1 || <= alpha^((n - 1) / 2) * || x || for every x in lattice
bool check_condition6(const std::vector<std::vector<Rational>>& lattice, const std::vector<std::vector<Rational>>& reduced_basis, const Rational& alpha)
{
	const int n = lattice.size();
	Rational reduced_basis_0_norm = calculate_norm(reduced_basis[0]);

	for(int i = 0 ; i < n ; ++i)
	{
		Rational right_side = std::pow(double(alpha), (n - 1) / 2.0) * calculate_norm(lattice[i]);
		if(reduced_basis_0_norm > right_side)
			return false;
	}

	return true;
}

bool are_properties_satisfied(const std::vector<std::vector<Rational>>& lattice, const std::vector<std::vector<Rational>>& reduced_basis, const Rational& delta)
{
	const Rational alpha = 1 / (delta - 1.0/4);

	std::vector<std::vector<Rational>> gram_schmidt_coefficients;
	std::vector<std::vector<Rational>> normalised_vectors;
	std::tie(gram_schmidt_coefficients, normalised_vectors, std::ignore) = calculate_gram_schmidt_coefficients_and_normalised_vector_norms(reduced_basis);

	bool condition1 = check_condition1(gram_schmidt_coefficients);
	if( !condition1 )
		return false;

	bool condition2 = check_condition2(gram_schmidt_coefficients, normalised_vectors, delta);
	if( !condition2 )
		return false;

	bool condition3 = check_condition3(lattice, reduced_basis, alpha);
	if( !condition3 )
		return false;

	bool condition4 = check_condition4(normalised_vectors, reduced_basis, alpha);
	if( !condition4 )
		return false;

	bool condition5 = check_condition5(lattice, reduced_basis, alpha);
	if( !condition5 )
		return false;

	bool condition6 = check_condition6(lattice, reduced_basis, alpha);
	if( !condition6 )
		return false;

	return true;
}

int main()
{
	std::vector<std::vector<Rational>> basis;
	// basis = {{201, 37}, {1648, 297}};
	// basis = {{1650, 3744}, {164534, 2342}};
	// basis = {{1, 1}, {1, 0}};

	// basis = std::vector<std::vector<Rational>>(10, std::vector<Rational>(10));
	// for(int i = 0 ; i < 10 ; ++i)
	// 	basis[i][i] = 1;
	// basis[2][2] = 2;

	// basis = {{2, 3, 4, 1}, {1, 1, -1, 1}, {-1, 1, 9, 5}, {12, 7, 3, 5}};

	basis = {{1121, 343, -4354, 9746}, {34342, 98341, -13224, -56942}, {1, 1, 1, 1}, {-53, 43, -932, -1}};

	Rational delta(99,100);

	const int m = basis.size();

	for(int i = 0 ; i < m ; ++i)
		assert(basis[i].size() == m);

	std::cout << "Input Basis:" << "\n";
	std::cout << basis << "\n";
	std::cout << "delta: " << delta << "\n";
	std::cout << std::endl;

	std::vector<std::vector<Rational>> reduced_basis = l3_reduction(basis, delta);

	std::cout << "L3 reduced Output Basis:" << "\n";
	std::cout << reduced_basis << "\n";
	std::cout << std::endl;

	bool check_validty = are_properties_satisfied(basis, reduced_basis, delta);
	std::cout << "L3 reduction is " << (check_validty ? "" : "not ") << "valid" << "\n";
	std::cout << std::endl;

	return 0;
}
