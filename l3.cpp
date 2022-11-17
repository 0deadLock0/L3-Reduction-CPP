
#include <iostream>
#include <cassert>
#include <vector>
#include <utility>
#include <tuple>
#include <cmath>

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

long long int custom_round(double val)
{
	return std::ceil(val - 1.0/2);
}

std::vector<double> subtract_vectors(const std::vector<double>& vector1, const std::vector<double>& vector2)
{
	const int m1 = vector1.size();
	const int m2 = vector2.size();
	assert(m1 == m2);

	std::vector<double> result(m1);
	for(int i = 0 ; i < m1 ; ++i)
		result[i] = vector1[i] - vector2[i];

	return result;
}

std::vector<double> add_vectors(const std::vector<double>& vector1, const std::vector<double>& vector2)
{
	const int m1 = vector1.size();
	const int m2 = vector2.size();
	assert(m1 == m2);

	std::vector<double> result(m1);
	for(int i = 0 ; i < m1 ; ++i)
		result[i] = vector1[i] + vector2[i];

	return result;
}

std::vector<double> multipy_scaler_to_vector(int scaler_value, const std::vector<double>& vector)
{
	const int m = vector.size();
	std::vector<double> scaled_vector(m);
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

double custom_pow(double num, int exponent)
{
	double result = 1;
	while(exponent > 0)
	{
		if(exponent & 1)
			result *= num;
		num *= num;
		exponent >>= 1;
	}
	return result;
}

template<typename T1, typename T2>
double calculate_inner_product(const std::vector<T1>& vector1, const std::vector<T2>& vector2)
{
	const int m1 = vector1.size();
	const int m2 = vector2.size();
	assert(m1 == m2);

	double inner_product = 0;
	for(int i = 0 ; i < m1 ; ++i)
	{
		double terms_product = vector1[i] * vector2[i];
		inner_product += terms_product;
	}

	return inner_product;
}

template<typename T>
double calculate_norm(const std::vector<T>& vector)
{
	double norm = std::sqrt(calculate_inner_product(vector, vector));
	return norm;
}

template<typename T1, typename T2>
double calculate_gram_schmidt_coefficient(const std::vector<T1>& vector, const std::vector<T2>& normalised_vector)
{
	double numerator = calculate_inner_product(vector, normalised_vector);
	double denominator = calculate_inner_product(normalised_vector, normalised_vector);
	double coefficient = numerator / (denominator + std::numeric_limits<double>::epsilon());

	return coefficient;
}

template<typename T>
std::tuple<std::vector<std::vector<double>>, std::vector<std::vector<double>>, std::vector<double>> calculate_gram_schmidt_coefficients_and_normalised_vector_norms(const std::vector<std::vector<T>>& vectors)
{
	const int m = vectors.size();

	std::vector<std::vector<double>> normalised_vectors(m, std::vector<double>(m));

	std::vector<std::vector<double>> gram_schmidt_coefficients(m, std::vector<double>(m));
	for(int i = 0 ; i < m ; ++i)
		gram_schmidt_coefficients[i][i] = 1;

	for(int i = 0 ; i < m ; ++i)
	{
		for(int j = 0 ; j < m ; ++j)
			normalised_vectors[i][j] = vectors[i][j];
		for(int j = 0 ; j < i ; ++j)
		{
			gram_schmidt_coefficients[i][j] = calculate_gram_schmidt_coefficient(vectors[i], normalised_vectors[j]);
			std::vector<double> projection_vector = multipy_scaler_to_vector(gram_schmidt_coefficients[i][j], normalised_vectors[j]);

			normalised_vectors[i] = subtract_vectors(normalised_vectors[i], projection_vector);
		}
	}

	std::vector<double> normalised_vector_norms_squared(m);
	for(int i = 0 ; i < m ; ++i)
		normalised_vector_norms_squared[i] = calculate_inner_product(normalised_vectors[i], normalised_vectors[i]);

	return std::make_tuple(gram_schmidt_coefficients, normalised_vectors, normalised_vector_norms_squared);
}

void size_reduce_basis_k_and_update_gram_schmidt_coefficients(int k, std::vector<std::vector<double>>* reduced_basis, std::vector<std::vector<double>>* gram_schmidt_coefficients)
{
	// std::cout << *gram_schmidt_coefficients << std::endl;

	const int m = reduced_basis->size();

	for(int j = k - 1 ; j >= 0 ; --j)
	{
		if(std::abs((*gram_schmidt_coefficients)[k][j]) > 1.0/2)
		{
			std::vector<double> scaled_vector = multipy_scaler_to_vector(custom_round((*gram_schmidt_coefficients)[k][j]), (*reduced_basis)[j]);
			(*reduced_basis)[k] = subtract_vectors((*reduced_basis)[k], scaled_vector);

			for(int i = 0 ; i < m ; ++i)
			{
				(*gram_schmidt_coefficients)[k][i] -= custom_round((*gram_schmidt_coefficients)[k][j]) * (*gram_schmidt_coefficients)[j][i];
			}
		}
	}

	// std::cout << *gram_schmidt_coefficients << std::endl;
}

std::vector<std::vector<double>> l3_reduction(const std::vector<std::vector<int>>& basis, double delta)
{
	const int m = basis.size();
	std::vector<std::vector<double>> reduced_basis(m, std::vector<double>(m));
	for(int i = 0 ; i < m ; ++i)
		for(int j = 0 ; j < m ; ++j)
			reduced_basis[i][j] = basis[i][j];

	std::vector<std::vector<double>> gram_schmidt_coefficients;
	std::vector<double> normalised_vector_norms_squared;
	std::tie(gram_schmidt_coefficients, std::ignore, normalised_vector_norms_squared) = calculate_gram_schmidt_coefficients_and_normalised_vector_norms(reduced_basis);

	// std::cout << gram_schmidt_coefficients << std::endl ;
	// std::cout << normalised_vector_norms_squared << std::endl;

	int k = 1; //1 as 0 base indexing
	while( k < m )
	{
		size_reduce_basis_k_and_update_gram_schmidt_coefficients(k, &reduced_basis, &gram_schmidt_coefficients);
		if( delta * normalised_vector_norms_squared[k - 1]  > normalised_vector_norms_squared[k] + std::pow(gram_schmidt_coefficients[k][k-1], 2) * normalised_vector_norms_squared[k - 1] )
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
bool check_condition1(const std::vector<std::vector<double>>& gram_schmidt_coefficients)
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
bool check_condition2(const std::vector<std::vector<double>>& gram_schmidt_coefficients, const std::vector<std::vector<double>>& normalised_vectors, const double& delta)
{
	const int n = normalised_vectors.size();
	for(int i = 0 ; i < n - 1; ++i)
	{
		double left_side_term = delta * calculate_inner_product(normalised_vectors[i], normalised_vectors[i]);
		std::vector<double> right_side_term1 = multipy_scaler_to_vector(gram_schmidt_coefficients[i + 1][i], normalised_vectors[i]);
		std::vector<double> right_side_term2 = add_vectors(right_side_term1, normalised_vectors[i + 1]);
		double right_side_term = calculate_inner_product(right_side_term2, right_side_term2);

		if( left_side_term > right_side_term )
			return false;
	}

	return true;
}

//det(lattice) <= PI_i_1_to_n || b_i(reduced_basis) || <= det(lattice) * (PI_i_1_to_n (1 + (aplha/4) * ((alpha^(i - 1) - 1 ) / (alpha -1)))) ^ 0.5
bool check_condition3(const std::vector<std::vector<int>>& lattice, const std::vector<std::vector<double>>& reduced_basis, const double& alpha)
{//Code needs to be looked for correctness

	const int n = lattice.size();
	long long int determinant = calculate_determinant<long long int, int>(lattice);
	long long int abs_determinant = std::abs(determinant);
	long long int norm_product = 1;
	for(int i = 0 ; i < n ; ++i)
		norm_product *= calculate_norm(reduced_basis[i]);

	if( abs_determinant > norm_product )
		return false;

	double product = 1;
	for(int i = 0 ; i < n ; ++i) //0-based indexing
	{
		double term = 1 + (alpha / 4.0) * ((custom_pow(alpha, i) - 1) / (alpha - 1));
		// std::cout << i << " -> " << alpha << " | " << alpha - 1 << " " << custom_pow(alpha, i) << " : " << term << "\n";
		product *= term;
	}

	double right_side_term = abs_determinant * std::sqrt(product);

	// std::cout << alpha << "\n";
	// std::cout << abs_determinant << " " << norm_product << " " << right_side_term << "\n";

	if( norm_product > right_side_term )
		return false;

	return true;
}

//|| b_j || <= alpha ^ ((i - j)/2) * (1 + (aplha/4) * ((alpha^(i - 1) - 1 ) / (alpha -1))) ^ 0.5 * || ~(b_i) || for all i > j
bool check_condition4(const std::vector<std::vector<double>>& normalised_vectors, const std::vector<std::vector<double>>& reduced_basis, const double& alpha)
{
	const int n = normalised_vectors.size();

	std::vector<double> normalised_vectors_norm(n);
	for(int i = 0 ; i < n ; ++i)
		normalised_vectors_norm[i] = calculate_norm(normalised_vectors[i]);

	std::vector<double> reduced_basis_norm(n);
		for(int i = 0 ; i < n ; ++i)
			reduced_basis_norm[i] = calculate_norm(reduced_basis[i]);

	for(int j = 0 ; j < n ; ++j)
	{
		for(int i = j + 1 ; i < n ; ++i)
		{
			double term1 = std::pow(alpha, (i - j) / 2.0);
			double term2 = 1 + (alpha / 4.0) * ((custom_pow(alpha, i) - 1) / (alpha - 1));
			double term3 = normalised_vectors_norm[i];

			double term = term1 * std::sqrt(term2) * term3;

			if(reduced_basis_norm[j] > term)
				return false;
		}
	}

	return true;
}

// || b_1 || <= alpha^((n - 1) / 4) * det(lattice) ^ 1 / n
bool check_condition5(const std::vector<std::vector<int>>& lattice, const std::vector<std::vector<double>>& reduced_basis, const double& alpha)
{
	const int n = lattice.size();
	long long int determinant = calculate_determinant<long long int, int>(lattice);

	double reduced_basis_0_norm = calculate_norm(reduced_basis[0]);
	double right_side = std::pow(alpha, (n - 1) / 4.0) * std::pow(determinant, n);

	if(reduced_basis_0_norm > right_side)
		return false;

	return true;
}

// || b_1 || <= alpha^((n - 1) / 2) * || x || for every x in lattice
bool check_condition6(const std::vector<std::vector<int>>& lattice, const std::vector<std::vector<double>>& reduced_basis, const double& alpha)
{
	const int n = lattice.size();
	double reduced_basis_0_norm = calculate_norm(reduced_basis[0]);

	for(int i = 0 ; i < n ; ++i)
	{
		double right_side = std::pow(alpha, (n - 1) / 2.0) * calculate_norm(lattice[i]);
		if(reduced_basis_0_norm > right_side)
			return false;
	}

	return true;
}

bool are_properties_satisfied(const std::vector<std::vector<int>>& lattice, const std::vector<std::vector<double>>& reduced_basis, const double& delta)
{
	const double alpha = 1 / (delta - 1.0/4);

	std::vector<std::vector<double>> gram_schmidt_coefficients;
	std::vector<std::vector<double>> normalised_vectors;
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
	std::vector<std::vector<int>> basis = {{201, 37}, {1648, 297}};
	// std::vector<std::vector<int>> basis = {{1650, 3744}, {164534, 2342}};
	// std::vector<std::vector<int>> basis = {{1, 1}, {1, 0}};
	double delta = 0.99;

	const int m = basis.size();

	for(int i = 0 ; i < m ; ++i)
		assert(basis[i].size() == m);

	std::cout << "Input Basis:" << "\n";
	std::cout << basis << "\n";
	std::cout << "delta: " << delta << "\n";
	std::cout << std::endl;

	std::vector<std::vector<double>> reduced_basis = l3_reduction(basis, delta);

	std::cout << "L3 reduced Output Basis:" << "\n";
	std::cout << reduced_basis << "\n";
	std::cout << std::endl;

	bool check_validty = are_properties_satisfied(basis, reduced_basis, delta);
	std::cout << "L3 reduction is " << (check_validty ? "" : "not ") << "valid" << "\n";
	std::cout << std::endl;

	return 0;
}
