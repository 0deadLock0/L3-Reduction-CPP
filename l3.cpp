
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

std::vector<double> multipy_scaler_to_vector(int scaler_value, const std::vector<double>& vector)
{
	const int m = vector.size();
	std::vector<double> scaled_vector(m);
	for(int i = 0 ; i < m ; ++i)
		scaled_vector[i] = scaler_value * vector[i];

	return scaled_vector;
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

template<typename T1, typename T2>
double calculate_gram_schmidt_coefficient(const std::vector<T1>& vector, const std::vector<T2>& normalised_vector)
{
	double numerator = calculate_inner_product(vector, normalised_vector);
	double denominator = calculate_inner_product(normalised_vector, normalised_vector);
	double coefficient = numerator / (denominator + std::numeric_limits<double>::epsilon());

	return coefficient;
}

template<typename T>
std::pair<std::vector<std::vector<double>>, std::vector<double>> calculate_gram_schmidt_coefficients_and_normalised_vector_norms(const std::vector<std::vector<T>>& vectors)
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

	return std::make_pair(gram_schmidt_coefficients, normalised_vector_norms_squared);
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
	std::tie(gram_schmidt_coefficients, normalised_vector_norms_squared) = calculate_gram_schmidt_coefficients_and_normalised_vector_norms(reduced_basis);

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

			std::tie(gram_schmidt_coefficients, normalised_vector_norms_squared) = calculate_gram_schmidt_coefficients_and_normalised_vector_norms(reduced_basis);
		}
		else
			++k;
		// std::cout << "k=" << k << std::endl;
		// std::cout << reduced_basis << std::endl;
		// std::cout << gram_schmidt_coefficients << std::endl;
	}

	return reduced_basis;
}

int main()
{
	// std::vector<std::vector<int>> basis = {{201, 37}, {1648, 297}};
	// std::vector<std::vector<int>> basis = {{1650, 3744}, {164534, 2342}};
	std::vector<std::vector<int>> basis = {{1, 1}, {1, 0}};
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

	return 0;
}
