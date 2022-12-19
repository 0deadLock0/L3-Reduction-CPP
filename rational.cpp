/*
Orginally uploaded at https://gist.github.com/0deadLock0/7b78eda0bce2173710e626e67db5f945
*/


#include <iostream>
#include <cmath>
#include <stdexcept>

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
			Rational(double r, double precision = 1e-5, int64_t cycles = 10)
			{//Reference: https://gist.github.com/mikeando/7073d62385a34a61a6f7

				int64_t Aprev[2] = {1, 0};
				int64_t Bprev[2] = {0, 1};

				double x = r;
				double dif = precision + 1;
				while(cycles > 0 && dif > precision)
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

					double approx = (double)B/(double)A;
					double dif = std::abs(approx - r);

					--cycles;
				}

				this->numerator = Bprev[1];
				this->denominator = Aprev[1];

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
using namespace std;

int main()
{
	Rational r1(2, 3);
	Rational r2(1, 3);

	cout << r1 << " " << r2 << "\n";

	Rational r3 = r2 - 1;
	// r3 -= r2;
	cout << r3 << "\n";

	Rational r4 = float(1.71717);
	cout << r4 << "\n";

	cout << boolalpha << (r4 < r1) << "\n";

	float n = double(r4);
	cout << n << "\n";

	// Rational r5;
	// cin >> r5;
	// cout << r5;

	Rational r6(2, 3);
	r6 = -1 / r6;
	Rational r7(3, 2);
	cout << r6 << " " << r7 << "\n";

	Rational r8 = r6;
	r8 *= r7;
	cout <<  r8 << "\n";

	cout << r8 / -1 << "\n";

	cout << r1 << " " << floor(r1) << " " << round(r1) << " " << ceil(r1) << "\n";
	cout << r2 << " " << floor(r2) << " " << round(r2) << " " << ceil(r2) << "\n";
	cout << r3 << " " << floor(r3) << " " << round(r3) << " " << ceil(r3) << "\n";
	// cout << r5 << " " << floor(r5) << " " << round(r5) << " " << ceil(r5) << "\n";
	cout << r6 << " " << floor(r6) << " " << round(r6) << " " << ceil(r6) << "\n";
	cout << r7 << " " << floor(r7) << " " << round(r7) << " " << ceil(r7) << "\n";

	cout << r6 << " " << abs(r6) << "\n";

	cout << r8 << " " << pow(r8, 4) << " " << pow(r8, -3) << "\n";

	return 0;
}
