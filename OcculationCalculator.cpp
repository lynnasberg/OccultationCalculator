// OcculationCalculator.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <xmmintrin.h>
#include <immintrin.h>
#include <time.h>
#include <ctime>
#include <iomanip>
#include <chrono>
#pragma warning(disable : 4996)

struct double4
{
	union {
		__m256d data4;
		double data[4];
	};
};

double4 _radii[9];
double4 _mass[9];
double4 _x[9];
double4 _v[9];
const char* names[9];

double G = 6.67408E-11; // from https://ssd.jpl.nasa.gov/?planet_phys_par
time_t currentTime;

void SetData()
{
	for ( int i = 0; i < 9; i++ )
	{
		_x[i].data[3] = 0;
		_v[i].data[3] = 0;
	}

	// masses from https://ssd.jpl.nasa.gov/?planet_phys_par

	// sun
	names[0] = "Sun";
	_mass[0].data4 = _mm256_set1_pd( G * 1.98847E30 );
	_radii[0].data4 = _mm256_set1_pd( 696340E3 );
	_x[0].data[0] = 5.064702305943604E8;
	_x[0].data[1] = 5.889102377519094E8;
	_x[0].data[2] = -2.335506963855098E7;

	_v[0].data[0] = -4.528770171532725f;
	_v[0].data[1] = 1.162145548819547E1f;
	_v[0].data[2] = 9.255801942556659E-2f;

	// mercury
	names[1] = "Mercury";
	_mass[1].data4 = _mm256_set1_pd( G * 0.330114E24 );
	_radii[1].data4 = _mm256_set1_pd( 2440E3 );
	_x[1].data[0] = -2.094166757815767E10;
	_x[1].data[1] = 4.303122811605530E10;
	_x[1].data[2] = 5.412389196360279E9;

	_v[1].data[0] = -5.327883840950205E4;
	_v[1].data[1] = -2.011402533603810E4;
	_v[1].data[2] = 3.243060694058319E3;

	// venus
	names[2] = "Venus";
	_mass[2].data4 = _mm256_set1_pd( G * 4.86747E24 );
	_radii[2].data4 = _mm256_set1_pd( 6052E3 );
	_x[2].data[0] = 7.041720536543189E10f;
	_x[2].data[1] = 8.299087742435262E10f;
	_x[2].data[2] = -2.927748819523938E9f;

	_v[2].data[0] = -2.681660121315567E4f;
	_v[2].data[1] = 2.251800819839024E4f;
	_v[2].data[2] = 1.855924201511814E3f;

	// earth
	names[3] = "Earth";
	_mass[3].data4 = _mm256_set1_pd( G * ( 5.97237E24 + 7.3476731E22 ) );
	_radii[3].data4 = _mm256_set1_pd( 6371E3 );
	_x[3].data[0] = -2.636334964502048E10;
	_x[3].data[1] = 1.452193717485576E11;
	_x[3].data[2] = -2.884341217565536E7;

	_v[3].data[0] = -2.978625055952658E4;
	_v[3].data[1] = -5.548408871854586E3;
	_v[3].data[2] = 1.465217303959765;

	// mars
	names[4] = "Mars";
	_mass[4].data4 = _mm256_set1_pd( G * 0.641712E24 );
	_radii[4].data4 = _mm256_set1_pd( 3390E3 );
	_x[4].data[0] = 2.031671350128849E11;
	_x[4].data[1] = 5.846457202241653E10;
	_x[4].data[2] = -3.784260025553051E9;

	_v[4].data[0] = -5.730507905139778E3;
	_v[4].data[1] = 2.538300685947983E4;
	_v[4].data[2] = 6.722935979870517E2;

	// jupiter
	names[5] = "Jupiter";
	_mass[5].data4 = _mm256_set1_pd( G * ( 1898.187E24 + 3.9310888E23 ) );
	_radii[5].data4 = _mm256_set1_pd( 69911E3 );
	_x[5].data[0] = -8.012982678594466E11;
	_x[5].data[1] = -1.509045181233322E11;
	_x[5].data[2] = 1.854711461009607E10;

	_v[5].data[0] = 2.266628102575662E3;
	_v[5].data[1] = -1.222160269052195E4;
	_v[5].data[2] = 9.014396886453113E-2;

	// saturn
	names[6] = "Saturn";
	_mass[6].data4 = _mm256_set1_pd( G * ( 568.3174E24 + 1.4054061E23 ) );
	_radii[6].data4 = _mm256_set1_pd( 58232E3 );
	_x[6].data[0] = -2.790165349370531E11;
	_x[6].data[1] = -1.475897874611040E12;
	_x[6].data[2] = 3.676627677393824E10;

	_v[6].data[0] = 8.960458302019051E3;
	_v[6].data[1] = -1.825209783284202E3;
	_v[6].data[2] = -3.246327508861824E2;

	// uranus
	names[7] = "Uranus";
	_mass[7].data4 = _mm256_set1_pd( G * ( 86.8127E24 + 8.8901806E21 ) );
	_radii[7].data4 = _mm256_set1_pd( 25362E3 );
	_x[7].data[0] = 2.744017451075625E12;
	_x[7].data[1] = 1.171374128846561E12;
	_x[7].data[2] = -3.119877803055137E10;

	_v[7].data[0] = -2.723695774691808E3;
	_v[7].data[1] = 5.945761125340464E3;
	_v[7].data[2] = 5.712391934823868E1;

	// neptune
	names[8] = "Neptune";
	_mass[8].data4 = _mm256_set1_pd( G * ( 102.4126E24 + 2.1487923E22 ) );
	_radii[8].data4 = _mm256_set1_pd( 24622 );
	_x[8].data[0] = 4.239524025371338E12;
	_x[8].data[1] = -1.449525865099807E12;
	_x[8].data[2] = -6.785384248451066E10;

	_v[8].data[0] = 1.723050820526497E3;
	_v[8].data[1] = 5.175614088548725E3;
	_v[8].data[2] = -1.465667215811188E2;
}

__m256d GetLength( __m256d a )
{
	double4 d4;
	d4.data4 = _mm256_mul_pd( a, a ); // squared length
	return _mm256_sqrt_pd( _mm256_set1_pd( d4.data[0] + d4.data[1] + d4.data[2] ) ); // length
}

__m256d GetCubedLength( __m256d a )
{
	double4 d4;
	d4.data4 = _mm256_mul_pd( a, a ); // square differences
	d4.data4 = _mm256_set1_pd( d4.data[0] + d4.data[1] + d4.data[2] ); // square length
	return _mm256_mul_pd( d4.data4, _mm256_sqrt_pd( d4.data4 ) );
}

float year = 1900;
bool occultation[81];
int main()
{
	SetData();

	currentTime = time_t( 1483225200 );

	float timestep = 60;
	__m256d timestep4 = _mm256_set1_pd( timestep );

	for ( int step = 0; step < 1E9f; step++ )
	{
		currentTime += time_t( timestep );

		//
		// Gravity simulation
		//

		// update velocities
		for ( int i = 0; i < 9; i++ )
		{
			// loop over planets, calculate gravitational force, update velocities
			for ( int j = 0; j < 9; j++ )
			{
				if ( i == j ) continue;
				__m256d distanceVector = _mm256_sub_pd( _x[j].data4, _x[i].data4 );
				__m256d distanceCubed = GetCubedLength( distanceVector );
				__m256d acceleration = _mm256_mul_pd( _mm256_div_pd( _mass[j].data4, distanceCubed ), distanceVector );
				_v[i].data4 = _mm256_add_pd( _v[i].data4, _mm256_mul_pd( acceleration, timestep4 ) );
			}
		}

		// update positions
		for ( int i = 0; i < 9; i++ )
		{
			_x[i].data4 = _mm256_add_pd( _x[i].data4, _mm256_mul_pd( _v[i].data4, timestep4 ) );
		}

		//
		// check if there is an occultation by constructing a vector from the earth to each planet, then checking the angles between them
		//

		__m256d vectors[9];
		double angularRadii[9];

		// construct vectors and angular radii
		for ( int i = 0; i < 9; i++ )
		{
			if ( i == 3 ) continue;
			vectors[i] = _mm256_sub_pd( _x[i].data4, _x[3].data4 );

			// normalize vector
			double4 distance;
			distance.data4 = GetLength( vectors[i] );
			vectors[i] = _mm256_div_pd( vectors[i], distance.data4 );
			angularRadii[i] = _radii[i].data[0] / distance.data[0]; // approximated for distance >> radius
		}

		// get angles
		double angles[81];
		for ( int i = 0; i < 9; i++ )
		{
			if ( i == 3 || i == 0 ) continue; // no occulation with earth or sun
			for ( int j = 0; j < 9; j++ )
			{
				if ( j == 3 || j == 0 || i == j ) continue; // no occultation with earth, sun or self
				int n = i * 9 + j;

				// calculate angle between the vectors to the two planets using a dot product
				double4 dot4;
				dot4.data4 = _mm256_mul_pd( vectors[i], vectors[j] );
				dot4.data[3] = dot4.data[0] + dot4.data[1] + dot4.data[2];
				dot4.data4 = _mm256_acos_pd( dot4.data4 );
				angles[n] = dot4.data[3];

				if ( angles[n] < angularRadii[i] + angularRadii[j] )
				{
					if ( !occultation[n] )
					{
						// occultation starts
						std::cout << "Occultation between " << names[i] << " and " << names[j] << " at " << std::put_time( std::localtime( &currentTime ), "%F %T" ) << std::endl;
						occultation[n] = true;
					}
				}
				else
				{
					occultation[n] = false;
				}
			}
		}

		float newYear = gmtime( &currentTime )->tm_year + 1900;
		if ( newYear > year )
		{
			year = newYear;
			std::cout << year << std::endl;
		}
	}
}
