#include <iostream>
#include <xmmintrin.h>
#include <immintrin.h>
#include <ctime>
#include <iomanip>
#include <chrono>
#pragma warning(disable : 4996)

#define COUNT 10

struct double4
{
	union {
		__m256d data4;
		double data[4];
	};
};

double4 _radii[COUNT];
double4 _mass[COUNT]; // mass of the planet + moons times the gravitational constant
double4 _x[COUNT];
double4 _v[COUNT];
const char* names[COUNT];

time_t currentTime;

void SetData()
{
	for ( int i = 0; i < COUNT; i++ )
	{
		_x[i].data[3] = 0;
		_v[i].data[3] = 0;
	}

	double G = 6.67408E-11; // from https://ssd.jpl.nasa.gov/?planet_phys_par

	// masses and radii from https://ssd.jpl.nasa.gov/?planet_phys_par
	// approximate cumulative mass of moons is added to the planets' masses
	// positions and velocities from https://ssd.jpl.nasa.gov/horizons.cgi

	// sun
	names[0] = "Sun";
	_mass[0].data4 = _mm256_set1_pd(G * 1.98847E30);
	_radii[0].data4 = _mm256_set1_pd(696340E3);
	_x[0].data[0] = -1.185107277176128E+09;
	_x[0].data[1] = 6.868809980765145E+08;
	_x[0].data[2] = 2.207810241420846E+07;

	_v[0].data[0] = -8.889457110342993E+00;
	_v[0].data[1] = -1.304250923538896E+01;
	_v[0].data[2] = 3.163853782834460E-01;

	// mercury
	names[1] = "Mercury";
	_mass[1].data4 = _mm256_set1_pd(G * 0.330114E24);
	_radii[1].data4 = _mm256_set1_pd(2440E3);
	_x[1].data[0] = -2.528619763481253E+10;
	_x[1].data[1] = 4.204832238691485E+10;
	_x[1].data[2] = 5.612793694435976E+09;

	_v[1].data[0] = -5.190738464527591E+04;
	_v[1].data[1] = -2.267367123092821E+04;
	_v[1].data[2] = 2.909180247640673E+03;

	// venus
	names[2] = "Venus";
	_mass[2].data4 = _mm256_set1_pd(G * 4.86747E24);
	_radii[2].data4 = _mm256_set1_pd(6052E3);
	_x[2].data[0] = -9.468973844425946E+10;
	_x[2].data[1] = -5.326287411690450E+10;
	_x[2].data[2] = 4.677363463577840E+09;

	_v[2].data[0] = 1.724869176950871E+04;
	_v[2].data[1] = -3.051159682535456E+04;
	_v[2].data[2] = -1.414101576603557E+03;

	// earth
	names[3] = "Earth";
	_mass[3].data4 = _mm256_set1_pd(G * 5.97237E24);
	_radii[3].data4 = _mm256_set1_pd(6371E3);
	_x[3].data[0] = 9.174344465476196E+10;
	_x[3].data[1] = -1.194191501816959E+11;
	_x[3].data[2] = 2.773489751581848E+07;

	_v[3].data[0] = 2.307173316736675E+04;
	_v[3].data[1] = 1.809496693541192E+04;
	_v[3].data[2] = -1.444899977126823E+00;

	// moon
	names[4] = "Moon";
	_mass[4].data4 = _mm256_set1_pd(G * 7.3476731E22);
	_radii[4].data4 = _mm256_set1_pd(1737E3);
	_x[4].data[0] = 9.208442853733037E+10;
	_x[4].data[1] = -1.192091206674976E+11;
	_x[4].data[2] = 6.181517798915505E+06;

	_v[4].data[0] = 2.259347030425684E+04;
	_v[4].data[1] = 1.894181681459600E+04;
	_v[4].data[2] = 6.600945431891603E+01;

	// mars
	names[5] = "Mars";
	_mass[5].data4 = _mm256_set1_pd(G * 0.641712E24);
	_radii[5].data4 = _mm256_set1_pd(3390E3);
	_x[5].data[0] = -2.404779300094771E+11;
	_x[5].data[1] = 6.908361964286420E+10;
	_x[5].data[2] = 7.325276039005000E+09;

	_v[5].data[0] = -5.764148340222327E+03;
	_v[5].data[1] = -2.124067255936238E+04;
	_v[5].data[2] = -3.033728502941466E+02;

	// jupiter
	names[6] = "Jupiter";
	_mass[6].data4 = _mm256_set1_pd(G * (1898.187E24 + 3.9310888E23));
	_radii[6].data4 = _mm256_set1_pd(69911E3);
	_x[6].data[0] = 6.162893437854434E+11;
	_x[6].data[1] = -4.292871209667132E+11;
	_x[6].data[2] = -1.200698000106323E+10;

	_v[6].data[0] = 7.309742699059670E+03;
	_v[6].data[1] = 1.133550722172970E+04;
	_v[6].data[2] = -2.105543749874936E+02;

	// saturn
	names[7] = "Saturn";
	_mass[7].data4 = _mm256_set1_pd(G * (568.3174E24 + 1.4054061E23));
	_radii[7].data4 = _mm256_set1_pd(58232E3);
	_x[7].data[0] = 9.518842035357896E+11;
	_x[7].data[1] = -1.142629610224875E+12;
	_x[7].data[2] = -1.803016850966245E+10;

	_v[7].data[0] = 6.880523794997542E+03;
	_v[7].data[1] = 6.162006013694501E+03;
	_v[7].data[2] = -3.814032289550697E+02;

	// uranus
	names[8] = "Uranus";
	_mass[8].data4 = _mm256_set1_pd(G * (86.8127E24 + 8.8901806E21));
	_radii[8].data4 = _mm256_set1_pd(25362E3);
	_x[8].data[0] = 2.214195292062304E+12;
	_x[8].data[1] = 1.954074207224361E+12;
	_x[8].data[2] = -2.142782863786554E+10;

	_v[8].data[0] = -4.556211357425507E+03;
	_v[8].data[1] = 4.788736453348525E+03;
	_v[8].data[2] = 7.687111845396810E+01;

	// neptune
	names[9] = "Neptune";
	_mass[9].data4 = _mm256_set1_pd(G * (102.4126E24 + 2.1487923E22));
	_radii[9].data4 = _mm256_set1_pd(24622E3);
	_x[9].data[0] = 4.421803545573483E+12;
	_x[9].data[1] = -6.834354150245397E+11;
	_x[9].data[2] = -8.783097535269541E+10;

	_v[9].data[0] = 7.934208025468410E+02;
	_v[9].data[1] = 5.403427551804345E+03;
	_v[9].data[2] = -1.292465944069530E+02;
}

__m256d GetLength( __m256d a )
{
	double4 d4;
	d4.data4 = _mm256_mul_pd(a, a);

	return _mm256_set1_pd(sqrt(d4.data[0] + d4.data[1] + d4.data[2]));
}

float year = 1900;
bool occultation[COUNT * COUNT];
int main()
{
	SetData();

	currentTime = time_t(1627682400);

	float timestep = 60;
	__m256d timestep4 = _mm256_set1_pd( timestep );

	std::cout << "==========================================" << std::endl;
	std::cout << "  Welcome to the Occultation Calculator.  " << std::endl;
	std::cout << "  Each year searched will be printed, as  " << std::endl;
	std::cout << "   well as each occultation between two   " << std::endl;
	std::cout << "                 planets.                 " << std::endl;
	std::cout << "==========================================" << std::endl;
	std::cout << std::endl;

	for ( int step = 0; step < 1E9f; step++ )
	{
		currentTime += time_t( timestep );

		//
		// Gravity simulation
		//

		// update velocities
		for ( int i = 0; i < COUNT; i++ )
		{
			// loop over planets, calculate gravitational acceleration, update velocities
			for ( int j = i + 1; j < COUNT; j++ )
			{
				//if (i == j) continue;
				__m256d distanceVector = _mm256_sub_pd( _x[j].data4, _x[i].data4 );
				__m256d distance = GetLength(distanceVector);
				__m256d distanceCubed = _mm256_mul_pd(distance, _mm256_mul_pd(distance, distance));
				__m256d d = _mm256_div_pd(distanceVector, distanceCubed);

				//_a[i].data4 = _mm256_mul_pd(_mm256_mul_pd(d, _mass[j].data4), timestep4);
				//_v[i].data4 = _mm256_add_pd( _v[i].data4, _mm256_mul_pd(_a[i].data4, timestep4 ) );
				_v[i].data4 = _mm256_add_pd( _v[i].data4, _mm256_mul_pd(_mm256_mul_pd(d, _mass[j].data4), timestep4));
				_v[j].data4 = _mm256_sub_pd( _v[j].data4, _mm256_mul_pd(_mm256_mul_pd(d, _mass[i].data4), timestep4 ));
			}
		}

		// update positions
		for ( int i = 0; i < COUNT; i++ )
		{
			_x[i].data4 = _mm256_add_pd( _x[i].data4, _mm256_mul_pd( _v[i].data4, timestep4 ) );
			_x[i].data[3] = 0;
		}

		//
		// check if there is an occultation by constructing a vector from the earth to each planet, then checking the angles between them
		//

		__m256d vectors[COUNT];
		double angularRadii[COUNT];

		// construct vectors and angular radii
		for ( int i = 0; i < COUNT; i++ )
		{
			if ( i == 3 ) continue;
			vectors[i] = _mm256_sub_pd( _x[i].data4, _x[3].data4 );

			// normalize vector
			double4 distance;
			distance.data4 = GetLength( vectors[i] );
			vectors[i] = _mm256_div_pd( vectors[i], distance.data4 );

			// calculate angular radius in radians
			angularRadii[i] = 2 * asin(_radii[i].data[0] / distance.data[0]);
		}

		// get angles
		for ( int i = 0; i < COUNT; i++ )
		{
			if ( i == 3 || i == 0 || i == 4 ) continue; // no occulation with earth, sun or moon
			for ( int j = i + 1; j < COUNT; j++ )
			{
				if ( j == 3 || j == 0 || j == 4) continue; // no occultation with earth, sun or moon
				int n = i * COUNT + j;

				// calculate angle between the vectors to the two planets using a dot product
				double4 dot4;
				dot4.data4 = _mm256_mul_pd(vectors[i], vectors[j]);
				double4 length;
				length.data4 = vectors[j];

				double angle = acos(dot4.data[0] + dot4.data[1] + dot4.data[2]);

				if ( angle < angularRadii[i] + angularRadii[j] )
				{
					if ( !occultation[n] )
					{
						// occultation starts
						std::cout << std::endl << std::endl;
						std::cout << "  > Occultation between " << names[i] << " and " << names[j]
							<< " at " << std::put_time( std::localtime( &currentTime ), "%F" )
							<< " from " << std::put_time(std::localtime(&currentTime), "%R");
						occultation[n] = true;
					}
				}
				else
				{
					if (occultation[n]) {
						std::cout << " until " << std::put_time(std::localtime(&currentTime), "%R");
						std::cout << std::endl << std::endl;
						occultation[n] = false;
					}
				}
			}
		}

		float newYear = gmtime( &currentTime )->tm_year + 1900;
		if ( newYear > year )
		{
			year = newYear;
			std::cout << year << "...";
			if (int(year) % 10 == 0) std::cout << std::endl;
		}
	}
}