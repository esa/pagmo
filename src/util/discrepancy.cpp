# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <fstream>
# include <cmath>
# include <ctime>
# include <cstring>

# include "discrepancy.h"

using namespace std;
namespace pagmo{ namespace util {namespace discrepancy {

base::~base() {}



/// Van Der Corput sequence
/**
 * Returns the n-th number in the Halton sequence
 *
 * @param[in] n selects which number of the Halton sequence to return
 * @param[in] base prime number to be used as a base of the sequence
 * @returns the n-th number in the van_der_corput sequence
 *
 * @see http://en.wikipedia.org/wiki/Van_der_Corput_sequence
**/
double van_der_corput(unsigned int n, unsigned int base) {
	double retval = 0;
	double f = 1.0 / base;
	unsigned int i = n;
	while (i > 0) {
		retval += f * (i % base);
		i = floor(i / base);
		f = f / base;
	}
	return retval;
}

# define PRIME_MAX 1601

/// Returns a prime number
/**
 * Returns any of the first 1600 prime numbers.

 * @param[in] n the index of the desired prime number (0 is 1)
 *
 * @throw value_error if n is larger than 1600 or smaller than 0
**/
unsigned int prime ( int n )
{
	unsigned int npvec[PRIME_MAX] = { 1,
		2,    3,    5,    7,   11,   13,   17,   19,   23,   29,
		31,   37,   41,   43,   47,   53,   59,   61,   67,   71,
		73,   79,   83,   89,   97,  101,  103,  107,  109,  113,
		127,  131,  137,  139,  149,  151,  157,  163,  167,  173,
		179,  181,  191,  193,  197,  199,  211,  223,  227,  229,
		233,  239,  241,  251,  257,  263,  269,  271,  277,  281,
		283,  293,  307,  311,  313,  317,  331,  337,  347,  349,
		353,  359,  367,  373,  379,  383,  389,  397,  401,  409,
		419,  421,  431,  433,  439,  443,  449,  457,  461,  463,
		467,  479,  487,  491,  499,  503,  509,  521,  523,  541,
		547,  557,  563,  569,  571,  577,  587,  593,  599,  601,
		607,  613,  617,  619,  631,  641,  643,  647,  653,  659,
		661,  673,  677,  683,  691,  701,  709,  719,  727,  733,
		739,  743,  751,  757,  761,  769,  773,  787,  797,  809,
		811,  821,  823,  827,  829,  839,  853,  857,  859,  863,
		877,  881,  883,  887,  907,  911,  919,  929,  937,  941,
		947,  953,  967,  971,  977,  983,  991,  997, 1009, 1013,
		1019, 1021, 1031, 1033, 1039, 1049, 1051, 1061, 1063, 1069,
		1087, 1091, 1093, 1097, 1103, 1109, 1117, 1123, 1129, 1151,
		1153, 1163, 1171, 1181, 1187, 1193, 1201, 1213, 1217, 1223,
		1229, 1231, 1237, 1249, 1259, 1277, 1279, 1283, 1289, 1291,
		1297, 1301, 1303, 1307, 1319, 1321, 1327, 1361, 1367, 1373,
		1381, 1399, 1409, 1423, 1427, 1429, 1433, 1439, 1447, 1451,
		1453, 1459, 1471, 1481, 1483, 1487, 1489, 1493, 1499, 1511,
		1523, 1531, 1543, 1549, 1553, 1559, 1567, 1571, 1579, 1583,
		1597, 1601, 1607, 1609, 1613, 1619, 1621, 1627, 1637, 1657,
		1663, 1667, 1669, 1693, 1697, 1699, 1709, 1721, 1723, 1733,
		1741, 1747, 1753, 1759, 1777, 1783, 1787, 1789, 1801, 1811,
		1823, 1831, 1847, 1861, 1867, 1871, 1873, 1877, 1879, 1889,
		1901, 1907, 1913, 1931, 1933, 1949, 1951, 1973, 1979, 1987,
		1993, 1997, 1999, 2003, 2011, 2017, 2027, 2029, 2039, 2053,
		2063, 2069, 2081, 2083, 2087, 2089, 2099, 2111, 2113, 2129,
		2131, 2137, 2141, 2143, 2153, 2161, 2179, 2203, 2207, 2213,
		2221, 2237, 2239, 2243, 2251, 2267, 2269, 2273, 2281, 2287,
		2293, 2297, 2309, 2311, 2333, 2339, 2341, 2347, 2351, 2357,
		2371, 2377, 2381, 2383, 2389, 2393, 2399, 2411, 2417, 2423,
		2437, 2441, 2447, 2459, 2467, 2473, 2477, 2503, 2521, 2531,
		2539, 2543, 2549, 2551, 2557, 2579, 2591, 2593, 2609, 2617,
		2621, 2633, 2647, 2657, 2659, 2663, 2671, 2677, 2683, 2687,
		2689, 2693, 2699, 2707, 2711, 2713, 2719, 2729, 2731, 2741,
		2749, 2753, 2767, 2777, 2789, 2791, 2797, 2801, 2803, 2819,
		2833, 2837, 2843, 2851, 2857, 2861, 2879, 2887, 2897, 2903,
		2909, 2917, 2927, 2939, 2953, 2957, 2963, 2969, 2971, 2999,
		3001, 3011, 3019, 3023, 3037, 3041, 3049, 3061, 3067, 3079,
		3083, 3089, 3109, 3119, 3121, 3137, 3163, 3167, 3169, 3181,
		3187, 3191, 3203, 3209, 3217, 3221, 3229, 3251, 3253, 3257,
		3259, 3271, 3299, 3301, 3307, 3313, 3319, 3323, 3329, 3331,
		3343, 3347, 3359, 3361, 3371, 3373, 3389, 3391, 3407, 3413,
		3433, 3449, 3457, 3461, 3463, 3467, 3469, 3491, 3499, 3511,
		3517, 3527, 3529, 3533, 3539, 3541, 3547, 3557, 3559, 3571,
		3581, 3583, 3593, 3607, 3613, 3617, 3623, 3631, 3637, 3643,
		3659, 3671, 3673, 3677, 3691, 3697, 3701, 3709, 3719, 3727,
		3733, 3739, 3761, 3767, 3769, 3779, 3793, 3797, 3803, 3821,
		3823, 3833, 3847, 3851, 3853, 3863, 3877, 3881, 3889, 3907,
		3911, 3917, 3919, 3923, 3929, 3931, 3943, 3947, 3967, 3989,
		4001, 4003, 4007, 4013, 4019, 4021, 4027, 4049, 4051, 4057,
		4073, 4079, 4091, 4093, 4099, 4111, 4127, 4129, 4133, 4139,
		4153, 4157, 4159, 4177, 4201, 4211, 4217, 4219, 4229, 4231,
		4241, 4243, 4253, 4259, 4261, 4271, 4273, 4283, 4289, 4297,
		4327, 4337, 4339, 4349, 4357, 4363, 4373, 4391, 4397, 4409,
		4421, 4423, 4441, 4447, 4451, 4457, 4463, 4481, 4483, 4493,
		4507, 4513, 4517, 4519, 4523, 4547, 4549, 4561, 4567, 4583,
		4591, 4597, 4603, 4621, 4637, 4639, 4643, 4649, 4651, 4657,
		4663, 4673, 4679, 4691, 4703, 4721, 4723, 4729, 4733, 4751,
		4759, 4783, 4787, 4789, 4793, 4799, 4801, 4813, 4817, 4831,
		4861, 4871, 4877, 4889, 4903, 4909, 4919, 4931, 4933, 4937,
		4943, 4951, 4957, 4967, 4969, 4973, 4987, 4993, 4999, 5003,
		5009, 5011, 5021, 5023, 5039, 5051, 5059, 5077, 5081, 5087,
		5099, 5101, 5107, 5113, 5119, 5147, 5153, 5167, 5171, 5179,
		5189, 5197, 5209, 5227, 5231, 5233, 5237, 5261, 5273, 5279,
		5281, 5297, 5303, 5309, 5323, 5333, 5347, 5351, 5381, 5387,
		5393, 5399, 5407, 5413, 5417, 5419, 5431, 5437, 5441, 5443,
		5449, 5471, 5477, 5479, 5483, 5501, 5503, 5507, 5519, 5521,
		5527, 5531, 5557, 5563, 5569, 5573, 5581, 5591, 5623, 5639,
		5641, 5647, 5651, 5653, 5657, 5659, 5669, 5683, 5689, 5693,
		5701, 5711, 5717, 5737, 5741, 5743, 5749, 5779, 5783, 5791,
		5801, 5807, 5813, 5821, 5827, 5839, 5843, 5849, 5851, 5857,
		5861, 5867, 5869, 5879, 5881, 5897, 5903, 5923, 5927, 5939,
		5953, 5981, 5987, 6007, 6011, 6029, 6037, 6043, 6047, 6053,
		6067, 6073, 6079, 6089, 6091, 6101, 6113, 6121, 6131, 6133,
		6143, 6151, 6163, 6173, 6197, 6199, 6203, 6211, 6217, 6221,
		6229, 6247, 6257, 6263, 6269, 6271, 6277, 6287, 6299, 6301,
		6311, 6317, 6323, 6329, 6337, 6343, 6353, 6359, 6361, 6367,
		6373, 6379, 6389, 6397, 6421, 6427, 6449, 6451, 6469, 6473,
		6481, 6491, 6521, 6529, 6547, 6551, 6553, 6563, 6569, 6571,
		6577, 6581, 6599, 6607, 6619, 6637, 6653, 6659, 6661, 6673,
		6679, 6689, 6691, 6701, 6703, 6709, 6719, 6733, 6737, 6761,
		6763, 6779, 6781, 6791, 6793, 6803, 6823, 6827, 6829, 6833,
		6841, 6857, 6863, 6869, 6871, 6883, 6899, 6907, 6911, 6917,
		6947, 6949, 6959, 6961, 6967, 6971, 6977, 6983, 6991, 6997,
		7001, 7013, 7019, 7027, 7039, 7043, 7057, 7069, 7079, 7103,
		7109, 7121, 7127, 7129, 7151, 7159, 7177, 7187, 7193, 7207,
		7211, 7213, 7219, 7229, 7237, 7243, 7247, 7253, 7283, 7297,
		7307, 7309, 7321, 7331, 7333, 7349, 7351, 7369, 7393, 7411,
		7417, 7433, 7451, 7457, 7459, 7477, 7481, 7487, 7489, 7499,
		7507, 7517, 7523, 7529, 7537, 7541, 7547, 7549, 7559, 7561,
		7573, 7577, 7583, 7589, 7591, 7603, 7607, 7621, 7639, 7643,
		7649, 7669, 7673, 7681, 7687, 7691, 7699, 7703, 7717, 7723,
		7727, 7741, 7753, 7757, 7759, 7789, 7793, 7817, 7823, 7829,
		7841, 7853, 7867, 7873, 7877, 7879, 7883, 7901, 7907, 7919,
		7927, 7933, 7937, 7949, 7951, 7963, 7993, 8009, 8011, 8017,
		8039, 8053, 8059, 8069, 8081, 8087, 8089, 8093, 8101, 8111,
		8117, 8123, 8147, 8161, 8167, 8171, 8179, 8191, 8209, 8219,
		8221, 8231, 8233, 8237, 8243, 8263, 8269, 8273, 8287, 8291,
		8293, 8297, 8311, 8317, 8329, 8353, 8363, 8369, 8377, 8387,
		8389, 8419, 8423, 8429, 8431, 8443, 8447, 8461, 8467, 8501,
		8513, 8521, 8527, 8537, 8539, 8543, 8563, 8573, 8581, 8597,
		8599, 8609, 8623, 8627, 8629, 8641, 8647, 8663, 8669, 8677,
		8681, 8689, 8693, 8699, 8707, 8713, 8719, 8731, 8737, 8741,
		8747, 8753, 8761, 8779, 8783, 8803, 8807, 8819, 8821, 8831,
		8837, 8839, 8849, 8861, 8863, 8867, 8887, 8893, 8923, 8929,
		8933, 8941, 8951, 8963, 8969, 8971, 8999, 9001, 9007, 9011,
		9013, 9029, 9041, 9043, 9049, 9059, 9067, 9091, 9103, 9109,
		9127, 9133, 9137, 9151, 9157, 9161, 9173, 9181, 9187, 9199,
		9203, 9209, 9221, 9227, 9239, 9241, 9257, 9277, 9281, 9283,
		9293, 9311, 9319, 9323, 9337, 9341, 9343, 9349, 9371, 9377,
		9391, 9397, 9403, 9413, 9419, 9421, 9431, 9433, 9437, 9439,
		9461, 9463, 9467, 9473, 9479, 9491, 9497, 9511, 9521, 9533,
		9539, 9547, 9551, 9587, 9601, 9613, 9619, 9623, 9629, 9631,
		9643, 9649, 9661, 9677, 9679, 9689, 9697, 9719, 9721, 9733,
		9739, 9743, 9749, 9767, 9769, 9781, 9787, 9791, 9803, 9811,
		9817, 9829, 9833, 9839, 9851, 9857, 9859, 9871, 9883, 9887,
		9901, 9907, 9923, 9929, 9931, 9941, 9949, 9967, 9973,10007,
		10009,10037,10039,10061,10067,10069,10079,10091,10093,10099,
		10103,10111,10133,10139,10141,10151,10159,10163,10169,10177,
		10181,10193,10211,10223,10243,10247,10253,10259,10267,10271,
		10273,10289,10301,10303,10313,10321,10331,10333,10337,10343,
		10357,10369,10391,10399,10427,10429,10433,10453,10457,10459,
		10463,10477,10487,10499,10501,10513,10529,10531,10559,10567,
		10589,10597,10601,10607,10613,10627,10631,10639,10651,10657,
		10663,10667,10687,10691,10709,10711,10723,10729,10733,10739,
		10753,10771,10781,10789,10799,10831,10837,10847,10853,10859,
		10861,10867,10883,10889,10891,10903,10909,10937,10939,10949,
		10957,10973,10979,10987,10993,11003,11027,11047,11057,11059,
		11069,11071,11083,11087,11093,11113,11117,11119,11131,11149,
		11159,11161,11171,11173,11177,11197,11213,11239,11243,11251,
		11257,11261,11273,11279,11287,11299,11311,11317,11321,11329,
		11351,11353,11369,11383,11393,11399,11411,11423,11437,11443,
		11447,11467,11471,11483,11489,11491,11497,11503,11519,11527,
		11549,11551,11579,11587,11593,11597,11617,11621,11633,11657,
		11677,11681,11689,11699,11701,11717,11719,11731,11743,11777,
		11779,11783,11789,11801,11807,11813,11821,11827,11831,11833,
		11839,11863,11867,11887,11897,11903,11909,11923,11927,11933,
		11939,11941,11953,11959,11969,11971,11981,11987,12007,12011,
		12037,12041,12043,12049,12071,12073,12097,12101,12107,12109,
		12113,12119,12143,12149,12157,12161,12163,12197,12203,12211,
		12227,12239,12241,12251,12253,12263,12269,12277,12281,12289,
		12301,12323,12329,12343,12347,12373,12377,12379,12391,12401,
		12409,12413,12421,12433,12437,12451,12457,12473,12479,12487,
		12491,12497,12503,12511,12517,12527,12539,12541,12547,12553,
		12569,12577,12583,12589,12601,12611,12613,12619,12637,12641,
		12647,12653,12659,12671,12689,12697,12703,12713,12721,12739,
		12743,12757,12763,12781,12791,12799,12809,12821,12823,12829,
		12841,12853,12889,12893,12899,12907,12911,12917,12919,12923,
		12941,12953,12959,12967,12973,12979,12983,13001,13003,13007,
		13009,13033,13037,13043,13049,13063,13093,13099,13103,13109,
		13121,13127,13147,13151,13159,13163,13171,13177,13183,13187,
		13217,13219,13229,13241,13249,13259,13267,13291,13297,13309,
		13313,13327,13331,13337,13339,13367,13381,13397,13399,13411,
		13417,13421,13441,13451,13457,13463,13469,13477,13487,13499 };
	if ( n == -1) {
		return npvec[PRIME_MAX-1];
	} else if ( n < PRIME_MAX ) {
		return npvec[n];
	}
	else
	{
		pagmo_throw(value_error,"n must be in [-1,1600]");
	}
}

/// Returns the smallest prime greater than or equal to n
/**
 * Returns the smallest prime greater than or equal to n

 * @param[in] n the number to be bounded.
 * @returns the smallest prime greater than or equal to n.
 *
 * @throw value_error if n is larger than 13499
 * @author dario.izzo@gmail.com
 *
**/
unsigned int prime_ge ( unsigned int n )
{
	unsigned int lb = 0;
	unsigned int ub = PRIME_MAX-1;
	unsigned int min_p = prime(lb);
	unsigned int max_p = prime(ub);
	if (n<min_p || n>max_p) {
		pagmo_throw(value_error,"n is out of range");
	}
	if (n == min_p) return min_p;
	if (n == max_p) return max_p;
	while(!(ub-lb==1)) {
		unsigned int tmp = (ub+lb)/2;
		unsigned int new_p = prime(tmp);
		if (new_p >= n) {
			ub = tmp;
			max_p = new_p;
		} else {
			lb = tmp;
			min_p = new_p;
		}
	}
	return max_p;
}
# undef PRIME_MAX


//****************************************************************************80

int *faure::binomial_table ( int qs, int m, int n )

//****************************************************************************80
//
//  Purpose:
//
//    BINOMIAL_TABLE computes a table of bionomial coefficients MOD QS.
//
//  Discussion:
//
//    Thanks to Michael Baudin for pointing out an error in a previous
//    version of this function, 07 December 2009.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 December 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int QS, the base for the MOD operation.
//
//    Input, int M, N, the limits of the binomial table.
//
//    Output, int BINOMIAL_TABLE[(M+1)*(N+1)], the table of binomial
//    coefficients modulo QS.
//
{
  int *coef;
  int i;
  int j;

  coef = new int[(m+1)*(n+1)];

  for ( j = 0; j <= n; j++ )
  {
	for ( i = 0; i <= m; i++ )
	{
	  coef[i+j*(m+1)] = 0;
	}
  }

  coef[0] = 1;

  j = 0;
  for ( i = 1; i <= m; i++ )
  {
	coef[i+j*(m+1)] = 1;
  }

  for ( i = 1; i <= i4_min ( m, n ); i++ )
  {
	j = i;
	coef[i+j*(m+1)] = 1;
  }

  for( j = 1; j <= n; j++ )
  {
	for ( i = j + 1; i <= m; i++ )
	{
	  coef[i+j*(m+1)] = ( coef[i-1+j*(m+1)] + coef[i-1+(j-1)*(m+1)] ) % qs;
	}
  }

  return coef;
}
//****************************************************************************80

void faure::faure_orig ( unsigned int dim_num, unsigned int *seed, double quasi[] )

//****************************************************************************80
//
//  Purpose:
//
//    FAURE generates a new quasirandom Faure vector with each call.
//
//  Discussion:
//
//    This routine implements the Faure method for computing
//    quasirandom numbers.  It is a merging and adaptation of
//    Bennett Fox's routines INFAUR and GOFAUR from ACM TOMS Algorithm 647.
//
//    Michael Baudin suggested a change so that whenever HISUM is altered,
//    the binomial table is recomputed, which allows a user to go back or
//    forth in the sequence.  08 December 2009.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 December 2009
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Henri Faure,
//    Discrepance de suites associees a un systeme de numeration
//    (en dimension s),
//    Acta Arithmetica,
//    Volume 41, 1982, pages 337-351.
//
//    Bennett Fox,
//    Algorithm 647:
//    Implementation and Relative Efficiency of Quasirandom
//    Sequence Generators,
//    ACM Transactions on Mathematical Software,
//    Volume 12, Number 4, December 1986, pages 362-376.
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension, which should be
//    at least 2.
//
//    Input/output, int *SEED, the seed, which can be used to index
//    the values.  On first call, set the input value of SEED to be 0
//    or negative.  The routine will automatically initialize data,
//    and set SEED to a new value.  Thereafter, to compute successive
//    entries of the sequence, simply call again without changing
//    SEED.  On the first call, if SEED is negative, it will be set
//    to a positive value that "skips over" an early part of the sequence
//    (This is recommended for better results).
//
//    Output, double QUASI[DIM_NUM], the next quasirandom vector.
//
{
  int hisum;
  int ktemp;
  int ltemp;
  int mtemp;
  double r;
  int ztemp;
//
//  Initialization required or requested?
//
  if ( m_qs <= 0 || *seed <= 0 )
  {
	m_qs = prime_ge ( dim_num );

	if ( m_qs < 1 )
	{
	  cout << "\n";
	  cout << "FAURE - Fatal error!\n";
	  cout << "  PRIME_GE failed.\n";
	  exit ( 1 );
	}
	m_hisum_save = -1;
  }

  if ( *seed == 0 )
  {
	hisum = 0;
  }
  else
  {
	hisum = i4_log_i4 ( *seed, m_qs );
  }
//
//  Is it necessary to recompute the coefficient table?
//
  if ( m_hisum_save != hisum )
  {
	if ( m_coef != NULL )
	{
	  delete [] m_coef;
	}

	if ( m_ytemp != NULL )
	{
	  delete [] m_ytemp;
	}

	m_hisum_save = hisum;

	m_coef = binomial_table ( m_qs, hisum, hisum );

	m_ytemp = new int[hisum+1];
  }
//
//  Find QUASI(1) using the method of Faure.
//
//  SEED has a representation in base QS of the form:
//
//    Sum ( 0 <= J <= HISUM ) YTEMP(J) * QS**J
//
//  We now compute the YTEMP(J)'s.
//
  ktemp = i4_power ( m_qs, hisum + 1 );
  ltemp = *seed;

  for ( int i = hisum; 0 <= i; i-- )
  {
	ktemp = ktemp / m_qs;
	mtemp = ltemp % ktemp;
	m_ytemp[i] = ( ltemp - mtemp ) / ktemp;
	ltemp = mtemp;
  }
//
//  QUASI(K) has the form
//
//    Sum ( 0 <= J <= HISUM ) YTEMP(J) / QS**(J+1)
//
//  Compute QUASI(1) using nested multiplication.
//
  r = ( ( double ) m_ytemp[hisum] );
  for ( int i = hisum-1; 0 <= i; i-- )
  {
	r = ( ( double ) m_ytemp[i] ) + r / ( ( double ) m_qs );
  }

  quasi[0] = r / ( ( double ) m_qs );
//
//  Find components QUASI(2:DIM_NUM) using the Faure method.
//
  for ( unsigned int k = 1; k < dim_num; k++ )
  {
	quasi[k] = 0.0;
	r = 1.0 / ( ( double ) m_qs );

	for ( int j = 0; j <= hisum; j++ )
	{
	  ztemp = 0;
	  for ( int i = j; i <= hisum; i++ )
	  {
		ztemp = ztemp + m_ytemp[i] * m_coef[i+j*(hisum+1)];
	  }
//
//  New YTEMP(J) is:
//
//    Sum ( J <= I <= HISUM ) ( old ytemp(i) * binom(i,j) ) mod QS.
//
	  m_ytemp[j] = ztemp % m_qs;
	  quasi[k] = quasi[k] + ( ( double ) m_ytemp[j] ) * r;
	  r = r / ( ( double ) m_qs );
	}
  }
//
//  Update SEED.
//
  *seed = *seed + 1;

  return;
}
//****************************************************************************80

int faure::i4_log_i4 ( int i4, int j4 )

//****************************************************************************80
//
//  Purpose:
//
//    I4_LOG_I4 returns the logarithm of an I4 to an I4 base.
//
//  Discussion:
//
//    Only the integer part of the logarithm is returned.
//
//    If
//
//      K4 = I4_LOG_J4 ( I4, J4 ),
//
//    then we ordinarily have
//
//      J4^(K4-1) < I4 <= J4^K4.
//
//    The base J4 should be positive, and at least 2.  If J4 is negative,
//    a computation is made using the absolute value of J4.  If J4 is
//    -1, 0, or 1, the logarithm is returned as 0.
//
//    The number I4 should be positive and at least 2.  If I4 is negative,
//    a computation is made using the absolute value of I4.  If I4 is
//    -1, 0, or 1, then the logarithm is returned as 0.
//
//    An I4 is an integer ( kind = 4 ) value.
//
//  Example:
//
//    I4  J4  K4
//
//     0   3   0
//     1   3   0
//     2   3   0
//     3   3   1
//     4   3   1
//     8   3   1
//     9   3   2
//    10   3   2
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    09 June 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I4, the number whose logarithm is desired.
//
//    Input, int J4, the base of the logarithms.
//
//    Output, int I4_LOG_I4, the integer part of the logarithm
//    base abs(J4) of abs(I4).
//
{
  int i4_abs;
  int j4_abs;
  int value;

  value = 0;

  i4_abs = abs ( i4 );

  if ( 2 <= i4_abs )
  {
	j4_abs = abs ( j4 );

	if ( 2 <= j4_abs )
	{
	  while ( j4_abs <= i4_abs )
	  {
		i4_abs = i4_abs / j4_abs;
		value = value + 1;
	  }
	}
  }
  return value;
}
//****************************************************************************80

int faure::i4_min ( int i1, int i2 )

//****************************************************************************80
//
//  Purpose:
//
//    I4_MIN returns the minimum of two I4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 October 1998
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I1, I2, two integers to be compared.
//
//    Output, int I4_MIN, the smaller of I1 and I2.
//
{
  int value;

  if ( i1 < i2 )
  {
	value = i1;
  }
  else
  {
	value = i2;
  }
  return value;
}
//****************************************************************************80

int faure::i4_power ( int i, int j )

//****************************************************************************80
//
//  Purpose:
//
//    I4_POWER returns the value of I^J.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 April 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I, J, the base and the power.  J should be nonnegative.
//
//    Output, int I4_POWER, the value of I^J.
//
{
  int k;
  int value;

  if ( j < 0 )
  {
	if ( i == 1 )
	{
	  value = 1;
	}
	else if ( i == 0 )
	{
	  cout << "\n";
	  cout << "I4_POWER - Fatal error!\n";
	  cout << "  I^J requested, with I = 0 and J negative.\n";
	  exit ( 1 );
	}
	else
	{
	  value = 0;
	}
  }
  else if ( j == 0 )
  {
	if ( i == 0 )
	{
	  cout << "\n";
	  cout << "I4_POWER - Fatal error!\n";
	  cout << "  I^J requested, with I = 0 and J = 0.\n";
	  exit ( 1 );
	}
	else
	{
	  value = 1;
	}
  }
  else if ( j == 1 )
  {
	value = i;
  }
  else
  {
	value = 1;
	for ( k = 1; k <= j; k++ )
	{
	  value = value * i;
	}
  }
  return value;
}

/// Constructor
/**
 * @param[in] dim dimension of the hypercube
 * @param[in] count starting point of the sequence (the first point is  [0.5,0.33333, ....])
 *
 * @throws value_error if dim not in [1,10]
*/
halton::halton(unsigned int dim, unsigned int count) : base(dim,count), m_primes() {
	if (dim >10 || dim==0) {
		pagmo_throw(value_error,"Halton sequences should not be used in dimension >10");
	}
	if (count == 0) {
		pagmo_throw(value_error,"The first element of the sequence has index 1");
	}
	for (size_t i=1; i<=dim; ++i) {
		m_primes.push_back(prime(i));
	}
}

/// Clone method.
base_ptr halton::clone() const
{
	return base_ptr(new halton(*this));
}
/// Operator ()
/**
 * Returns the next point in the sequence
 *
 * @return an std::vector<double> containing the next point
 */
std::vector<double> halton::operator()() {
	std::vector<double> retval;
	for (size_t i=0; i<m_dim; ++i) {
		retval.push_back(van_der_corput(m_count,m_primes.at(i)));
	}
	m_count++;
	return retval;
}
/// Operator (size_t n)
/**
 * Returns the n-th point in the sequence
 *
 * @param[in] n the point along the sequence to be returned
 * @return an std::vector<double> containing the n-th point
 */
std::vector<double> halton::operator()(unsigned int n) {
	if (n == 0) {
		pagmo_throw(value_error,"Halton sequence first point id is 1");
	}
	std::vector<double> retval;
	for (size_t i=0; i<m_primes.size(); ++i) {
		retval.push_back(van_der_corput(n,m_primes.at(i)));
	}
	m_count = n+1;
	return retval;
}



/// Constructor
/**
 * @param[in] dim dimension of the hypercube
 * @param[in] count starting point of the sequence
 *
 * @throws value_error if dim not in [2,23]
*/
faure::faure(unsigned int dim, unsigned int count) : base(dim, count), m_coef(NULL), m_hisum_save(-1), m_qs(-1), m_ytemp(NULL) {
		if (dim >23 || dim <2) {
			pagmo_throw(value_error,"Faure sequences can have dimension [2,23]");
		}
	}
/// Clone method.
base_ptr faure::clone() const
{
	return base_ptr(new faure(*this));
}
/// Operator ()
/**
 * Returns the next point in the sequence
 *
 * @return an std::vector<double> containing the next point
 */
std::vector<double> faure::operator()() {
	std::vector<double> retval(m_dim,0.0);
	faure_orig(m_dim, &m_count, &retval[0]);
	return retval;
}
/// Operator (unsigned int n)
/**
 * Returns the n-th point in the sequence
 *
 * @param[in] n the point along the sequence to be returned
 * @return an std::vector<double> containing the n-th point
 */
std::vector<double> faure::operator()(unsigned int n) {
	m_count = n;
	std::vector<double> retval(m_dim,0.0);
	faure_orig(m_dim, &m_count, &retval[0]);
	return retval;
}

/// Constructor
/**
 * @param[in] dim dimension of the hypercube
 * @param[in] count starting point of the sequence
 *
 * @throws value_error if dim not in [1..10] or count ==0
*/
simplex::simplex(unsigned int dim, unsigned int count) : base(dim,count), m_generator(dim-1), m_projector(dim) {
	if (dim >10 || dim==0) {
		pagmo_throw(value_error,"Halton sequences should not be used in dimension >10");
	}
	if (count == 0) {
		pagmo_throw(value_error,"The first element of the sequence has index 1");
	}
}

/// Clone method.
base_ptr simplex::clone() const
{
	return base_ptr(new simplex(*this));
}
/// Operator ()
/**
 * Returns the next point in the sequence
 *
 * @return an std::vector<double> containing the next point
 */
std::vector<double> simplex::operator()() {
	std::vector<double> tmp = m_generator();
//std::cout << tmp[0] << " " << tmp[1] << " " << std::endl;
	std::vector<double> retval = m_projector(tmp);
	return retval;
}
/// Operator (unsigned int n)
/**
 * Returns the n-th point in the sequence
 *
 * @param[in] n the point along the sequence to be returned
 * @return an std::vector<double> containing the n-th point
 */
std::vector<double> simplex::operator()(unsigned int n) {
	std::vector<double> retval = m_projector(m_generator(n));
	return retval;
}

}}} //namespaces
