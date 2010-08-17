/*
 * orbital.c
 * Created by Serhii Tsyba (sertsy@gmail.com) on 29.07.10.
 */

#include "orbital.h"
#include <math.h>

/*
 * The following definitions and a pair of eight arrays define coefficients of ecliptic variables given in
 * semi-analytic planetary theory VSOP82 by P. Bretagnon which are used to compute the following mean orbital planetary
 * elements:
 *      semi-major axis of planetary orbit (a)
 *      mean ecliptic longitude of the planet (λ)
 *      eccentricity of the planetary orbit (e)
 *      longitude of the perhelion of the planetary orbit (ϖ)
 *      inclination of the planetary orbit from the ecliptic (i)
 *      longitude of the ascending node (Ω)
 *
 * Each orbital element can be computed as a sum of a serie using the following expression:
 *
 *                                      p = Σ pᵢtⁱ, i = 1..n
 *
 * where p is the orbital elemet to be computed, pᵢ are coefficients provided by VSOP82 and t is the time instant
 * instant expressed in Julian millenia since the beginning of the epoch J2000 to a date at which orbital coefficient
 * is to be computed. n the maximum power of t that coefficients are given for.
 *
 * Planetary theory VSOP82 provides coefficients for orbital elements a and λ directly. Elements e, ϖ, i and Ω can be
 * found solving the following parametric expressions:
 *
 *                                      k = ecosϖ, h = esinϖ
 *                                      q = sin(i/2)cosΩ, p = sin(i/2)sinΩ
 *
 * where coefficients for parameters k, h, q and p are given.
 *
 * Mean orbital elemetns a and e are not effected by precession of Earth's rotational axis. However, the rest of the
 * elements is, so planetary theory VSOP82 provides coefficients for these elements referred to the standard equinox of
 * J2000 for the first group and mean dynamical ecliptic and equinox of date for the second one. Last part of the
 * variable name refers to the reference equinox of coefficients.
 *
 * Each array holds the coefficients for orbital elements of the major planet it names, which are given in the
 * following order: a, λ, k, h, q, p.
 *
 * Source: P. Bretagnon. Théorie du mouvement de l'ensemble des planetès. Solution VSOP82. Astronomy and Astrophisics,
 *         vol. 114, 1982, pp. 278-287.
 */

#define TOTAL_ORBITAL_COEFFICIENTS_OF_J2000 5   // total number of coefficients used to compute mean orbital elements
                                                // referred to the standard equinox of J2000
#define TOTAL_ORBITAL_COEFFICIENTS_OF_DATE 8    // total number of coefficients used to compute mean orbital elements
                                                // referred to the equinox of date

// coefficients used to compute mean orbital arguments referred to the standard equinox of J2000
static double mercury_orbital_coefficients_of_J2000[TOTAL_ORBITAL_ELEMENTS * TOTAL_ORBITAL_COEFFICIENTS_OF_J2000] = {
    0.38709830982, 0.0, 0.0, 0.0, 0.0,
    4.40260884240, 26087.90314157420, -0.00000934290, 0.00000003100, 0.0,
    0.04466059760, -0.00552114624, -0.00001860397, 0.00000063362, 0.0,
    0.20072331368, 0.00143750118, -0.00007974689, -0.00000026309, 0.0,
    0.04061563384, 0.00065433117, -0.00001071215, 0.00000021149, 0.0,
    0.04563550461, -0.00127633657, -0.00000913350, 0.00000018004, 0.0
};
static double venus_orbital_coefficients_of_J2000[TOTAL_ORBITAL_ELEMENTS * TOTAL_ORBITAL_COEFFICIENTS_OF_J2000] = {
    0.72332981996, 0.0, 0.0, 0.0, 0.0,
    3.17614669689, 10213.28554621100, 0.00000287555, -0.00000003038, 0.0,
    -0.00449282133, 0.00031259019, 0.00000605913, -0.00000069239, 0.0,
    0.00506684726, -0.00036121239, 0.00001839627, -0.00000000971, 0.0,
    0.00682410142, 0.00138133826, -0.00001090942, -0.00000185920, 0.0,
    0.02882285775, -0.00040384791, -0.00006232891, 0.00000025137, 0.0
};
static double earth_orbital_coefficients_of_J2000[TOTAL_ORBITAL_ELEMENTS * TOTAL_ORBITAL_COEFFICIENTS_OF_J2000] = {
    1.00000101778, 0.0, 0.0, 0.0, 0.0,
    1.75347031435, 6283.07584918000, -0.00000991890, 0.00000000073, 0.0,
    -0.00374081650, -0.00082266699, 0.00002748939, 0.00000104217, 0.0,
    0.01628447663, -0.00062030259, -0.00003353888, 0.00000071185, 0.0,
    0.0, -0.00113469002, 0.00001237314, 0.00000127050, 0.0,
    0.0, 0.00010180391, 0.00004701998, -0.00000053829, 0.0
};
static double mars_orbital_coefficients_of_J2000[TOTAL_ORBITAL_ELEMENTS * TOTAL_ORBITAL_COEFFICIENTS_OF_J2000] = {
    1.52367934191, 0.00000000031, 0.0, 0.0, 0.0,
    6.20348091341, 3340.61243149230, 0.00000454761, -0.00000005057, 0.0,
    0.08536560252, 0.00376330152, -0.00024657416, -0.00000395241, 0.0,
    -0.03789973236, 0.00624657465, 0.00015527232, -0.00000671940, 0.0,
    0.01047042574, 0.00017138526, -0.00004077591, -0.00000138600, 0.0,
    0.01228449307, -0.00108020083, -0.00001922195, 0.00000088373, 0.0
};
static double jupiter_orbital_coefficients_of_J2000[TOTAL_ORBITAL_ELEMENTS * TOTAL_ORBITAL_COEFFICIENTS_OF_J2000] = {
    5.20260319132, 0.00000191323, 0.0, 0.0, 0.0,
    0.59954649739, 529.69096509460, -0.00014837133, 0.00000007482, 0.0,
    0.04698572124, 0.00113010377, -0.00010930126, -0.00000428748, 0.00000020539,
    0.01200385748, 0.00217149360, 0.00009858539, -0.00000513109, -0.00000009007,
    -0.00206561098, -0.00031340156, -0.00001667392, 0.00000076926, 0.0,
    0.01118377157, -0.00023427562, 0.00002086760, 0.00000050721, 0.0
};
static double saturn_orbital_coefficients_of_J2000[TOTAL_ORBITAL_ELEMENTS * TOTAL_ORBITAL_COEFFICIENTS_OF_J2000] = {
    9.55490959574, -0.00002138917, 0.0, 0.0, 0.0,
    0.87401675650, 213.299095438, 0.00036659741, -0.00000033330, 0.00000000217,
    -0.00296003595, -0.00529602626, 0.00030928405, 0.00001296215, -0.00000059959,
    0.05542964254, -0.00375593887, -0.00031990236, 0.00001598633, 0.00000032451,
    -0.00871747436, 0.00080171499, 0.00004142282, -0.00000196049, -0.00000009439,
    0.01989147301, 0.00059439766, -0.00005235717, -0.00000127219, 0.00000008295
};
static double uranus_orbital_coefficients_of_J2000[TOTAL_ORBITAL_ELEMENTS * TOTAL_ORBITAL_COEFFICIENTS_OF_J2000] = {
    19.21844606178, -0.00000037163, 0.00000009791, 0.0, 0.0,
    5.48129387159, 74.78159856730, -0.00000848828, 0.00000010450, 0.0,
    -0.04595132376, 0.00018344050, -0.00000080849, -0.00000045396, 0.00000002185,
    0.00563791307, -0.00074964350, 0.00001210200, -0.00000042088, -0.00000001714,
    0.00185915075, -0.00012449382, -0.00000207373, 0.00000007621, 0.0,
    0.00648617008, -0.00011744733, 0.00000317799, 0.00000007317, 0.0
};
static double neptune_orbital_coefficients_of_J2000[TOTAL_ORBITAL_ELEMENTS * TOTAL_ORBITAL_COEFFICIENTS_OF_J2000] = {
    30.11038686942, 0.00000166346, 0.00000006857, 0.0, 0.0,
    5.31188628676, 38.13303563780, 0.00000102311, -0.00000004340, 0.0,
    0.00599977571, 0.00000871279, -0.00000119902, -0.00000004034, 0.0,
    0.00669242413, 0.000007824336, 0.00000080801, -0.00000003955, 0.0,
    -0.01029147819, -0.00000072727, -0.00000006568, 0.00000001668, 0.0,
    0.01151683985, 0.00002575536, 0.00000019377, 0.00000001331, 0.0
};

// coefficients used to compute mean orbital arguments referred to equinox of date
static double mercury_orbital_coefficients_of_date[TOTAL_ORBITAL_ELEMENTS * TOTAL_ORBITAL_COEFFICIENTS_OF_DATE] = {
    0.38709830982, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    4.40260884240, 26088.1470711010, 0.0005305220, 0.0000003097, 0.0000000022, 0.0, 0.0, 0.0,
    0.04466059760, -0.0544834872, -0.0018063052, 0.0006631851, 0.0000146668, -0.0000023759, -0.0000000510, 0.0000000038,
    0.20072331368, 0.0123315378, -0.0073740873, 0.0001850021, 0.0000445323, 0.0000009447, -0.0000001047, -0.0000000022,
    0.04061563384, -0.0093423887, -0.0009194385, 0.0000651855, 0.0000036846, -0.0000001337, -0.0000000055, 0.0,
    0.04563550461, 0.0085271267, -0.0009554616, -0.0000671422, 0.0000033240, 0.0000001594, -0.0000000043, 0.0
};
static double venus_orbital_coefficients_of_date[TOTAL_ORBITAL_ELEMENTS * TOTAL_ORBITAL_COEFFICIENTS_OF_DATE] = {
    0.72332981996, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    3.17614669689, 10213.5294305233, 0.0005421057, 0.0000002537, 0.0, 0.0, 0.0, 0.0,
    -0.00449282133, -0.0009231346, 0.0002250365, -0.0000014408, -0.0000016821, 0.0000000623, 0.0000000051, 0.0,
    0.00506684726, -0.0014569408, -0.0000584775, 0.0000225733, -0.0000006030, -0.0000001007, 0.0000000042, 0.0,
    0.00682410142, -0.0045129495, -0.0001184303, 0.0000177664, 0.0000004721, -0.0000000036, 0.0000000002, 0.0,
    0.02882285775, 0.0011584562, -0.0003492012, -0.0000087791, 0.0000005940, 0.0000000171, 0.0000000013, 0.0
};
static double earth_orbital_coefficients_of_date[TOTAL_ORBITAL_ELEMENTS * TOTAL_ORBITAL_COEFFICIENTS_OF_DATE] = {
    1.00000101778, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    1.75347031435, 6283.3199666635, 0.0005300181, 0.0000003692, 0.0, 0.0, 0.0, 0.0,
    -0.00374081650, -0.0047931064, 0.0002811275, 0.0000738309, -0.0000026511, -0.0000003676, 0.0000000082, 0.0,
    0.01628447663, -0.0015323789, -0.0007201706, 0.0000322989, 0.0000057885, -0.0000001659, -0.0000000196, 0.0,
    0.0, -0.00113469002, 0.00001237314, 0.00000127050, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.00010180391, 0.00004701998, -0.00000053829, 0.0, 0.0, 0.0, 0.0
};
static double mars_orbital_coefficients_of_date[TOTAL_ORBITAL_ELEMENTS * TOTAL_ORBITAL_COEFFICIENTS_OF_DATE] = {
    1.52367934191, 0.00000000031, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    6.20348091341, 3340.8562789898, 0.0005427485, 0.0000002663, 0.0, 0.0, 0.0, 0.0,
    0.08536560252, 0.0130050525, -0.0042873758, -0.0002598371, 0.0000354359, 0.0000015619, -0.0000001107, -0.0000000041,
    -0.03789973236, 0.0270627603, 0.0022456771, -0.0004518251, -0.0000225665, 0.0000021906, 0.00000000898, -0.0000000045,
    0.01047042574, -0.0016894313, -0.0000828113, 0.0000036129, -0.0000000219, 0.0000000241, 0.0000000011, 0.0,
    0.01228449307, 0.0013710386, -0.0001073560, -0.0000026035, -0.0000000679, -0.0000000106, 0.0000000024, 0.0
};
static double jupiter_orbital_coefficients_of_date[TOTAL_ORBITAL_ELEMENTS * TOTAL_ORBITAL_COEFFICIENTS_OF_DATE] = {
    5.20260319132, 0.00000191323, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.59954649739, 529.9348075394, 0.0003904990, 0.0000004374, 0.0, 0.0, 0.0, 0.0,
    0.04698572124, -0.0017969488, -0.0020421373, -0.0000402623, 0.0000168535, 0.0000005800, -0.0000000608, 0.0000000022,
    0.01200385748, 0.0136286045, 0.0000426023, -0.0002108274, -0.0000061282, 0.0000011026, 0.0000000412, -0.0000000027,
    -0.00206561098, -0.0019057263, 0.0001082728, 0.0000089342, -0.0000004218, -0.0000000122, 0.0000000011, 0.0,
    0.01118377157, -0.008397307, -0.0001594888, 0.0000079163, 0.0000003776, -0.0000000205, -0.0000000003, 0.0
};
static double saturn_orbital_coefficients_of_date[TOTAL_ORBITAL_ELEMENTS * TOTAL_ORBITAL_COEFFICIENTS_OF_DATE] = {
    9.55490959574, -0.00002138917, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.87401675650, 213.5429562981, 0.0009067343, -0.0000000557, 0.0000000017, 0.0, 0.0, 0.0,
    -0.00296003595, -0.0188131383, 0.0012832847, 0.0003848105, -0.0000214460, -0.0000025225, 0.0000001112, 0.0000000076,
    0.05542964254, -0.0044777706, -0.0032611427, 0.0002000724, 0.0000346612, -0.0000017338, -0.0000001542, 0.0000000055,
    -0.00871747436, -0.00029141827, 0.0001573509, 0.0000123816, -0.0000007530, -0.0000000283, 0.0000000040, 0.0,
    0.01989147301, -0.0016330436, -0.0002233231, 0.0000111925, 0.0000005881, -0.0000000548, -0.0000000016, 0.0
};
static double uranus_orbital_coefficients_of_date[TOTAL_ORBITAL_ELEMENTS * TOTAL_ORBITAL_COEFFICIENTS_OF_DATE] = {
    19.21844606178, -0.00000037163, 0.00000009791, 0.0, 0.0, 0.0, 0.0, 0.0,
    5.48129387159, 75.0254311493, 0.0005311712, 0.0000004525, 0.0, 0.0, 0.0, 0.0,
    -0.04595132376, -0.0011912664, 0.0015449390, 0.0000112133, -0.0000083600, -0.0000000380, 0.0000000168, 0.0,
    0.00563791307, -0.0119540732, -0.0001355665, 0.0001320329, 0.0000007337, -0.0000004157, -0.0000000015, 0.0,
    0.00185915075, -0.0005713203, -0.0000197534, -0.0000049897, 0.0000000189, 0.0000000323, 0.0000000006, 0.0,
    0.00648617008, 0.0002340586, 0.0000106597, 0.0000011931, -0.0000004831, -0.0000000048, 0.0000000018, 0.0
};
static double neptune_orbital_coefficients_of_date[TOTAL_ORBITAL_ELEMENTS * TOTAL_ORBITAL_COEFFICIENTS_OF_DATE] = {
    30.11038686942, 0.00000166346, 0.00000006857, 0.0, 0.0, 0.0, 0.0, 0.0,
    5.31188628676, 38.3768771649, 0.0005397652, 0.0000003066, 0.0, 0.0, 0.0, 0.0,
    0.00599977571, -0.0016231781, -0.0002022529, 0.0000148426, 0.0000012224, -0.0000000341, -0.0000000031, 0.0,
    0.00669242413, 0.0015412378, -0.0001927965, -0.0000180283, 0.0000008224, 0.0000000668, -0.0000000010, 0.0,
    -0.01029147819, -0.0016743180, 0.0003058262, 0.0000056754, -0.0000014012, -0.0000000046, 0.0000000030, 0.0,
    0.01151683985, -0.0025854023, -0.0001182724, 0.0000237388, 0.0000002091, -0.0000000686, 0.0000000001, 0.0
};

/*
 * An enumeration provided for access to the elements of the arrays of coefficients for mean orbital planetary
 * arguments.
 */
enum {
    A = 0,                      // index of semi-major axis of planetary orbit (a)
    L = 1,                      // index of mean ecliptic longitude of the planet (λ)
    K,                          // index of parametric variable k
    H,                          // index of parametric variable h
    Q,                          // index of parametric variable q
    P = 5                       // index of parametric variable p
};

/*
 * Computes mean orbital arguments of a major planet according to planetary theory VSOP82 given time instant (t) at
 * which elements are to be computed expressed in Julian millenia since the beginning of the epoch J2000 and
 * coefficients (coefficients) for the ecliptic variables provided by planetary theory VSOP82.
 * Output is written into a given array (elements).
 * Result orbital elements are expressed in radians (except for the mean eccentricity of the orbit).
 *
 * Source: P. Bretagnon. Théorie du mouvement de l'ensemble des planetès. Solution VSOP82. Astronomy and Astrophisics,
 *         vol. 114, 1982, p. 285-286.
 */
static void compute_orbital_elements_series(double t, double coefficients[], int m, double elements[])
{
    double tn;                              // accumulative variable holding n-th power of t
    int i, j;                               // loop index variables
    double ve[TOTAL_ORBITAL_ELEMENTS];      // array holding intermediate computations of VSOP82 ecliptic variables
                                            // a, λ, k, h, q, p
    
    // computing ecliptic variables a, λ, k, h, q, p of VSOP82 theory
    for (i = 0; i < TOTAL_ORBITAL_ELEMENTS; i++){
        for (j = 0, ve[i] = 0.0, tn = 1.0; j < m; j++, tn *= t)
            ve[i] += coefficients[i * m + j] * tn;
    }
    
    // mean orbital elements λ and a are computed directly from given coefficients
    elements[PLANET_MEAN_LONGITUDE] = ve[L];
    elements[SEMI_MAJOR_AXIS] = ve[A];
    
    // the rest of the mean orbital elements is deduced from the expressions for k, h, q and p:
    //
    //                              k = ecosϖ, h = esinϖ
    //                              q = sin(i/2)cosΩ, p = sin(i/2)sinΩ
    //
    // solving above systems results in the following expressions for ϖ, e, Ω and i:
    //
    //                              ϖ = arctan(h/k), e = k/cosϖ (e = h/sinϖ)
    //                              Ω = arctanp/q, i = 2arcsin(q/cosΩ) (i = 2arcsin(p/sinΩ))
    //
    elements[PERHELION_LONGITUDE] = atan2(ve[H], ve[K]);
    elements[ECCENTRICITY_OF_ORBIT] = ve[K] / cos(elements[PERHELION_LONGITUDE]);
    elements[ASCENDING_NODE_LONGITUDE] = atan2(ve[P], ve[Q]);
    elements[ECLIPTIC_INCLINATION] = 2.0 * asin(ve[Q] / cos(elements[ASCENDING_NODE_LONGITUDE]));
}

void compute_orbital_elements_of_J2000(date d, Planet p, double elements[])
{
    double t;                   // time instant in measured in Julian millenia since the beginning of the epoch J2000
                                // till a given date
    
    // computing time instant in Julian centuries since the beginning of the epoch J2000
    t = (julian_ephemeris_date(d) - J2000) / DAYS_IN_JULIAN_MILLENIUM;
    
    // computing mean orbital elements for a given planet
    switch (p) {
        case MERCURY:
            compute_orbital_elements_series(t, mercury_orbital_coefficients_of_J2000, TOTAL_ORBITAL_COEFFICIENTS_OF_J2000, elements);
            break;
        case VENUS:
            compute_orbital_elements_series(t, venus_orbital_coefficients_of_J2000, TOTAL_ORBITAL_COEFFICIENTS_OF_J2000, elements);
            break;
        case EARTH:
            compute_orbital_elements_series(t, earth_orbital_coefficients_of_J2000, TOTAL_ORBITAL_COEFFICIENTS_OF_J2000, elements);
            break;
        case MARS:
            compute_orbital_elements_series(t, mars_orbital_coefficients_of_J2000, TOTAL_ORBITAL_COEFFICIENTS_OF_J2000, elements);
            break;
        case JUPITER:
            compute_orbital_elements_series(t, jupiter_orbital_coefficients_of_J2000, TOTAL_ORBITAL_COEFFICIENTS_OF_J2000, elements);
            break;
        case SATURN:
            compute_orbital_elements_series(t, saturn_orbital_coefficients_of_J2000, TOTAL_ORBITAL_COEFFICIENTS_OF_J2000, elements);
            break;
        case URANUS:
            compute_orbital_elements_series(t, uranus_orbital_coefficients_of_J2000, TOTAL_ORBITAL_COEFFICIENTS_OF_J2000, elements);
            break;
        case NEPTUNE:
            compute_orbital_elements_series(t, neptune_orbital_coefficients_of_J2000, TOTAL_ORBITAL_COEFFICIENTS_OF_J2000, elements);
            break;
        default:
            break;
    }
}

void compute_orbital_elements_of_date(date d, Planet p, double elements[])
{
    double t;                   // time instant in measured in Julian millenia since the beginning of the epoch J2000
                                // till a given date
    
    // computing time instant in Julian centuries since the beginning of the epoch J2000
    t = (julian_ephemeris_date(d) - J2000) / DAYS_IN_JULIAN_MILLENIUM;
    
    // computing mean orbital elements for a given planet
    switch (p) {
        case MERCURY:
            compute_orbital_elements_series(t, mercury_orbital_coefficients_of_date, TOTAL_ORBITAL_COEFFICIENTS_OF_DATE, elements);
            break;
        case VENUS:
            compute_orbital_elements_series(t, venus_orbital_coefficients_of_date, TOTAL_ORBITAL_COEFFICIENTS_OF_DATE, elements);
            break;
        case EARTH:
            compute_orbital_elements_series(t, earth_orbital_coefficients_of_date, TOTAL_ORBITAL_COEFFICIENTS_OF_DATE, elements);
            break;
        case MARS:
            compute_orbital_elements_series(t, mars_orbital_coefficients_of_date, TOTAL_ORBITAL_COEFFICIENTS_OF_DATE, elements);
            break;
        case JUPITER:
            compute_orbital_elements_series(t, jupiter_orbital_coefficients_of_date, TOTAL_ORBITAL_COEFFICIENTS_OF_DATE, elements);
            break;
        case SATURN:
            compute_orbital_elements_series(t, saturn_orbital_coefficients_of_date, TOTAL_ORBITAL_COEFFICIENTS_OF_DATE, elements);
            break;
        case URANUS:
            compute_orbital_elements_series(t, uranus_orbital_coefficients_of_date, TOTAL_ORBITAL_COEFFICIENTS_OF_DATE, elements);
            break;
        case NEPTUNE:
            compute_orbital_elements_series(t, neptune_orbital_coefficients_of_date, TOTAL_ORBITAL_COEFFICIENTS_OF_DATE, elements);
            break;
        default:
            break;
    }
}