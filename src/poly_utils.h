#ifndef POLY_UTILS_H
#define POLY_UTILS_H

#include "types.h"
#include <math.h>
#include <stdint.h>

double positive_fmod(double x, double m);

int64_t poly_degree(const Poly * __restrict p);

double get_coeff(Poly p, int64_t degree);

void set_coeff(Poly *p, int64_t degree, double value);

Poly coeff_mod(const Poly * __restrict p, double modulus);

Poly poly_add(Poly a, Poly b);

Poly poly_mul_scalar(Poly p, double scalar);

Poly poly_mul(const Poly * __restrict a, const Poly * __restrict b);

void poly_divmod(const Poly * __restrict numerator, const Poly * __restrict denominator, Poly *quotient,
                 Poly *remainder);

Poly create_poly(void);

#endif