#include "poly_utils.h"
#include <assert.h>
#include <stdlib.h>
#include <string.h>

Poly create_poly(void) {
  Poly p;
  for (int i = 0; i < MAX_POLY_DEGREE; i++) {
    p.coeffs[i] = 0.0;
  }
  return p;
}

inline double positive_fmod(double x, double m) {
  assert(m > 0.0);
  double r = fmod(x, m);
  if (r < 0.0)
    r += m;
  return r;
}

int64_t poly_degree(const Poly * __restrict p) {
  for (int64_t i = MAX_POLY_DEGREE - 1; i >= 0; i--) {
    if (fabs(p->coeffs[i]) > 1e-9) {
      return i;
    }
  }
  return 0;
}

double get_coeff(Poly p, int64_t degree) {
  if (degree >= MAX_POLY_DEGREE || degree < 0) {
    return 0.0;
  }
  return p.coeffs[degree];
}

void set_coeff(Poly *p, int64_t degree, double value) {
  if (degree >= MAX_POLY_DEGREE || degree < 0) {
    return;
  }
  p->coeffs[degree] = value;
}

Poly coeff_mod(const Poly * __restrict p, double modulus) {
  Poly out = create_poly();
  const int deg = poly_degree(p);
  for (int i = 0; i <= deg; i++) {
    double c = p->coeffs[i];
    if(fabs(c) < 1e-9) {
      continue;
    }
    double rounded = nearbyint(c);
    out.coeffs[i] = positive_fmod(c, modulus);
  }
  return out;
}

Poly poly_add(Poly a, Poly b) {
  Poly sum = create_poly();
  for (int i = 0; i < MAX_POLY_DEGREE; i++) {
    sum.coeffs[i] = a.coeffs[i] + b.coeffs[i];
  }
  return sum;
}

Poly poly_mul_scalar(Poly p, double scalar) {
  Poly res = create_poly();
  for (int i = 0; i < MAX_POLY_DEGREE; i++) {
    res.coeffs[i] = p.coeffs[i] * scalar;
  }
  return res;
}

Poly poly_mul(const Poly * __restrict a, const Poly * __restrict b) {
  Poly res = create_poly();
  const int deg_a = poly_degree(a);
  const int deg_b = poly_degree(b);
  assert(deg_a + deg_b < MAX_POLY_DEGREE);

  for (int i = 0; i <= deg_a; i++) {
    const double a_coeff = a->coeffs[i];
    if (fabs(a_coeff) < 1e-9) 
      continue;
    double * __restrict res_ptr = &res.coeffs[i];
    int j = 0;
    for (; j <= deg_b - 3; j += 4) {
      res_ptr[j]     += a_coeff * b->coeffs[j];
      res_ptr[j + 1] += a_coeff * b->coeffs[j + 1];
      res_ptr[j + 2] += a_coeff * b->coeffs[j + 2];
      res_ptr[j + 3] += a_coeff * b->coeffs[j + 3];
    }
    for (; j <= deg_b; j++) {
      res_ptr[j] += a_coeff * b->coeffs[j];
    }
  }
  return res;
}

void poly_divmod(const Poly * __restrict num, const Poly * __restrict den, Poly *quot, Poly *rem) {
  // In our case `den` should always be (x^n + 1)
  size_t ndeg = poly_degree(num);
  size_t ddeg = poly_degree(den);

  *quot = create_poly();
  *rem = *num;

  if (ndeg < ddeg) {
    return;
  }

  for (int64_t k = ndeg; k >= ddeg; --k) {
    double coeff = rem->coeffs[k];
    if (fabs(coeff) < 1e-9) {
      continue;
    }
    quot->coeffs[k - ddeg] += coeff;
    rem->coeffs[k] -= coeff;
    rem->coeffs[k - ddeg]  -= coeff;
  }
}

Poly poly_round_div_scalar(Poly x, double divisor) {
  Poly out = create_poly();
  assert(fabs(divisor) > 1e-9);

  for (int i = 0; i < MAX_POLY_DEGREE; i++) {
    double v = x.coeffs[i];
    out.coeffs[i] = round(v / divisor);
  }
  return out;
}
