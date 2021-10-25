#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>

#define COMPOSITE 0
#define PRIME 1

#ifdef USE_PARALLEL
#include <omp.h>
#else
#include <time.h>
#endif

typedef struct {
  mpz_t table;
  unsigned int size;		/* Bit set 1: composite */
} Sieve;



void initialize_sieve (Sieve* p_sieve)
{
  p_sieve->size = 2;
  mpz_init(p_sieve->table);
}



void destroy_sieve (Sieve* p_sieve)
{
  mpz_clear(p_sieve->table);
  p_sieve->size = 0;
}



int sieve_primality_test (unsigned int n, Sieve* p_sieve)
{
  if (n <= p_sieve->size) {
    return !mpz_tstbit(p_sieve->table, n);
  }
  unsigned int prev_size = p_sieve->size;
  p_sieve->size *= 2;
  unsigned int i;
  for (i = 2; i <= prev_size; i++) {
    if (mpz_tstbit(p_sieve->table, i) == 0) {
      unsigned int j;
      for (j = i * 2; j <= p_sieve->size; j += i) {
	mpz_setbit(p_sieve->table, j);
      }
    }
  }
  return sieve_primality_test(n, p_sieve);
}



typedef struct {
  mpz_t* coef;
  unsigned int deg;
} Polynomial;



void initialize_polynomial (Polynomial** pp_poly, unsigned int deg)
{
  (*pp_poly) = (Polynomial*)malloc(sizeof(Polynomial));
  (*pp_poly)->coef = (mpz_t*)malloc((deg + 1) * sizeof(mpz_t));
  (*pp_poly)->deg = deg;
  unsigned int i;
  for (i = 0; i <= deg; i++) {
    mpz_init_set_ui((*pp_poly)->coef[i], 0);
  }
}



void destroy_polynomial (Polynomial** pp_poly)
{
  unsigned int i;
  for (i = 0; i <= (*pp_poly)->deg; i++) {
    mpz_clear((*pp_poly)->coef[i]);
  }
  free((*pp_poly)->coef);
  (*pp_poly)->deg = 0;
  free(*pp_poly);
}



void clone_polynomial(Polynomial** pp_poly_cloned, Polynomial* p_poly)
{
  (*pp_poly_cloned) = (Polynomial*)malloc(sizeof(Polynomial));
  (*pp_poly_cloned)->coef = (mpz_t*)malloc((p_poly->deg + 1) * sizeof(mpz_t));
  (*pp_poly_cloned)->deg = p_poly->deg;
  unsigned int i;
  for (i = 0; i <= p_poly->deg; i++) {
    mpz_init_set((*pp_poly_cloned)->coef[i], p_poly->coef[i]);
  }
}



void compact_polynomial (Polynomial* p_poly)
{
  unsigned int i;
  for (i = p_poly->deg; i > 0; i--) {
    if (mpz_cmp_ui(p_poly->coef[i], 0) != 0) {
      break;
    }
  }
  if (i < p_poly->deg) {
    unsigned int j;
    for (j = i + 1; j <= p_poly->deg; j++) {
      mpz_clear(p_poly->coef[j]);
    }
    p_poly->coef = (mpz_t*)realloc(p_poly->coef, sizeof(mpz_t) * (i + 1));
    p_poly->deg = i;
  }
}



/* Return 1 if two polynomials are equal */
int is_equal_polynomial (Polynomial* p_poly0, Polynomial* p_poly1)
{
  if (p_poly0->deg != p_poly1->deg) {
    return 0;
  }
  unsigned int i;
  for (i = 0; i <= p_poly0->deg; i++) {
    if (mpz_cmp(p_poly0->coef[i], p_poly1->coef[i]) != 0) {
      return 0;
    }
  }  
  return 1;
}



static inline void get_polynomial_coef (mpz_t* p_coef, 
				 Polynomial* p_poly, unsigned int order)
{
  if (order > p_poly->deg) {
    mpz_init_set_ui(*p_coef, 0);
    return;
  }
  mpz_set(*p_coef, p_poly->coef[order]);
}



void set_polynomial_coef (Polynomial* p_poly, unsigned int order, 
			  const mpz_t* p_coef)
{
  if (order <= p_poly->deg) {
    mpz_set(p_poly->coef[order], *p_coef);
    return;
  }
  p_poly->coef = (mpz_t*)realloc(p_poly->coef, sizeof(mpz_t) * (order + 1));
  unsigned int i;
  for (i = p_poly->deg + 1; i < order; i++) {
    mpz_init_set_ui(p_poly->coef[i], 0);
  }
  mpz_init_set(p_poly->coef[order], *p_coef);
  p_poly->deg = order;
}



void set_polynomial_coef_si (Polynomial* p_poly, unsigned int order, 
			     int coef_si)
{
  if (order <= p_poly->deg) {
    mpz_set_si(p_poly->coef[order], coef_si);
    return;
  }
  p_poly->coef = (mpz_t*)realloc(p_poly->coef, sizeof(mpz_t) * (order + 1));
  unsigned int i;
  for (i = p_poly->deg + 1; i < order; i++) {
    mpz_init_set_ui(p_poly->coef[i], 0);
  }
  mpz_init_set_si(p_poly->coef[order], coef_si);
  p_poly->deg = order;
}



void polynomial_modular_multiplication (Polynomial** pp_poly_res, 
					Polynomial* p_poly0, Polynomial* p_poly1, 
					mpz_t n, unsigned int r)
{
  unsigned int max_output_deg = p_poly0->deg + p_poly1->deg < r? 
    p_poly0->deg + p_poly1->deg: r - 1;
  initialize_polynomial(pp_poly_res, max_output_deg);
  mpz_t coef;
  mpz_init(coef);
  unsigned int i;
  for (i = 0; i < r; i++) {	/* For every result coefficient */
    mpz_set_ui(coef, 0);
    mpz_t c0, c1;
    mpz_init(c0);
    mpz_init(c1);
    unsigned int jmin = i > p_poly1->deg? i - p_poly1->deg: 0;
    unsigned int jmax = i < p_poly0->deg? i: p_poly0->deg;
    unsigned int j;
    for (j = jmin; j <= jmax; j++) { /* For all a(j) * b(i - j) */
      get_polynomial_coef(&c0, p_poly0, j);
      get_polynomial_coef(&c1, p_poly1, i - j);
      mpz_mul(c0, c0, c1);
      mpz_add(coef, c0, coef);	/* coef += c0 * c1 */
    }
    jmin = i + r - p_poly1->deg;
    jmax = p_poly0->deg;
    for (j = jmin; j <= jmax; j++) { /* For all a(j) * b(i - j + r) */
      get_polynomial_coef(&c0, p_poly0, j);
      get_polynomial_coef(&c1, p_poly1, i + r - j);
      mpz_mul(c0, c0, c1);
      mpz_add(coef, c0, coef);	/* coef += c0 * c1 */
    }
    mpz_clear(c0);
    mpz_clear(c1);
    mpz_mod(coef, coef, n);	/* coef = coef % n */
    if (mpz_cmp_ui(coef, 0) != 0) {
      set_polynomial_coef((*pp_poly_res), i, &coef);
    }
  }
  compact_polynomial((*pp_poly_res));
  mpz_clear(coef);
}



/* Compute ((*p_poly_base) ^ n) % (X ^ r - 1) */
void polynomial_modular_power (Polynomial** pp_poly_res, Polynomial* p_poly_base, 
			       mpz_t n, unsigned int r)
{
  initialize_polynomial(pp_poly_res, 0);
  set_polynomial_coef_si((*pp_poly_res), 0, 1);
  Polynomial* p_poly_temp = NULL;
  unsigned int i;
  for (i = mpz_sizeinbase(n, 2) + 1; i> 0; ) {
    clone_polynomial(&p_poly_temp, (*pp_poly_res));
    destroy_polynomial(pp_poly_res);
    polynomial_modular_multiplication(pp_poly_res, p_poly_temp, p_poly_temp, n, r);
    destroy_polynomial(&p_poly_temp);
    if (mpz_tstbit(n, --i) == 1) {
      clone_polynomial(&p_poly_temp, (*pp_poly_res));
      destroy_polynomial(pp_poly_res);
      polynomial_modular_multiplication(pp_poly_res, p_poly_temp, p_poly_base, n, r);
      destroy_polynomial(&p_poly_temp);
    }
  }
}



int aks (mpz_t n)
{
  /* Step 1: perfect power */
  if (mpz_perfect_power_p(n)) {
    return COMPOSITE;
  }
  /* Step 2: witness search */
  mpz_t r;
  mpz_init_set_ui(r, 2);
  unsigned int r_ui = 2;
  unsigned int logn_ui = mpz_sizeinbase(n, 2);
  mpz_t logn;
  mpz_init_set_ui(logn, logn_ui);
  mpz_t imax;			/* Upper bound of i = 4 * logn ^ 2*/
  mpz_init(imax);
  mpz_mul(imax, logn, logn);
  mpz_mul_ui(imax, imax, 4);
  Sieve sieve;			/* Sieve of Eratosthenes */
  initialize_sieve(&sieve);
  while (mpz_cmp(r, n) < 0) {
    if (mpz_divisible_p(n, r)) {
      mpz_clear(r);
      mpz_clear(logn);
      mpz_clear(imax);
      destroy_sieve(&sieve);
      return COMPOSITE;
    }
    if (sieve_primality_test(r_ui, &sieve) == PRIME) {
      int is_break = 0;
      mpz_t pwm;
      mpz_init(pwm);
      mpz_t i;
      mpz_init_set_ui(i, 1);
      while (mpz_cmp(i, imax) <= 0) {
	mpz_powm(pwm, n, i, r);
	if (mpz_cmp_ui(pwm, 1) != 0) {
	  is_break = 1;
	  break;
	}
	mpz_add_ui(i, i, 1);
      }
      mpz_clear(pwm);
      mpz_clear(i);
      if (is_break) {
	break;
      }
    }
    mpz_add_ui(r, r, 1);
    r_ui++;
  }
  mpz_clear(imax);
  destroy_sieve(&sieve);
  if (mpz_cmp(r, n) == 0) {
    mpz_clear(r);
    mpz_clear(logn);
    return PRIME;
  }
  /* Step 3: polynomial check */
  mpz_t amax;			/* Upper bound of a = 2 * sqrt(r) * logn */
  mpz_init_set_ui(amax, 2);
  mpz_t sqrtr;
  mpz_init(sqrtr);
  mpz_sqrt(sqrtr, r);
  mpz_mul(amax, amax, sqrtr);
  mpz_mul(amax, amax, logn);
  mpz_clear(sqrtr);
  mpz_clear(logn);
  mpz_t power_right;
  mpz_init(power_right);
  mpz_mod(power_right, n, r);
  unsigned int power_right_ui = mpz_get_ui(power_right);
  mpz_clear(r);
  mpz_clear(power_right);
  unsigned int is_returns = 0;
#ifndef USE_PARALLEL
  Polynomial* p_poly_right;
  Polynomial* p_poly_left; 
  Polynomial* p_poly_left_base;
  initialize_polynomial(&p_poly_right, power_right_ui);
  set_polynomial_coef_si(p_poly_right, power_right_ui, 1); /* X ^ (n % r) */
  initialize_polynomial(&p_poly_left_base, 1);
  set_polynomial_coef_si(p_poly_left_base, 1, 1); /* X */
  mpz_t a;
  mpz_init_set_ui(a, 1);
  mpz_t a_mod_n;
  mpz_init(a_mod_n);
  while (mpz_cmp(a, amax) <= 0) {
    mpz_mod(a_mod_n, a, n);
    set_polynomial_coef(p_poly_right, 0, &a_mod_n); /* X ^ (n % r) + a % n */
    set_polynomial_coef(p_poly_left_base, 0, &a); /* X + a */
    polynomial_modular_power(&p_poly_left, p_poly_left_base, n, r_ui);
    if (!is_equal_polynomial(p_poly_left, p_poly_right)) {
      destroy_polynomial(&p_poly_left);
      destroy_polynomial(&p_poly_right);
      destroy_polynomial(&p_poly_left_base);
      mpz_clear(amax);
      mpz_clear(a);
      mpz_clear(a_mod_n);
      return COMPOSITE;
    }
    destroy_polynomial(&p_poly_left);
    mpz_add_ui(a, a, 1);
  }
  destroy_polynomial(&p_poly_right);
  destroy_polynomial(&p_poly_left_base);
  mpz_clear(amax);
  mpz_clear(a);
  mpz_clear(a_mod_n);
#else
  mpz_t a;
#pragma omp parallel private(a)
{
  mpz_init_set_ui(a, omp_get_thread_num());

  for (; is_returns == 0 && mpz_cmp(a, amax) < 0; mpz_add_ui(a, a, omp_get_num_threads())) {
    mpz_t a_mod_n;
    mpz_init(a_mod_n);
    mpz_mod(a_mod_n, a, n);
    Polynomial* p_poly_right;
    Polynomial* p_poly_left;
    Polynomial* p_poly_left_base;
    initialize_polynomial(&p_poly_right, power_right_ui);
    set_polynomial_coef_si(p_poly_right, power_right_ui, 1);
    set_polynomial_coef(p_poly_right, 0, &a_mod_n);
    initialize_polynomial(&p_poly_left_base, 1);
    set_polynomial_coef_si(p_poly_left_base, 1, 1);
    set_polynomial_coef(p_poly_left_base, 0, &a);
    polynomial_modular_power(&p_poly_left, p_poly_left_base, n, r_ui);
    if (!is_equal_polynomial(p_poly_left, p_poly_right)) {
      #pragma omp atomic
      is_returns++;
    }
    mpz_clear(a_mod_n);
    destroy_polynomial(&p_poly_right);
    destroy_polynomial(&p_poly_left);
    destroy_polynomial(&p_poly_left_base);
  } /* end while a */
} /* end omp parallel block */
  mpz_clear(amax);
#endif

  /* Step 4: after all... */
  return is_returns > 0 ? COMPOSITE : PRIME;
}



/* int main (int argc, char* argv[]) */
/* { */
/*   mpz_t n; */
/*   mpz_init_set_str(n, argv[1], 10); */
/*   printf("%d\n", aks(n)); */
/*   mpz_clear(n); */
/*   return 0; */
/* } */



/* int main (int argc, char* argv[]) */
/* { */
/*   char n_str[1000]; */
/*   scanf("%s", n_str); */
/*   mpz_t n; */
/*   mpz_init_set_str(n, n_str, 10); */
/*   gmp_printf("%Zd\n", n); */
/*   printf("%d\n", aks(n)); */
/*   mpz_clear(n); */
/*   return 0; */
/* } */



int main (int argc __attribute__((unused)), char* argv[])
{
  int stats[2];
  int res;
  char n_str[1000];
  stats[COMPOSITE] = stats[PRIME] = 0;
  FILE* fp = fopen(argv[1], "r");
  mpz_t n;
  mpz_init(n);
  while (fscanf(fp, "%s", n_str) != EOF) {
    #ifdef USE_PARALLEL
    double elapsed_time = omp_get_wtime();
    #else
    clock_t start = clock();
    #endif
    mpz_set_str(n, n_str, 10);
    gmp_printf("%Zd: ", n);
    res = aks(n);
    printf("%d\n", res);
    stats[res]++;
    printf("Time: %f\n",
    #ifdef USE_PARALLEL
          omp_get_wtime() - elapsed_time
    #else
          (double)(clock() - start) / (double)CLOCKS_PER_SEC
    #endif
          );
  }
  mpz_clear(n);
  fclose(fp);
  printf("%-10s: %4d\n%-10s: %4d\n", "PRIME", stats[PRIME], "COMPOSITE", stats[COMPOSITE]);
  return 0;
}
