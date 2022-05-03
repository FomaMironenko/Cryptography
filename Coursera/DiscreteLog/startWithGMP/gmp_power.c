#include <gmp.h>
#include <stdio.h>

#define MAX_DIGITS 4096

int main() {
    mpz_t base, pow, mod;
    mpz_inits(base, pow, mod, NULL);
    char base_c[MAX_DIGITS], pow_c[MAX_DIGITS], mod_c[MAX_DIGITS];
    printf("Enter the base:   > ");
    scanf("%s", base_c);
    printf("Enter the power:  > ");
    scanf("%s", pow_c);
    printf("Emter the modulo: > ");
    scanf("%s", mod_c);
    mpz_set_str(base, base_c, 10);
    mpz_set_str(pow, pow_c, 10);
    mpz_set_str(mod, mod_c, 10);
    
    mpz_t result;
    mpz_init(result);
    mpz_powm(result, base, pow, mod);
    gmp_printf("%Zd ^ %Zd (mod %Zd) == %Zd\n", base, pow, mod, result);
    
    mpz_clears(base, pow, mod, result, NULL);
    return 0;
}
