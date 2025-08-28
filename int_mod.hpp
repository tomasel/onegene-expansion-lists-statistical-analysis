#ifndef TESI_INT_MOD_HPP
#define TESI_INT_MOD_HPP
#include <cstdlib>
#include <cmath>

// Calcolo dell'integrale col metodo dei trapezi, dx costante
double Integral(double* arr, size_t length, double dx);

// Calcolo del modulo di un vettore
double Module(double *arr, size_t length);

#endif
