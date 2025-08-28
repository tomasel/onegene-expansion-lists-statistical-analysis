#include "int_mod.hpp"

double Integral(const double* arr, size_t length, double dx) {
	if (arr == nullptr || length == 0 || dx == 0)
		return 0;
	double result = 0;
	for (int i = 1; i < length; ++i)
		result += arr[i-1] + arr[i];
	return result * dx / 2;
}

double Module(const double *arr, size_t length) {
	if (arr == nullptr || length == 0)
		return 0;
	double result = 0;
	for (int i = 0; i < length; ++i)
		result += pow(arr[i], 2);
	return dsqrtl(result);
}