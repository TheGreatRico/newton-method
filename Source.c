#include <stdio.h>
#include <stdlib.h>
#define _USE_MATH_DEFINES
#include <math.h>
const double EPSILON = 10e-9;

double f1(double x1, double x2) {
	return atan(x1 - 1);
}
double f2(double x1, double x2) {
	return atan(-x1 + x2 * x2);
}
double df1_dx1(double x1, double x2) {
	return (f1(x1 + EPSILON, x2) - f1(x1, x2)) / EPSILON;
}
double df1_dx2(double x1, double x2) {
	return (f1(x1, x2 + EPSILON) - f1(x1, x2)) / EPSILON;
}
double df2_dx1(double x1, double x2) {
	return (f2(x1 + EPSILON, x2) - f2(x1, x2)) / EPSILON;
}
double df2_dx2(double x1, double x2) {
	return (f2(x1, x2 + EPSILON) - f2(x1, x2)) / EPSILON;
}

void multinewton(double x_arr[]) {
	double x1_0 = x_arr[0], x2_0 = x_arr[1];
	int t = 0, k = 0;
	double f_0[2], df1_dx1_0, df1_dx2_0, df2_dx1_0, df2_dx2_0, jacobian[2][2], coef, jacobian_inverse[2][2], j_inv_mult_f_0s[2], x_new[2], error[2] = { 10e10, 10e10 };
	//while (fabs(error[0]) > EPSILON || fabs(error[1]) > EPSILON) {
	while (sqrt(f1(x1_0, x2_0) * f1(x1_0, x2_0) + f2(x1_0, x2_0) * f2(x1_0, x2_0)) > EPSILON) {
		if (t > 10e3) {
			printf("\n\n\t\titeration limit exceeded allotted value, exiting program...\n\n");
			exit(0);
		}
		f_0[0] = f1(x1_0, x2_0);
		f_0[1] = f2(x1_0, x2_0);
		df1_dx1_0 = df1_dx1(x1_0, x2_0), df1_dx2_0 = df1_dx2(x1_0, x2_0), df2_dx1_0 = df2_dx1(x1_0, x2_0), df2_dx2_0 = df2_dx2(x1_0, x2_0);
		jacobian[0][0] = df1_dx1_0;
		jacobian[1][0] = df1_dx2_0;
		jacobian[0][1] = df2_dx1_0;
		jacobian[1][1] = df2_dx2_0;
		coef = 1 / (df1_dx1_0 * df2_dx2_0 - df1_dx2_0 * df2_dx1_0);
		jacobian_inverse[0][0] = df2_dx2_0 * coef;
		jacobian_inverse[1][0] = -df1_dx2_0 * coef;
		jacobian_inverse[0][1] = -df2_dx1_0 * coef;
		jacobian_inverse[1][1] = df1_dx1_0 * coef;
		j_inv_mult_f_0s[0] = jacobian_inverse[0][0] * f_0[0] + jacobian_inverse[1][0] * f_0[1];
		j_inv_mult_f_0s[1] = jacobian_inverse[0][1] * f_0[0] + jacobian_inverse[1][1] * f_0[1];
		for (k = 0; k < 60; k++) {
			x_new[0] = x1_0 - pow(2., -k) * j_inv_mult_f_0s[0];
			x_new[1] = x2_0 - pow(2., -k) * j_inv_mult_f_0s[1];
			t += 1;
			printf("\tx1 = %.2f,\tx2 = %.2f,\titeration: %i\n", x_new[0], x_new[1], t);
			if (fabs(f1(x_new[0], x_new[1])) <= fabs(f1(x1_0, x2_0)) && fabs(f2(x_new[0], x_new[1])) <= fabs(f2(x1_0, x2_0))) {
				//if (sqrt(f1(x_new[0], x_new[1])*f1(x_new[0], x_new[1]) + f2(x_new[0], x_new[1])*f2(x_new[0], x_new[1])) < sqrt(f1(x1_0, x2_0)*f1(x1_0, x2_0) + f2(x1_0, x2_0)*f2(x1_0, x2_0))) {
				x1_0 = x_new[0];
				x2_0 = x_new[1];
				printf("\n\t\tBREAKING! on k = %i\n\n", k);
				break;
			}
		}
		/*
		error[0] = x_new[0] - x1_0;
		error[1] = x_new[1] - x2_0;
		x1_0 = x_new[0];
		x2_0 = x_new[1];
		*/
	}
	x_arr[0] = x1_0;
	x_arr[1] = x2_0;
}

int main() {
	double x[2] = { 50, 50 }; // initial values
	multinewton(x);
	printf("\n\n\tx1 = %.2f,\t\tx2 = %.2f\n\tf1(x1,2) = %.2f,\tf2(x1,2) = %.2f", x[0], x[1], f1(x[0], x[1]), f2(x[0], x[1]));
	return 1;
}