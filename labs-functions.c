#include <stdio.h>
#include <math.h>
#define sign(x) (x==0 ? 0 : (x<0 ? -1 : 1))

int factorial(int n) {
	if (n<2) return 1;
	else return n*factorial(n-1);
}

void plot(double (*f)(double x), FILE* output, double N, double a, double b, double infinity) {
	double st = (b-a)/N, x, y;
	for (x=a; x<=b; x+=st) {
		if (fabs(x)<1.E-15) {
			y=0.;
		}
		else {
			y=f(x);
		}
		if (fabs(y)>infinity) y=INFINITY;
		if (isnan(y)) fprintf(output, "\n");
		else fprintf(output, "%lf %lf\n", x, y);
	}
}

double ms(double (*f)(double x), double x1, double x2, double eps) {
	double x0, e0, la, y1, y2;

	int i=0;

	y1=f(x1); y2=f(x2);
	la = (x2-x1)/(y2-y1);
	x0 = (y2*x1-y1*x2)/(y2-y1);
	do {
		i++;
		e0=la*f(x0);
		x0-=e0;
	} while(fabs(e0)>eps);
	return x0;
}

double fr(double (*f)(double x), double x1, double x2, double eps) {
	int i=0;
	while(fabs(x2-x1)>eps) {
		i++;
		x1 = x2 - (x2-x1)*f(x2)/(f(x2)-f(x1));
		x2 = x1 - (x1-x2)*f(x1)/(f(x1)-f(x2));
	}
	
	return x2;
}

double bisect(double (*f)(double x), double a, double b, double eps) {
	int i=0;
	double x0;
	double fa = 0, fb = 0;
	while (fabs(b-a) > eps || fabs(fb-fa) > eps) {
		x0 = (b+a)/2.;
		fa = f(a); fb = f(x0);
		if (sign(fa)==sign(fb)) a=x0;
		else b=x0;
		i++;
	}
	return x0;
}

// -----------------------------------

double rect(double (*f)(double x), double a, double b, int N) {
	double st = (b-a)/N, sum = 0, x = a - st/2.;
	int i;

	for (i=0; i<N; i++) {
		x+=st; sum += f(x);
	}
	return sum * st;
}

double trap(double (*f)(double x), double a, double b, int N) {
	double st = (b-a)/N, sum = (f(a) + f(b))/2., x = a;
	int i;

	for (i=0; i<N-1; i++) {
		x+=st; sum += f(x);
	}
	return sum * st;
}

double parab(double (*f)(double x), double a, double b, int N) {
	double st = (b-a)/N, s1=0, s2=0, x = a;
	int i;

	for (i=0; i<N-1; i++) {
		x+=st; s1 += f(x); s2 += f(x-st/2.);
	}

	s2 += f(b-st/2.);
	return st/6.*(f(a)+f(b)+2*s1+4*s2);
}

//-------------------------

double dfdx(double (*f)(double x), double x, double h) {
	return (f(x+h)-f(x-h))/(2*h);
}

double d2fdx2(double (*f)(double x), double x, double h) {
	return (f(x+h)-2*f(x)+f(x-h))/(h*h);
}

//--------------------------

double ljx(double* mx, int n, int j, double x) {
	double v=1.0; int k;
	for (k=0; k<=n; k++) {
		if (k!=j) v*=(x-mx[k])/(mx[j]-mx[k]);
	}
	return v;
}

double lx(double* mx, double* my, int n, double x) {
	double v=0.0; int j;
	for (j=0; j<=n; j++) {
		v+=my[j]*ljx(mx, n, j, x);
	}
	return v;
}

double dfdx_recur(double (*f)(double x), int i, double x, double h) {
	if (i==0) return f(x);
	else return (
		 dfdx_recur(f, i-1, x+h, h)
		-dfdx_recur(f, i-1, x-h, h)
		)/(2*h);
}
