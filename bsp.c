#include "../stuff/labs-functions.c"

double dw(int i, int m, int k, double *xs){ 
	int j;
	double w=1, q=0;
		for(j=k-1; j<=k+m; j++) {
			q = xs[i] - xs[j];
			if (q==0) q=1.;
			w *= q;
		}
	return w;
}

double BasisSplain(double *xs, int m, double x) {
	int i, k=1;
	double B=0, p;

	for(i=k-1; i<=k+m; i++) {
		p = xs[i] - x;
		if (p >= 0)
			B+=((m+1)*pow(p,m))/dw(i,m,k,xs);
	}
	return B;
}

int main(){
	double x0[2]={-0.5,0.5},
	       x1[3]={-1.,0.,1.},
	       x2[4]={-1.5,-0.5,0.5,1.5},
	       x3[5]={-2.,-1.,0.,1.,2.};

	double myspline(double x) {
		// 2 параметр - степень сплайна
		return BasisSplain(x2, 2, x);
	}

	int N=333;
	double a=-2, b=2;

	plot_s(myspline, "Bsplain.dat", N, a, b, 100500);

	FILE* myplot = fopen("gnuplot.gp", "w");
	fprintf(myplot, "#!/usr/bin/gnuplot\n\nset term png\nset output \"1.png\"\n");
	fprintf(myplot, "plot 'Bsplain.dat' u 1:2 title 'функция' w l lt 9");
	fprintf(myplot, "\nexit");
	fclose(myplot);

	return 0;
}
