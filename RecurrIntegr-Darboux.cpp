#include<bits/stdc++.h>

double straight_recurrint(void)
{
	double I0=log(7./6), I1, I2;
	int n,N=31;
	
	for(n=1,I1=I0; n<=N; n++)
	{
		I2 = 1./n - 6*I1;
		I1=I2;
	}
	return I2;
}

double reverse_recurrint(void)
{
	double I62=0., I1, I2;
	int n,N=62;
	
	for(n=N,I2=I62; n>=31; n--)
	{
		I1 = 1./(6*n) - I2/6;
		I2=I1;
	}
	return I1;
}

double f(double x)
{
	double ans;
	ans=pow(x,31)/(x+6);
	return ans;
}

double darboux_sum(void)
{
	double x=0,s=0;
	int i,k=1e3;
	for(i=0;i<k;i++)
	{
		s+=f(x)/k;
		x+=1./k;
	}
	return s;
}

//lower or upper Darboux sum - if x is added before s or vice versa

double Darboux_sum(void)
{
	double x=0,s=0;
	int i,k=1e3;
	for(i=0;i<k;i++)
	{
		x+=1./k;
		s+=f(x)/k;
	}
	return s;
}

int main(void)
{
	printf("\nstraight method: I_31 = %g\n", straight_recurrint());
	printf("\nreverse method: I_31 = %.10lf\n", reverse_recurrint());
	printf("\nlower Darboux sum: %.14lf\n", darboux_sum());
	printf("\nupper Darboux sum: %.14lf\n\n", Darboux_sum());

	return 0;
}