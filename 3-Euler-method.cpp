#include<bits/stdc++.h>
using namespace std;

double g(double x) //y'=g
{
	return 7*(2*x+1)/(2*sqrt(pow(x,2)+x+1))+9*exp(-x)*(cos(x)-sin(x))+2/(x*pow(log(x),3));
}

double y(double x) //y=f
{
	return 7*sqrt(pow(x,2)+x+1)+9*exp(-x)*sin(x)-1/pow(log(x),2);
}

void euler(string name, double h) // x \in [x0, X] = [0.5, 0.9]
{
	fstream f_out;
	f_out.open(name, ios_base::out);
	if (f_out.is_open())
	{
		double y0,y1,y2,x0,x,X;

		x0=0.5; X=0.9; y0=y(x0);
		
		for(y1=y0,x=x0;x<X;x+=h)
		{
			y2=y1+h*g(x);
			f_out<<"\n"<<x<<"\t"<<y2;
			y1=y2;
		}
		printf("\nh = %g\tdelta = %g\n",h,abs(y2-y(X)));
		printf("\nData for h=%g saved in a file.\n",h);
	}
	else cout<<"\nCannot open file!\n\n"<<endl;
	f_out.close();
	return;
}

int main(void)
{
	euler("fh1e-1", 1e-1);
	euler("fh1e-2", 1e-2);
	euler("fh1e-3", 1e-3);
	printf("\n");
	return 0;
}
