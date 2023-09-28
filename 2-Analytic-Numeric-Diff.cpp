#include<bits/stdc++.h>
using namespace std;

double f(double x)
{
	return 7*sqrt(pow(x,2)+x+1)+9*exp(-x)*sin(x)-1/pow(log(x),2);
}

double deriv(double x)
{
	return 7*(2*x+1)/(2*sqrt(pow(x,2)+x+1))+9*exp(-x)*(cos(x)-sin(x))+2/(x*pow(log(x),3));
}

double R1(double x0, double h)
{
	return abs(deriv(x0)-(f(x0+h)-f(x0))/h);
}

void table(string name, double x0)
{
	fstream f_out;
	f_out.open(name, ios_base::out);
	if (f_out.is_open())
	{
		double h;
		for(h=1;h>1e-20;h*=0.1)
		{
			f_out<<"\n"<<h<<"\t"<<R1(x0,h);
		}
	}
	else cout<<"\nCannot open file!\n\n"<<endl;
	f_out.close();
	return;
}

int main(void)
{
	table("output",M_PI); //f'(pi)
	return 0;
}
