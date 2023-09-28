#include<bits/stdc++.h>
using namespace std;

double fx(double x, double z) //x'=fx(y)
{
	return z;
}

double fz(double x, double z) //z'=fz(y)
{
	return -x;
}

double Fx(double t) //x=sin(t)
{
	return sin(t);
}

double Fz(double t) //z=cos(t)
{
	return cos(t);
}

//T=2pi, 1e1pi, 1e2pi, ... , 1e5pi; h=1e-1, 1e-2, 1e-3

void euler(fstream& f_out_plot, fstream& f_out_table, double h, double T)
{
	if (f_out_plot.is_open() && f_out_table.is_open())
	{
		double t0,t;
		double x0,x1,x2;
		double z0,z1,z2;

		t0=0; x0=0; z0=1;


		for(x1=x0, z1=z0, t=t0; t<=T; t+=h)
		{
			x2=x1+h*fx(x1,z1);
			z2=z1+h*fz(x1,z1);

			f_out_plot<<T/M_PI<<"pi\t"<<h<<'\t'<<x2<<'\t'<<z2<<endl;

			x1=x2; z1=z2;
		}
		
		f_out_plot<<"\n\n";
		f_out_table<<T/M_PI<<"pi\t"<<h<<'\t'<<abs(x2-Fx(T))<<'\t'<<abs(z2-Fz(T))<<endl;
	}
	else
	{
		cout<<"\nCannot open file!\n\n"<<endl;
	}
	return;
}

void Euler(string name1, string name2)
{
	fstream f_out_plot, f_out_table;
	double h, T;

	f_out_plot.open(name1, ios_base::out);
	f_out_table.open(name2, ios_base::out);

	f_out_table<<"T\th\tdelta_x\tdelta_z"<<endl;

	if (f_out_plot.is_open() && f_out_table.is_open())
	{
		for(h=1e-1;h>=1e-3;h*=1e-1)
		{
			for(T=2*M_PI;T<=1e6*M_PI;)
			{
				euler(f_out_plot, f_out_table, h, T);
				if (T==2*M_PI) T+=8*M_PI;
				else T*=10;
			}
		}
	}
	else cout<<"\nCannot open file!\n\n"<<endl;
	f_out_plot.close();
	f_out_table.close();
	return;
}


int main(void)
{
	Euler("data_new.txt","table");
	return 0;
}
