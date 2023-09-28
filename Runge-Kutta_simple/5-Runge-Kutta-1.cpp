#include<bits/stdc++.h>
using namespace std;

double g(double x) //y'=g(x)
{
	return 7*(2*x+1)/(2*sqrt(pow(x,2)+x+1))+9*exp(-x)*(cos(x)-sin(x))+2/(x*pow(log(x),3));
}

double y(double x) //y=f
{
	return 7*sqrt(pow(x,2)+x+1)+9*exp(-x)*sin(x)-1/pow(log(x),2);
}

double Rungy(double yhh, double y2h, int s)
{
	return fabs((yhh-y2h)/(1-pow(2,s)));
}

void runge_kutta(fstream& f_out, double h, int s) // x \in [x0, X] = [0.5, 0.9]
{
		double x0, x, X;
		double y0, y1, yhh, y2h;
		double bk, bk2, k[s], k2[s];
		
		double a[7][7]={
		{0.,0.,0.,0.,0.,0.,0.},
		{1./5, 0.,0.,0.,0.,0.,0.},
		{3./40, 9./40, 0.,0.,0.,0.,0.},
		{44./45, -56./15, 32./9, 0.,0.,0.,0.},
		{19372./6561, -25360./2187, 64448./6561, -212./729, 0.,0.,0.},
		{9017./3168, -355./33, 46732./5247, 49./176, -5103./18656, 0.,0.},
		{35./384, 0., 500./1113, 125./192, -2187./6784, 11./84, 0.}
		};
		
		double b[7]={35./384, 0., 500./1113, 125./192, -2187./6784, 11./84, 0.};
		double c[7]={0., 1./5, 3./10, 4./5, 8./9, 1., 1.};

		x0=0.5; X=0.9; y0=y(x0);

		for(y1=yhh=y2h=y0, x=x0; x<=X; x+=h)
		{
			bk=bk2=0;
			
			for(int i=0; i<s; i++)
			{
					k[i]=g(x+c[i]*h);
					k2[i]=g(x+c[i]*2*h);
					
					bk+=b[i]*k[i];
					bk2+=b[i]*k2[i];
			}

			yhh=y1;
			y1+=h*bk;
			yhh+=h*bk;
			y2h+=2*h*bk2;
			
			f_out<<h<<'\t'<<fabs(y1-y(x))<<'\t'<<Rungy(yhh,y2h,s)<<endl;
		}
		f_out<<"\n\n";
	return;
}

void Runge_Kutta(string name)
{
	fstream f_out;
	f_out.open(name, ios_base::out);
	if (f_out.is_open())
	{
		int s = 7;
		
		f_out<<"h\tdelta_y\tdelta_Runge"<<endl;
		
		for(double h=1e-1; h>=1e-3; h*=1e-1) runge_kutta(f_out, h, s);
	}
	
	else cout<<"\nCannot open file!\n\n"<<endl;
	
	f_out.close();
	
	return;
}

int main(void)
{
	Runge_Kutta("task-1.6a_out");
	return 0;
}
