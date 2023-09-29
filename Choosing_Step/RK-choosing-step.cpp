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

double Error(int p, double x1, double z1, double x2, double z2)
{
	double error;
	error = sqrt((x1-x2)*(x1-x2)+(z1-z2)*(z1-z2))/(pow(2,p)-1);
	return error;
}

double Factor(double facmax, double facmin, double fac, double tol, double err, int p)
{
	double factor;
	factor = fmin(facmax, fmax(facmin, fac*pow(tol/err, 1./(p+1))));
	return factor;
}

void runge_kutta_step(fstream& f_out_table, fstream& f_out_glob, double T, int p, double tol)
{
	double t0,t;
	double x0,x,sx,sxcap,akx,bkx,bcapkx;
	double z0,z,sz,szcap,akz,bkz,bcapkz;
		
	double fac=0.98, facmax=1.5, facmin=0.7, fctr;
		
	int N_steps, flag;
		
	double err, h, h1, delta_glob;
		
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
	double bcap[7]={5179./57600, 0., 7571./16695, 393./640, -92097./339200, 187./2100, 1./40};
		
	double c[7]={0., 1./5, 3./10, 4./5, 8./9, 1., 1.};

	double kx[7]={0.,0.,0.,0.,0.,0.,0.};
	double kz[7]={0.,0.,0.,0.,0.,0.,0.};


	//t0=0; x0=0; z0=1; h=0.1; delta_glob=0; N_steps=0;

	for(x=x0=0, z=z0=1, t=t0=0, h=0.1, delta_glob=0, N_steps=0, flag=1; flag!=0;)
	{
		bkx=bkz=0; bcapkx=bcapkz=0;

		for(int i=0; i<7; i++)
		{
			akx=akz=0;
			sx=sz=sxcap=szcap=0;
				
				for(int j=0; j<7; j++)
				{
					akx+=a[i][j]*kx[j];
					akz+=a[i][j]*kz[j];
				}
					
				kx[i]=fx(x+h*akx,z+h*akz);
				kz[i]=fz(x+h*akx,z+h*akz);
					
				bkx+=b[i]*kx[i];
				bkz+=b[i]*kz[i];
					
				bcapkx+=bcap[i]*kx[i];
				bcapkz+=bcap[i]*kz[i];
		}

		sx+=h*bkx;
		sz+=h*bkz;

		sxcap+=h*bcapkx;
		szcap+=h*bcapkz;

		err=Error(p,sx,sz,sxcap,szcap);
		delta_glob+=err;
			
		fctr=Factor(facmax,facmin,fac,tol,err,p);
			
		h1=h;
		h*=fctr;
		
		if(t+h1>T)
		{
			cout<<"\nWe are here\n"<<endl;
			h1=T-t;
			flag=0;
		}

		if(err<=tol)
		{
			x+=sx;
			z+=sz;
			t+=h1;
			N_steps++;
		}
	}
	
	printf("\nt = %.10e\tT = %.10e\tt-T = %.10e\n",t,T,t-T);
	
	//f_out_table<<T<<'\t'<<t<<'\t'<<T-t<<'\t'<<tol<<'\t'<<N_steps <<'\t'<<h1<<'\t'<<fabs(x-Fx(T))<<'\t'<<fabs(z-Fz(T))<<'\t'<<x<<'\t'<<z<<endl;
	
	f_out_table<<T/M_PI<<"pi\t"<<tol<<'\t'<<N_steps<<'\t'<<h1<<'\t'<<fabs(x-Fx(T))<<'\t'<<fabs(z-Fz(T))<<'\t'<<x<<'\t'<<z<<endl;
		
	f_out_glob<<T/M_PI<<"pi\t"<<tol<<'\t'<<delta_glob<<endl;

	return;
}

double Runge_number(double x0, double x1, double x2)
{
	double ans;
	ans = fabs((x0-x1)/(x1-x2));
	return ans;
}

void Runge_Kutta_step(string name1, string name2, string name3)
{
	fstream f_out_table, f_out_numbers, f_out_glob;
	
	double T, tol;
	
	double x0,x1,x2;
	double z0,z1,z2;
	string a;
	double b,c,d,e,f,g;
	
	int p;

	f_out_table.open(name1, ios_base::out);
	f_out_glob.open(name3, ios_base::out);

	if (f_out_table.is_open() && f_out_glob.is_open())
	{
		f_out_table <<"T\ttol\tN_steps\tstep\tdelta_x\tdelta_z\tx\tz"<<endl;
		
		f_out_glob<<"T\ttol\tdelta_glob"<<endl;

		for(T=1e2*M_PI, p=5; T<=1e6*M_PI; T*=10)
		{
			for(tol=1e-7;tol>=1e-11;tol*=1e-2)
			runge_kutta_step(f_out_table,f_out_glob,T,p,tol);
		}
	}
	else cout<<"\nCannot open file!\n\n"<<endl;
	
	f_out_table.close();
	f_out_glob.close();
	
	f_out_table.open(name1, ios_base::in);
	f_out_numbers.open(name2, ios_base::out);
	
	if (f_out_table.is_open() && f_out_numbers.is_open())
	{
		f_out_numbers<<"T\tR_x\tR_z"<<endl;
		
		for(int i=0; i<4; i++)
		{
			getline(f_out_table,a);
			f_out_table>>a>>b>>c>>d>>e>>f>>x0>>z0;
			f_out_table>>a>>b>>c>>d>>e>>f>>x1>>z1;
			f_out_table>>a>>b>>c>>d>>e>>f>>x2>>z2;
			f_out_numbers<<a<<'\t'<<Runge_number(x0,x1,x2)<<'\t'<< Runge_number(z0,z1,z2)<<endl;
		}
	}
	else cout<<"\nCannot open file!\n\n"<<endl;
	
	f_out_table.close();
	f_out_numbers.close();
	
	return;
}

int main(void)
{
	Runge_Kutta_step("task-1.7_table.csv","task-1.7_numbers.csv","task-1.7_delta-glob.csv");
	return 0;
}
