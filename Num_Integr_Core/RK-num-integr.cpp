#include<bits/stdc++.h>
using namespace std;

double fx(double t, double x, double z) //x'=fx(y)
{
	return z;
}

double fz(double t, double x, double z) //z'=fz(y)
{
	return -x*(1-0.2*x*x)+cos(t);
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

double intersect(double x1, double x2, double z1, double z2)
{
	double k,ans;
	
	k = fabs(z1)/fabs(z2);
	ans = (x1+k*x2)/(1+k);
	
	return ans;
}

void finding_solution(fstream& f_out_table, int p, double tol, double start)
{
	double t0,t,tprev,T,Tprev;
	double x0,x1,x2,x,xprev,delta_x,delta_x_prev;
	double sx,sxcap,akx,bkx,bcapkx;
	double z0,z,zprev,sz,szcap,akz,bkz,bcapkz;
		
	double fac=0.98, facmax=1.5, facmin=0.7, fctr;
		
	int N_steps, flag;
		
	double err,h,hprev;
		
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

	for(T=Tprev=t=tprev=t0=0, x=xprev=x0=x1=x2=start, delta_x=delta_x_prev=0, z=zprev=z0=0, h=hprev=0.1, N_steps=0, flag=1; N_steps<15000;)
	{
		bkx=bkz=bcapkx=bcapkz=0;

		for(int i=0; i<7; i++)
		{
			akx=akz=0;
			sx=sz=sxcap=szcap=0;
				
				for(int j=0; j<7; j++)
				{
					akx+=a[i][j]*kx[j];
					akz+=a[i][j]*kz[j];
				}
					
				kx[i]=fx(t+h*c[i],x+h*akx,z+h*akz);
				kz[i]=fz(t+h*c[i],x+h*akx,z+h*akz);
					
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
			
		fctr=Factor(facmax,facmin,fac,tol,err,p);
			
		hprev=h;
		h*=fctr;

		if (err<=tol)
		{
			f_out_table<<t<<'\t'<<x<<'\t'<<z<<'\t'<<N_steps<<endl;
			
			xprev=x; zprev=z; tprev=t;
			x+=sx; z+=sz; t+=hprev;
			
			N_steps++;

			if(zprev*z<0 && x>0)
			{
				t=intersect(tprev,t,zprev,z);
				T=fabs(t-Tprev);
				Tprev=t;
				
				x1=x2;
				x2=intersect(xprev,x,zprev,z);
				
				delta_x_prev=delta_x;
				delta_x=fabs(x1-x2);
				
				//cout<<"\nIntersect: delta_x = "<<delta_x<<"\tT = "<<T<<endl;
				
				//cout<<"\nIntersect: x = "<<x2<<"\tdelta_t = "<<t-T<<endl;
				
				if(delta_x<1e-2)
				{
					//T=fabs(t-Tprev);
					cout<<"\nPeriod: "<<T<<endl;
					return;
				}
				
				if(delta_x>delta_x_prev && flag)
				{
					cout<<"\nWe are flying away\n"<<endl;
					h*=-1;
					flag=0;
				}
				
				x=xprev=x2; z=zprev=z0; N_steps++;
			}
		}
	}
	return;
}


void Solving_Equations(string name)
{
	fstream f_out_table;
	
	double tol=1e-11;
	
	int p=5;

	f_out_table.open(name, ios_base::out);

	if (f_out_table.is_open())
	{
		f_out_table<<"t\tx\tz\tN_steps"<<endl;
		finding_solution(f_out_table,p,tol,-0.1);
		f_out_table<<"\n\n";
		finding_solution(f_out_table,p,tol,-0.2);
		f_out_table<<"\n\n";
		finding_solution(f_out_table,p,tol,-0.5);
		f_out_table<<"\n\n";
		finding_solution(f_out_table,p,tol,-1.5);
		f_out_table<<"\n\n";
		//finding_solution(f_out_table,p,tol,-2.5);
		f_out_table<<"\n\n";
		//finding_solution(f_out_table,p,tol,5.5);
		f_out_table<<"\n\n";
		//finding_solution(f_out_table,p,tol,-5.5);
		f_out_table<<"\n\n";
	}
	else cout<<"\nCannot open file!\n\n"<<endl;
	
	f_out_table.close();
	
	return;
}

int main(void)
{
	Solving_Equations("task-1.8_table.txt");
	return 0;
}
