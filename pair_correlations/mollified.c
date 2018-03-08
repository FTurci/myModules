// Radial distribution functions for ctypes
// 
// Francesco Turci (2018)
// 
// gcc -shared -fPIC mollified.c -o mollified.so
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

// total radial distribution function
void radial_distribution(double *gr, double *r, double *x, double* y, double* z, double *box , int nbins, int N){
	
	double Lx,Ly,Lz, Lmax;
	double dx, dy,dz,dr,binwidth;
	double volume,Ideal;
	double rho,V;
	int i,j, index;

	Lx = box[0]; Ly = box[1]; Lz = box[2];
	V = Lx*Ly*Lz;
	rho = N/V;

	Lmax = Lx; if (Ly<Lmax) Lmax=Ly; if (Lmax<Lz) Lmax=Lz;
	Lmax/=2.;

	binwidth = Lmax/nbins;

	for ( i = 0; i < N-1; ++i)
	{
		for ( j = i+1; j < N; ++j)
		{
			dx=x[i]-x[j]; dy=y[i]-y[j]; dz=z[i]-z[j];

			if(dx>Lx/2) dx-=Lx; if(dx<-Lx/2) dx+=Lx;
			if(dy>Ly/2) dy-=Ly; if(dy<-Ly/2) dy+=Ly;
			if(dz>Lz/2) dz-=Lz; if(dz<-Lz/2) dz+=Lz;

			dr=sqrt(dx*dx+dy*dy+dz*dz);
			if(dr<Lmax){
				index=(int)(dr/binwidth);
				gr[index]+=2.;
			}
		}
	}
	
	for (i = 0; i < nbins; ++i)
	{
		r[i]=binwidth*(i+0.5);
		volume=((i+1)*(i+1)*(i+1)-i*i*i)*binwidth*binwidth*binwidth;
		Ideal=(4./3.)*M_PI*volume*rho;
		gr[i]/=(double)(N*Ideal);
	}
	
}

// mollified radial distribution for each particle
void  mollified_radial_distribution(double *gr, double *r, double *x, double* y, double* z, double *box , int *int_params, double *double_params){

	double Lx,Ly,Lz, Lmax;
	double dx, dy,dz,dr,binwidth;
	double volume,Ideal;
	double rho,V;

	// unpack parameters
	double stdev = double_params[0];
	double rcut = double_params[1];

	int i = (int)int_params[0];
	int N = (int)int_params[1];
	int nbins = (int)int_params[2];
	// printf("i %d nbins %d N %d\n",i,nbins,N);
	// printf("stdev %g rcut %g\n",stdev,rcut);

	double Sigma = 2*stdev*stdev;
	double Norm = sqrt(2*M_PI*stdev*stdev);


	Lx = box[0]; Ly = box[1]; Lz = box[2];
	V = Lx*Ly*Lz;
	rho = N/V;

	binwidth = rcut/nbins;
	// prepare g(r) and r
	for ( int k = 0; k < nbins; ++k){
		r[k] = binwidth*(k+0.5);
		gr[k] = 0;
	}



	for (int j = 0; j < N; ++j)
		if (j!=i)
		{
			dx=x[i]-x[j]; dy=y[i]-y[j]; dz=z[i]-z[j];

			if(dx>Lx/2) dx-=Lx; if(dx<-Lx/2) dx+=Lx;
			if(dy>Ly/2) dy-=Ly; if(dy<-Ly/2) dy+=Ly;
			if(dz>Lz/2) dz-=Lz; if(dz<-Lz/2) dz+=Lz;

			dr=sqrt(dx*dx+dy*dy+dz*dz);
			
			for (int k = 0; k < nbins; ++k)
			{	
				gr[k]+=exp(-(r[k]-dr)*(r[k]-dr)/Sigma);
			}
		}

	for (int k = 0; k < nbins; ++k)
	{
		gr[k]/=4*M_PI*rho*r[k]*r[k]*Norm;
		
	}
}

double trapz(double *y,double*x,int nbins){
	double sum=0;
	for (int k = 0; k < nbins-1; ++k)
	{
		sum += (y[k+1]+y[k])*(x[k+1]-x[k]);
	}
	return 0.5*sum;

}
void s2 (double *s2,double *r,double* gr, double *x, double* y, double* z, double *box , int *int_params, double *double_params){
	
	double epsilon = 1e-8;
	double Lx,Ly,Lz, Lmax;
	double dx, dy,dz,dr,binwidth;
	double volume,Ideal;
	double rho,V;
	// unpack parameters
	double stdev = double_params[0];
	double rcut = double_params[1];
	

	int N = (int)int_params[0];
	int nbins = (int)int_params[1];

	double *integrand = malloc(sizeof(double)*nbins);

	Lx = box[0]; Ly = box[1]; Lz = box[2];
	V = Lx*Ly*Lz;
	rho = N/V;

	Lmax = Lx; if (Ly<Lmax) Lmax=Ly; if (Lmax<Lz) Lmax=Lz;
	Lmax/=2.;

	
	double Sigma = 2*stdev*stdev;
	double Norm = sqrt(2*M_PI*stdev*stdev);


	Lx = box[0]; Ly = box[1]; Lz = box[2];
	V = Lx*Ly*Lz;
	rho = N/V;

	binwidth = rcut/nbins;

	// prepare g(r) and r
	for ( int k = 0; k < nbins; ++k){
		r[k] = binwidth*(k+0.5);
		integrand[k]=0;
	}

	for (int i = 0; i < N-1; ++i)
	{
	
		// if (iprintf("\bIn progress %d% ", i*100./N);
		// fflush(stdout);
		for (int j = i+1; j < N; ++j)
		{
			dx=x[i]-x[j]; dy=y[i]-y[j]; dz=z[i]-z[j];

			if(dx>Lx/2) dx-=Lx; if(dx<-Lx/2) dx+=Lx;
			if(dy>Ly/2) dy-=Ly; if(dy<-Ly/2) dy+=Ly;
			if(dz>Lz/2) dz-=Lz; if(dz<-Lz/2) dz+=Lz;

			dr=sqrt(dx*dx+dy*dy+dz*dz);
			
			for (int k = 0; k < nbins; ++k)
			{	
				gr[i*nbins+k]+=exp(-(r[k]-dr)*(r[k]-dr)/Sigma);				
				gr[j*nbins+k]+=exp(-(r[k]-dr)*(r[k]-dr)/Sigma);
			}
		}
			// printf("\r");
	}

	for (int i = 0; i < N; ++i) {
		for (int k = 0; k < nbins; ++k) gr[i*nbins+k]/=4*M_PI*rho*r[k]*r[k]*Norm;
	}
	printf("\n");
	
	// integrate (trapezoidal rule)
	for (int i = 0; i < N; ++i)
	{	
		for (int k = 0; k < nbins; ++k)
		{	
			// integrand[k]=0;
			if (gr[i*nbins+k]<epsilon) integrand[k]=0;
			else{
				integrand[k] = gr[i*nbins+k]*log(gr[i*nbins+k])-gr[i*nbins+k]+1.;
			}
		}
		s2[i] =-2.*M_PI *rho*trapz(integrand,r,nbins);
	}
	free(integrand);
	
}
// compute locally averaged value with a cutoff radius for a scalar field defined for every particle 
void local_average(double *output, double *input, double *x, double* y, double* z, double *box , int N, double cutoff){

	double Lx,Ly,Lz;
	double dx, dy,dz,dr2;
	double sum; int count;
	double cutoff2= cutoff*cutoff;

	printf("%d\n",N );
	Lx = box[0]; Ly = box[1]; Lz = box[2];
	

	for (int i = 0; i < N; ++i)
	{
		sum = 0;
		count = 0;
		 // find the neighbours
		for (int j = 0; j < N; ++j){
			dx=x[i]-x[j]; dy=y[i]-y[j]; dz=z[i]-z[j];

			if(dx>Lx/2) dx-=Lx; if(dx<-Lx/2) dx+=Lx;
			if(dy>Ly/2) dy-=Ly; if(dy<-Ly/2) dy+=Ly;
			if(dz>Lz/2) dz-=Lz; if(dz<-Lz/2) dz+=Lz;

			dr2=dx*dx+dy*dy+dz*dz;
			
			if (dr2<cutoff2){
				sum += input[j];
				count++;
			}


		}
		output[i] = sum/count;
	}

}


void check(double my){
	printf("%g\n",my );
}
