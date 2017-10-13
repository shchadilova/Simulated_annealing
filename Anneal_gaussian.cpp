#include "Anneal_headers.cpp"
#include <iomanip>
#include <math.h>

int const NK=20, Ntheta=10, Nsweeps1=10, Nsweeps2=30;//20 10 100
double const Kmax=2.; 
double vr=1.;
double M=20.; 
double q=vr*M;
double SUM=0;
double dg=0.01;
double GammaMAX=1.;//1
double SumComplete=0.;
double Eth=-1000.;


ofstream debug("debug.dat");

double kT=0.001, E=0;

complex ** alpha, ** alpha_aux;
complex ** gamma0;
double ** beta;


double V(double k)
{
	return (1.+1./M)*sqrt(sqrt(k*k/(2+k*k)))/sqrt(Pi); 
;}

double omega(double k)
{
	return k*sqrt(1+k*k/2.);
;}

double W(double k, double th)
{
	return omega(k)+k*k/(2.*M)-k*cos(th)*(q-SUM)/M; 
;}

double W2(double k, double th)
{
	return omega(k)+k*k/(2.*M)-k*cos(th)*(q-.5*SUM)/M; 
;}


void Ini()
{
// alpha- displacement
// gamma- diagonal part of mu12

	{
	alpha=new complex *[NK];
	gamma0=new complex *[NK];
    beta=new double *[NK];
	for (int k=0; k<NK; k++)
	{
		alpha[k]=new complex [Ntheta];
		gamma0[k]=new complex [Ntheta];
		beta[k]=new double [Ntheta];
		for (int th=0; th<Ntheta; th++){ 
		alpha[k][th]=0.; gamma0[k][th]=0.;beta[k][th]=0.;}
	;}
	;}
;}


double MU11(int k, int th, complex dgamma)
{
    double Delta=abs(gamma0[k][th]+dgamma)*abs(gamma0[k][th]+dgamma)
    -beta[k][th]*beta[k][th];
     return  sinhc(Delta)*sinhc(Delta)
     *abs(gamma0[k][th]+dgamma)*abs(gamma0[k][th]+dgamma);
;}

complex MU12(int k, int th, complex dgamma)
{
    double Delta=abs(gamma0[k][th]+dgamma)*abs(gamma0[k][th]+dgamma)
    -beta[k][th]*beta[k][th];
     return  -I*(gamma0[k][th]+dgamma)*(sinhc(Delta)*sinhc(Delta)*beta[k][th]
            +I*sinhc(2.*Delta));
;}

double calc_SumComplete()
{
    double sum=0.;
	{
		for (int k=0; k<NK; k++)
		for (int th=0; th<Ntheta; th++)
		{
			double k_=k*Kmax/NK, th_=th*Pi/Ntheta;      
			            
            double mu11=MU11(k,th,0.); // because it is mu22
			
			sum+=k_*cos(th_)*2.*Pi*k_*k_*sin(th_)*(Kmax/NK)*(Pi/Ntheta)*
			(norm2(alpha[k][th])+mu11);                    

		;}
	;}
	return sum;
;}


double Htdcomplete(double & sum)
{
    SUM=calc_SumComplete();

	complex s=0.;
	{
	    double mu11=0.;
	    complex mu12=0.;
	    complex alph=0.;
	    complex gamm=0.;
	    
		for (int k=0; k<NK; k++)
		for (int th=0; th<Ntheta; th++)
		{
			double k_=k*Kmax/NK, th_=th*Pi/Ntheta;
			
			alph=alpha[k][th];
			
            mu11=MU11(k, th, 0.);
            mu12=MU12(k, th, 0.);
            
            double dV=2.*Pi*k_*k_*sin(th_)*(Kmax/NK)*(Pi/Ntheta);
            double kx=k_*cos(th_);
            
			s+=(2.*real(V(k_)*alph)+norm2(alph)*W2(k_, th_)+W2(k_, th_)*mu11)
			*dV;
			s+=kx*kx/(2.*M)*
			(conj(alph)*mu11*alph+
			conj(alph)*mu12*conj(alph)+
			alph*conj(mu12)*alph+
			alph*mu11*conj(alph)+
			mu11*mu11+mu12*conj(mu12))
			*dV;
			
	}
	;}
	
	
	return real(s);
;}


double dH_alpha(int k, int th, complex da)
{
	double k_=k*Kmax/NK, th_=th*Pi/Ntheta;
	
    double mu11=MU11(k, th, 0.);
    complex mu12=MU12(k, th, 0.);
	
	double da2=2.*real(conj(da)*alpha[k][th])+norm2(da);
	complex daa=2.*alpha[k][th]*da+da*da;
	double dV=2.*Pi*k_*k_*sin(th_)*(Kmax/NK)*(Pi/Ntheta);
	double kx=k_*cos(th_);
	
	double s=2.*V(k_)*real(da)+da2*W(k_,th_)+da2*da2*kx*kx*dV/(2.*M);
	    s+=kx*kx/(2.*M)*
	    (2*mu11*da2+
		2*real(daa*conj(mu12))
		);
	return s*dV;
	
	
;}


double dH_gamma(int k, int th, complex dgamma)
{
     SUM=calc_SumComplete();
     
     complex alph=0.;
     
     double k_=k*Kmax/NK, th_=th*Pi/Ntheta;
     
     double mu11=MU11(k, th, 0.);
     complex mu12=MU12(k, th, 0.);
     
     double dmu11=MU11(k, th, dgamma)-mu11;
     complex dmu12=MU12(k, th, dgamma)-mu12;
     
     double dV=2.*Pi*k_*k_*sin(th_)*(Kmax/NK)*(Pi/Ntheta);
     double kx=k_*cos(th_);
     
     alph=alpha[k][th];
     
     complex s=(W(k_, th_)*dmu11)+kx*kx*dmu11*dmu11*dV/(2.*M);
     s+=kx*kx/(2.*M)*
     (conj(alph)*dmu11*alph+
     conj(alph)*dmu12*conj(alph)+
     alph*conj(dmu12)*alph+
     alph*dmu11*conj(alph)+
     2.*dmu11*mu11+dmu11*dmu11+
     dmu12*conj(mu12)+mu12*conj(dmu12)+dmu12*conj(dmu12))
     ;

	return real(s)*dV;
;}

void step_alpha()
{
	int k=rnd(NK), th=rnd(Ntheta);
	
	SUM=calc_SumComplete();
	
	complex da(rnd()-0.5, rnd()-0.5);
	
	double dE=dH_alpha(k, th, da);
	
	if (exp(-dE/kT)>rnd()) 
	{
		alpha[k][th]+=da; E+=dE;
		SUM=calc_SumComplete();
	;}
;}


void step_gamma()
{
	int k=rnd(NK), th=rnd(Ntheta);
	
	SUM=calc_SumComplete();
	
	complex dgamma((rnd()-0.5)*dg, (rnd()-0.5)*dg);
	
	double dE=dH_gamma(k, th, dgamma);
	
	if (exp(-dE/kT)>rnd() && abs(gamma0[k][th])<GammaMAX) 
	{
        gamma0[k][th]+=dgamma; E+=dE;
	    SUM=calc_SumComplete();
	;}
;}

void sweep_alpha()
{
	for (int i=0; i<NK*Ntheta; i++) step_alpha();
;}

void sweep_gamma()
{
	for (int i=0; i<NK*Ntheta; i++) step_gamma();
;}


int main(int argc, char **argv)
{
	if (argc>1) {vr=double(atof(argv[1])); M=double(atof(argv[2]));
	cout<<"vr="<<vr<<"\n"<<"M="<<M<<"\n"<<flush;}
	Ini();

    // read stuff from file
    if (1==0) 
    {
        std::ifstream ou_file;
        ou_file.open("alpha.dat");
        double tmp;
		for (int k=1; k<NK; k++)
		{
            for (int th=1; th<Ntheta; th++) 
            {
              ou_file >> tmp; // k*Kmax/NK
              ou_file >> tmp; // th*Pi/Ntheta
              ou_file >> tmp; // abs(alpha[k][th])
              alpha[k][th]=tmp;
              }
		};
        ou_file.close();
			
		ou_file.open("gamma.dat");
		for (int k=1; k<NK; k++)
		{
			for (int th=1; th<Ntheta; th++) 
			    {
                ou_file >> tmp; // k*Kmax/NK
                ou_file >> tmp; // th*Pi/Ntheta
                ou_file >> tmp; // abs(gamma0[k][th])
                gamma0[k][th]=tmp;
                };
		};
        ou_file.close();
    };
    
	q=vr*M;
	for (kT=0.1; kT>1e-12; kT/=1.2)
	for (int i=0; i<10; i++) //100
	{   
//	    cout <<  calc_SumComplete() <<"\n"<<endl;
	    
		for (int i=0; i<Nsweeps1; i++) 
		{
		    sweep_alpha();
		    double SUM=calc_SumComplete();
//		    double hc=Htdcomplete(SUM);
//		    debug<<E-hc <<"\t"<<SUM<< "\n"<<flush;
		}
		
		debug<<"===================================="<<"\n"<<flush;
		
		for (int i=0; i<Nsweeps2; i++) 
		{
		    sweep_gamma();
		    double SUM=calc_SumComplete();
//		    double hc=Htdcomplete(SUM);
//		    debug<<E-hc<<"\t"<<SUM<<"\n"<<flush;
		}
		
		debug<<"===================================="<<"\n"<<flush;
		
        SUM=calc_SumComplete();
        double hc1=Htdcomplete(SUM);
		cout<< std::setprecision(9)<<E<<"\t"<<SUM<<"\t"<<E-hc1<<"\n"<<flush;
		E=hc1;
		if(E<Eth) { std::cerr << "EXIT" << std::endl; exit(1); }
	;}
	
	{
		ofstream ou_file("alpha.dat");
		for (int k=1; k<NK; k++)
		{
			for (int th=1; th<Ntheta; th++)
			ou_file<<k*Kmax/NK<<"\t"<<th*Pi/Ntheta<<"\t"<<abs(alpha[k][th])<<"\n";
			ou_file<<"\n"<<flush;
		;}
		;}	
		
	{
		ofstream ou_file("gamma.dat");
		for (int k=1; k<NK; k++)
		{
			for (int th=1; th<Ntheta; th++)
			ou_file<<k*Kmax/NK<<"\t"<<th*Pi/Ntheta<<"\t"
			<<abs(gamma0[k][th])<<"\n";
			ou_file<<"\n"<<flush;
		;}
	;}	
	
	
	return 0;
;}
