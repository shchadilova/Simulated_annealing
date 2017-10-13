#include "Anneal_headers.cpp"

int const NK=2, Ntheta=2, NsweepsVector=10, NsweepsMatrix=100;//50 30 100000
double const Kmax=5; 
double M=20; double q=10;
double SUM=0;
double SumComplete=0.;


ofstream debug("debug.dat");

double kT=0.001, E=0;

complex ** alpha, ** alpha_aux;
complex **** mu11, **** mu12;


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
	{
	alpha=new complex *[NK];
	alpha_aux=new complex *[NK];
	for (int k=0; k<NK; k++)
	{
		alpha[k]=new complex [Ntheta];
		alpha_aux[k]=new complex [Ntheta];
		for (int th=0; th<Ntheta; th++) alpha[k][th]=0.
	;}
	;}

// mu11 - Hermitian
// mu12 - Arbitrary
	
	mu11=new complex ***[NK];
	mu12=new complex ***[NK];
	for (int k=0; k<NK; k++)
	{
		mu11[k]=new complex **[Ntheta];
		mu12[k]=new complex **[Ntheta];
		for (int th=0; th<Ntheta; th++) 
		{
		    mu11[k][th]=new complex *[NK];
		    mu12[k][th]=new complex *[NK];
		    for (int k1=0; k1<NK; k1++)
	        {
	            mu11[k][th][k1]=new complex [Ntheta];
	            mu12[k][th][k1]=new complex [Ntheta];
	            for (int th1=0; th1<Ntheta; th1++) 
		        {
		            mu11[k][th][k1][th1]=0.;
		            mu12[k][th][k1][th1]=0.;
		        ;}
	        ;}
		;}
	;}
;}


double calc_SumComplete()
{
    double sum=0.;
	{
		for (int k=0; k<NK; k++)
		for (int th=0; th<Ntheta; th++)
		{
			double k_=k*Kmax/NK, th_=th*Pi/Ntheta;      
			sum+=k_*cos(th_)*2.*Pi*k_*k_*sin(th_)*(Kmax/NK)*(Pi/Ntheta)*(norm2(alpha[k][th])+real(mu11[k][th][k][th]));                    

		;}
	;}
	return sum;
;}

double Htensor(double & sum)
{
    sum=calc_SumComplete();
    complex s=0.;
    {
    for (int k=0; k<NK; k++)
	for (int th=0; th<Ntheta; th++)
	{
	double k_=k*Kmax/NK, th_=th*Pi/Ntheta;
	
	s+=(2.*real(V(k_)*alpha[k][th])+norm2(alpha[k][th])*W2(k_, th_))
	*2.*Pi*k_*k_*sin(th_)*(Kmax/NK)*(Pi/Ntheta);
	
	for (int k1=0; k1<NK; k1++)
	for (int th1=0; th1<Ntheta; th1++)
	{
	double k1_=k1*Kmax/NK, th1_=th1*Pi/Ntheta;
    s+=1./(2.*M)*conj(alpha[k][th])*(k_*cos(th_)*k1_*cos(th1_))*
    (mu11[k][th][k1][th1]*alpha[k1][th1]+mu12[k][th][k1][th1]*conj(alpha[k1][th1]))   
    *2.*Pi*k_*k_*sin(th_)*(Kmax/NK)*(Pi/Ntheta)
    *2.*Pi*k1_*k1_*sin(th1_)*(Kmax/NK)*(Pi/Ntheta);
    s+=1./(2.*M)*alpha[k][th]*(k_*cos(th_)*k1_*cos(th1_))*
    (conj(mu12[k1][th1][k][th])*alpha[k1][th1]+mu11[k][th][k1][th1]*conj(alpha[k1][th1]))    
    *2.*Pi*k_*k_*sin(th_)*(Kmax/NK)*(Pi/Ntheta)
    *2.*Pi*k1_*k1_*sin(th1_)*(Kmax/NK)*(Pi/Ntheta);
    ;}
    ;}
    ;}
    return real(s);

;}

double Hperm(double & sum)
{
    sum=calc_SumComplete();
    complex perm=0.;
    {
    for (int k=0; k<NK; k++)
	for (int th=0; th<Ntheta; th++)
	{
	double k_=k*Kmax/NK, th_=th*Pi/Ntheta;
	
	perm+=W2(k_, th_)*mu11[k][th][k][th]
	*2.*Pi*k_*k_*sin(th_)*(Kmax/NK)*(Pi/Ntheta);
	
	for (int k1=0; k1<NK; k1++)
	for (int th1=0; th1<Ntheta; th1++)
	{
	double k1_=k1*Kmax/NK, th1_=th1*Pi/Ntheta;
    
    perm+=1./(2.*M)*(k_*cos(th_)*k1_*cos(th1_))
    *(mu11[k][th][k1][th1]*mu11[k][th][k1][th1]
    +mu12[k][th][k1][th1]*conj(mu12[k1][th1][k][th]))
    	*2.*Pi*k_*k_*sin(th_)*(Kmax/NK)*(Pi/Ntheta)
    	*2.*Pi*k1_*k1_*sin(th1_)*(Kmax/NK)*(Pi/Ntheta);

    ;}
    ;}
    ;}
    
    return  real(perm);
}


double Hcomplete(double & sum)
{
     sum=calc_SumComplete();

	double s=0;
	{
		for (int k=0; k<NK; k++)
		for (int th=0; th<Ntheta; th++)
		{
			double k_=k*Kmax/NK, th_=th*Pi/Ntheta;
			s+=(2*real(V(k_)*alpha[k][th])+norm2(alpha[k][th])*W2(k_, th_))
			*2.*Pi*k_*k_*sin(th_)*(Kmax/NK)*(Pi/Ntheta);
		;}
	;}
	return s;
;}

double dH(int k, int th, complex da)
{
	double k_=k*Kmax/NK, th_=th*Pi/Ntheta;
	double dbb=2.*real(conj(da)*alpha[k][th])+norm2(da);
	double dV=2.*Pi*k_*k_*sin(th_)*(Kmax/NK)*(Pi/Ntheta);
	double s=2*V(k_)*real(da)+dbb*W(k_,th_)+sqr(dbb*k_*cos(th_))*dV/(2.*M);
	return s*2.*Pi*k_*k_*sin(th_)*(Kmax/NK)*(Pi/Ntheta);
;}

double dH1(int k, int th, int k1, int th1, complex dmu11)
{
    
    mu11[k][th][k1][th1]+=dmu11;
    if(!(k==k1 && th==th1)) mu11[k1][th1][k][th]+=conj(dmu11);
    double E1=Htensor(SUM)+Hperm(SUM);
    mu11[k][th][k1][th1]-=dmu11;
    if(!(k==k1 && th==th1)) mu11[k1][th1][k][th]-=conj(dmu11);
    E1-=Htensor(SUM)+Hperm(SUM);
	return E1;
;}

double dH2(int k, int th, int k1, int th1, complex dmu12)
{
    
    mu12[k][th][k1][th1]+=dmu12;
    double E2=Htensor(SUM)+Hperm(SUM);
    mu12[k][th][k1][th1]-=dmu12;
    E2-=Htensor(SUM)+Hperm(SUM);
	return E2;
;}



void step_alpha()
{
	int k=rnd(NK), th=rnd(Ntheta);
	complex da(rnd()-0.5, rnd()-0.5);
	double dE=dH(k, th, da);
	if (exp(-dE/kT)>rnd()) 
	{
		alpha[k][th]+=da; E+=dE;
		double k_=k*Kmax/NK, th_=th*Pi/Ntheta;
		double s=k_*cos(th_)*(norm2(alpha[k][th])-norm2(alpha[k][th]-da));
		SUM+=s*2.*Pi*k_*k_*sin(th_)*(Kmax/NK)*(Pi/Ntheta);
//		debug<<"accepted";
	;}
//	debug<<"\n";
;}

void step_mu11()
{
	int k=rnd(NK), th=rnd(Ntheta), k1=rnd(NK), th1=rnd(Ntheta);
	complex dmu11(rnd()-0.5, rnd()-0.5);
	if(k==k1 && th==th1) dmu11=rnd()-0.5;
	double dE=dH1(k, th, k1, th1, dmu11);
	if (exp(-dE/kT)>rnd()) 
	{
		mu11[k][th][k1][th1]+=dmu11;
        if(!(k==k1 && th==th1)) mu11[k1][th1][k][th]+=conj(dmu11);
		E+=dE;
		SUM=calc_SumComplete();
//		debug<<"accepted";
	;}
//	debug<<"\n";
;}

void step_mu11_diag()
{
	int k=rnd(NK), th=rnd(Ntheta);
	complex dmu11=rnd()-0.5;
	double dE=dH1(k, th, k, th, dmu11);
	if (exp(-dE/kT)>rnd()) 
	{
		mu11[k][th][k][th]+=dmu11+conj(dmu11);
		E+=dE;
		SUM=calc_SumComplete();
//		debug<<"accepted";
	;}
//	debug<<"\n";
;}

void step_mu12()
{
	int k=rnd(NK), th=rnd(Ntheta), k1=rnd(NK), th1=rnd(Ntheta);
	complex dmu12(rnd()-0.5, rnd()-0.5);
	double dE=dH2(k, th, k1, th1, dmu12);
	if (exp(-dE/kT)>rnd()) 
	{
		mu12[k][th][k1][th1]+=dmu12;
		E+=dE;
		SUM=calc_SumComplete();
//		debug<<"accepted";
	;}
//	debug<<"\n";
;}

void sweep_alpha()
{
	for (int i=0; i<NK*Ntheta; i++) step_alpha();
;}

void sweep_mu11()
{
	for (int i=0; i<NK*NK*Ntheta*Ntheta; i++) step_mu11_diag();
;}

void sweep_mu12()
{
	for (int i=0; i<NK*Ntheta*NK*Ntheta; i++) step_mu12();
;}


int main(int argc, char **argv)
{
	if (argc>1) {q=atof(argv[1]);cout<<"q="<<q<<"\n"<<flush;}
	Ini();
//	for (kT=0.1; kT>1e-12; kT/=1.2)
	//for (int i=0; i<100; i++)
	{
		for (int i=0; i<NsweepsVector; i++) 
		{
		    sweep_alpha();
		    double hc1=Hperm(SUM)+Htensor(SUM);
		    double hc2=Hcomplete(SUM);
		    debug<<hc1<<"=="<<hc2<<"\n"<<flush;
		}
		for (int i=0; i<NsweepsMatrix; i++) {
		sweep_mu11();
			double hc1=Hperm(SUM)+Htensor(SUM);
		    double hc2=Hcomplete(SUM);
		    debug<<hc1<<"=="<<hc2<<"\n"<<flush;
		}
//		for (int i=0; i<NsweepsMatrix; i++) sweep_mu12();
		
        
        double hc1=Hperm(SUM)+Htensor(SUM);
		cout<<E<<"\t"<<hc1<<"\n"<<flush;
		{
		ofstream ou_file("alpha.dat");
		for (int k=1; k<NK; k++)
		{
			for (int th=1; th<Ntheta; th++)
			ou_file<<k*Kmax/NK<<"\t"<<th*Pi/Ntheta<<"\t"<<abs(alpha[k][th])<<"\n";
			ou_file<<"\n"<<flush;
		;}
		;}	
	;}
	
//	ofstream ou(alpha)
	return 0;
;}
