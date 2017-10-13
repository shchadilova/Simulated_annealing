#include <iostream>
#include <fstream>
#include <complex>
#include <cstdlib>
//#include <mpi.h>

//int MPI_size, MPI_rank;

#define complex complex<double>
using namespace std;

#define do_once  {static int flag_once=0; if (flag_once==0) {flag_once=1;
#define end_do_once ;};}

double sinh(double x)
{
return 0.5*(exp(x)-exp(-x)); 
}

double sinhc(double x)
{
    double out=1.;
    if(abs(x)>1e-8) out=0.5*(exp(x)-exp(-x))/x;
    return out; 
}

int sgn(double x) {if (x>0) return 1; if (x<0) return -1; return 0;}
double sqr(double x) {return x*x;}
complex sqr(const complex & x) {return x*x;}

complex I(0,1); double Pi=4.*atan(1.);
inline double norm2(complex x) {double a=real(x), b=imag(x); return a*a+b*b;}


int initial_random=time(NULL);
int INT_RANDOM=initial_random;//  //from 1 to 2^31-1

int int_rnd()
{
//for (int j=1; j<3; j++)
{
	 int k=INT_RANDOM/127773;
	 INT_RANDOM=16807*(INT_RANDOM-k*127773)-2836*k;
	 if (INT_RANDOM<0) INT_RANDOM+=2147483647;
;}
	 return INT_RANDOM;
;}

int INT_RANDOM_AUX_VARIABLE=int_rnd()+int_rnd()+int_rnd()+int_rnd();  //randomising...

double rnd ()
{
	return int_rnd()/2147483647.0
;}

int rnd (int k)
{
	int d=2147483647%k, d1=2147483647-d, r=int_rnd();
	if (r>=d1) return rnd(k);
	return r/(d1/k)
;}

void swing(int & x, int & y) {int a=x; x=y; y=a;}
void swing(double & x, double & y) {double a=x; x=y; y=a;}

int min(int & x, int &y) {if (x<y) return x; else return y;}
int max(int & x, int &y) {if (x>y) return x; else return y;}

double max_abs(const double & x, const double &y) {if (abs(x)>abs(y)) return x; else return y;}


double * file_to_vector(const char* filename, int& n, int use_old_memory=0)
{
	static double * table=NULL; static int n_=0;
	if (use_old_memory==1) delete [] table;
	ifstream inp(filename); n_=0;
	double d; while(!inp.eof()) {inp>>d; n_++;};
	inp.close();
	inp.clear();
	table=new double [n_];
	inp.open(filename); 
	{for (int i=0; i<n_; i++) inp>>table[i];}
	n=n_;
	//if (MPI_rank==0) 
	cout<<filename<<" imported";
	return table;
;}

double ** file_to_matrix(const char* filename, int& n1, int& n2, int use_old_memory=0)//n1:number of rows, n2: number of columns
{
	static double ** table=NULL; static int n1_=0, n2_=0;
	if (use_old_memory==1) {for (int i=0; i<n1_; i++) delete [] table[i]; delete [] table;} 
	
	int N=0; char buf[4096];
	//{
	ifstream inp(filename); 
	double d; while(!inp.eof()) {inp>>buf;N++;}//{inp>>d; N++;};
	//;}
	
	inp.clear();inp.seekg(0, ios::beg);
	
	//{	ifstream inp(filename);
	n1_=0; 

	while(!inp.eof()) {inp.getline(buf, 4096); n1_++;}// {inp.getline(buf, 255); n1_++;}
	if (char_traits<char>::length(buf)<2) {n1_--; N--;}
	//;}
	
	inp.clear();inp.seekg(0, ios::beg);

	//ifstream inp(filename); 	
	n2_=N/n1_;
//	cout<<filename<<":  "<<n1_<<"  "<<n2_<<"  "<<N<<"\n"<<flush;

	{table=new double * [n1_]; for (int i=0; i<n1_; i++) table[i]=new double [n2_];}
	{for (int j=0; j<n1_; j++) for (int i=0; i<n2_; i++)  inp>>table[j][i];}
//	{for (int j=0; j<6; j++) {for (int i=0; i<n1_; i++)  cout<<table[i][j]<<"\t"; cout<<"\n";}}

	
	//if (MPI_rank==0) 
	cout<<filename<<" imported\n";
	n1=n1_; n2=n2_;  return table;
;}  
