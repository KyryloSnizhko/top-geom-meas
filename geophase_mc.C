#include <cmath>
#include <iostream>
#include <fstream>

#define pi 3.1415926535898

using namespace std;

double ran2(long *);
double phase(double cr,double ci);

long idum;

int main(void)
{
	double psiur,psiui,psidr,psidi,sqrteta,sqrtmeta,msqrtmeta,cosphi,sinphi,pup,eta,norm,geophase,geophaseold;
	double theta,phi,psi2ur,psi2ui,psi2dr,psi2di,c,cost,sint,phaser,phasei,averagephaser,averagephasei,puplast;
	double cstart,cend,dc,thetastart,thetaend,dtheta;
	int n,nmax,nn,nnmax,protocol;
	char name[100];
				
	cout << "Number of measurements:";cin >> nmax;
	cout << "Number of realizations:"; cin >> nnmax;
	cout << "Vary parameter c starting from ... until ... with stepsize ...:"; cin >> cstart; cin >> cend; cin >> dc;
	cout << "Vary parameter theta/pi starting from ... until ... with stepsize ...:"; cin >> thetastart; cin >> thetaend; cin >> dtheta;
	cout << "Seed for random number generator (must be a negative integer number):"; cin >> idum;
	cout << "Protocol 1 ( < exp(i phi) > ) or protocol 2 ( < exp(2 i phi) >):"; cin >> protocol;
	cout << "Filename for output data:"; cin >> name;
	ofstream Data(name,ios::out);
	
	geophaseold=0.0;
	for (theta=thetastart;theta<=thetaend+0.5*dtheta;theta+=dtheta) for (c=cstart;c<=cend+0.5*dc;c+=dc) 
	{
		eta=4.0*c/(double)(nmax);
		sqrteta=sqrt(eta);
		sqrtmeta=sqrt(1.0-eta);
		msqrtmeta=1.0-sqrtmeta;
		cost=cos(theta*pi/2.0);
		sint=sin(theta*pi/2.0);
		averagephaser=0.0;
		averagephasei=0.0;
		puplast=0.0;
		
		for (nn=0;nn<nnmax;nn++)
		{
	
			psiur=cost;// Initial state |\psi> = (psiur + I psiui) |u> + (psidr + I psidi) |d>
			psiui=0.0;
			psidr=sint;
			psidi=0.0;
	
			for (n=1;n<nmax;n++)
			{
				phi=2.0*pi*double(n)/(double)(nmax);
				cosphi=cos(phi);
				sinphi=sin(phi);
		
				psi2ur=(cost*cost*msqrtmeta+sqrtmeta)*psiur+cost*sint*msqrtmeta*(cosphi*psidr+sinphi*psidi); // New state after detection "up"
				psi2ui=(cost*cost*msqrtmeta+sqrtmeta)*psiui+cost*sint*msqrtmeta*(cosphi*psidi-sinphi*psidr);
				psi2dr=cost*sint*msqrtmeta*(cosphi*psiur-sinphi*psiui)+(1.0-cost*cost*msqrtmeta)*psidr;
				psi2di=cost*sint*msqrtmeta*(cosphi*psiui+sinphi*psiur)+(1.0-cost*cost*msqrtmeta)*psidi;
		
				pup=psi2ur*psi2ur+psi2ui*psi2ui+psi2dr*psi2dr+psi2di*psi2di; // Probability for detection "up"
		
				if (ran2(&idum)<pup) // Detection "up"
				{
					psiur=psi2ur/sqrt(pup);
					psiui=psi2ui/sqrt(pup);
					psidr=psi2dr/sqrt(pup);
					psidi=psi2di/sqrt(pup);
				}
				else // Detection "down"
				{
					psi2ur=(sint*sint*psiur-cost*sint*(cosphi*psidr+sinphi*psidi))*sqrteta;
					psi2ui=(sint*sint*psiui-cost*sint*(cosphi*psidi-sinphi*psidr))*sqrteta;
					psi2dr=(-cost*sint*(cosphi*psiur-sinphi*psiui)+cost*cost*psidr)*sqrteta;
					psi2di=(-cost*sint*(cosphi*psiui+sinphi*psiur)+cost*cost*psidi)*sqrteta;
			
					psiur=psi2ur/sqrt(1.0-pup);
					psiui=psi2ui/sqrt(1.0-pup);
					psidr=psi2dr/sqrt(1.0-pup);
					psidi=psi2di/sqrt(1.0-pup);
				}		
			}
		
			phaser=cost*psiur+sint*psidr;
			phasei=cost*psiui+sint*psidi;
	
			norm=phaser*phaser+phasei*phasei; // Probability for last (strong) measurement "up"
	
			phaser/=sqrt(norm);
			phasei/=sqrt(norm);	

			if (protocol==1)
			{
				averagephaser+=norm*phaser;
				averagephasei+=norm*phasei;
			}
			else
			{
				averagephaser+=norm*(phaser*phaser-phasei*phasei);
				averagephasei+=norm*2.0*phaser*phasei;
			}
			
			puplast+=norm;
		}
		averagephaser/=puplast;
		averagephasei/=puplast;
				
		geophase=phase(averagephaser,averagephasei);
		
		while (geophase-geophaseold>pi) geophase-=2.0*pi; // adjust phase to make curve continous
		while (geophase-geophaseold<-pi) geophase+=2.0*pi;
		geophaseold=geophase;
		
		if (protocol==2) 
		{
			geophase/=2.0;
			cout << "c=" << c << ", theta=" << theta << ", <exp(2 i phi)> = " << averagephaser << " + I " << averagephasei << ", phi = " << geophase/pi << " pi, e^(-alpha) = " << sqrt(averagephaser*averagephaser+averagephasei*averagephasei) << endl;
		}
		else cout << "c=" << c << ", theta=" << theta << ", <exp(i phi)> = " << averagephaser << " + I " << averagephasei << ", phi = " << geophase/pi << " pi, e^(-alpha) = " << sqrt(averagephaser*averagephaser+averagephasei*averagephasei) << endl;
			
		if ((cend-cstart>dc) && (thetaend-thetastart>dtheta)) Data << theta << " " << c << " " << geophase/pi  << " " << sqrt(averagephaser*averagephaser+averagephasei*averagephasei) << endl;
		if ((cend-cstart>dc) && (thetaend-thetastart<dtheta)) Data << c << " " << geophase/pi  << " " << sqrt(averagephaser*averagephaser+averagephasei*averagephasei) << endl;
		if ((cend-cstart<dc) && (thetaend-thetastart>dtheta)) Data << theta << " " << geophase/pi  << " " << sqrt(averagephaser*averagephaser+averagephasei*averagephasei) << endl;
		
	}
}  

/* note #undef's at end of file */
#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

double ran2(long *idum)
{
	int j;
	long k;
	static long idum2=123456789;
	static long iy=0;
	static long iv[NTAB];
	double temp;

	if (*idum <= 0) {
		if (-(*idum) < 1) *idum=1;
		else *idum = -(*idum);
		idum2=(*idum);
		for (j=NTAB+7;j>=0;j--) {
			k=(*idum)/IQ1;
			*idum=IA1*(*idum-k*IQ1)-k*IR1;
			if (*idum < 0) *idum += IM1;
			if (j < NTAB) iv[j] = *idum;
		}
		iy=iv[0];
	}
	k=(*idum)/IQ1;
	*idum=IA1*(*idum-k*IQ1)-k*IR1;
	if (*idum < 0) *idum += IM1;
	k=idum2/IQ2;
	idum2=IA2*(idum2-k*IQ2)-k*IR2;
	if (idum2 < 0) idum2 += IM2;
	j=iy/NDIV;
	iy=iv[j]-idum2;
	iv[j] = *idum;
	if (iy < 1) iy += IMM1;
	if ((temp=AM*iy) > RNMX) return RNMX;
	else return temp;
}
#undef IM1
#undef IM2
#undef AM
#undef IMM1
#undef IA1
#undef IA2
#undef IQ1
#undef IQ2
#undef IR1
#undef IR2
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX

double phase(double z1r,double z1i)
{
	double phi;
	if (fabs(z1r)>1.0e-10) 
	{
		phi=atan(z1i/z1r);
		if (z1r<0)
		{
			phi+=pi;
			if (phi>pi) phi-=2.0*pi;
		}	
	}
	else 
	{
		if (z1i>=0.0) phi=0.5*pi;
		else phi=-0.5*pi;
	}
	return(phi);	
}