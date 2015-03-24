
#include <iostream>
#include <iomanip>

using namespace std;

void Regress(double x[], double y[], double n, double &a1, double &a0, double &syx, double &r2){
	double sumx = 0, sumxy = 0, st = 0, sumy = 0, sumx2 = 0, sr = 0;

	for(int i = 0; i<n; i++){
		sumx = sumx + x[i];
		sumy = sumy + y[i];
		sumxy = sumxy + x[i]*y[i];
		sumx2 = sumx2 + x[i]*x[i];
	}
	double xm = sumx/n;
	double ym = sumy/n;
	a1 = (n*sumxy - sumx*sumy)/(n*sumx2 - sumx*sumx);
	a0 = ym - a1*xm;
	for(int i = 0; i < n; i++){
		st = st + pow((y[i] - ym), 2);
		sr = sr + pow((y[i] - a1*x[i] - a0), 2);
	}
	syx = sqrt(sr/(n-2));
	r2 = (st - sr)/st;
}

double FitData(double a0, double a1, double x){
	return a0+a1*x;
}


int main(){
	int m, n;
	double *x, *y;
	double a1, a0, syx, r2;

	cout<<"Enter number of x elements (n): "; 
	cin>>n;
	cout<<"Enter number of y elements (m): ";
	cin>>m;

	if(n != m){
		cout<<"No regression possible."<<endl;
		return 0;
	}
	
	x = new double [n];
	y = new double [m]; 
	
	for(int i = 0; i < n; i++){
		cout<<"Enter x value " <<i+1<< ": ";
		cin>>x[i];
	}
	for(int i = 0; i < m; i++){
		cout<<"Enter y value " <<i+1<< ": ";
		cin>>y[i];
	}
	cout<<endl;

	Regress(x, y, n, a1, a0, syx, r2);

	cout<<"n:   "<<n<<endl;
	cout<<"a1:  "<<a1<<endl;
	cout<<"a0:  "<<a0<<endl;
	cout<<"syx: "<<syx<<endl;
	cout<<"r2:  "<<r2<<endl;
	cout<<"r:   "<<sqrt(r2)<<endl<<endl;

	for(int i = 0; i < n; i++){
		cout<<setw(10)<<x[i]<<setw(10)<<y[i]<<setw(10)<<FitData(a0, a1, x[i])<<endl;
	}
	cout<<endl;

	return 0;
}

