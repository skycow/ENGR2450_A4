#include <iostream>
#include <iomanip>

using namespace std;

const int nn = 4, mm = 10;

void NewtInt(double x[], double y[], int n, double xi, double yint[], double ea[]){
	double fdd[nn][nn];
	for(int i = 0; i < n; i++){
		fdd[i][0] = y[i];
	}
	for(int j = 1; j < n; j++){
		for(int i = 0; i < n-j; i++){
			fdd[i][j] = (fdd[i+1][j-1] - fdd[i][j-1])/(x[i+j] - x[i]);
		}
	}
	double xterm = 1;
	yint[0] = fdd[0][0];
	for(int order = 1; order < n; order++){
		xterm = xterm * (xi - x[order-1]);
		double yint2 = yint[order-1] + fdd[0][order] * xterm;
		ea[order-1] = yint2 - yint[order-1];
		yint[order] = yint2;
	}
}

void printvector(double u[], int n)
{//Prints out vector u of size n
	int i;
	for (i = 0; i < n; i++){
		cout << setw(10) << u[i]; cout << endl;
	};
}

int main(){

	//double x[mm] = {1, 2, 3, 5, 7, 8};
	//double y[mm] = {3, 6, 19, 99, 291, 444};
	//double ea[mm], xi = 4.0, yint[mm];
	//double x[5] = {2.5,3.2,3.9,4.2,4.5,5.2,6.3,7.1,8.2};
	//double y[5] = {34,70,130,170,205,331,618,900,1400};
	double x[5] = {1, 2, 3, 5, 7};
	double y[5] = {3, 6, 19, 99, 291};
	double ea[4], xi = 4.0, yint[4];
	int n = 4;

	NewtInt(x, y, n, xi, yint, ea); 

	printvector(yint, n);
	cout<<endl;
	printvector(ea, n);

	return 0;
}