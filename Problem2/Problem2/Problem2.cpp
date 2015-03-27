# include <iostream>
# include <iomanip>

using namespace std;
const int nn = 5, mm = 4;

void mtranspose(double M[nn][mm], double MT[mm][nn], int n, int m)
{//Calculates transpose of matrix M(nxm) as MT(mxn)
	int i; int j;
	for (i = 0; i < n; i++){
		for (j = 0; j < m; j++){
			MT[j][i] = M[i][j];
		};
	};
}

void multiply_matrices(double A[mm][nn], double B[nn][mm], double C[mm][mm], \
	int n, int m, int p)
{//calculates C(nxm) = A(nxp)*B(pxm)
	int i; int j; int k;
	for (i = 0; i < n; i++){
		for (j = 0; j < m; j++){
			C[i][j] = 0;
			for (k = 0; k < p; k++){
				C[i][j] = C[i][j] + A[i][k] * B[k][j];
			};
		};
	};
}

void Pivot(double A[mm][mm], int O[], double s[], int n, int k)
{//Pivot routine for LU decomposition 
 //Use it together with the MatrixInverse subroutine
	int p; double big; double dummyB; int ii; int dummyP;
	p = k; big = abs(A[O[k]][k] / s[O[k]]);
	for (ii = k + 1; ii<n; ii++)
	{
		dummyB = abs(A[O[ii]][k] / s[O[ii]]);
		if (dummyB>big)
		{
			big = dummyB; p = ii;
		};
	};
	dummyP = O[p]; O[p] = O[k]; O[k] = dummyP;
};

void Substitute(double A[mm][mm], int O[], int n, double b[], double x[])
{//Substitution routine for LU decomposition
 //Use it together with the MatrixInverse subroutine
	double sum; int i; int j;
	for (i = 1; i<n; i++)
	{
		sum = b[O[i]];
		for (j = 0; j <= i - 1; j++)
		{
			sum = sum - A[O[i]][j] * b[O[j]];
		};
		b[O[i]] = sum;
	};
	x[n - 1] = b[O[n - 1]] / A[O[n - 1]][n - 1];
	for (i = n - 2; i >= 0; i--)
	{
		sum = 0.0;
		for (j = i + 1; j<n; j++)
		{
			sum = sum + A[O[i]][j] * x[j];
		};
		x[i] = (b[O[i]] - sum) / A[O[i]][i];
	};
};

void Decompose(double A[mm][mm], int n, double tol, int O[], double s[], int er)
{//Decomposition step in the LU decomposition algorithm
 //Use it together with the MatrixInverse subroutine
	int i; int j; int k; double factor;
	for (i = 0; i<n; i++)
	{
		O[i] = i; s[i] = abs(A[i][0]);
		for (j = 1; j<n; j++)
		{
			if (abs(A[i][j])>s[i])
			{
				s[i] = abs(A[i][j]);
			};
		};
	};
	for (k = 0; k<n - 1; k++)
	{
		Pivot(A, O, s, n, k);
		if (abs(A[O[k]][k] / s[O[k]])<tol)
		{
			er = -1; cout << "trouble1" << endl;
			break;
		};
		for (i = k + 1; i<n; i++)
		{
			factor = A[O[i]][k] / A[O[k]][k];
			A[O[i]][k] = factor;
			for (j = k + 1; j<n; j++)
			{
				A[O[i]][j] = A[O[i]][j] - factor*A[O[k]][j];
			}
		};
	};
	if (abs(A[O[k]][k] / s[O[k]]) < tol)
	{
		er = -1; cout << "trouble2" << endl;
		//er = -1; cout << endl;
	}
};

void MatrixInverse(double A[mm][mm], double AI[mm][mm], int n, double tol, int er)
{//Calculates the inverse of matrix A(nxn), i.e., AI(nxn)
 //This subroutine is the main driver for Matrix Inverse with LU decomposition algorithm
	int O[mm]; double s[mm]; double b[mm]; double x[mm];
	int i; int j;
	Decompose(A, n, tol, O, s, er);
	if (er == 0)
	{
		for (i = 0; i < n; i++)
		{
			for (j = 0; j < n; j++)
			{
				if (i == j)
				{
					b[j] = 1;
				}
				else
				{
					b[j] = 0;
				}
			}
			Substitute(A, O, n, b, x);
			for (j = 0; j < n; j++)
			{
				AI[j][i] = x[j];
			}
		}
	}
	else
	{
		cout << "Matrix is singular" << endl;
	}
}

void multiply_matrix_to_vector(double A[mm][nn], double u[], double v[], \
	int n, int m)
{//calculates v(n) = A(nxm)*u(m)
	int i; int k;
	for (i = 0; i < n; i++){
		v[i] = 0;
		for (k = 0; k < m; k++){
			v[i] = v[i] + A[i][k] * u[k];
		};
	};
}

void multiply_matrix_to_vector(double A[mm][mm], double u[], double v[], \
	int n, int m)
{//calculates v(n) = A(nxm)*u(m)
	int i; int k;
	for (i = 0; i < n; i++){
		v[i] = 0;
		for (k = 0; k < m; k++){
			v[i] = v[i] + A[i][k] * u[k];
		};
	};
}

void NLRegress(double Z[nn][mm], double y[], double a[], int n, int m){
	int er = 0;
	double tol = 0.0001;
	double ZT[mm][nn], ZTZ[mm][mm], ZTZI[mm][mm];
	double ZTY[mm];

	mtranspose(Z, ZT, nn, mm);
	multiply_matrices(ZT, Z, ZTZ, mm, mm, nn);
	MatrixInverse(ZTZ, ZTZI, mm, tol, er);
	multiply_matrix_to_vector(ZT, y, ZTY, mm, nn);
	multiply_matrix_to_vector(ZTZI, ZTY, a, mm, mm);
	
}

void BuildZP(double Z[nn][mm], double x[],int n,int m){
	for(int i=0; i<n; i++){
		for(int j=0; j<m; j++){
			Z[i][j] = pow(x[i],j);
		}
	}
}

/*void BuildZP(double Z[nn][mm], double x[],int n,int m){
	for(int i=0; i<n; i++){
		for(int j=0; j<m-1; j++){
			Z[i][j] = pow(x[i],j);
		}
		Z[i][m-1] = sqrt(x[i]);
	}
}*/

void BuildZM(double Z[nn][mm], double x[nn][mm], int n, int m){
	for (int i = 0; i < n; i++){
		Z[i][0] = 1.0;
		for (int j = 1; j < m; j++){
				Z[i][j] = x[i][j - 1];
			}
	}
}

void printvector(double u[], int n)
{//Prints out vector u of size n
	int i;
	for (i = 0; i < n; i++){
		cout << setw(10) << u[i]; cout << endl;
	};
}

void printmatrix(double A[nn][mm], int n, int m)
{//Prints out matrix A(nxm)
	int i; int j;
	for (i = 0; i < n; i++){
		for (j = 0; j < m; j++){
			cout << setw(10) << A[i][j] << " ";
		}; cout << endl;
	};
}

int main(){

	//double x[nn] = {0.2, 0.5, 0.8, 1.2, 1.7, 2, 2.3};
	//double y[nn] = {500, 700, 1000, 1200, 2200, 2650, 3750};
	//double x[nn] = {2.5, 3.5, 4.5, 5.5, 6.5, 7.5};
	//double y[nn] = {10, 9, 8, 6, 2, -2};
	double xx[nn][mm] = { { 2.0, 1.0, 4.5 }, { 3.2, 1.2, 5.5 }, { 4.0, 1.5, 2.5 }, { 3.0, 1.8, 4.0 }, { 2.5, 2.0, 5.0 }};
	double y[nn] = { 30.1, 31.2, 24.0, 32.1, 35.1 };

	double Z[nn][mm], a[mm], yf[nn];

	int m = mm - 1;

	//BuildZP(Z, x, nn, mm);
	
	BuildZM(Z, xx, nn, mm);
	printmatrix(xx, nn, mm);
	cout << endl;
	printmatrix(Z, nn, mm);
	cout << endl;
	NLRegress(Z, y, a, nn, mm);
	cout << endl;

	multiply_matrix_to_vector(Z, a, yf, nn, mm);
	for(int i = 0; i < nn; i++){
		cout<<setw(10)<<xx[i]<<setw(10)<<y[i]<<setw(10)<<yf[i]<<endl;
	}
	cout<<endl;

	printvector(a, mm);

	return 0;
}

