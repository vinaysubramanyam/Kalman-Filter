#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <cstring>
#include <fstream>
#include <string>
#include <cstdint>
#include <sstream>
#include <vector>
#include <iterator>
using namespace std;


inline double closed_interval_rand(double x0, double x1)
{
	return x0 + (x1 - x0) * rand() / ((double) RAND_MAX);
}

int main(){
	
	int StateDim=2;
	int ObsDim=1;
	int N=100;

	float Var_PNoise=0.1;
	float Mu_PNoise=0;
	float Std_PNoise=pow(Var_PNoise,0.5);

	float Var_ONoise=2;
	float Mu_ONoise=0;
	float Std_ONoise=pow(Var_ONoise,0.5);
	
	vector<vector<float>> A;
	A.assign(2,vector<float>(2,0));
	vector<vector<float>> A1;
	vector<vector<float>> A2;
	vector<vector<float>> PNoise;
	vector<vector<float>> B1;
	vector<vector<float>> B2;
	vector<vector<float>> ONoise;
	vector<vector<float>> X;
	vector<vector<float>> C;
	vector<vector<float>> y;


	X.assign(StateDim,vector<float>(N,0));
	C.assign(ObsDim,vector<float>(StateDim,0));
	y.assign(ObsDim,vector<float>(N,0));
	C[0][0]=1;
	A[0][0]=1.9223;
	A[0][1]=-0.9604;
	A[1][0]=1.0000;
	A[1][1]=0;
	//X[0][0]=1;
	
	A1.assign(StateDim,vector<float>(N,1));
	A2.assign(StateDim,vector<float>(N,1));
	PNoise.assign(StateDim,vector<float>(N,0));

	B1.assign(ObsDim,vector<float>(N,1));
	B2.assign(ObsDim,vector<float>(N,1));
	ONoise.assign(ObsDim,vector<float>(N,0));


	/* PNoise */
	for (int i = 0; i < A1.size(); i++)
	{
	    for (int j = 0; j < A1[i].size(); j++)
	    {
	    	
			float r=closed_interval_rand(0,1);
			A1[i][j]=r*Std_PNoise;
	        
	    }
	}

	for (int i = 0; i < A2.size(); i++)
	{
	    for (int j = 0; j < A2[i].size(); j++)
	    {
	    	
			A2[i][j]*=Mu_PNoise;
	       
	    }
	}

	for (int i = 0; i < PNoise.size(); i++)
	{
	    for (int j = 0; j < PNoise[i].size(); j++)
	    {
	    	
			PNoise[i][j]=A1[i][j]+A2[i][j];
	        
	    }
	}

	/* ONoise */
	for (int i = 0; i < B1.size(); i++)
	{
	    for (int j = 0; j < B1[i].size(); j++)
	    {
	    	
			float r=closed_interval_rand(0,1);
			B1[i][j]=r*Std_ONoise;
	        
	    }
	}

	for (int i = 0; i < B2.size(); i++)
	{
	    for (int j = 0; j < B2[i].size(); j++)
	    {
	    	
			B2[i][j]*=Mu_PNoise;
	       
	    }
	}

	for (int i = 0; i < ONoise.size(); i++)
	{
	    for (int j = 0; j < ONoise[i].size(); j++)
	    {
	    	
			ONoise[i][j]=B1[i][j]+B2[i][j];
			//printf("%6.3f\n",ONoise[i][j]);
	        
	    }
	}


	vector<vector<float>> PNoise_transpose(PNoise[0].size(),vector<float>(PNoise.size()));
	for (int i=0;i<PNoise.size();i++){
		for (int j=0;j<PNoise[0].size();j++){
			PNoise_transpose[j][i]=PNoise[i][j];
		}
	}

	// Covariance calculation
	vector<vector<float>> PNoise_mul(PNoise_transpose.size(),vector<float>(PNoise_transpose[0].size(),0));

	float PNoise_transpose_col1sum=0;
	float PNoise_transpose_col2sum=0;

	// 100*2 matrix
	for (int i=0;i<PNoise_transpose.size();i++){
		for (int j=0;j<PNoise_transpose[0].size();j++){
			PNoise_transpose_col1sum+=PNoise_transpose[i][0];
			PNoise_transpose_col2sum+=PNoise_transpose[i][1];

		}
	}

	PNoise_transpose_col1sum/=N;
	PNoise_transpose_col2sum/=N;
	//printf("%f\n%f",PNoise_transpose_col1sum,PNoise_transpose_col2sum);
	
	for (int i=0;i<PNoise_transpose.size();i++){
		for (int j=0;j<PNoise_transpose[0].size();j++){
			PNoise_mul[i][0]=PNoise_transpose_col1sum;
			PNoise_mul[i][1]=PNoise_transpose_col2sum;
		}
	}

	//Difference
	vector<vector<float>> PNoise_mul_aa(PNoise_transpose.size(),vector<float>(PNoise_transpose[0].size(),0));
	for (int i=0;i<PNoise_transpose.size();i++){
		for (int j=0;j<PNoise_transpose[0].size();j++){
			PNoise_mul_aa[i][j]=PNoise_transpose[i][j]-PNoise_mul[i][j];
			//printf("%f\n",PNoise_mul_aa[i][j]);
		}
	}

	vector<vector<float>> PNoise_mul_aa_trans(PNoise_mul_aa[0].size(),vector<float>(PNoise_mul_aa.size(),0));
	for (int i=0;i<PNoise_mul_aa.size();i++){
		for (int j=0;j<PNoise_mul_aa[0].size();j++){
			PNoise_mul_aa_trans[j][i]=PNoise_mul_aa[i][j];
		}
	}
	
	vector<vector<float>> Q(PNoise_mul_aa_trans.size(),vector<float>(PNoise_mul_aa_trans.size(),0));
	
	// A*A` calculation
	for (int i=0;i<PNoise_mul_aa_trans.size();i++){
		for (int j=0;j<PNoise_mul_aa_trans.size();j++){
			Q[i][j]=0;
			for (int k=0;k<PNoise_mul_aa.size();k++){
				Q[i][j]+=PNoise_mul_aa_trans[i][k]*PNoise_mul_aa[k][j];
				//printf("%f\n",Q[i][j]);
			}
		}
	}

	
	/* Important */
	for (int i=0;i<Q.size();i++){
		for (int j=0;j<Q[0].size();j++){
			Q[i][j]/=N; 
			//printf("%f\n",Q[i][j]);
			
		}
	}


////////////////////////////////////////////
	// ONoise covariance
	float ONoise_colsum=0;
	vector<vector<float>> ONoise_transpose(ONoise[0].size(),vector<float>(ONoise.size()));
	for (int i=0;i<ONoise.size();i++){
		for (int j=0;j<ONoise[0].size();j++){
			ONoise_colsum+=ONoise[i][j];
			ONoise_transpose[j][i]=ONoise[i][j];
		}
	}

	ONoise_colsum/=N;
	vector<vector<float>> ONoise_a(ONoise_transpose.size(),vector<float>(ONoise_transpose[0].size()));
	for (int i=0;i<ONoise_a.size();i++){
		for (int j=0;j<ONoise_a[0].size();j++){
			ONoise_a[i][j]=ONoise_transpose[i][j]-ONoise_colsum;
		}
	}

	/* Important */	
	float R=0;
	vector<vector<float>> ONoise_cov(ONoise_transpose[0].size(),vector<float>(ONoise_transpose[0].size()));
	for (int i=0;i<ONoise_cov.size();i++){
		for (int j=0;j<ONoise_cov[0].size();j++){
			ONoise_cov[i][j]=(ONoise_a[i][j]*ONoise_a[i][j])/N;
			R=ONoise_cov[i][j]; // ONoise Covariance
		}
	}

	/* Initial values for model */
	X[0][0]=1;
	X[1][0]=0;

	/* Initial state */
	y[0][0]=C[0][0]*X[0][0]+C[0][1]*X[1][0]+ONoise[0][0];


	/* Simulating states and real observations */
	for (int state=1;state<N;state++){

		X[0][state]=A[0][0]*X[0][state-1]+A[0][1]*X[1][state-1]+PNoise[0][state];
		X[1][state]=A[1][0]*X[0][state-1]+A[1][1]*X[1][state-1]+PNoise[1][state];
		y[0][state]=C[0][0]*X[0][state]+C[0][1]*X[1][state]+ONoise[0][state];

	}

	// for (int i=0;i<y.size();i++){
	// 	for (int j=0;j<y[0].size();j++){
	// 		printf("%f\t",y[i][j] );
	// 	}
	// }
	vector<vector<float>> Px(StateDim,vector<float>(StateDim,0));
	Px[0][0]=1;
	Px[0][1]=0;
	Px[1][0]=0;
	Px[1][1]=1;

	vector<vector<float>> Px_(StateDim,vector<float>(StateDim,0));

	vector<vector<float>> xh(StateDim,vector<float>(N,0));
	for (int i=0;i<xh.size();i++){
		for (int j=0;j<xh[0].size();j++){
			xh[i][j]=0.01*closed_interval_rand(1,-1);
			//printf("%f\t",xh[i][j]);
		}
	}
	vector<vector<float>> xh_(StateDim,vector<float>(N,0));
	vector<vector<float>> yh_(ObsDim,vector<float>(N,0));
	vector<vector<float>> inov(ObsDim,vector<float>(N,0));
	vector<vector<float>> K(StateDim,vector<float>(ObsDim,0));
	// R=2.2891;
	// Q[0][0]=0.10542;
	// Q[0][1]=0.0090144;
	// Q[1][0]=0.0090144;
	// Q[1][1]=0.1123068;
	for (int i=0;i<N;i++){

		/* priori estimate */
		xh_[0][i]=A[0][0]*xh[0][i]+A[0][1]*xh[1][i];
		xh_[1][i]=A[1][0]*xh[0][i]+A[1][1]*xh[1][i];

		/* priori estimate of state covariance matrix */
		Px_[0][0]=((1.9223)*(((1.9223)*(Px[0][0]))+((-0.9604)*(Px[1][0]))))+((-0.9604)*(((1.9223)*(Px[0][1]))+((-0.9604)*(Px[1][1])))) + Q[0][0];
		Px_[0][1]=(1*(((1.9223)*(Px[0][0]))+((-0.9604)*(Px[1][0]))))+((0)*(((1.9223)*(Px[0][1]))+((-0.9604)*(Px[1][1])))) + Q[0][1];
		Px_[1][0]=((1.9223)*(((1)*(Px[0][0]))+((0)*Px[1][0])))+((-0.9604)*(((1)*Px[0][1])+((0)*Px[1][1]))) + Q[1][0];
		Px_[1][1]=(1*(((1)*(Px[0][0]))+((0)*Px[1][0])))+((0)*(((1)*(Px[0][1]))+((0)*Px[1][1]))) + Q[1][1];

		/* Kalman filter coefficient */
		K[0][0]=Px_[0][0]*(1/((Px_[0][0])+R));
		K[1][0]=Px_[1][0]*(1/((Px_[0][0])+R));

		/* Estimated observation */
		yh_[0][i]=xh_[0][i]+R;

		/* Innovation error */
		inov[0][i]=y[0][i]-yh_[0][i];
		//printf("%f\t",inov[0][i]);

		/* Updated estimate of current state */
		xh[0][i+1]=xh_[0][i]+K[0][0]*inov[0][i];
		xh[1][i+1]=xh_[1][i]+K[1][0]*inov[0][i];

		/* Updated state covariance matrix */
		Px[0][0]=Px_[0][0]-K[0][0]*Px_[0][0];
		Px[0][1]=Px_[0][1]-K[0][0]*Px_[0][1];
		Px[1][0]=Px_[1][0]-K[1][0]*Px_[0][0];
		Px[1][1]=Px_[1][1]-K[1][0]*Px_[0][1];

		//printf("y:%f\n",y[0][i]);
		//printf("yh:%f\n",yh_[0][i]);
		printf("xh:%f\t%f\n",xh[0][i],xh[1][i]);
		cout << "#######" << endl;
		printf("xh_:%f\t%f\n",xh_[0][i],xh_[1][i]);
		
	}

	return 0;
}
