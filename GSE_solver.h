#ifndef _GSE_SOLVER_H_
#define _GSE_SOLVER_H_

/*
	GSE_KRYLOV定义及一些内部函数

*/
enum {	ARE_RET_NORMAL=0x0,ARE_AI_CONVERGE=10	};

#define ARE_trace printf

/*
	将dimH改为ldh

	注意：
		1、Q,H采用列优先存储
*/
typedef struct {
	int info,nRestart,MAX_RESTART,checkOffBase,init_V0;
	int method,nEV,nGV,nRV,ldh,nCurEV,V_m,V_p,ldq;
	int dim,nnz,*ptr,*ind,*A_ptr,*A_ind,*i_temp;		
	int A_type,val_type,sft_mode,ivt_mode,stop_mode,restart_mode;
	int balance_mode,adapt_Vp;		//internal parameter

	GRUS_FLOATING *V,*AV,*H,*Q,*rv_A,*ev_H,*work,shift,delta,sft_old;
	double *d_temp,*D,a_norm,tole,ortho_thresh;
	GRUS_FLOATING *w_H,*Base;

	char enviro_str[200];

/*#ifdef ZINT
#else
	double *wr_H,*wi_H;
#endif*/
	void *val,*A;
	void *hLU;
}ARE_R_SOLVER;		//arnoldi eigenvalue solver based on rational transformation
typedef ARE_R_SOLVER GSE_KRYLOV; 
typedef ARE_R_SOLVER GSE_RATIOK; 

//************	GSE_KRYLOV的标志变量	***********//
enum	{	
	GSE_KRYLOV_REPOS=0x0,GSE_KRYLOV_ONAV=0x10000	
};



double Eigen_Defect( ARE_R_SOLVER *hAR_R,GRUS_FLOATING lenda,GRUS_FLOATING *ev,GRUS_FLOATING *r );
void GSE_Krylov_update_( GSE_KRYLOV* hSolver,int j,int flag );
void GSE_Krylov_expand_( GSE_KRYLOV* hSolver,int p,int m  );
double GSE_Krylov_Resid( GSE_KRYLOV *hSolver,int ev_0,int ev_1,GRUS_FLOATING *res,int flag );


extern double mch_eps;
#endif