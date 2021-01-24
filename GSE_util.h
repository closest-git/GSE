#ifndef _GSE_UTIL_H_
#define _GSE_UTIL_H_

enum{
	GSE_OUTPUT_SHORT=0x10
};

#ifdef __cplusplus
extern "C" {
#endif
	int GSE_Balance_1( int dim,int loop,double *D,double *d_temp,GRUS_FLOATING *work,int flag );

	int GSE_Select( int m,DoubleComplex *w,int K,int *select,int *order,DoubleComplex off,double thresh,int flag );
	int GSE_MS_sort( int n,int K,int t,GRUS_FLOATING *T,int ldt,GRUS_FLOATING *Q,int ldq,GRUS_FLOATING *w,int *select,int flag );
	int GSE_MG_evs( int n,int t,GRUS_FLOATING *M,int ldm,GRUS_FLOATING * w,GRUS_FLOATING *Q,int ldq,GRUS_FLOATING *ev,int ldev,GRUS_FLOATING * work,double* rwork,int *flag );
	void	GSE_random_1( int dim,GRUS_FLOATING *V );

	double Optimal_Ellipse_0( int m_0,int m_1,GRUS_FLOATING *w_H,GRUS_FLOATING *e,GRUS_FLOATING *c );

	//operation on V
	int ZV_load( int dim,DoubleComplex *z,char* sPath )	;
	void ZV_order_x( int n,DoubleComplex *w,DoubleComplex x,int flag );
//	void ZV_select( int m,DoubleComplex *w,int K,int *select,int *temp,int flag );

	//verify function
	int ZU_verify( int m,int n,DoubleComplex* mat,int ldm,int flag );
	int ZV_select_verify( int m,DoubleComplex *w,int K,int flag );

#ifdef __cplusplus
  }
#endif

#define GRUS_SE_INFO_ITEM 50
extern double GRUS_SE_INFO[GRUS_SE_INFO_ITEM];

// ‰≥ˆ–≈œ¢ 
#define GRUS_SE_STATUS			0	
#define GRUS_SE_NROW			1		
#define GRUS_SE_NCOL			2		
#define GRUS_SE_ENTRY			3		/*entries in A */
#define GRUS_SE_SYMMETRY		4		
#define GRUS_SE_COND_1			10		
#define GRUS_SE_COND_F			11		

#define GRUS_SE_MAX_OFF_ORTHONORMAL			16
#define GRUS_SE_MAX_INVERT_RES				17
#define GRUS_SE_MAX_OFF_H					18
#define GRUS_SE_MAX_OFF_KRYLOV				19

#define GRUS_SE_TOTAL_TIME		20
#define GRUS_SE_MV_TIME			21
#define GRUS_SE_MV_COUNT		22
#define GRUS_SE_ARNORDI_TIME	23		    
#define GRUS_SE_INVERT_TIME		24	
#define GRUS_SE_UPDATE_TIME		25	
#define GRUS_SE_VERIFY_TIME		26	
#define GRUS_SE_RESTART_TIME	27	
#define GRUS_SE_RESTART_COUNT	28	
#define GRUS_SE_REORTHO_TIME	29	
#define GRUS_SE_RESTART_RATIO	31

#define GRUS_SE_PEAK_MEMORY		41
#define GRUS_SE_USERSPACE		42

#endif
