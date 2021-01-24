#ifndef _GSE_H_
#define _GSE_H_

#ifdef _GRUS_STATIC_
	#define GRUS_GSE_DLL 
#else
	#ifndef GRUS_GSE_DLL
	#define GRUS_GSE_DLL __declspec(dllimport)
	#endif
#endif
/*
#define PI 3.1415926535897932384626433832795

//错误代码
#define GSE_OK					0x0
#define GSE_LOST_ORTHOGONAL		-10
#define GSE_FAIL_CLUSTER		-11
#define GSE_MISS_ENGENVALUE		-12
#define GSE_INVERT_FAIL			-20
#define GSE_FOPEN_ERR			-100

#define GSE_JDQR_CONVERGE	0x100
#define GSE_JDQR_LOOP		0x101
//参数
#define GSE_METHOD					0x0
#define GSE_PARAM_NEV				0x2
#define GSE_PARAM_NGV				0x3
#define GSE_DATA_TYPE				0x4
#define GSE_MAT_KIND				0x5
#define GSE_PARAM_MAX_RESTART		0x6
#define GSE_INVERT_MODE				0x7
#define GSE_SHIFT_MODE				0x8
#define GSE_STOP_MODE				0x9
#define GSE_RESTART_MODE			0xa
//缺省值30		常用的值[45,30,5,0]	
#define GSE_REORTHO_THRESH			0xb			
//缺省值0.05		常用的值[0.01,0.05,0.1]	
#define GSE_CLUSTER_THRESH			0xc
//1.0		常用的值[10,1,0.1,0.0]
#define GSE_TOLERANCE_MULTIPLE		0xd
#define GSE_PARAM_MOST				0x15

enum {	GSE_METHOD_ARE=0x1,GSE_METHOD_KSCHUR,GSE_METHOD_JDQR	};
enum {	GSE_DATA_DOUBLE=0x1,GSE_DATA_COMPLEX	};
enum {	GSE_MAT_HESSENBERG=0x1,GSE_MAT_HERMITIAN,GSE_MAT_GENERAL };
enum {	GSE_INVERT_NO=0x0,GSE_INVERT_GSS=0x01,GSE_INVERT_MF_1=0x02	};		//invert mode
enum {	GSE_STOP_DIRECT=0x0,GSE_STOP_KRYLOV=0x01,GSE_STOP_BLOCK=0x02	};		//invert mode
enum {	GSE_SHIFT_NO=0x0,GSE_SHIFT_MA=0x01,GSE_SHIFT_NORMAL=0x02	};		//invert mode
enum {	GSE_RESTART_CHEBYCHEV_0=-1,GSE_RESTART_JACOBI=-2,GSE_RESTART_SORENSEN=-4,GSE_RESTART_RANDOM=-13,					//restart mode
		GSE_RESTART_V1=0,
		GSE_RESTART_MOGAN=1	};		
#define GRUS_SE_INFO_ITEM 50
extern GRUS_GSE_DLL double GRUS_SE_INFO[GRUS_SE_INFO_ITEM];
//输出信息 
#define GRUS_SE_STATUS			0	
#define GRUS_SE_NROW			1		
#define GRUS_SE_NCOL			2		
#define GRUS_SE_ENTRY			3		
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
#define GRUS_SE_X_TIME			29	
#define GRUS_SE_RESTART_RATIO	31

#define GRUS_SE_PEAK_MEMORY		41
#define GRUS_SE_USERSPACE		42		
*/

#ifdef __cplusplus
extern "C" {
#endif
	GRUS_GSE_DLL void *GSE_Krylov_Init( int dim,int *ptr,int *ind,void *val,int *param,double *shift );
	GRUS_GSE_DLL void GSE_Krylov_Clear( void* hSolver );
	GRUS_GSE_DLL int GSE_Krylov_Post( void* hSolver );

	GRUS_GSE_DLL int ARE_main( void* hSolver,double *w );
	GRUS_GSE_DLL int KSchur_core( void* hSolver,double *w );
	GRUS_GSE_DLL int JDQR_engine( void* hSolver );

	GRUS_GSE_DLL int ARE_get_engine( void* hSolver );
	GRUS_GSE_DLL void ARE_dump( int );

	GRUS_GSE_DLL void GSE_QR_D( int nCol,int *ptr,int *ind,double *val,int flag,double sft[2] );
	GRUS_GSE_DLL void GSE_QR_Z( int nCol,int *ptr,int *ind,DoubleComplex *val,int flag,DoubleComplex sft );

#ifdef __cplusplus
  }
#endif

#endif