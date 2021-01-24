#ifndef _GSE_H_
#define _GSE_H_

#ifdef _GRUS_STATIC_
	#define GRUS_GSE_DLL 
#else
	#ifndef GRUS_GSE_DLL
	#define GRUS_GSE_DLL __declspec(dllimport)
	#endif
#endif

#include "GSE_def.h"
#include "GSE_solver.h"

#ifdef __cplusplus
extern "C" {
#endif
	GRUS_GSE_DLL void GSE_driver_core( void* hSolver,void **u_0,void **u_1,int *code,int flag  );

	GRUS_GSE_DLL void *GSE_Krylov_Init( int dim,double *param,double *shift );
	GRUS_GSE_DLL void *GSE_Krylov_Init_ccs( int dim,int *ptr,int *ind,void *val,double *param,double *shift );
	GRUS_GSE_DLL void GSE_Krylov_Clear( void* hSolver );
	GRUS_GSE_DLL int GSE_Krylov_Post( void* hSolver,GRUS_FLOATING *w,GRUS_FLOATING *z,int *nEV,int *flag );
	GRUS_GSE_DLL void GSE_Krylov_Shift( void *hSolver,double *shift,int flag );

	GRUS_GSE_DLL void GSE_Krylov_Dump( void* hSolver,int );
	int GSE_Krylov_Balance( void *hSolver,int loop,int flag );
	int GSE_Krylov_Balance_2( void *hSolver,int loop,int flag );
	void GSE_Krylov_DynShift( GSE_KRYLOV* hSolver );

	GRUS_GSE_DLL int ARE_main( void* hSolver,double *w );
	GRUS_GSE_DLL int KSchur_core( void* hSolver,GRUS_FLOATING **u_0,GRUS_FLOATING **u_1,int *code,int flag );
//	GRUS_GSE_DLL int KSchur_core_1( void* hSolver,GRUS_FLOATING **u_0,GRUS_FLOATING **u_1,int *code,int flag );
	GRUS_GSE_DLL int JDQR_engine( void* hSolver );

	unsigned __stdcall KSchur_ThreadProc(  void* pArguments  );
	unsigned __stdcall KSchur_ThreadProc_1(  void* pArguments  );
	unsigned __stdcall ResidK_ThreadProc(  void* pArguments  );

	GRUS_GSE_DLL int ARE_get_engine( void* hSolver );

	GRUS_GSE_DLL void *GSE_invert( int dim,int*ptr,int *ind,GRUS_FLOATING *val,int ivt_mode,int val_type,int *ret );
//	GRUS_GSE_DLL void GSE_solve( void *hLU,int dim,int ivt_mode,GRUS_FLOATING *v,GRUS_FLOATING *z,int code,int *ret );
	GRUS_GSE_DLL void GSE_solve( void *hLU,int dim,int *ptr,int *ind,GRUS_FLOATING *val,GRUS_FLOATING *z,int ivt_mode,int code,int *ret );
	GRUS_GSE_DLL void GSE_clear( void *hLU,int ivt_mode );

	GRUS_GSE_DLL void GSE_QR_D( int nCol,int *ptr,int *ind,double *val,int flag,double sft[2] );
	GRUS_GSE_DLL void GSE_QR_Z( int nCol,int *ptr,int *ind,DoubleComplex *val,int flag,DoubleComplex sft );

	GRUS_GSE_DLL double CCS_DLA_F( int dim,int *ptr,int *ind,double *val );
	GRUS_GSE_DLL double CCS_DLA_1( int dim,int *ptr,int *ind,double *val );

#ifdef __cplusplus
  }
#endif

#endif