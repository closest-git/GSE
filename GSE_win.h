#ifndef _GSE_WIN_H_
#define _GSE_WIN_H_


#ifdef __cplusplus
extern "C" {
#endif

	int GSE_driver_init( unsigned ( __stdcall *ThreadProc )( void * ),void* pArguments );
	int GSE_driver_clear( );

	int GSE_driver_switch( void *u_0,void *u_1,int code,int noThread );

#ifdef __cplusplus
  }
#endif



#endif
