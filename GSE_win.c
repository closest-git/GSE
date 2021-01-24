#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <windows.h>
#include <process.h>    
#include <stddef.h>
#include "GRUS_debug.h"
#include "GSE_win.h"
#include "GSE_def.h"

#define GSE_MAX_THREADS 8
static int nThread=1;
static void *GSE_u0,*GSE_u1;
static int GSE_code;
static HANDLE event_A[GSE_MAX_THREADS],event_B[GSE_MAX_THREADS],threads[GSE_MAX_THREADS];

/*
	v0.1	cys
		8/13/2007
*/
void __declspec(dllexport) GSE_driver_core( void* hSolver,void **u_0,void **u_1,int *code,int flag  )	{
	int i;
	if( *code==GSE_DRIVER_START )	{
		for( i = 0; i < nThread; i++ )	
			ResumeThread( threads[i] );
	}else	{
		for( i = 0; i < nThread; i++ )	{
			SetEvent( event_B[i] );
		}
	}
	WaitForMultipleObjects( nThread,event_A,TRUE,INFINITE );

	*u_0 = GSE_u0,	*u_1=GSE_u1;
	*code = GSE_code;
	if( *code==GSE_DRIVER_FINISH )	{

		for( i = 0; i < nThread; i++ )	{
			SetEvent( event_B[i] );
		}
	}

	return;
}

/*
	v0.1	cys
		8/13/2007
*/
int GSE_driver_switch( void *u_0,void *u_1,int code,int noThread )	{
	int ret=GSE_OK,wait;
	ASSERT( nThread==1 && noThread==0 );
	GSE_u0 = u_0,	GSE_u1 = u_1;
	GSE_code = code;

	SetEvent( event_A[noThread] );
	wait = WaitForSingleObject( event_B[noThread],INFINITE );
	ASSERT( wait != WAIT_FAILED );

	return ret;
}

/*
	v0.1	cys
		8/13/2007
*/
int GSE_driver_init( unsigned ( __stdcall *ThreadProc )( void * ),void* pArguments )	{
	int i,ret=GSE_DRIVER_INIT_FAIL;
	unsigned threadID;

	for( i = 0; i < nThread; i++ )	{	
		event_A[i] = CreateEvent( NULL,FALSE,FALSE,NULL );
		if( event_A[i]==NULL )
			goto END;
		event_B[i] = CreateEvent( NULL,FALSE,FALSE,NULL );
		if( event_B[i]==NULL )
			goto END;

		threads[i] = (HANDLE)_beginthreadex(NULL, 0x0, ThreadProc,pArguments, CREATE_SUSPENDED, &threadID ); 
		if( threads[i]==NULL )
			goto END;
	}
	ret=GSE_OK;
END:
	return ret;
}

/*
	v0.1	cys
		8/13/2007
*/
int GSE_driver_clear( )	{
	int i,ret=GSE_DRIVER_INIT_FAIL,bRet;
	unsigned threadID;

	for( i = 0; i < nThread; i++ )	{	
		bRet = CloseHandle( event_A[i] );			ASSERT( bRet != 0 );
		bRet = CloseHandle( event_B[i] );			ASSERT( bRet != 0 );
		bRet = CloseHandle( threads[i] );			ASSERT( bRet != 0 );
	}
	ret=GSE_OK;
END:
	return ret;
}