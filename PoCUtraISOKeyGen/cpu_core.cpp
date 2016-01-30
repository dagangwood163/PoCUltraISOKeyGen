#if defined(WIN32) || defined(WIN64)
#include <windows.h>
#elif defined(LINUX) || defined(SOLARIS) || defined(AIX) 
#include <sys/sysinfo.h>
#else 
#error  not support system
#endif 

//in logical
int get_CPU_core_num()
{
#if defined(WIN32) || defined(WIN64)
    SYSTEM_INFO info;
    GetNativeSystemInfo(&info);
    return info.dwNumberOfProcessors;
#elif defined(LINUX) || defined(SOLARIS) || defined(AIX) 
    return get_nprocs();   //GNU fuction 
#else 
#error  not support system
#endif 
}