
typedef signed char __int8_t;
typedef unsigned char __uint8_t;
typedef short __int16_t;
typedef unsigned short __uint16_t;
typedef int __int32_t;
typedef unsigned int __uint32_t;
typedef long long __int64_t;
typedef unsigned long long __uint64_t;
typedef long __darwin_intptr_t;
typedef unsigned int __darwin_natural_t;
typedef int __darwin_ct_rune_t;
typedef union {
 char __mbstate8[128];
 long long _mbstateL;
} __mbstate_t;
typedef __mbstate_t __darwin_mbstate_t;
typedef long int __darwin_ptrdiff_t;
typedef long unsigned int __darwin_size_t;
typedef __builtin_va_list __darwin_va_list;
typedef int __darwin_wchar_t;
typedef __darwin_wchar_t __darwin_rune_t;
typedef int __darwin_wint_t;
typedef unsigned long __darwin_clock_t;
typedef __uint32_t __darwin_socklen_t;
typedef long __darwin_ssize_t;
typedef long __darwin_time_t;
typedef __int64_t __darwin_blkcnt_t;
typedef __int32_t __darwin_blksize_t;
typedef __int32_t __darwin_dev_t;
typedef unsigned int __darwin_fsblkcnt_t;
typedef unsigned int __darwin_fsfilcnt_t;
typedef __uint32_t __darwin_gid_t;
typedef __uint32_t __darwin_id_t;
typedef __uint64_t __darwin_ino64_t;
typedef __darwin_ino64_t __darwin_ino_t;
typedef __darwin_natural_t __darwin_mach_port_name_t;
typedef __darwin_mach_port_name_t __darwin_mach_port_t;
typedef __uint16_t __darwin_mode_t;
typedef __int64_t __darwin_off_t;
typedef __int32_t __darwin_pid_t;
typedef __uint32_t __darwin_sigset_t;
typedef __int32_t __darwin_suseconds_t;
typedef __uint32_t __darwin_uid_t;
typedef __uint32_t __darwin_useconds_t;
typedef unsigned char __darwin_uuid_t[16];
typedef char __darwin_uuid_string_t[37];
struct __darwin_pthread_handler_rec {
 void (*__routine)(void *);
 void *__arg;
 struct __darwin_pthread_handler_rec *__next;
};
struct _opaque_pthread_attr_t {
 long __sig;
 char __opaque[56];
};
struct _opaque_pthread_cond_t {
 long __sig;
 char __opaque[40];
};
struct _opaque_pthread_condattr_t {
 long __sig;
 char __opaque[8];
};
struct _opaque_pthread_mutex_t {
 long __sig;
 char __opaque[56];
};
struct _opaque_pthread_mutexattr_t {
 long __sig;
 char __opaque[8];
};
struct _opaque_pthread_once_t {
 long __sig;
 char __opaque[8];
};
struct _opaque_pthread_rwlock_t {
 long __sig;
 char __opaque[192];
};
struct _opaque_pthread_rwlockattr_t {
 long __sig;
 char __opaque[16];
};
struct _opaque_pthread_t {
 long __sig;
 struct __darwin_pthread_handler_rec *__cleanup_stack;
 char __opaque[8176];
};
typedef struct _opaque_pthread_attr_t __darwin_pthread_attr_t;
typedef struct _opaque_pthread_cond_t __darwin_pthread_cond_t;
typedef struct _opaque_pthread_condattr_t __darwin_pthread_condattr_t;
typedef unsigned long __darwin_pthread_key_t;
typedef struct _opaque_pthread_mutex_t __darwin_pthread_mutex_t;
typedef struct _opaque_pthread_mutexattr_t __darwin_pthread_mutexattr_t;
typedef struct _opaque_pthread_once_t __darwin_pthread_once_t;
typedef struct _opaque_pthread_rwlock_t __darwin_pthread_rwlock_t;
typedef struct _opaque_pthread_rwlockattr_t __darwin_pthread_rwlockattr_t;
typedef struct _opaque_pthread_t *__darwin_pthread_t;
typedef int __darwin_nl_item;
typedef int __darwin_wctrans_t;
typedef __uint32_t __darwin_wctype_t;

typedef enum {
 P_ALL,
 P_PID,
 P_PGID
} idtype_t;
typedef __darwin_pid_t pid_t;
typedef __darwin_id_t id_t;
typedef int sig_atomic_t;
typedef signed char int8_t;
typedef short int16_t;
typedef int int32_t;
typedef long long int64_t;

typedef unsigned char u_int8_t;
typedef unsigned short u_int16_t;
typedef unsigned int u_int32_t;
typedef unsigned long long u_int64_t;
typedef int64_t register_t;
typedef __darwin_intptr_t intptr_t;
typedef unsigned long uintptr_t;
typedef u_int64_t user_addr_t;
typedef u_int64_t user_size_t;
typedef int64_t user_ssize_t;
typedef int64_t user_long_t;
typedef u_int64_t user_ulong_t;
typedef int64_t user_time_t;
typedef int64_t user_off_t;
typedef u_int64_t syscall_arg_t;
struct __darwin_arm_exception_state
{
 __uint32_t __exception;
 __uint32_t __fsr;
 __uint32_t __far;
};
struct __darwin_arm_exception_state64
{
 __uint64_t __far;
 __uint32_t __esr;
 __uint32_t __exception;
};
struct __darwin_arm_exception_state64_v2
{
 __uint64_t __far;
 __uint64_t __esr;
};
struct __darwin_arm_thread_state
{
 __uint32_t __r[13];
 __uint32_t __sp;
 __uint32_t __lr;
 __uint32_t __pc;
 __uint32_t __cpsr;
};
struct __darwin_arm_thread_state64
{
 __uint64_t __x[29];
 __uint64_t __fp;
 __uint64_t __lr;
 __uint64_t __sp;
 __uint64_t __pc;
 __uint32_t __cpsr;
 __uint32_t __pad;
};
struct __darwin_arm_vfp_state
{
 __uint32_t __r[64];
 __uint32_t __fpscr;
};
struct __darwin_arm_neon_state64
{
 __uint128_t __v[32];
 __uint32_t __fpsr;
 __uint32_t __fpcr;
};
struct __darwin_arm_neon_state
{
 __uint128_t __v[16];
 __uint32_t __fpsr;
 __uint32_t __fpcr;
};
struct __arm_pagein_state
{
 int __pagein_error;
};
struct __darwin_arm_sme_state
{
 __uint64_t __svcr;
 __uint64_t __tpidr2_el0;
 __uint16_t __svl_b;
};
struct __darwin_arm_sve_z_state
{
 char __z[16][256];
} __attribute__((aligned(4)));
struct __darwin_arm_sve_p_state
{
 char __p[16][256 / 8];
} __attribute__((aligned(4)));
struct __darwin_arm_sme_za_state
{
 char __za[4096];
} __attribute__((aligned(4)));
struct __darwin_arm_sme2_state
{
 char __zt0[64];
} __attribute__((aligned(4)));
struct __arm_legacy_debug_state
{
 __uint32_t __bvr[16];
 __uint32_t __bcr[16];
 __uint32_t __wvr[16];
 __uint32_t __wcr[16];
};
struct __darwin_arm_debug_state32
{
 __uint32_t __bvr[16];
 __uint32_t __bcr[16];
 __uint32_t __wvr[16];
 __uint32_t __wcr[16];
 __uint64_t __mdscr_el1;
};
struct __darwin_arm_debug_state64
{
 __uint64_t __bvr[16];
 __uint64_t __bcr[16];
 __uint64_t __wvr[16];
 __uint64_t __wcr[16];
 __uint64_t __mdscr_el1;
};
struct __darwin_arm_cpmu_state64
{
 __uint64_t __ctrs[16];
};
struct __darwin_mcontext32
{
 struct __darwin_arm_exception_state __es;
 struct __darwin_arm_thread_state __ss;
 struct __darwin_arm_vfp_state __fs;
};
struct __darwin_mcontext64
{
 struct __darwin_arm_exception_state64 __es;
 struct __darwin_arm_thread_state64 __ss;
 struct __darwin_arm_neon_state64 __ns;
};
typedef struct __darwin_mcontext64 *mcontext_t;

typedef __darwin_pthread_attr_t pthread_attr_t;

struct __darwin_sigaltstack
{
 void *ss_sp;
 __darwin_size_t ss_size;
 int ss_flags;
};
typedef struct __darwin_sigaltstack stack_t;
struct __darwin_ucontext
{
 int uc_onstack;
 __darwin_sigset_t uc_sigmask;
 struct __darwin_sigaltstack uc_stack;
 struct __darwin_ucontext *uc_link;
 __darwin_size_t uc_mcsize;
 struct __darwin_mcontext64 *uc_mcontext;
};
typedef struct __darwin_ucontext ucontext_t;
typedef __darwin_sigset_t sigset_t;
typedef __darwin_size_t size_t;
typedef __darwin_uid_t uid_t;

union sigval {
 int sival_int;
 void *sival_ptr;
};
struct sigevent {
 int sigev_notify;
 int sigev_signo;
 union sigval sigev_value;
 void (*sigev_notify_function)(union sigval);
 pthread_attr_t *sigev_notify_attributes;
};
typedef struct __siginfo {
 int si_signo;
 int si_errno;
 int si_code;
 pid_t si_pid;
 uid_t si_uid;
 int si_status;
 void *si_addr;
 union sigval si_value;
 long si_band;
 unsigned long __pad[7];
} siginfo_t;
union __sigaction_u {
 void (*__sa_handler)(int);
 void (*__sa_sigaction)(int, struct __siginfo *,
     void *);
};
struct __sigaction {
 union __sigaction_u __sigaction_u;
 void (*sa_tramp)(void *, int, int, siginfo_t *, void *);
 sigset_t sa_mask;
 int sa_flags;
};
struct sigaction {
 union __sigaction_u __sigaction_u;
 sigset_t sa_mask;
 int sa_flags;
};
typedef void (*sig_t)(int);
struct sigvec {
 void (*sv_handler)(int);
 int sv_mask;
 int sv_flags;
};
struct sigstack {
 char *ss_sp;
 int ss_onstack;
};
void(*signal(int, void (*)(int)))(int);
typedef unsigned char uint8_t;
typedef unsigned short uint16_t;
typedef unsigned int uint32_t;
typedef unsigned long long uint64_t;
typedef int8_t int_least8_t;
typedef int16_t int_least16_t;
typedef int32_t int_least32_t;
typedef int64_t int_least64_t;
typedef uint8_t uint_least8_t;
typedef uint16_t uint_least16_t;
typedef uint32_t uint_least32_t;
typedef uint64_t uint_least64_t;
typedef int8_t int_fast8_t;
typedef int16_t int_fast16_t;
typedef int32_t int_fast32_t;
typedef int64_t int_fast64_t;
typedef uint8_t uint_fast8_t;
typedef uint16_t uint_fast16_t;
typedef uint32_t uint_fast32_t;
typedef uint64_t uint_fast64_t;
typedef long int intmax_t;
typedef long unsigned int uintmax_t;
struct timeval
{
 __darwin_time_t tv_sec;
 __darwin_suseconds_t tv_usec;
};
typedef __uint64_t rlim_t;
struct rusage {
 struct timeval ru_utime;
 struct timeval ru_stime;
 long ru_maxrss;
 long ru_ixrss;
 long ru_idrss;
 long ru_isrss;
 long ru_minflt;
 long ru_majflt;
 long ru_nswap;
 long ru_inblock;
 long ru_oublock;
 long ru_msgsnd;
 long ru_msgrcv;
 long ru_nsignals;
 long ru_nvcsw;
 long ru_nivcsw;
};
typedef void *rusage_info_t;
struct rusage_info_v0 {
 uint8_t ri_uuid[16];
 uint64_t ri_user_time;
 uint64_t ri_system_time;
 uint64_t ri_pkg_idle_wkups;
 uint64_t ri_interrupt_wkups;
 uint64_t ri_pageins;
 uint64_t ri_wired_size;
 uint64_t ri_resident_size;
 uint64_t ri_phys_footprint;
 uint64_t ri_proc_start_abstime;
 uint64_t ri_proc_exit_abstime;
};
struct rusage_info_v1 {
 uint8_t ri_uuid[16];
 uint64_t ri_user_time;
 uint64_t ri_system_time;
 uint64_t ri_pkg_idle_wkups;
 uint64_t ri_interrupt_wkups;
 uint64_t ri_pageins;
 uint64_t ri_wired_size;
 uint64_t ri_resident_size;
 uint64_t ri_phys_footprint;
 uint64_t ri_proc_start_abstime;
 uint64_t ri_proc_exit_abstime;
 uint64_t ri_child_user_time;
 uint64_t ri_child_system_time;
 uint64_t ri_child_pkg_idle_wkups;
 uint64_t ri_child_interrupt_wkups;
 uint64_t ri_child_pageins;
 uint64_t ri_child_elapsed_abstime;
};
struct rusage_info_v2 {
 uint8_t ri_uuid[16];
 uint64_t ri_user_time;
 uint64_t ri_system_time;
 uint64_t ri_pkg_idle_wkups;
 uint64_t ri_interrupt_wkups;
 uint64_t ri_pageins;
 uint64_t ri_wired_size;
 uint64_t ri_resident_size;
 uint64_t ri_phys_footprint;
 uint64_t ri_proc_start_abstime;
 uint64_t ri_proc_exit_abstime;
 uint64_t ri_child_user_time;
 uint64_t ri_child_system_time;
 uint64_t ri_child_pkg_idle_wkups;
 uint64_t ri_child_interrupt_wkups;
 uint64_t ri_child_pageins;
 uint64_t ri_child_elapsed_abstime;
 uint64_t ri_diskio_bytesread;
 uint64_t ri_diskio_byteswritten;
};
struct rusage_info_v3 {
 uint8_t ri_uuid[16];
 uint64_t ri_user_time;
 uint64_t ri_system_time;
 uint64_t ri_pkg_idle_wkups;
 uint64_t ri_interrupt_wkups;
 uint64_t ri_pageins;
 uint64_t ri_wired_size;
 uint64_t ri_resident_size;
 uint64_t ri_phys_footprint;
 uint64_t ri_proc_start_abstime;
 uint64_t ri_proc_exit_abstime;
 uint64_t ri_child_user_time;
 uint64_t ri_child_system_time;
 uint64_t ri_child_pkg_idle_wkups;
 uint64_t ri_child_interrupt_wkups;
 uint64_t ri_child_pageins;
 uint64_t ri_child_elapsed_abstime;
 uint64_t ri_diskio_bytesread;
 uint64_t ri_diskio_byteswritten;
 uint64_t ri_cpu_time_qos_default;
 uint64_t ri_cpu_time_qos_maintenance;
 uint64_t ri_cpu_time_qos_background;
 uint64_t ri_cpu_time_qos_utility;
 uint64_t ri_cpu_time_qos_legacy;
 uint64_t ri_cpu_time_qos_user_initiated;
 uint64_t ri_cpu_time_qos_user_interactive;
 uint64_t ri_billed_system_time;
 uint64_t ri_serviced_system_time;
};
struct rusage_info_v4 {
 uint8_t ri_uuid[16];
 uint64_t ri_user_time;
 uint64_t ri_system_time;
 uint64_t ri_pkg_idle_wkups;
 uint64_t ri_interrupt_wkups;
 uint64_t ri_pageins;
 uint64_t ri_wired_size;
 uint64_t ri_resident_size;
 uint64_t ri_phys_footprint;
 uint64_t ri_proc_start_abstime;
 uint64_t ri_proc_exit_abstime;
 uint64_t ri_child_user_time;
 uint64_t ri_child_system_time;
 uint64_t ri_child_pkg_idle_wkups;
 uint64_t ri_child_interrupt_wkups;
 uint64_t ri_child_pageins;
 uint64_t ri_child_elapsed_abstime;
 uint64_t ri_diskio_bytesread;
 uint64_t ri_diskio_byteswritten;
 uint64_t ri_cpu_time_qos_default;
 uint64_t ri_cpu_time_qos_maintenance;
 uint64_t ri_cpu_time_qos_background;
 uint64_t ri_cpu_time_qos_utility;
 uint64_t ri_cpu_time_qos_legacy;
 uint64_t ri_cpu_time_qos_user_initiated;
 uint64_t ri_cpu_time_qos_user_interactive;
 uint64_t ri_billed_system_time;
 uint64_t ri_serviced_system_time;
 uint64_t ri_logical_writes;
 uint64_t ri_lifetime_max_phys_footprint;
 uint64_t ri_instructions;
 uint64_t ri_cycles;
 uint64_t ri_billed_energy;
 uint64_t ri_serviced_energy;
 uint64_t ri_interval_max_phys_footprint;
 uint64_t ri_runnable_time;
};
struct rusage_info_v5 {
 uint8_t ri_uuid[16];
 uint64_t ri_user_time;
 uint64_t ri_system_time;
 uint64_t ri_pkg_idle_wkups;
 uint64_t ri_interrupt_wkups;
 uint64_t ri_pageins;
 uint64_t ri_wired_size;
 uint64_t ri_resident_size;
 uint64_t ri_phys_footprint;
 uint64_t ri_proc_start_abstime;
 uint64_t ri_proc_exit_abstime;
 uint64_t ri_child_user_time;
 uint64_t ri_child_system_time;
 uint64_t ri_child_pkg_idle_wkups;
 uint64_t ri_child_interrupt_wkups;
 uint64_t ri_child_pageins;
 uint64_t ri_child_elapsed_abstime;
 uint64_t ri_diskio_bytesread;
 uint64_t ri_diskio_byteswritten;
 uint64_t ri_cpu_time_qos_default;
 uint64_t ri_cpu_time_qos_maintenance;
 uint64_t ri_cpu_time_qos_background;
 uint64_t ri_cpu_time_qos_utility;
 uint64_t ri_cpu_time_qos_legacy;
 uint64_t ri_cpu_time_qos_user_initiated;
 uint64_t ri_cpu_time_qos_user_interactive;
 uint64_t ri_billed_system_time;
 uint64_t ri_serviced_system_time;
 uint64_t ri_logical_writes;
 uint64_t ri_lifetime_max_phys_footprint;
 uint64_t ri_instructions;
 uint64_t ri_cycles;
 uint64_t ri_billed_energy;
 uint64_t ri_serviced_energy;
 uint64_t ri_interval_max_phys_footprint;
 uint64_t ri_runnable_time;
 uint64_t ri_flags;
};
struct rusage_info_v6 {
 uint8_t ri_uuid[16];
 uint64_t ri_user_time;
 uint64_t ri_system_time;
 uint64_t ri_pkg_idle_wkups;
 uint64_t ri_interrupt_wkups;
 uint64_t ri_pageins;
 uint64_t ri_wired_size;
 uint64_t ri_resident_size;
 uint64_t ri_phys_footprint;
 uint64_t ri_proc_start_abstime;
 uint64_t ri_proc_exit_abstime;
 uint64_t ri_child_user_time;
 uint64_t ri_child_system_time;
 uint64_t ri_child_pkg_idle_wkups;
 uint64_t ri_child_interrupt_wkups;
 uint64_t ri_child_pageins;
 uint64_t ri_child_elapsed_abstime;
 uint64_t ri_diskio_bytesread;
 uint64_t ri_diskio_byteswritten;
 uint64_t ri_cpu_time_qos_default;
 uint64_t ri_cpu_time_qos_maintenance;
 uint64_t ri_cpu_time_qos_background;
 uint64_t ri_cpu_time_qos_utility;
 uint64_t ri_cpu_time_qos_legacy;
 uint64_t ri_cpu_time_qos_user_initiated;
 uint64_t ri_cpu_time_qos_user_interactive;
 uint64_t ri_billed_system_time;
 uint64_t ri_serviced_system_time;
 uint64_t ri_logical_writes;
 uint64_t ri_lifetime_max_phys_footprint;
 uint64_t ri_instructions;
 uint64_t ri_cycles;
 uint64_t ri_billed_energy;
 uint64_t ri_serviced_energy;
 uint64_t ri_interval_max_phys_footprint;
 uint64_t ri_runnable_time;
 uint64_t ri_flags;
 uint64_t ri_user_ptime;
 uint64_t ri_system_ptime;
 uint64_t ri_pinstructions;
 uint64_t ri_pcycles;
 uint64_t ri_energy_nj;
 uint64_t ri_penergy_nj;
 uint64_t ri_secure_time_in_system;
 uint64_t ri_secure_ptime_in_system;
 uint64_t ri_neural_footprint;
 uint64_t ri_lifetime_max_neural_footprint;
 uint64_t ri_interval_max_neural_footprint;
 uint64_t ri_reserved[9];
};
typedef struct rusage_info_v6 rusage_info_current;
struct rlimit {
 rlim_t rlim_cur;
 rlim_t rlim_max;
};
struct proc_rlimit_control_wakeupmon {
 uint32_t wm_flags;
 int32_t wm_rate;
};
int getpriority(int, id_t);
int getiopolicy_np(int, int) __attribute__((availability(macosx,introduced=10.5)));
int getrlimit(int, struct rlimit *) __asm("_" "getrlimit" );
int getrusage(int, struct rusage *);
int setpriority(int, id_t, int);
int setiopolicy_np(int, int, int) __attribute__((availability(macosx,introduced=10.5)));
int setrlimit(int, const struct rlimit *) __asm("_" "setrlimit" );
static inline
__uint16_t
_OSSwapInt16(
 __uint16_t _data
 )
{
 return (__uint16_t)(_data << 8 | _data >> 8);
}
static inline
__uint32_t
_OSSwapInt32(
 __uint32_t _data
 )
{
 _data = __builtin_bswap32(_data);
 return _data;
}
static inline
__uint64_t
_OSSwapInt64(
 __uint64_t _data
 )
{
 return __builtin_bswap64(_data);
}
union wait {
 int w_status;
 struct {
  unsigned int w_Termsig:7,
      w_Coredump:1,
      w_Retcode:8,
      w_Filler:16;
 } w_T;
 struct {
  unsigned int w_Stopval:8,
      w_Stopsig:8,
      w_Filler:16;
 } w_S;
};
pid_t wait(int *) __asm("_" "wait" );
pid_t waitpid(pid_t, int *, int) __asm("_" "waitpid" );
int waitid(idtype_t, id_t, siginfo_t *, int) __asm("_" "waitid" );
pid_t wait3(int *, int, struct rusage *);
pid_t wait4(pid_t, int *, int, struct rusage *);

void * alloca(size_t __size);
typedef __darwin_ct_rune_t ct_rune_t;
typedef __darwin_rune_t rune_t;
typedef __darwin_wchar_t wchar_t;
typedef struct {
 int quot;
 int rem;
} div_t;
typedef struct {
 long quot;
 long rem;
} ldiv_t;
typedef struct {
 long long quot;
 long long rem;
} lldiv_t;
extern int __mb_cur_max;
typedef unsigned long long malloc_type_id_t;

__attribute__((availability(macos,introduced=14.0))) __attribute__((availability(ios,introduced=17.0))) __attribute__((availability(tvos,introduced=17.0))) __attribute__((availability(watchos,introduced=10.0))) __attribute__((availability(visionos,introduced=1.0))) __attribute__((availability(driverkit,introduced=23.0)))
void * malloc_type_malloc(size_t size, malloc_type_id_t type_id) __attribute__((__warn_unused_result__)) __attribute__((alloc_size(1)));
__attribute__((availability(macos,introduced=14.0))) __attribute__((availability(ios,introduced=17.0))) __attribute__((availability(tvos,introduced=17.0))) __attribute__((availability(watchos,introduced=10.0))) __attribute__((availability(visionos,introduced=1.0))) __attribute__((availability(driverkit,introduced=23.0)))
void * malloc_type_calloc(size_t count, size_t size, malloc_type_id_t type_id) __attribute__((__warn_unused_result__)) __attribute__((alloc_size(1,2)));
__attribute__((availability(macos,introduced=14.0))) __attribute__((availability(ios,introduced=17.0))) __attribute__((availability(tvos,introduced=17.0))) __attribute__((availability(watchos,introduced=10.0))) __attribute__((availability(visionos,introduced=1.0))) __attribute__((availability(driverkit,introduced=23.0)))
void malloc_type_free(void * ptr, malloc_type_id_t type_id);
__attribute__((availability(macos,introduced=14.0))) __attribute__((availability(ios,introduced=17.0))) __attribute__((availability(tvos,introduced=17.0))) __attribute__((availability(watchos,introduced=10.0))) __attribute__((availability(visionos,introduced=1.0))) __attribute__((availability(driverkit,introduced=23.0)))
void * malloc_type_realloc(void * ptr, size_t size, malloc_type_id_t type_id) __attribute__((__warn_unused_result__)) __attribute__((alloc_size(2)));
__attribute__((availability(macos,introduced=14.0))) __attribute__((availability(ios,introduced=17.0))) __attribute__((availability(tvos,introduced=17.0))) __attribute__((availability(watchos,introduced=10.0))) __attribute__((availability(visionos,introduced=1.0))) __attribute__((availability(driverkit,introduced=23.0)))
void * malloc_type_valloc(size_t size, malloc_type_id_t type_id) __attribute__((__warn_unused_result__)) __attribute__((alloc_size(1)));
__attribute__((availability(macos,introduced=14.0))) __attribute__((availability(ios,introduced=17.0))) __attribute__((availability(tvos,introduced=17.0))) __attribute__((availability(watchos,introduced=10.0))) __attribute__((availability(visionos,introduced=1.0))) __attribute__((availability(driverkit,introduced=23.0)))
void * malloc_type_aligned_alloc(size_t alignment, size_t size, malloc_type_id_t type_id) __attribute__((__warn_unused_result__)) __attribute__((alloc_align(1))) __attribute__((alloc_size(2)));
__attribute__((availability(macos,introduced=14.0))) __attribute__((availability(ios,introduced=17.0))) __attribute__((availability(tvos,introduced=17.0))) __attribute__((availability(watchos,introduced=10.0))) __attribute__((availability(visionos,introduced=1.0))) __attribute__((availability(driverkit,introduced=23.0)))
int malloc_type_posix_memalign(void * *memptr, size_t alignment, size_t size, malloc_type_id_t type_id) ;
typedef struct _malloc_zone_t malloc_zone_t;
__attribute__((availability(macos,introduced=14.0))) __attribute__((availability(ios,introduced=17.0))) __attribute__((availability(tvos,introduced=17.0))) __attribute__((availability(watchos,introduced=10.0))) __attribute__((availability(visionos,introduced=1.0))) __attribute__((availability(driverkit,introduced=23.0)))
void * malloc_type_zone_malloc(malloc_zone_t *zone, size_t size, malloc_type_id_t type_id) __attribute__((__warn_unused_result__)) __attribute__((alloc_size(2)));
__attribute__((availability(macos,introduced=14.0))) __attribute__((availability(ios,introduced=17.0))) __attribute__((availability(tvos,introduced=17.0))) __attribute__((availability(watchos,introduced=10.0))) __attribute__((availability(visionos,introduced=1.0))) __attribute__((availability(driverkit,introduced=23.0)))
void * malloc_type_zone_calloc(malloc_zone_t *zone, size_t count, size_t size, malloc_type_id_t type_id) __attribute__((__warn_unused_result__)) __attribute__((alloc_size(2,3)));
__attribute__((availability(macos,introduced=14.0))) __attribute__((availability(ios,introduced=17.0))) __attribute__((availability(tvos,introduced=17.0))) __attribute__((availability(watchos,introduced=10.0))) __attribute__((availability(visionos,introduced=1.0))) __attribute__((availability(driverkit,introduced=23.0)))
void malloc_type_zone_free(malloc_zone_t *zone, void * ptr, malloc_type_id_t type_id);
__attribute__((availability(macos,introduced=14.0))) __attribute__((availability(ios,introduced=17.0))) __attribute__((availability(tvos,introduced=17.0))) __attribute__((availability(watchos,introduced=10.0))) __attribute__((availability(visionos,introduced=1.0))) __attribute__((availability(driverkit,introduced=23.0)))
void * malloc_type_zone_realloc(malloc_zone_t *zone, void * ptr, size_t size, malloc_type_id_t type_id) __attribute__((__warn_unused_result__)) __attribute__((alloc_size(3)));
__attribute__((availability(macos,introduced=14.0))) __attribute__((availability(ios,introduced=17.0))) __attribute__((availability(tvos,introduced=17.0))) __attribute__((availability(watchos,introduced=10.0))) __attribute__((availability(visionos,introduced=1.0))) __attribute__((availability(driverkit,introduced=23.0)))
void * malloc_type_zone_valloc(malloc_zone_t *zone, size_t size, malloc_type_id_t type_id) __attribute__((__warn_unused_result__)) __attribute__((alloc_size(2)));
__attribute__((availability(macos,introduced=14.0))) __attribute__((availability(ios,introduced=17.0))) __attribute__((availability(tvos,introduced=17.0))) __attribute__((availability(watchos,introduced=10.0))) __attribute__((availability(visionos,introduced=1.0))) __attribute__((availability(driverkit,introduced=23.0)))
void * malloc_type_zone_memalign(malloc_zone_t *zone, size_t alignment, size_t size, malloc_type_id_t type_id) __attribute__((__warn_unused_result__)) __attribute__((alloc_align(2))) __attribute__((alloc_size(3)));
void * malloc(size_t __size) __attribute__((__warn_unused_result__)) __attribute__((alloc_size(1))) ;
void * calloc(size_t __count, size_t __size) __attribute__((__warn_unused_result__)) __attribute__((alloc_size(1,2))) ;
void free(void * );
void * realloc(void * __ptr, size_t __size) __attribute__((__warn_unused_result__)) __attribute__((alloc_size(2))) ;
void * reallocf(void * __ptr, size_t __size) __attribute__((__warn_unused_result__)) __attribute__((alloc_size(2)));
void * valloc(size_t __size) __attribute__((__warn_unused_result__)) __attribute__((alloc_size(1))) ;
void * aligned_alloc(size_t __alignment, size_t __size) __attribute__((__warn_unused_result__)) __attribute__((alloc_align(1))) __attribute__((alloc_size(2))) __attribute__((availability(macosx,introduced=10.15))) __attribute__((availability(ios,introduced=13.0))) __attribute__((availability(tvos,introduced=13.0))) __attribute__((availability(watchos,introduced=6.0)));
int posix_memalign(void * *__memptr, size_t __alignment, size_t __size) __attribute__((availability(macosx,introduced=10.6)));
void abort(void) __attribute__((__cold__)) __attribute__((__noreturn__));
int abs(int) __attribute__((__const__));
int atexit(void (* _Nonnull)(void));
int at_quick_exit(void (*)(void));
double atof(const char *);
int atoi(const char *);
long atol(const char *);
long long
  atoll(const char *);
void *bsearch(const void * __key, const void * __base, size_t __nel,
     size_t __width, int (* _Nonnull __compar)(const void *, const void *));
div_t div(int, int) __attribute__((__const__));
void exit(int) __attribute__((__noreturn__));
char * getenv(const char *);
long labs(long) __attribute__((__const__));
ldiv_t ldiv(long, long) __attribute__((__const__));
long long
  llabs(long long);
lldiv_t lldiv(long long, long long);
int mblen(const char * __s, size_t __n);
size_t mbstowcs(wchar_t * restrict , const char * restrict, size_t __n);
int mbtowc(wchar_t * restrict, const char * restrict , size_t __n);
void qsort(void * __base, size_t __nel, size_t __width,
     int (* _Nonnull __compar)(const void *, const void *));
void quick_exit(int) __attribute__((__noreturn__));
int rand(void) __attribute__((__availability__(swift, unavailable, message="Use arc4random instead.")));
void srand(unsigned) __attribute__((__availability__(swift, unavailable, message="Use arc4random instead.")));
double strtod(const char *, char * *) __asm("_" "strtod" );
float strtof(const char *, char * *) __asm("_" "strtof" );
long strtol(const char *__str, char * *__endptr, int __base);
long double
  strtold(const char *, char * *);
long long
  strtoll(const char *__str, char * *__endptr, int __base);
unsigned long
  strtoul(const char *__str, char * *__endptr, int __base);
unsigned long long
  strtoull(const char *__str, char * *__endptr, int __base);
__attribute__((__availability__(swift, unavailable, message="Use posix_spawn APIs or NSTask instead. (On iOS, process spawning is unavailable.)")))
__attribute__((availability(macos,introduced=10.0))) __attribute__((availability(ios,unavailable)))
__attribute__((availability(watchos,unavailable))) __attribute__((availability(tvos,unavailable)))
int system(const char *) __asm("_" "system" );
size_t wcstombs(char * restrict , const wchar_t * restrict, size_t __n);
int wctomb(char *, wchar_t);
void _Exit(int) __attribute__((__noreturn__));
long a64l(const char *);
double drand48(void);
char * ecvt(double, int, int *restrict, int *restrict);
double erand48(unsigned short[3]);
char * fcvt(double, int, int *restrict, int *restrict);
char * gcvt(double, int, char *) ;
int getsubopt(char * *, char * const *, char * *);
int grantpt(int);
char *
  initstate(unsigned, char *, size_t __size);
long jrand48(unsigned short[3]) __attribute__((__availability__(swift, unavailable, message="Use arc4random instead.")));
char *l64a(long);
void lcong48(unsigned short[7]);
long lrand48(void) __attribute__((__availability__(swift, unavailable, message="Use arc4random instead.")));
__attribute__((__deprecated__("This function is provided for compatibility reasons only.  Due to security concerns inherent in the design of mktemp(3), it is highly recommended that you use mkstemp(3) instead.")))
char * mktemp(char *);
int mkstemp(char *);
long mrand48(void) __attribute__((__availability__(swift, unavailable, message="Use arc4random instead.")));
long nrand48(unsigned short[3]) __attribute__((__availability__(swift, unavailable, message="Use arc4random instead.")));
int posix_openpt(int);
char * ptsname(int);
int ptsname_r(int fildes, char * buffer, size_t buflen) __attribute__((availability(macos,introduced=10.13.4))) __attribute__((availability(ios,introduced=11.3))) __attribute__((availability(tvos,introduced=11.3))) __attribute__((availability(watchos,introduced=4.3)));
int putenv(char *) __asm("_" "putenv" );
long random(void) __attribute__((__availability__(swift, unavailable, message="Use arc4random instead.")));
int rand_r(unsigned *) __attribute__((__availability__(swift, unavailable, message="Use arc4random instead.")));
char * realpath(const char * restrict, char * restrict ) __asm("_" "realpath" "$DARWIN_EXTSN");
unsigned short * seed48(unsigned short[3]);
int setenv(const char * __name, const char * __value, int __overwrite) __asm("_" "setenv" );
void setkey(const char *) __asm("_" "setkey" );
char * setstate(const char *);
void srand48(long);
void srandom(unsigned);
int unlockpt(int);
int unsetenv(const char *) __asm("_" "unsetenv" );
typedef __darwin_dev_t dev_t;
typedef __darwin_mode_t mode_t;
uint32_t arc4random(void);
void arc4random_addrandom(unsigned char * , int __datlen)
    __attribute__((availability(macosx,introduced=10.0))) __attribute__((availability(macosx,deprecated=10.12,message="use arc4random_stir")))
    __attribute__((availability(ios,introduced=2.0))) __attribute__((availability(ios,deprecated=10.0,message="use arc4random_stir")))
    __attribute__((availability(tvos,introduced=2.0))) __attribute__((availability(tvos,deprecated=10.0,message="use arc4random_stir")))
    __attribute__((availability(watchos,introduced=1.0))) __attribute__((availability(watchos,deprecated=3.0,message="use arc4random_stir")));
void arc4random_buf(void * __buf, size_t __nbytes) __attribute__((availability(macosx,introduced=10.7)));
void arc4random_stir(void);
uint32_t
  arc4random_uniform(uint32_t __upper_bound) __attribute__((availability(macosx,introduced=10.7)));
int atexit_b(void (^ _Nonnull)(void)) __attribute__((availability(macosx,introduced=10.6)));
void *bsearch_b(const void * __key, const void * __base, size_t __nel,
     size_t __width, int (^ _Nonnull __compar)(const void *, const void *) __attribute__((__noescape__)))
     __attribute__((availability(macosx,introduced=10.6)));
char * cgetcap(char *, const char *, int);
int cgetclose(void);
int cgetent(char * *, char * *, const char *);
int cgetfirst(char * *, char * *);
int cgetmatch(const char *, const char *);
int cgetnext(char * *, char * *);
int cgetnum(char *, const char *, long *);
int cgetset(const char *);
int cgetstr(char *, const char *, char * *);
int cgetustr(char *, const char *, char * *);
int daemon(int, int) __asm("_" "daemon" ) __attribute__((availability(macosx,introduced=10.0,deprecated=10.5,message="Use posix_spawn APIs instead."))) __attribute__((availability(watchos,unavailable))) __attribute__((availability(tvos,unavailable)));
char * devname(dev_t, mode_t);
char * devname_r(dev_t, mode_t, char * buf, int len);
char * getbsize(int *, long *);
int getloadavg(double [], int __nelem);
const char
 *getprogname(void);
void setprogname(const char *);
int heapsort(void * __base, size_t __nel, size_t __width,
     int (* _Nonnull __compar)(const void *, const void *));
int heapsort_b(void * __base, size_t __nel, size_t __width,
     int (^ _Nonnull __compar)(const void *, const void *) __attribute__((__noescape__)))
     __attribute__((availability(macosx,introduced=10.6)));
int mergesort(void * __base, size_t __nel, size_t __width,
     int (* _Nonnull __compar)(const void *, const void *));
int mergesort_b(void * __base, size_t __nel, size_t __width,
     int (^ _Nonnull __compar)(const void *, const void *) __attribute__((__noescape__)))
     __attribute__((availability(macosx,introduced=10.6)));
void psort(void * __base, size_t __nel, size_t __width,
     int (* _Nonnull __compar)(const void *, const void *))
     __attribute__((availability(macosx,introduced=10.6)));
void psort_b(void * __base, size_t __nel, size_t __width,
     int (^ _Nonnull __compar)(const void *, const void *) __attribute__((__noescape__)))
     __attribute__((availability(macosx,introduced=10.6)));
void psort_r(void * __base, size_t __nel, size_t __width, void *,
     int (* _Nonnull __compar)(void *, const void *, const void *))
     __attribute__((availability(macosx,introduced=10.6)));
void qsort_b(void * __base, size_t __nel, size_t __width,
     int (^ _Nonnull __compar)(const void *, const void *) __attribute__((__noescape__)))
     __attribute__((availability(macosx,introduced=10.6)));
void qsort_r(void * __base, size_t __nel, size_t __width, void *,
     int (* _Nonnull __compar)(void *, const void *, const void *));
int radixsort(const unsigned char * * __base, int __nel, const unsigned char * __table,
     unsigned __endbyte);
int rpmatch(const char *)
 __attribute__((availability(macos,introduced=10.15))) __attribute__((availability(ios,introduced=13.0))) __attribute__((availability(tvos,introduced=13.0))) __attribute__((availability(watchos,introduced=6.0)));
int sradixsort(const unsigned char * * __base, int __nel, const unsigned char * __table,
     unsigned __endbyte);
void sranddev(void);
void srandomdev(void);
long long
 strtonum(const char *__numstr, long long __minval, long long __maxval, const char * *__errstrp)
 __attribute__((availability(macos,introduced=11.0))) __attribute__((availability(ios,introduced=14.0))) __attribute__((availability(tvos,introduced=14.0))) __attribute__((availability(watchos,introduced=7.0)));
long long
  strtoq(const char *__str, char * *__endptr, int __base);
unsigned long long
  strtouq(const char *__str, char * *__endptr, int __base);
extern char * suboptarg;
typedef __darwin_va_list va_list;

int renameat(int, const char *, int, const char *) __attribute__((availability(macosx,introduced=10.10)));
int renamex_np(const char *, const char *, unsigned int) __attribute__((availability(macosx,introduced=10.12))) __attribute__((availability(ios,introduced=10.0))) __attribute__((availability(tvos,introduced=10.0))) __attribute__((availability(watchos,introduced=3.0)));
int renameatx_np(int, const char *, int, const char *, unsigned int) __attribute__((availability(macosx,introduced=10.12))) __attribute__((availability(ios,introduced=10.0))) __attribute__((availability(tvos,introduced=10.0))) __attribute__((availability(watchos,introduced=3.0)));
int printf(const char * restrict, ...) __attribute__((__format__ (__printf__, 1, 2)));
typedef __darwin_off_t fpos_t;
struct __sbuf {
 unsigned char * _base;
 int _size;
};
struct __sFILEX;
typedef struct __sFILE {
 unsigned char * _p;
 int _r;
 int _w;
 short _flags;
 short _file;
 struct __sbuf _bf;
 int _lbfsize;
 void *_cookie;
 int (* _Nullable _close)(void *);
 int (* _Nullable _read) (void *, char *, int __n);
 fpos_t (* _Nullable _seek) (void *, fpos_t, int);
 int (* _Nullable _write)(void *, const char *, int __n);
 struct __sbuf _ub;
 struct __sFILEX *_extra;
 int _ur;
 unsigned char _ubuf[3];
 unsigned char _nbuf[1];
 struct __sbuf _lb;
 int _blksize;
 fpos_t _offset;
} FILE;
extern FILE *__stdinp __attribute__((__swift_attr__("nonisolated(unsafe)")));
extern FILE *__stdoutp __attribute__((__swift_attr__("nonisolated(unsafe)")));
extern FILE *__stderrp __attribute__((__swift_attr__("nonisolated(unsafe)")));
void clearerr(FILE *);
int fclose(FILE *);
int feof(FILE *);
int ferror(FILE *);
int fflush(FILE *);
int fgetc(FILE *);
int fgetpos(FILE * restrict, fpos_t *);
char * fgets(char * restrict , int __size, FILE *);
FILE *fopen(const char * restrict __filename, const char * restrict __mode) __asm("_" "fopen" );
int fprintf(FILE * restrict, const char * restrict, ...) __attribute__((__format__ (__printf__, 2, 3)));
int fputc(int, FILE *);
int fputs(const char * restrict, FILE * restrict) __asm("_" "fputs" );
size_t fread(void * restrict __ptr, size_t __size, size_t __nitems, FILE * restrict __stream);
FILE *freopen(const char * restrict, const char * restrict,
     FILE * restrict) __asm("_" "freopen" );
int fscanf(FILE * restrict, const char * restrict, ...) __attribute__((__format__ (__scanf__, 2, 3)));
int fseek(FILE *, long, int);
int fsetpos(FILE *, const fpos_t *);
long ftell(FILE *);
size_t fwrite(const void * restrict __ptr, size_t __size, size_t __nitems, FILE * restrict __stream) __asm("_" "fwrite" );
int getc(FILE *);
int getchar(void);
__attribute__((__deprecated__("This function is provided for compatibility reasons only.  Due to security concerns inherent in the design of gets(3), it is highly recommended that you use fgets(3) instead.")))
char * gets(char *) ;
void perror(const char *) __attribute__((__cold__));
int putc(int, FILE *);
int putchar(int);
int puts(const char *);
int remove(const char *);
int rename (const char *__old, const char *__new);
void rewind(FILE *);
int scanf(const char * restrict, ...) __attribute__((__format__ (__scanf__, 1, 2)));
void setbuf(FILE * restrict, char * restrict );
int setvbuf(FILE * restrict, char * restrict , int, size_t __size);
__attribute__((__availability__(swift, unavailable, message="Use snprintf instead.")))
__attribute__((__deprecated__("This function is provided for compatibility reasons only.  Due to security concerns inherent in the design of sprintf(3), it is highly recommended that you use snprintf(3) instead.")))
int sprintf(char * restrict , const char * restrict, ...) __attribute__((__format__ (__printf__, 2, 3))) ;
int sscanf(const char * restrict, const char * restrict, ...) __attribute__((__format__ (__scanf__, 2, 3)));
FILE *tmpfile(void);
__attribute__((__availability__(swift, unavailable, message="Use mkstemp(3) instead.")))
__attribute__((__deprecated__("This function is provided for compatibility reasons only.  Due to security concerns inherent in the design of tmpnam(3), it is highly recommended that you use mkstemp(3) instead.")))
char * tmpnam(char *);
int ungetc(int, FILE *);
int vfprintf(FILE * restrict, const char * restrict, va_list) __attribute__((__format__ (__printf__, 2, 0)));
int vprintf(const char * restrict, va_list) __attribute__((__format__ (__printf__, 1, 0)));
__attribute__((__availability__(swift, unavailable, message="Use vsnprintf instead.")))
__attribute__((__deprecated__("This function is provided for compatibility reasons only.  Due to security concerns inherent in the design of sprintf(3), it is highly recommended that you use vsnprintf(3) instead.")))
int vsprintf(char * restrict , const char * restrict, va_list) __attribute__((__format__ (__printf__, 2, 0))) ;
char * ctermid(char *);
FILE *fdopen(int, const char *) __asm("_" "fdopen" );
int fileno(FILE *);
int pclose(FILE *) __attribute__((__availability__(swift, unavailable, message="Use posix_spawn APIs or NSTask instead. (On iOS, process spawning is unavailable.)")));
FILE *popen(const char *, const char *) __asm("_" "popen" ) __attribute__((__availability__(swift, unavailable, message="Use posix_spawn APIs or NSTask instead. (On iOS, process spawning is unavailable.)")));
int __srget(FILE *);
int __svfscanf(FILE *, const char *, va_list) __attribute__((__format__ (__scanf__, 2, 0)));
int __swbuf(int, FILE *);
inline __attribute__ ((__always_inline__)) int __sputc(int _c, FILE *_p) {
 if (--_p->_w >= 0 || (_p->_w >= _p->_lbfsize && (char)_c != '\n'))
  return (*_p->_p++ = _c);
 else
  return (__swbuf(_c, _p));
}
void flockfile(FILE *);
int ftrylockfile(FILE *);
void funlockfile(FILE *);
int getc_unlocked(FILE *);
int getchar_unlocked(void);
int putc_unlocked(int, FILE *);
int putchar_unlocked(int);
int getw(FILE *);
int putw(int, FILE *);
__attribute__((__availability__(swift, unavailable, message="Use mkstemp(3) instead.")))
__attribute__((__deprecated__("This function is provided for compatibility reasons only.  Due to security concerns inherent in the design of tempnam(3), it is highly recommended that you use mkstemp(3) instead.")))
char * tempnam(const char *__dir, const char *__prefix) __asm("_" "tempnam" );
typedef __darwin_off_t off_t;
int fseeko(FILE * __stream, off_t __offset, int __whence);
off_t ftello(FILE * __stream);
int snprintf(char * restrict __str, size_t __size, const char * restrict __format, ...) __attribute__((__format__ (__printf__, 3, 4)));
int vfscanf(FILE * restrict __stream, const char * restrict __format, va_list) __attribute__((__format__ (__scanf__, 2, 0)));
int vscanf(const char * restrict __format, va_list) __attribute__((__format__ (__scanf__, 1, 0)));
int vsnprintf(char * restrict __str, size_t __size, const char * restrict __format, va_list) __attribute__((__format__ (__printf__, 3, 0)));
int vsscanf(const char * restrict __str, const char * restrict __format, va_list) __attribute__((__format__ (__scanf__, 2, 0)));
typedef __darwin_ssize_t ssize_t;
int dprintf(int, const char * restrict, ...) __attribute__((__format__ (__printf__, 2, 3))) __attribute__((availability(macosx,introduced=10.7)));
int vdprintf(int, const char * restrict, va_list) __attribute__((__format__ (__printf__, 2, 0))) __attribute__((availability(macosx,introduced=10.7)));
ssize_t getdelim(char * *restrict __linep, size_t * restrict __linecapp, int __delimiter, FILE * restrict __stream) __attribute__((availability(macosx,introduced=10.7)));
ssize_t getline(char * *restrict __linep, size_t * restrict __linecapp, FILE * restrict __stream) __attribute__((availability(macosx,introduced=10.7)));
FILE *fmemopen(void * restrict __buf , size_t __size, const char * restrict __mode) __attribute__((availability(macos,introduced=10.13))) __attribute__((availability(ios,introduced=11.0))) __attribute__((availability(tvos,introduced=11.0))) __attribute__((availability(watchos,introduced=4.0)));
FILE *open_memstream(char * *__bufp, size_t *__sizep) __attribute__((availability(macos,introduced=10.13))) __attribute__((availability(ios,introduced=11.0))) __attribute__((availability(tvos,introduced=11.0))) __attribute__((availability(watchos,introduced=4.0)));
extern const int sys_nerr;
extern const char *const sys_errlist[];
int asprintf(char * *restrict, const char * restrict, ...) __attribute__((__format__ (__printf__, 2, 3)));
char * ctermid_r(char *);
char * fgetln(FILE *, size_t *__len);
const char *fmtcheck(const char *, const char *) __attribute__((format_arg(2)));
int fpurge(FILE *);
void setbuffer(FILE *, char *, int __size);
int setlinebuf(FILE *);
int vasprintf(char * *restrict, const char * restrict, va_list) __attribute__((__format__ (__printf__, 2, 0)));
FILE *funopen(const void *,
     int (* _Nullable)(void *, char *, int __n),
     int (* _Nullable)(void *, const char *, int __n),
     fpos_t (* _Nullable)(void *, fpos_t, int),
     int (* _Nullable)(void *));
extern int __snprintf_chk (char * restrict , size_t __maxlen, int, size_t,
     const char * restrict, ...);
extern int __vsnprintf_chk (char * restrict , size_t __maxlen, int, size_t,
     const char * restrict, va_list);
extern int __sprintf_chk (char * restrict , int, size_t,
     const char * restrict, ...);
extern int __vsprintf_chk (char * restrict , int, size_t,
     const char * restrict, va_list);
void *
  memchr(const void * __s, int __c, size_t __n);
int memcmp(const void * __s1, const void * __s2,
  size_t __n);
void *
  memcpy(void * __dst, const void * __src,
  size_t __n);
void *
  memmove(void * __dst,
  const void * __src, size_t __len);
void *
  memset(void * __b, int __c, size_t __len);
char *
  strcat(char * __s1, const char *__s2)
                                  ;
char * strchr(const char *__s, int __c);
int strcmp(const char *__s1, const char *__s2);
int strcoll(const char *__s1, const char *__s2);
char *
  strcpy(char * __dst, const char *__src)
                                  ;
size_t strcspn(const char *__s, const char *__charset);
char * strerror(int __errnum) __asm("_" "strerror" );
size_t strlen(const char *__s);
char *
  strncat(char * __s1,
  const char * __s2, size_t __n)
                                  ;
int strncmp(const char * __s1,
  const char * __s2, size_t __n);
char *
  strncpy(char * __dst,
        const char * __src, size_t __n)
                                        ;
char * strpbrk(const char *__s, const char *__charset);
char * strrchr(const char *__s, int __c);
size_t strspn(const char *__s, const char *__charset);
char * strstr(const char *__big, const char *__little);
char * strtok(char * __str, const char *__sep);
size_t strxfrm(char * __s1, const char *__s2, size_t __n);
char *
        strtok_r(char * __str, const char *__sep,
        char * *__lasts);
int strerror_r(int __errnum, char * __strerrbuf,
        size_t __buflen);
char * strdup(const char *__s1);
void *
        memccpy(void * __dst, const void * __src,
        int __c, size_t __n);
char *
        stpcpy(char * __dst, const char *__src) ;
char *
        stpncpy(char * __dst,
        const char * __src, size_t __n)
        __attribute__((availability(macosx,introduced=10.7)))
                                        ;
char * strndup(const char * __s1, size_t __n) __attribute__((availability(macosx,introduced=10.7)));
size_t strnlen(const char * __s1, size_t __n) __attribute__((availability(macosx,introduced=10.7)));
char * strsignal(int __sig);
typedef __darwin_size_t rsize_t;
typedef int errno_t;
errno_t memset_s(void * __s, rsize_t __smax, int __c, rsize_t __n) __attribute__((availability(macosx,introduced=10.9)));
void *
        memmem(const void * __big, size_t __big_len,
        const void * __little, size_t __little_len) __attribute__((availability(macosx,introduced=10.7)));
void memset_pattern4(void * __b, const void * __pattern4, size_t __len) __attribute__((availability(macosx,introduced=10.5)));
void memset_pattern8(void * __b, const void * __pattern8, size_t __len) __attribute__((availability(macosx,introduced=10.5)));
void memset_pattern16(void * __b, const void * __pattern16, size_t __len) __attribute__((availability(macosx,introduced=10.5)));
char *
        strcasestr(const char *__big, const char *__little);
__attribute__((availability(macosx,introduced=15.4))) __attribute__((availability(ios,introduced=18.4)))
__attribute__((availability(tvos,introduced=18.4))) __attribute__((availability(watchos,introduced=11.4)))
char *
        strchrnul(const char *__s, int __c);
char *
        strnstr(const char * __big, const char *__little, size_t __len);
size_t strlcat(char * __dst, const char *__source, size_t __size);
size_t strlcpy(char * __dst, const char *__source, size_t __size);
void strmode(int __mode, char * __bp);
char *
        strsep(char * *__stringp, const char *__delim);
void swab(const void * restrict, void * restrict, ssize_t __len);
__attribute__((availability(macosx,introduced=10.12.1))) __attribute__((availability(ios,introduced=10.1)))
__attribute__((availability(tvos,introduced=10.0.1))) __attribute__((availability(watchos,introduced=3.1)))
int timingsafe_bcmp(const void * __b1, const void * __b2, size_t __len);
__attribute__((availability(macosx,introduced=11.0))) __attribute__((availability(ios,introduced=14.0)))
__attribute__((availability(tvos,introduced=14.0))) __attribute__((availability(watchos,introduced=7.0)))
int strsignal_r(int __sig, char * __strsignalbuf, size_t __buflen);
int bcmp(const void *, const void *, size_t __n) ;
void bcopy(const void *, void *, size_t __n) ;
void bzero(void *, size_t __n) ;
char * index(const char *, int) ;
char * rindex(const char *, int) ;
int ffs(int);
int strcasecmp(const char *, const char *);
int strncasecmp(const char *, const char *, size_t);
int ffsl(long) __attribute__((availability(macosx,introduced=10.5)));
int ffsll(long long) __attribute__((availability(macosx,introduced=10.9)));
int fls(int) __attribute__((availability(macosx,introduced=10.5)));
int flsl(long) __attribute__((availability(macosx,introduced=10.5)));
int flsll(long long) __attribute__((availability(macosx,introduced=10.9)));
 typedef float float_t;
    typedef double double_t;
extern int __math_errhandling(void);
extern int __fpclassifyf(float);
extern int __fpclassifyd(double);
extern int __fpclassifyl(long double);
inline __attribute__ ((__always_inline__)) int __inline_isfinitef(float);
inline __attribute__ ((__always_inline__)) int __inline_isfinited(double);
inline __attribute__ ((__always_inline__)) int __inline_isfinitel(long double);
inline __attribute__ ((__always_inline__)) int __inline_isinff(float);
inline __attribute__ ((__always_inline__)) int __inline_isinfd(double);
inline __attribute__ ((__always_inline__)) int __inline_isinfl(long double);
inline __attribute__ ((__always_inline__)) int __inline_isnanf(float);
inline __attribute__ ((__always_inline__)) int __inline_isnand(double);
inline __attribute__ ((__always_inline__)) int __inline_isnanl(long double);
inline __attribute__ ((__always_inline__)) int __inline_isnormalf(float);
inline __attribute__ ((__always_inline__)) int __inline_isnormald(double);
inline __attribute__ ((__always_inline__)) int __inline_isnormall(long double);
inline __attribute__ ((__always_inline__)) int __inline_signbitf(float);
inline __attribute__ ((__always_inline__)) int __inline_signbitd(double);
inline __attribute__ ((__always_inline__)) int __inline_signbitl(long double);
inline __attribute__ ((__always_inline__)) int __inline_isfinitef(float __x) {
    return __x == __x && __builtin_fabsf(__x) != __builtin_inff();
}
inline __attribute__ ((__always_inline__)) int __inline_isfinited(double __x) {
    return __x == __x && __builtin_fabs(__x) != __builtin_inf();
}
inline __attribute__ ((__always_inline__)) int __inline_isfinitel(long double __x) {
    return __x == __x && __builtin_fabsl(__x) != __builtin_infl();
}
inline __attribute__ ((__always_inline__)) int __inline_isinff(float __x) {
    return __builtin_fabsf(__x) == __builtin_inff();
}
inline __attribute__ ((__always_inline__)) int __inline_isinfd(double __x) {
    return __builtin_fabs(__x) == __builtin_inf();
}
inline __attribute__ ((__always_inline__)) int __inline_isinfl(long double __x) {
    return __builtin_fabsl(__x) == __builtin_infl();
}
inline __attribute__ ((__always_inline__)) int __inline_isnanf(float __x) {
    return __x != __x;
}
inline __attribute__ ((__always_inline__)) int __inline_isnand(double __x) {
    return __x != __x;
}
inline __attribute__ ((__always_inline__)) int __inline_isnanl(long double __x) {
    return __x != __x;
}
inline __attribute__ ((__always_inline__)) int __inline_signbitf(float __x) {
    union { float __f; unsigned int __u; } __u;
    __u.__f = __x;
    return (int)(__u.__u >> 31);
}
inline __attribute__ ((__always_inline__)) int __inline_signbitd(double __x) {
    union { double __f; unsigned long long __u; } __u;
    __u.__f = __x;
    return (int)(__u.__u >> 63);
}
inline __attribute__ ((__always_inline__)) int __inline_signbitl(long double __x) {
    union { long double __f; unsigned long long __u;} __u;
    __u.__f = __x;
    return (int)(__u.__u >> 63);
}
inline __attribute__ ((__always_inline__)) int __inline_isnormalf(float __x) {
    return __inline_isfinitef(__x) && __builtin_fabsf(__x) >= 1.17549435e-38F;
}
inline __attribute__ ((__always_inline__)) int __inline_isnormald(double __x) {
    return __inline_isfinited(__x) && __builtin_fabs(__x) >= 2.2250738585072014e-308;
}
inline __attribute__ ((__always_inline__)) int __inline_isnormall(long double __x) {
    return __inline_isfinitel(__x) && __builtin_fabsl(__x) >= 2.2250738585072014e-308L;
}
extern float acosf(float);
extern double acos(double);
extern long double acosl(long double);
extern float asinf(float);
extern double asin(double);
extern long double asinl(long double);
extern float atanf(float);
extern double atan(double);
extern long double atanl(long double);
extern float atan2f(float, float);
extern double atan2(double, double);
extern long double atan2l(long double, long double);
extern float cosf(float);
extern double cos(double);
extern long double cosl(long double);
extern float sinf(float);
extern double sin(double);
extern long double sinl(long double);
extern float tanf(float);
extern double tan(double);
extern long double tanl(long double);
extern float acoshf(float);
extern double acosh(double);
extern long double acoshl(long double);
extern float asinhf(float);
extern double asinh(double);
extern long double asinhl(long double);
extern float atanhf(float);
extern double atanh(double);
extern long double atanhl(long double);
extern float coshf(float);
extern double cosh(double);
extern long double coshl(long double);
extern float sinhf(float);
extern double sinh(double);
extern long double sinhl(long double);
extern float tanhf(float);
extern double tanh(double);
extern long double tanhl(long double);
extern float expf(float);
extern double exp(double);
extern long double expl(long double);
extern float exp2f(float);
extern double exp2(double);
extern long double exp2l(long double);
extern float expm1f(float);
extern double expm1(double);
extern long double expm1l(long double);
extern float logf(float);
extern double log(double);
extern long double logl(long double);
extern float log10f(float);
extern double log10(double);
extern long double log10l(long double);
extern float log2f(float);
extern double log2(double);
extern long double log2l(long double);
extern float log1pf(float);
extern double log1p(double);
extern long double log1pl(long double);
extern float logbf(float);
extern double logb(double);
extern long double logbl(long double);
extern float modff(float, float *);
extern double modf(double, double *);
extern long double modfl(long double, long double *);
extern float ldexpf(float, int);
extern double ldexp(double, int);
extern long double ldexpl(long double, int);
extern float frexpf(float, int *);
extern double frexp(double, int *);
extern long double frexpl(long double, int *);
extern int ilogbf(float);
extern int ilogb(double);
extern int ilogbl(long double);
extern float scalbnf(float, int);
extern double scalbn(double, int);
extern long double scalbnl(long double, int);
extern float scalblnf(float, long int);
extern double scalbln(double, long int);
extern long double scalblnl(long double, long int);
extern float fabsf(float);
extern double fabs(double);
extern long double fabsl(long double);
extern float cbrtf(float);
extern double cbrt(double);
extern long double cbrtl(long double);
extern float hypotf(float, float);
extern double hypot(double, double);
extern long double hypotl(long double, long double);
extern float powf(float, float);
extern double pow(double, double);
extern long double powl(long double, long double);
extern float sqrtf(float);
extern double sqrt(double);
extern long double sqrtl(long double);
extern float erff(float);
extern double erf(double);
extern long double erfl(long double);
extern float erfcf(float);
extern double erfc(double);
extern long double erfcl(long double);
extern float lgammaf(float);
extern double lgamma(double);
extern long double lgammal(long double);
extern float tgammaf(float);
extern double tgamma(double);
extern long double tgammal(long double);
extern float ceilf(float);
extern double ceil(double);
extern long double ceill(long double);
extern float floorf(float);
extern double floor(double);
extern long double floorl(long double);
extern float nearbyintf(float);
extern double nearbyint(double);
extern long double nearbyintl(long double);
extern float rintf(float);
extern double rint(double);
extern long double rintl(long double);
extern long int lrintf(float);
extern long int lrint(double);
extern long int lrintl(long double);
extern float roundf(float);
extern double round(double);
extern long double roundl(long double);
extern long int lroundf(float);
extern long int lround(double);
extern long int lroundl(long double);
extern long long int llrintf(float);
extern long long int llrint(double);
extern long long int llrintl(long double);
extern long long int llroundf(float);
extern long long int llround(double);
extern long long int llroundl(long double);
extern float truncf(float);
extern double trunc(double);
extern long double truncl(long double);
extern float fmodf(float, float);
extern double fmod(double, double);
extern long double fmodl(long double, long double);
extern float remainderf(float, float);
extern double remainder(double, double);
extern long double remainderl(long double, long double);
extern float remquof(float, float, int *);
extern double remquo(double, double, int *);
extern long double remquol(long double, long double, int *);
extern float copysignf(float, float);
extern double copysign(double, double);
extern long double copysignl(long double, long double);
extern float nanf(const char *);
extern double nan(const char *);
extern long double nanl(const char *);
extern float nextafterf(float, float);
extern double nextafter(double, double);
extern long double nextafterl(long double, long double);
extern double nexttoward(double, long double);
extern float nexttowardf(float, long double);
extern long double nexttowardl(long double, long double);
extern float fdimf(float, float);
extern double fdim(double, double);
extern long double fdiml(long double, long double);
extern float fmaxf(float, float);
extern double fmax(double, double);
extern long double fmaxl(long double, long double);
extern float fminf(float, float);
extern double fmin(double, double);
extern long double fminl(long double, long double);
extern float fmaf(float, float, float);
extern double fma(double, double, double);
extern long double fmal(long double, long double, long double);
extern float __exp10f(float) __attribute__((availability(macos,introduced=10.9))) __attribute__((availability(ios,introduced=7.0)));
extern double __exp10(double) __attribute__((availability(macos,introduced=10.9))) __attribute__((availability(ios,introduced=7.0)));
inline __attribute__ ((__always_inline__)) void __sincosf(float __x, float *__sinp, float *__cosp);
inline __attribute__ ((__always_inline__)) void __sincos(double __x, double *__sinp, double *__cosp);
extern float __cospif(float) __attribute__((availability(macos,introduced=10.9))) __attribute__((availability(ios,introduced=7.0)));
extern double __cospi(double) __attribute__((availability(macos,introduced=10.9))) __attribute__((availability(ios,introduced=7.0)));
extern float __sinpif(float) __attribute__((availability(macos,introduced=10.9))) __attribute__((availability(ios,introduced=7.0)));
extern double __sinpi(double) __attribute__((availability(macos,introduced=10.9))) __attribute__((availability(ios,introduced=7.0)));
extern float __tanpif(float) __attribute__((availability(macos,introduced=10.9))) __attribute__((availability(ios,introduced=7.0)));
extern double __tanpi(double) __attribute__((availability(macos,introduced=10.9))) __attribute__((availability(ios,introduced=7.0)));
extern _Float16 __fabsf16(_Float16) __attribute__((availability(macos,introduced=15.0))) __attribute__((availability(ios,introduced=18.0))) __attribute__((availability(watchos,introduced=11.0))) __attribute__((availability(tvos,introduced=18.0)));
extern _Float16 __hypotf16(_Float16, _Float16) __attribute__((availability(macos,introduced=15.0))) __attribute__((availability(ios,introduced=18.0))) __attribute__((availability(watchos,introduced=11.0))) __attribute__((availability(tvos,introduced=18.0)));
extern _Float16 __sqrtf16(_Float16) __attribute__((availability(macos,introduced=15.0))) __attribute__((availability(ios,introduced=18.0))) __attribute__((availability(watchos,introduced=11.0))) __attribute__((availability(tvos,introduced=18.0)));
extern _Float16 __ceilf16(_Float16) __attribute__((availability(macos,introduced=15.0))) __attribute__((availability(ios,introduced=18.0))) __attribute__((availability(watchos,introduced=11.0))) __attribute__((availability(tvos,introduced=18.0)));
extern _Float16 __floorf16(_Float16) __attribute__((availability(macos,introduced=15.0))) __attribute__((availability(ios,introduced=18.0))) __attribute__((availability(watchos,introduced=11.0))) __attribute__((availability(tvos,introduced=18.0)));
extern _Float16 __rintf16(_Float16) __attribute__((availability(macos,introduced=15.0))) __attribute__((availability(ios,introduced=18.0))) __attribute__((availability(watchos,introduced=11.0))) __attribute__((availability(tvos,introduced=18.0)));
extern _Float16 __roundf16(_Float16) __attribute__((availability(macos,introduced=15.0))) __attribute__((availability(ios,introduced=18.0))) __attribute__((availability(watchos,introduced=11.0))) __attribute__((availability(tvos,introduced=18.0)));
extern _Float16 __truncf16(_Float16) __attribute__((availability(macos,introduced=15.0))) __attribute__((availability(ios,introduced=18.0))) __attribute__((availability(watchos,introduced=11.0))) __attribute__((availability(tvos,introduced=18.0)));
extern _Float16 __copysignf16(_Float16, _Float16) __attribute__((availability(macos,introduced=15.0))) __attribute__((availability(ios,introduced=18.0))) __attribute__((availability(watchos,introduced=11.0))) __attribute__((availability(tvos,introduced=18.0)));
extern _Float16 __nextafterf16(_Float16, _Float16) __attribute__((availability(macos,introduced=15.0))) __attribute__((availability(ios,introduced=18.0))) __attribute__((availability(watchos,introduced=11.0))) __attribute__((availability(tvos,introduced=18.0)));
extern _Float16 __fmaxf16(_Float16, _Float16) __attribute__((availability(macos,introduced=15.0))) __attribute__((availability(ios,introduced=18.0))) __attribute__((availability(watchos,introduced=11.0))) __attribute__((availability(tvos,introduced=18.0)));
extern _Float16 __fminf16(_Float16, _Float16) __attribute__((availability(macos,introduced=15.0))) __attribute__((availability(ios,introduced=18.0))) __attribute__((availability(watchos,introduced=11.0))) __attribute__((availability(tvos,introduced=18.0)));
extern _Float16 __fmaf16(_Float16, _Float16, _Float16) __attribute__((availability(macos,introduced=15.0))) __attribute__((availability(ios,introduced=18.0))) __attribute__((availability(watchos,introduced=11.0))) __attribute__((availability(tvos,introduced=18.0)));
inline __attribute__ ((__always_inline__)) void __sincospif(float __x, float *__sinp, float *__cosp);
inline __attribute__ ((__always_inline__)) void __sincospi(double __x, double *__sinp, double *__cosp);
struct __float2 { float __sinval; float __cosval; };
struct __double2 { double __sinval; double __cosval; };
extern struct __float2 __sincosf_stret(float);
extern struct __double2 __sincos_stret(double);
extern struct __float2 __sincospif_stret(float);
extern struct __double2 __sincospi_stret(double);
inline __attribute__ ((__always_inline__)) void __sincosf(float __x, float *__sinp, float *__cosp) {
    const struct __float2 __stret = __sincosf_stret(__x);
    *__sinp = __stret.__sinval; *__cosp = __stret.__cosval;
}
inline __attribute__ ((__always_inline__)) void __sincos(double __x, double *__sinp, double *__cosp) {
    const struct __double2 __stret = __sincos_stret(__x);
    *__sinp = __stret.__sinval; *__cosp = __stret.__cosval;
}
inline __attribute__ ((__always_inline__)) void __sincospif(float __x, float *__sinp, float *__cosp) {
    const struct __float2 __stret = __sincospif_stret(__x);
    *__sinp = __stret.__sinval; *__cosp = __stret.__cosval;
}
inline __attribute__ ((__always_inline__)) void __sincospi(double __x, double *__sinp, double *__cosp) {
    const struct __double2 __stret = __sincospi_stret(__x);
    *__sinp = __stret.__sinval; *__cosp = __stret.__cosval;
}
extern double j0(double) __attribute__((availability(macos,introduced=10.0))) __attribute__((availability(ios,introduced=3.2)));
extern double j1(double) __attribute__((availability(macos,introduced=10.0))) __attribute__((availability(ios,introduced=3.2)));
extern double jn(int, double) __attribute__((availability(macos,introduced=10.0))) __attribute__((availability(ios,introduced=3.2)));
extern double y0(double) __attribute__((availability(macos,introduced=10.0))) __attribute__((availability(ios,introduced=3.2)));
extern double y1(double) __attribute__((availability(macos,introduced=10.0))) __attribute__((availability(ios,introduced=3.2)));
extern double yn(int, double) __attribute__((availability(macos,introduced=10.0))) __attribute__((availability(ios,introduced=3.2)));
extern double scalb(double, double);
extern int signgam;
struct exception {
    int type;
    char *name;
    double arg1;
    double arg2;
    double retval;
};
typedef int MPI_Datatype;
typedef int MPI_Comm;
typedef int MPI_Group;
typedef int MPI_Win;
typedef int MPI_Session;
typedef struct ADIOI_FileD *MPI_File;
typedef int MPI_Op;
typedef int MPI_Errhandler;
typedef int MPI_Request;
typedef int MPI_Message;
typedef int MPIX_Grequest_class;
typedef int MPIX_Stream;
typedef int MPI_Info;
typedef long MPI_Aint;
typedef int MPI_Fint;
typedef long long MPI_Count;
typedef long long MPI_Offset;
typedef struct MPI_Status {
    int count_lo;
    int count_hi_and_cancelled;
    int MPI_SOURCE;
    int MPI_TAG;
    int MPI_ERROR;
} MPI_Status;
extern MPI_Fint * MPI_F_STATUS_IGNORE ;
extern MPI_Fint * MPI_F_STATUSES_IGNORE ;
typedef struct {
    MPI_Fint count_lo;
    MPI_Fint count_hi_and_cancelled;
    MPI_Fint MPI_SOURCE;
    MPI_Fint MPI_TAG;
    MPI_Fint MPI_ERROR;
} MPI_F08_status;
extern MPI_F08_status *MPI_F08_STATUS_IGNORE ;
extern MPI_F08_status *MPI_F08_STATUSES_IGNORE ;
typedef enum MPIR_Win_flavor {
    MPI_WIN_FLAVOR_CREATE = 1,
    MPI_WIN_FLAVOR_ALLOCATE = 2,
    MPI_WIN_FLAVOR_DYNAMIC = 3,
    MPI_WIN_FLAVOR_SHARED = 4
} MPIR_Win_flavor_t;
typedef enum MPIR_Win_model {
    MPI_WIN_SEPARATE = 1,
    MPI_WIN_UNIFIED = 2
} MPIR_Win_model_t;
typedef enum MPIR_Topo_type {
    MPI_GRAPH = 1,
    MPI_CART = 2,
    MPI_DIST_GRAPH = 3
} MPIR_Topo_type;
extern int * const MPI_UNWEIGHTED ;
extern int * const MPI_WEIGHTS_EMPTY ;
enum MPIR_Combiner_enum {
    MPI_COMBINER_NAMED = 1,
    MPI_COMBINER_DUP = 2,
    MPI_COMBINER_CONTIGUOUS = 3,
    MPI_COMBINER_VECTOR = 4,
    MPI_COMBINER_HVECTOR_INTEGER = 5,
    MPI_COMBINER_HVECTOR = 6,
    MPI_COMBINER_INDEXED = 7,
    MPI_COMBINER_HINDEXED_INTEGER = 8,
    MPI_COMBINER_HINDEXED = 9,
    MPI_COMBINER_INDEXED_BLOCK = 10,
    MPI_COMBINER_STRUCT_INTEGER = 11,
    MPI_COMBINER_STRUCT = 12,
    MPI_COMBINER_SUBARRAY = 13,
    MPI_COMBINER_DARRAY = 14,
    MPI_COMBINER_F90_REAL = 15,
    MPI_COMBINER_F90_COMPLEX = 16,
    MPI_COMBINER_F90_INTEGER = 17,
    MPI_COMBINER_RESIZED = 18,
    MPI_COMBINER_HINDEXED_BLOCK = 19,
    MPI_COMBINER_VALUE_INDEX = 20
};
enum {
    MPIX_ASYNC_NOPROGRESS = 0,
    MPIX_ASYNC_DONE = 1,
};
typedef struct MPIR_Async_thing * MPIX_Async_thing;
typedef void (MPI_Handler_function) ( MPI_Comm *, int *, ... );
typedef int (MPI_Comm_copy_attr_function)(MPI_Comm, int, void *, void *,
       void *, int *);
typedef int (MPI_Comm_delete_attr_function)(MPI_Comm, int, void *, void *);
typedef int (MPI_Type_copy_attr_function)(MPI_Datatype, int, void *, void *,
       void *, int *);
typedef int (MPI_Type_delete_attr_function)(MPI_Datatype, int, void *, void *);
typedef int (MPI_Win_copy_attr_function)(MPI_Win, int, void *, void *, void *,
      int *);
typedef int (MPI_Win_delete_attr_function)(MPI_Win, int, void *, void *);
typedef void (MPI_Comm_errhandler_function)(MPI_Comm *, int *, ...);
typedef void (MPI_File_errhandler_function)(MPI_File *, int *, ...);
typedef void (MPI_Win_errhandler_function)(MPI_Win *, int *, ...);
typedef void (MPI_Session_errhandler_function)(MPI_Session *, int *, ...);
typedef void (MPIX_Comm_errhandler_function_x)(MPI_Comm, int, void *);
typedef void (MPIX_File_errhandler_function_x)(MPI_File, int, void *);
typedef void (MPIX_Win_errhandler_function_x)(MPI_Win, int, void *);
typedef void (MPIX_Session_errhandler_function_x)(MPI_Session, int, void *);
typedef MPI_Comm_errhandler_function MPI_Comm_errhandler_fn;
typedef MPI_File_errhandler_function MPI_File_errhandler_fn;
typedef MPI_Win_errhandler_function MPI_Win_errhandler_fn;
typedef MPI_Session_errhandler_function MPI_Session_errhandler_fn;
typedef int (MPI_Copy_function) ( MPI_Comm, int, void *, void *, void *, int * );
typedef int (MPI_Delete_function) ( MPI_Comm, int, void *, void * );
typedef void (MPI_User_function) ( void *, void *, int *, MPI_Datatype * );
typedef void (MPI_User_function_c) ( void *, void *, MPI_Count *, MPI_Datatype * );
typedef void (MPIX_User_function_x) ( void *, void *, MPI_Count, MPI_Datatype, void *);
typedef void (MPIX_Destructor_function) (void *);
typedef int (MPI_Grequest_cancel_function)(void *, int);
typedef int (MPI_Grequest_free_function)(void *);
typedef int (MPI_Grequest_query_function)(void *, MPI_Status *);
typedef int (MPIX_Grequest_poll_function)(void *, MPI_Status *);
typedef int (MPIX_Grequest_wait_function)(int, void **, double, MPI_Status *);
typedef int (MPIX_Async_poll_function)(MPIX_Async_thing);
typedef int (MPI_Datarep_conversion_function)(void *, MPI_Datatype, int,
             void *, MPI_Offset, void *);
typedef int (MPI_Datarep_extent_function)(MPI_Datatype datatype, MPI_Aint *,
                      void *);
typedef int (MPI_Datarep_conversion_function_c)(void *, MPI_Datatype, MPI_Count,
             void *, MPI_Offset, void *);
struct MPIR_T_enum_s;
struct MPIR_T_cvar_handle_s;
struct MPIR_T_pvar_handle_s;
struct MPIR_T_pvar_session_s;
struct MPIR_T_event_registration_s;
struct MPIR_T_event_instance_s;
typedef struct MPIR_T_enum_s * MPI_T_enum;
typedef struct MPIR_T_cvar_handle_s * MPI_T_cvar_handle;
typedef struct MPIR_T_pvar_handle_s * MPI_T_pvar_handle;
typedef struct MPIR_T_pvar_session_s * MPI_T_pvar_session;
typedef struct MPIR_T_event_registration_s * MPI_T_event_registration;
typedef struct MPIR_T_event_instance_s * MPI_T_event_instance;
extern struct MPIR_T_pvar_handle_s * const MPI_T_PVAR_ALL_HANDLES ;
typedef enum MPIR_T_verbosity_t {
    MPI_T_VERBOSITY_INVALID = 0,
    MPI_T_VERBOSITY_USER_BASIC = 221,
    MPI_T_VERBOSITY_USER_DETAIL,
    MPI_T_VERBOSITY_USER_ALL,
    MPI_T_VERBOSITY_TUNER_BASIC,
    MPI_T_VERBOSITY_TUNER_DETAIL,
    MPI_T_VERBOSITY_TUNER_ALL,
    MPI_T_VERBOSITY_MPIDEV_BASIC,
    MPI_T_VERBOSITY_MPIDEV_DETAIL,
    MPI_T_VERBOSITY_MPIDEV_ALL
} MPIR_T_verbosity_t;
typedef enum MPIR_T_bind_t {
    MPI_T_BIND_INVALID = 0,
    MPI_T_BIND_NO_OBJECT = 9700,
    MPI_T_BIND_MPI_COMM,
    MPI_T_BIND_MPI_DATATYPE,
    MPI_T_BIND_MPI_ERRHANDLER,
    MPI_T_BIND_MPI_FILE,
    MPI_T_BIND_MPI_GROUP,
    MPI_T_BIND_MPI_OP,
    MPI_T_BIND_MPI_REQUEST,
    MPI_T_BIND_MPI_WIN,
    MPI_T_BIND_MPI_MESSAGE,
    MPI_T_BIND_MPI_INFO,
    MPI_T_BIND_MPI_SESSION
} MPIR_T_bind_t;
typedef enum MPIR_T_scope_t {
    MPI_T_SCOPE_INVALID = 0,
    MPI_T_SCOPE_CONSTANT = 60438,
    MPI_T_SCOPE_READONLY,
    MPI_T_SCOPE_LOCAL,
    MPI_T_SCOPE_GROUP,
    MPI_T_SCOPE_GROUP_EQ,
    MPI_T_SCOPE_ALL,
    MPI_T_SCOPE_ALL_EQ
} MPIR_T_scope_t;
typedef enum MPIR_T_pvar_class_t {
    MPI_T_PVAR_CLASS_INVALID = 0,
    MPIR_T_PVAR_CLASS_FIRST = 240,
    MPI_T_PVAR_CLASS_STATE = MPIR_T_PVAR_CLASS_FIRST,
    MPI_T_PVAR_CLASS_LEVEL,
    MPI_T_PVAR_CLASS_SIZE,
    MPI_T_PVAR_CLASS_PERCENTAGE,
    MPI_T_PVAR_CLASS_HIGHWATERMARK,
    MPI_T_PVAR_CLASS_LOWWATERMARK,
    MPI_T_PVAR_CLASS_COUNTER,
    MPI_T_PVAR_CLASS_AGGREGATE,
    MPI_T_PVAR_CLASS_TIMER,
    MPI_T_PVAR_CLASS_GENERIC,
    MPIR_T_PVAR_CLASS_LAST,
    MPIR_T_PVAR_CLASS_NUMBER = MPIR_T_PVAR_CLASS_LAST - MPIR_T_PVAR_CLASS_FIRST
} MPIR_T_pvar_class_t;
typedef enum MPI_T_cb_safety {
    MPI_T_CB_REQUIRE_NONE = 0,
    MPI_T_CB_REQUIRE_MPI_RESTRICTED,
    MPI_T_CB_REQUIRE_THREAD_SAFE,
    MPI_T_CB_REQUIRE_ASYNC_SIGNAL_SAFE
} MPI_T_cb_safety;
typedef enum MPI_T_source_order {
    MPI_T_SOURCE_ORDERED = 0,
    MPI_T_SOURCE_UNORDERED
} MPI_T_source_order;
typedef void (MPI_T_event_cb_function)(MPI_T_event_instance event_instance, MPI_T_event_registration event_registration, MPI_T_cb_safety cb_safety, void *user_data);
typedef void (MPI_T_event_free_cb_function)(MPI_T_event_registration event_registration, MPI_T_cb_safety cb_safety, void *user_data);
typedef void (MPI_T_event_dropped_cb_function)(MPI_Count count, MPI_T_event_registration event_registration, int source_index, MPI_T_cb_safety cb_safety, void *user_data);
typedef struct {
    void **storage_stack;
} QMPI_Context;
int QMPI_Register_tool_name(const char *tool_name,
                            void (*init_function_ptr) (int tool_id)) ;
int QMPI_Register_tool_storage(int tool_id, void *tool_storage) ;
int QMPI_Register_function(int calling_tool_id, int function_enum,
                           void (*function_ptr) (void)) ;
int QMPI_Get_function(int calling_tool_id, int function_enum,
                      void (**function_ptr) (void), int *next_tool_id) ;
int QMPI_Get_tool_storage(QMPI_Context context, int tool_id, void **storage) ;
int QMPI_Get_calling_address(QMPI_Context context, void **address) ;
typedef struct MPIX_Iov {
    void *iov_base;
    MPI_Aint iov_len;
} MPIX_Iov;
int MPIR_Dup_fn(MPI_Comm oldcomm, int keyval, void *extra_state, void *attribute_val_in,
               void *attribute_val_out, int *flag) ;
int MPI_Abi_get_fortran_booleans(int logical_size, void *logical_true, void *logical_false,
                                 int *is_set) ;
int MPI_Abi_get_fortran_info(MPI_Info *info) ;
int MPI_Abi_get_info(MPI_Info *info) ;
int MPI_Abi_get_version(int *abi_major, int *abi_minor) ;
int MPI_Abi_set_fortran_booleans(int logical_size, void *logical_true, void *logical_false)
                    ;
int MPI_Abi_set_fortran_info(MPI_Info info) ;
int MPI_Comm_toint(MPI_Comm comm) ;
MPI_Comm MPI_Comm_fromint(int comm) ;
int MPI_Errhandler_toint(MPI_Errhandler errhandler) ;
MPI_Errhandler MPI_Errhandler_fromint(int errhandler) ;
int MPI_Group_toint(MPI_Group group) ;
MPI_Group MPI_Group_fromint(int group) ;
int MPI_Info_toint(MPI_Info info) ;
MPI_Info MPI_Info_fromint(int info) ;
int MPI_Message_toint(MPI_Message message) ;
MPI_Message MPI_Message_fromint(int message) ;
int MPI_Op_toint(MPI_Op op) ;
MPI_Op MPI_Op_fromint(int op) ;
int MPI_Request_toint(MPI_Request request) ;
MPI_Request MPI_Request_fromint(int request) ;
int MPI_Session_toint(MPI_Session session) ;
MPI_Session MPI_Session_fromint(int session) ;
int MPI_Type_toint(MPI_Datatype datatype) ;
MPI_Datatype MPI_Type_fromint(int datatype) ;
int MPI_Win_toint(MPI_Win win) ;
MPI_Win MPI_Win_fromint(int win) ;
int MPIX_Async_start(MPIX_Async_poll_function *poll_fn, void *extra_state, MPIX_Stream stream)
                    ;
void * MPIX_Async_get_state(MPIX_Async_thing async_thing) ;
int MPIX_Async_spawn(MPIX_Async_thing async_thing, MPIX_Async_poll_function *poll_fn,
                     void *extra_state, MPIX_Stream stream) ;
int MPI_Comm_create_keyval(MPI_Comm_copy_attr_function *comm_copy_attr_fn,
                           MPI_Comm_delete_attr_function *comm_delete_attr_fn, int *comm_keyval,
                           void *extra_state) ;
int MPI_Keyval_create(MPI_Copy_function *copy_fn, MPI_Delete_function *delete_fn, int *keyval,
                      void *extra_state) ;
int MPI_Comm_delete_attr(MPI_Comm comm, int comm_keyval) ;
int MPI_Attr_delete(MPI_Comm comm, int keyval) ;
int MPI_Comm_free_keyval(int *comm_keyval) ;
int MPI_Keyval_free(int *keyval) ;
int MPI_Comm_get_attr(MPI_Comm comm, int comm_keyval, void *attribute_val, int *flag)
                    ;
int MPI_Attr_get(MPI_Comm comm, int keyval, void *attribute_val, int *flag) ;
int MPI_Comm_set_attr(MPI_Comm comm, int comm_keyval, void *attribute_val) ;
int MPI_Attr_put(MPI_Comm comm, int keyval, void *attribute_val) ;
int MPI_Type_create_keyval(MPI_Type_copy_attr_function *type_copy_attr_fn,
                           MPI_Type_delete_attr_function *type_delete_attr_fn, int *type_keyval,
                           void *extra_state) ;
int MPI_Type_delete_attr(MPI_Datatype datatype, int type_keyval) ;
int MPI_Type_free_keyval(int *type_keyval) ;
int MPI_Type_get_attr(MPI_Datatype datatype, int type_keyval, void *attribute_val, int *flag)
                    ;
int MPI_Type_set_attr(MPI_Datatype datatype, int type_keyval, void *attribute_val)
                    ;
int MPI_Win_create_keyval(MPI_Win_copy_attr_function *win_copy_attr_fn,
                          MPI_Win_delete_attr_function *win_delete_attr_fn, int *win_keyval,
                          void *extra_state) ;
int MPI_Win_delete_attr(MPI_Win win, int win_keyval) ;
int MPI_Win_free_keyval(int *win_keyval) ;
int MPI_Win_get_attr(MPI_Win win, int win_keyval, void *attribute_val, int *flag) ;
int MPI_Win_set_attr(MPI_Win win, int win_keyval, void *attribute_val) ;
int MPIX_Op_create_x(MPIX_User_function_x *user_fn_x, MPIX_Destructor_function *destructor_fn,
                     int commute, void *extra_state, MPI_Op *op) ;
int MPIX_Comm_create_errhandler_x(MPIX_Comm_errhandler_function_x *comm_errhandler_fn_x,
                                  MPIX_Destructor_function *destructor_fn, void *extra_state,
                                  MPI_Errhandler *errhandler) ;
int MPIX_Win_create_errhandler_x(MPIX_Win_errhandler_function_x *comm_errhandler_fn_x,
                                 MPIX_Destructor_function *destructor_fn, void *extra_state,
                                 MPI_Errhandler *errhandler) ;
int MPIX_File_create_errhandler_x(MPIX_File_errhandler_function_x *comm_errhandler_fn_x,
                                  MPIX_Destructor_function *destructor_fn, void *extra_state,
                                  MPI_Errhandler *errhandler) ;
int MPIX_Session_create_errhandler_x(MPIX_Session_errhandler_function_x *comm_errhandler_fn_x,
                                     MPIX_Destructor_function *destructor_fn, void *extra_state,
                                     MPI_Errhandler *errhandler) ;
int MPI_Allgather(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
                  int recvcount, MPI_Datatype recvtype, MPI_Comm comm)
                                                                                                              ;
int MPI_Allgather_init(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
                       int recvcount, MPI_Datatype recvtype, MPI_Comm comm, MPI_Info info,
                       MPI_Request *request)
                                                                                                                   ;
int MPI_Allgatherv(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
                   const int recvcounts[], const int displs[], MPI_Datatype recvtype,
                   MPI_Comm comm)
                                                                                                               ;
int MPI_Allgatherv_init(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
                        const int recvcounts[], const int displs[], MPI_Datatype recvtype,
                        MPI_Comm comm, MPI_Info info, MPI_Request *request)
                                                                                                                    ;
int MPI_Allreduce(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op,
                  MPI_Comm comm)
                                                                                                              ;
int MPI_Allreduce_init(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype,
                       MPI_Op op, MPI_Comm comm, MPI_Info info, MPI_Request *request)
                                                                                                                   ;
int MPI_Alltoall(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
                 int recvcount, MPI_Datatype recvtype, MPI_Comm comm)
                                                                                                             ;
int MPI_Alltoall_init(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
                      int recvcount, MPI_Datatype recvtype, MPI_Comm comm, MPI_Info info,
                      MPI_Request *request)
                                                                                                                  ;
int MPI_Alltoallv(const void *sendbuf, const int sendcounts[], const int sdispls[],
                  MPI_Datatype sendtype, void *recvbuf, const int recvcounts[], const int rdispls[],
                  MPI_Datatype recvtype, MPI_Comm comm)
                                                                                                              ;
int MPI_Alltoallv_init(const void *sendbuf, const int sendcounts[], const int sdispls[],
                       MPI_Datatype sendtype, void *recvbuf, const int recvcounts[],
                       const int rdispls[], MPI_Datatype recvtype, MPI_Comm comm, MPI_Info info,
                       MPI_Request *request)
                                                                                                                   ;
int MPI_Alltoallw(const void *sendbuf, const int sendcounts[], const int sdispls[],
                  const MPI_Datatype sendtypes[], void *recvbuf, const int recvcounts[],
                  const int rdispls[], const MPI_Datatype recvtypes[], MPI_Comm comm)
                                  ;
int MPI_Alltoallw_init(const void *sendbuf, const int sendcounts[], const int sdispls[],
                       const MPI_Datatype sendtypes[], void *recvbuf, const int recvcounts[],
                       const int rdispls[], const MPI_Datatype recvtypes[], MPI_Comm comm,
                       MPI_Info info, MPI_Request *request) ;
int MPI_Barrier(MPI_Comm comm) ;
int MPI_Barrier_init(MPI_Comm comm, MPI_Info info, MPI_Request *request) ;
int MPI_Bcast(void *buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm)
                                                          ;
int MPI_Bcast_init(void *buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm,
                   MPI_Info info, MPI_Request *request)
                                                                         ;
int MPI_Exscan(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op,
               MPI_Comm comm)
                                                                                                           ;
int MPI_Exscan_init(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op,
                    MPI_Comm comm, MPI_Info info, MPI_Request *request)
                                                                                                                ;
int MPI_Gather(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
               int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm)
                                                                                                           ;
int MPI_Gather_init(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
                    int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm, MPI_Info info,
                    MPI_Request *request)
                                                                                                                ;
int MPI_Gatherv(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
                const int recvcounts[], const int displs[], MPI_Datatype recvtype, int root,
                MPI_Comm comm)
                                                                                                            ;
int MPI_Gatherv_init(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
                     const int recvcounts[], const int displs[], MPI_Datatype recvtype, int root,
                     MPI_Comm comm, MPI_Info info, MPI_Request *request)
                                                                                                                 ;
int MPI_Iallgather(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
                   int recvcount, MPI_Datatype recvtype, MPI_Comm comm, MPI_Request *request)
                                                                                                               ;
int MPI_Iallgatherv(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
                    const int recvcounts[], const int displs[], MPI_Datatype recvtype,
                    MPI_Comm comm, MPI_Request *request)
                                                                                                                ;
int MPI_Iallreduce(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op,
                   MPI_Comm comm, MPI_Request *request)
                                                                                                               ;
int MPI_Ialltoall(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
                  int recvcount, MPI_Datatype recvtype, MPI_Comm comm, MPI_Request *request)
                                                                                                              ;
int MPI_Ialltoallv(const void *sendbuf, const int sendcounts[], const int sdispls[],
                   MPI_Datatype sendtype, void *recvbuf, const int recvcounts[],
                   const int rdispls[], MPI_Datatype recvtype, MPI_Comm comm, MPI_Request *request)
                                                                                                               ;
int MPI_Ialltoallw(const void *sendbuf, const int sendcounts[], const int sdispls[],
                   const MPI_Datatype sendtypes[], void *recvbuf, const int recvcounts[],
                   const int rdispls[], const MPI_Datatype recvtypes[], MPI_Comm comm,
                   MPI_Request *request) ;
int MPI_Ibarrier(MPI_Comm comm, MPI_Request *request) ;
int MPI_Ibcast(void *buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm,
               MPI_Request *request) ;
int MPI_Iexscan(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op,
                MPI_Comm comm, MPI_Request *request)
                                                                                                            ;
int MPI_Igather(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
                int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm,
                MPI_Request *request)
                                                                                                            ;
int MPI_Igatherv(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
                 const int recvcounts[], const int displs[], MPI_Datatype recvtype, int root,
                 MPI_Comm comm, MPI_Request *request)
                                                                                                             ;
int MPI_Ineighbor_allgather(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
                            void *recvbuf, int recvcount, MPI_Datatype recvtype, MPI_Comm comm,
                            MPI_Request *request)
                                                                                                                        ;
int MPI_Ineighbor_allgatherv(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
                             void *recvbuf, const int recvcounts[], const int displs[],
                             MPI_Datatype recvtype, MPI_Comm comm, MPI_Request *request)
                                                                                                                         ;
int MPI_Ineighbor_alltoall(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
                           int recvcount, MPI_Datatype recvtype, MPI_Comm comm,
                           MPI_Request *request)
                                                                                                                       ;
int MPI_Ineighbor_alltoallv(const void *sendbuf, const int sendcounts[], const int sdispls[],
                            MPI_Datatype sendtype, void *recvbuf, const int recvcounts[],
                            const int rdispls[], MPI_Datatype recvtype, MPI_Comm comm,
                            MPI_Request *request)
                                                                                                                        ;
int MPI_Ineighbor_alltoallw(const void *sendbuf, const int sendcounts[], const MPI_Aint sdispls[],
                            const MPI_Datatype sendtypes[], void *recvbuf, const int recvcounts[],
                            const MPI_Aint rdispls[], const MPI_Datatype recvtypes[], MPI_Comm comm,
                            MPI_Request *request) ;
int MPI_Ireduce(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op,
                int root, MPI_Comm comm, MPI_Request *request)
                                                                                                            ;
int MPI_Ireduce_scatter(const void *sendbuf, void *recvbuf, const int recvcounts[],
                        MPI_Datatype datatype, MPI_Op op, MPI_Comm comm, MPI_Request *request)
                                                                                                                    ;
int MPI_Ireduce_scatter_block(const void *sendbuf, void *recvbuf, int recvcount,
                              MPI_Datatype datatype, MPI_Op op, MPI_Comm comm,
                              MPI_Request *request)
                                                                                                                          ;
int MPI_Iscan(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op,
              MPI_Comm comm, MPI_Request *request)
                                                                                                          ;
int MPI_Iscatter(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
                 int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm,
                 MPI_Request *request)
                                                                                                             ;
int MPI_Iscatterv(const void *sendbuf, const int sendcounts[], const int displs[],
                  MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype,
                  int root, MPI_Comm comm, MPI_Request *request)
                                                                                                              ;
int MPI_Neighbor_allgather(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
                           int recvcount, MPI_Datatype recvtype, MPI_Comm comm)
                                                                                                                       ;
int MPI_Neighbor_allgather_init(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
                                void *recvbuf, int recvcount, MPI_Datatype recvtype, MPI_Comm comm,
                                MPI_Info info, MPI_Request *request)
                                                                                                                            ;
int MPI_Neighbor_allgatherv(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
                            void *recvbuf, const int recvcounts[], const int displs[],
                            MPI_Datatype recvtype, MPI_Comm comm)
                                                                                                                        ;
int MPI_Neighbor_allgatherv_init(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
                                 void *recvbuf, const int recvcounts[], const int displs[],
                                 MPI_Datatype recvtype, MPI_Comm comm, MPI_Info info,
                                 MPI_Request *request)
                                                                                                                             ;
int MPI_Neighbor_alltoall(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
                          int recvcount, MPI_Datatype recvtype, MPI_Comm comm)
                                                                                                                      ;
int MPI_Neighbor_alltoall_init(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
                               void *recvbuf, int recvcount, MPI_Datatype recvtype, MPI_Comm comm,
                               MPI_Info info, MPI_Request *request)
                                                                                                                           ;
int MPI_Neighbor_alltoallv(const void *sendbuf, const int sendcounts[], const int sdispls[],
                           MPI_Datatype sendtype, void *recvbuf, const int recvcounts[],
                           const int rdispls[], MPI_Datatype recvtype, MPI_Comm comm)
                                                                                                                       ;
int MPI_Neighbor_alltoallv_init(const void *sendbuf, const int sendcounts[], const int sdispls[],
                                MPI_Datatype sendtype, void *recvbuf, const int recvcounts[],
                                const int rdispls[], MPI_Datatype recvtype, MPI_Comm comm,
                                MPI_Info info, MPI_Request *request)
                                                                                                                            ;
int MPI_Neighbor_alltoallw(const void *sendbuf, const int sendcounts[], const MPI_Aint sdispls[],
                           const MPI_Datatype sendtypes[], void *recvbuf, const int recvcounts[],
                           const MPI_Aint rdispls[], const MPI_Datatype recvtypes[], MPI_Comm comm)
                                           ;
int MPI_Neighbor_alltoallw_init(const void *sendbuf, const int sendcounts[],
                                const MPI_Aint sdispls[], const MPI_Datatype sendtypes[],
                                void *recvbuf, const int recvcounts[], const MPI_Aint rdispls[],
                                const MPI_Datatype recvtypes[], MPI_Comm comm, MPI_Info info,
                                MPI_Request *request) ;
int MPI_Reduce(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op,
               int root, MPI_Comm comm)
                                                                                                           ;
int MPI_Reduce_init(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op,
                    int root, MPI_Comm comm, MPI_Info info, MPI_Request *request)
                                                                                                                ;
int MPI_Reduce_local(const void *inbuf, void *inoutbuf, int count, MPI_Datatype datatype,
                     MPI_Op op)
                                                                                                                 ;
int MPI_Reduce_scatter(const void *sendbuf, void *recvbuf, const int recvcounts[],
                       MPI_Datatype datatype, MPI_Op op, MPI_Comm comm)
                                                                                                                   ;
int MPI_Reduce_scatter_block(const void *sendbuf, void *recvbuf, int recvcount,
                             MPI_Datatype datatype, MPI_Op op, MPI_Comm comm)
                                                                                                                         ;
int MPI_Reduce_scatter_block_init(const void *sendbuf, void *recvbuf, int recvcount,
                                  MPI_Datatype datatype, MPI_Op op, MPI_Comm comm, MPI_Info info,
                                  MPI_Request *request)
                                                                                                                              ;
int MPI_Reduce_scatter_init(const void *sendbuf, void *recvbuf, const int recvcounts[],
                            MPI_Datatype datatype, MPI_Op op, MPI_Comm comm, MPI_Info info,
                            MPI_Request *request)
                                                                                                                        ;
int MPI_Scan(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op,
             MPI_Comm comm)
                                                                                                         ;
int MPI_Scan_init(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op,
                  MPI_Comm comm, MPI_Info info, MPI_Request *request)
                                                                                                              ;
int MPI_Scatter(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
                int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm)
                                                                                                            ;
int MPI_Scatter_init(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
                     int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm, MPI_Info info,
                     MPI_Request *request)
                                                                                                                 ;
int MPI_Scatterv(const void *sendbuf, const int sendcounts[], const int displs[],
                 MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype,
                 int root, MPI_Comm comm)
                                                                                                             ;
int MPI_Scatterv_init(const void *sendbuf, const int sendcounts[], const int displs[],
                      MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype,
                      int root, MPI_Comm comm, MPI_Info info, MPI_Request *request)
                                                                                                                  ;
int MPI_Comm_compare(MPI_Comm comm1, MPI_Comm comm2, int *result) ;
int MPI_Comm_create(MPI_Comm comm, MPI_Group group, MPI_Comm *newcomm) ;
int MPI_Comm_create_group(MPI_Comm comm, MPI_Group group, int tag, MPI_Comm *newcomm)
                    ;
int MPI_Comm_dup(MPI_Comm comm, MPI_Comm *newcomm) ;
int MPI_Comm_dup_with_info(MPI_Comm comm, MPI_Info info, MPI_Comm *newcomm) ;
int MPI_Comm_free(MPI_Comm *comm) ;
int MPI_Comm_get_info(MPI_Comm comm, MPI_Info *info_used) ;
int MPI_Comm_get_name(MPI_Comm comm, char *comm_name, int *resultlen) ;
int MPI_Comm_group(MPI_Comm comm, MPI_Group *group) ;
int MPI_Comm_idup(MPI_Comm comm, MPI_Comm *newcomm, MPI_Request *request) ;
int MPI_Comm_idup_with_info(MPI_Comm comm, MPI_Info info, MPI_Comm *newcomm, MPI_Request *request)
                    ;
int MPI_Comm_rank(MPI_Comm comm, int *rank) ;
int MPI_Comm_remote_group(MPI_Comm comm, MPI_Group *group) ;
int MPI_Comm_remote_size(MPI_Comm comm, int *size) ;
int MPI_Comm_set_info(MPI_Comm comm, MPI_Info info) ;
int MPI_Comm_set_name(MPI_Comm comm, const char *comm_name) ;
int MPI_Comm_size(MPI_Comm comm, int *size) ;
int MPI_Comm_split(MPI_Comm comm, int color, int key, MPI_Comm *newcomm) ;
int MPI_Comm_split_type(MPI_Comm comm, int split_type, int key, MPI_Info info, MPI_Comm *newcomm)
                    ;
int MPI_Comm_test_inter(MPI_Comm comm, int *flag) ;
int MPI_Intercomm_create(MPI_Comm local_comm, int local_leader, MPI_Comm peer_comm,
                         int remote_leader, int tag, MPI_Comm *newintercomm) ;
int MPI_Intercomm_create_from_groups(MPI_Group local_group, int local_leader,
                                     MPI_Group remote_group, int remote_leader,
                                     const char *stringtag, MPI_Info info,
                                     MPI_Errhandler errhandler, MPI_Comm *newintercomm)
                                                     ;
int MPI_Intercomm_merge(MPI_Comm intercomm, int high, MPI_Comm *newintracomm) ;
int MPIX_Comm_test_threadcomm(MPI_Comm comm, int *flag) ;
int MPIX_Comm_revoke(MPI_Comm comm) ;
int MPIX_Comm_shrink(MPI_Comm comm, MPI_Comm *newcomm) ;
int MPIX_Comm_failure_ack(MPI_Comm comm) ;
int MPIX_Comm_failure_get_acked(MPI_Comm comm, MPI_Group *failedgrp) ;
int MPIX_Comm_agree(MPI_Comm comm, int *flag) ;
int MPIX_Comm_get_failed(MPI_Comm comm, MPI_Group *failedgrp) ;
int MPI_Get_address(const void *location, MPI_Aint *address) ;
int MPI_Address(void *location, MPI_Aint *address) ;
int MPI_Get_count(const MPI_Status *status, MPI_Datatype datatype, int *count) ;
int MPI_Get_elements(const MPI_Status *status, MPI_Datatype datatype, int *count) ;
int MPI_Get_elements_x(const MPI_Status *status, MPI_Datatype datatype, MPI_Count *count)
                    ;
int MPI_Pack(const void *inbuf, int incount, MPI_Datatype datatype, void *outbuf, int outsize,
             int *position, MPI_Comm comm) ;
int MPI_Pack_external(const char *datarep, const void *inbuf, int incount, MPI_Datatype datatype,
                      void *outbuf, MPI_Aint outsize, MPI_Aint *position) ;
int MPI_Pack_external_size(const char *datarep, int incount, MPI_Datatype datatype, MPI_Aint *size)
                    ;
int MPI_Pack_size(int incount, MPI_Datatype datatype, MPI_Comm comm, int *size) ;
int MPI_Status_set_elements(MPI_Status *status, MPI_Datatype datatype, int count) ;
int MPI_Status_set_elements_x(MPI_Status *status, MPI_Datatype datatype, MPI_Count count)
                    ;
int MPI_Type_commit(MPI_Datatype *datatype) ;
int MPI_Type_contiguous(int count, MPI_Datatype oldtype, MPI_Datatype *newtype) ;
int MPI_Type_create_darray(int size, int rank, int ndims, const int array_of_gsizes[],
                           const int array_of_distribs[], const int array_of_dargs[],
                           const int array_of_psizes[], int order, MPI_Datatype oldtype,
                           MPI_Datatype *newtype) ;
int MPI_Type_create_f90_complex(int p, int r, MPI_Datatype *newtype) ;
int MPI_Type_create_f90_integer(int r, MPI_Datatype *newtype) ;
int MPI_Type_create_f90_real(int p, int r, MPI_Datatype *newtype) ;
int MPI_Type_create_hindexed(int count, const int array_of_blocklengths[],
                             const MPI_Aint array_of_displacements[], MPI_Datatype oldtype,
                             MPI_Datatype *newtype) ;
int MPI_Type_hindexed(int count, int array_of_blocklengths[], MPI_Aint array_of_displacements[],
                      MPI_Datatype oldtype, MPI_Datatype *newtype) ;
int MPI_Type_create_hindexed_block(int count, int blocklength,
                                   const MPI_Aint array_of_displacements[], MPI_Datatype oldtype,
                                   MPI_Datatype *newtype) ;
int MPI_Type_create_hvector(int count, int blocklength, MPI_Aint stride, MPI_Datatype oldtype,
                            MPI_Datatype *newtype) ;
int MPI_Type_hvector(int count, int blocklength, MPI_Aint stride, MPI_Datatype oldtype,
                     MPI_Datatype *newtype) ;
int MPI_Type_create_indexed_block(int count, int blocklength, const int array_of_displacements[],
                                  MPI_Datatype oldtype, MPI_Datatype *newtype) ;
int MPI_Type_create_resized(MPI_Datatype oldtype, MPI_Aint lb, MPI_Aint extent,
                            MPI_Datatype *newtype) ;
int MPI_Type_create_struct(int count, const int array_of_blocklengths[],
                           const MPI_Aint array_of_displacements[],
                           const MPI_Datatype array_of_types[], MPI_Datatype *newtype)
                                           ;
int MPI_Type_struct(int count, int array_of_blocklengths[], MPI_Aint array_of_displacements[],
                    MPI_Datatype array_of_types[], MPI_Datatype *newtype) ;
int MPI_Type_create_subarray(int ndims, const int array_of_sizes[], const int array_of_subsizes[],
                             const int array_of_starts[], int order, MPI_Datatype oldtype,
                             MPI_Datatype *newtype) ;
int MPI_Type_dup(MPI_Datatype oldtype, MPI_Datatype *newtype) ;
int MPI_Type_free(MPI_Datatype *datatype) ;
int MPI_Type_get_contents(MPI_Datatype datatype, int max_integers, int max_addresses,
                          int max_datatypes, int array_of_integers[], MPI_Aint array_of_addresses[],
                          MPI_Datatype array_of_datatypes[]) ;
int MPI_Type_get_envelope(MPI_Datatype datatype, int *num_integers, int *num_addresses,
                          int *num_datatypes, int *combiner) ;
int MPI_Type_get_extent(MPI_Datatype datatype, MPI_Aint *lb, MPI_Aint *extent) ;
int MPI_Type_get_extent_x(MPI_Datatype datatype, MPI_Count *lb, MPI_Count *extent)
                    ;
int MPI_Type_get_name(MPI_Datatype datatype, char *type_name, int *resultlen) ;
int MPI_Type_get_true_extent(MPI_Datatype datatype, MPI_Aint *true_lb, MPI_Aint *true_extent)
                    ;
int MPI_Type_get_true_extent_x(MPI_Datatype datatype, MPI_Count *true_lb, MPI_Count *true_extent)
                    ;
int MPI_Type_get_value_index(MPI_Datatype value_type, MPI_Datatype index_type,
                             MPI_Datatype *pair_type) ;
int MPI_Type_indexed(int count, const int array_of_blocklengths[],
                     const int array_of_displacements[], MPI_Datatype oldtype,
                     MPI_Datatype *newtype) ;
int MPI_Type_match_size(int typeclass, int size, MPI_Datatype *datatype) ;
int MPI_Type_set_name(MPI_Datatype datatype, const char *type_name) ;
int MPI_Type_size(MPI_Datatype datatype, int *size) ;
int MPI_Type_size_x(MPI_Datatype datatype, MPI_Count *size) ;
int MPI_Type_vector(int count, int blocklength, int stride, MPI_Datatype oldtype,
                    MPI_Datatype *newtype) ;
int MPI_Unpack(const void *inbuf, int insize, int *position, void *outbuf, int outcount,
               MPI_Datatype datatype, MPI_Comm comm) ;
int MPI_Unpack_external(const char datarep[], const void *inbuf, MPI_Aint insize,
                        MPI_Aint *position, void *outbuf, int outcount, MPI_Datatype datatype)
                                        ;
int MPI_Type_extent(MPI_Datatype datatype, MPI_Aint *extent) ;
int MPI_Type_lb(MPI_Datatype datatype, MPI_Aint *displacement) ;
int MPI_Type_ub(MPI_Datatype datatype, MPI_Aint *displacement) ;
int MPIX_Type_iov_len(MPI_Datatype datatype, MPI_Count max_iov_bytes, MPI_Count *iov_len,
                      MPI_Count *actual_iov_bytes) ;
int MPIX_Type_iov(MPI_Datatype datatype, MPI_Count iov_offset, MPIX_Iov *iov, MPI_Count max_iov_len,
                  MPI_Count *actual_iov_len) ;
int MPI_Add_error_class(int *errorclass) ;
int MPI_Add_error_code(int errorclass, int *errorcode) ;
int MPI_Add_error_string(int errorcode, const char *string) ;
int MPI_Comm_call_errhandler(MPI_Comm comm, int errorcode) ;
int MPI_Comm_create_errhandler(MPI_Comm_errhandler_function *comm_errhandler_fn,
                               MPI_Errhandler *errhandler) ;
int MPI_Errhandler_create(MPI_Comm_errhandler_function *comm_errhandler_fn,
                          MPI_Errhandler *errhandler) ;
int MPI_Comm_get_errhandler(MPI_Comm comm, MPI_Errhandler *errhandler) ;
int MPI_Errhandler_get(MPI_Comm comm, MPI_Errhandler *errhandler) ;
int MPI_Comm_set_errhandler(MPI_Comm comm, MPI_Errhandler errhandler) ;
int MPI_Errhandler_set(MPI_Comm comm, MPI_Errhandler errhandler) ;
int MPI_Errhandler_free(MPI_Errhandler *errhandler) ;
int MPI_Error_class(int errorcode, int *errorclass) ;
int MPI_Error_string(int errorcode, char *string, int *resultlen) ;
int MPI_File_call_errhandler(MPI_File fh, int errorcode) ;
int MPI_File_create_errhandler(MPI_File_errhandler_function *file_errhandler_fn,
                               MPI_Errhandler *errhandler) ;
int MPI_File_get_errhandler(MPI_File file, MPI_Errhandler *errhandler) ;
int MPI_File_set_errhandler(MPI_File file, MPI_Errhandler errhandler) ;
int MPI_Remove_error_class(int errorclass) ;
int MPI_Remove_error_code(int errorcode) ;
int MPI_Remove_error_string(int errorcode) ;
int MPI_Session_call_errhandler(MPI_Session session, int errorcode) ;
int MPI_Session_create_errhandler(MPI_Session_errhandler_function *session_errhandler_fn,
                                  MPI_Errhandler *errhandler) ;
int MPI_Session_get_errhandler(MPI_Session session, MPI_Errhandler *errhandler) ;
int MPI_Session_set_errhandler(MPI_Session session, MPI_Errhandler errhandler) ;
int MPI_Win_call_errhandler(MPI_Win win, int errorcode) ;
int MPI_Win_create_errhandler(MPI_Win_errhandler_function *win_errhandler_fn,
                              MPI_Errhandler *errhandler) ;
int MPI_Win_get_errhandler(MPI_Win win, MPI_Errhandler *errhandler) ;
int MPI_Win_set_errhandler(MPI_Win win, MPI_Errhandler errhandler) ;
MPI_Fint MPI_Comm_c2f(MPI_Comm comm) ;
MPI_Comm MPI_Comm_f2c(MPI_Fint comm) ;
MPI_Fint MPI_Errhandler_c2f(MPI_Errhandler errhandler) ;
MPI_Errhandler MPI_Errhandler_f2c(MPI_Fint errhandler) ;
MPI_Fint MPI_Group_c2f(MPI_Group group) ;
MPI_Group MPI_Group_f2c(MPI_Fint group) ;
MPI_Fint MPI_Info_c2f(MPI_Info info) ;
MPI_Info MPI_Info_f2c(MPI_Fint info) ;
MPI_Fint MPI_Message_c2f(MPI_Message message) ;
MPI_Message MPI_Message_f2c(MPI_Fint message) ;
MPI_Fint MPI_Op_c2f(MPI_Op op) ;
MPI_Op MPI_Op_f2c(MPI_Fint op) ;
MPI_Fint MPI_Request_c2f(MPI_Request request) ;
MPI_Request MPI_Request_f2c(MPI_Fint request) ;
MPI_Fint MPI_Session_c2f(MPI_Session session) ;
MPI_Session MPI_Session_f2c(MPI_Fint session) ;
int MPI_Status_c2f(const MPI_Status *c_status, MPI_Fint *f_status) ;
int MPI_Status_c2f08(const MPI_Status *c_status, MPI_F08_status *f08_status) ;
int MPI_Status_f082c(const MPI_F08_status *f08_status, MPI_Status *c_status) ;
int MPI_Status_f082f(const MPI_F08_status *f08_status, MPI_Fint *f_status) ;
int MPI_Status_f2c(const MPI_Fint *f_status, MPI_Status *c_status) ;
int MPI_Status_f2f08(const MPI_Fint *f_status, MPI_F08_status *f08_status) ;
MPI_Fint MPI_Type_c2f(MPI_Datatype datatype) ;
MPI_Datatype MPI_Type_f2c(MPI_Fint datatype) ;
MPI_Fint MPI_Win_c2f(MPI_Win win) ;
MPI_Win MPI_Win_f2c(MPI_Fint win) ;
int MPI_Group_compare(MPI_Group group1, MPI_Group group2, int *result) ;
int MPI_Group_difference(MPI_Group group1, MPI_Group group2, MPI_Group *newgroup) ;
int MPI_Group_excl(MPI_Group group, int n, const int ranks[], MPI_Group *newgroup)
                    ;
int MPI_Group_free(MPI_Group *group) ;
int MPI_Group_incl(MPI_Group group, int n, const int ranks[], MPI_Group *newgroup)
                    ;
int MPI_Group_intersection(MPI_Group group1, MPI_Group group2, MPI_Group *newgroup)
                    ;
int MPI_Group_range_excl(MPI_Group group, int n, int ranges[][3], MPI_Group *newgroup)
                    ;
int MPI_Group_range_incl(MPI_Group group, int n, int ranges[][3], MPI_Group *newgroup)
                    ;
int MPI_Group_rank(MPI_Group group, int *rank) ;
int MPI_Group_size(MPI_Group group, int *size) ;
int MPI_Group_translate_ranks(MPI_Group group1, int n, const int ranks1[], MPI_Group group2,
                              int ranks2[]) ;
int MPI_Group_union(MPI_Group group1, MPI_Group group2, MPI_Group *newgroup) ;
int MPI_Info_create(MPI_Info *info) ;
int MPI_Info_create_env(int argc, char *argv[], MPI_Info *info) ;
int MPI_Info_delete(MPI_Info info, const char *key) ;
int MPI_Info_dup(MPI_Info info, MPI_Info *newinfo) ;
int MPI_Info_free(MPI_Info *info) ;
int MPI_Info_get(MPI_Info info, const char *key, int valuelen, char *value, int *flag)
                    ;
int MPI_Info_get_nkeys(MPI_Info info, int *nkeys) ;
int MPI_Info_get_nthkey(MPI_Info info, int n, char *key) ;
int MPI_Info_get_string(MPI_Info info, const char *key, int *buflen, char *value, int *flag)
                    ;
int MPI_Info_get_valuelen(MPI_Info info, const char *key, int *valuelen, int *flag)
                    ;
int MPI_Info_set(MPI_Info info, const char *key, const char *value) ;
int MPIX_Info_set_hex(MPI_Info info, const char *key, const void *value, int value_size)
                    ;
int MPI_Abort(MPI_Comm comm, int errorcode) ;
int MPI_Comm_create_from_group(MPI_Group group, const char *stringtag, MPI_Info info,
                               MPI_Errhandler errhandler, MPI_Comm *newcomm) ;
int MPI_Finalize(void) ;
int MPI_Finalized(int *flag) ;
int MPI_Group_from_session_pset(MPI_Session session, const char *pset_name, MPI_Group *newgroup)
                    ;
int MPI_Init(int *argc, char ***argv) ;
int MPI_Init_thread(int *argc, char ***argv, int required, int *provided) ;
int MPI_Initialized(int *flag) ;
int MPI_Is_thread_main(int *flag) ;
int MPI_Query_thread(int *provided) ;
int MPI_Session_finalize(MPI_Session *session) ;
int MPI_Session_get_info(MPI_Session session, MPI_Info *info_used) ;
int MPI_Session_get_nth_pset(MPI_Session session, MPI_Info info, int n, int *pset_len,
                             char *pset_name) ;
int MPI_Session_get_num_psets(MPI_Session session, MPI_Info info, int *npset_names)
                    ;
int MPI_Session_get_pset_info(MPI_Session session, const char *pset_name, MPI_Info *info)
                    ;
int MPI_Session_init(MPI_Info info, MPI_Errhandler errhandler, MPI_Session *session)
                    ;
MPI_Aint MPI_Aint_add(MPI_Aint base, MPI_Aint disp) ;
MPI_Aint MPI_Aint_diff(MPI_Aint addr1, MPI_Aint addr2) ;
int MPI_Get_library_version(char *version, int *resultlen) ;
int MPI_Get_processor_name(char *name, int *resultlen) ;
int MPI_Get_version(int *version, int *subversion) ;
int MPI_Pcontrol(const int level, ...) ;
int MPIX_GPU_query_support(int gpu_type, int *is_supported) ;
int MPIX_Query_cuda_support(void) ;
int MPIX_Query_ze_support(void) ;
int MPIX_Query_hip_support(void) ;
int MPI_Op_commutative(MPI_Op op, int *commute) ;
int MPI_Op_create(MPI_User_function *user_fn, int commute, MPI_Op *op) ;
int MPI_Op_free(MPI_Op *op) ;
int MPI_Parrived(MPI_Request request, int partition, int *flag) ;
int MPI_Pready(int partition, MPI_Request request) ;
int MPI_Pready_list(int length, const int array_of_partitions[], MPI_Request request)
                    ;
int MPI_Pready_range(int partition_low, int partition_high, MPI_Request request) ;
int MPI_Precv_init(void *buf, int partitions, MPI_Count count, MPI_Datatype datatype, int dest,
                   int tag, MPI_Comm comm, MPI_Info info, MPI_Request *request)
                                                                         ;
int MPI_Psend_init(const void *buf, int partitions, MPI_Count count, MPI_Datatype datatype,
                   int dest, int tag, MPI_Comm comm, MPI_Info info, MPI_Request *request)
                                                                         ;
int MPI_Bsend(const void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm)
                                                          ;
int MPI_Bsend_init(const void *buf, int count, MPI_Datatype datatype, int dest, int tag,
                   MPI_Comm comm, MPI_Request *request)
                                                                         ;
int MPI_Buffer_attach(void *buffer, int size) ;
int MPI_Buffer_detach(void *buffer_addr, int *size) ;
int MPI_Buffer_flush(void) ;
int MPI_Buffer_iflush(MPI_Request *request) ;
int MPI_Comm_attach_buffer(MPI_Comm comm, void *buffer, int size) ;
int MPI_Comm_detach_buffer(MPI_Comm comm, void *buffer_addr, int *size) ;
int MPI_Comm_flush_buffer(MPI_Comm comm) ;
int MPI_Comm_iflush_buffer(MPI_Comm comm, MPI_Request *request) ;
int MPI_Ibsend(const void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm,
               MPI_Request *request) ;
int MPI_Improbe(int source, int tag, MPI_Comm comm, int *flag, MPI_Message *message,
                MPI_Status *status) ;
int MPI_Imrecv(void *buf, int count, MPI_Datatype datatype, MPI_Message *message,
               MPI_Request *request) ;
int MPI_Iprobe(int source, int tag, MPI_Comm comm, int *flag, MPI_Status *status) ;
int MPI_Irecv(void *buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm,
              MPI_Request *request) ;
int MPI_Irsend(const void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm,
               MPI_Request *request) ;
int MPI_Isend(const void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm,
              MPI_Request *request) ;
int MPI_Isendrecv(const void *sendbuf, int sendcount, MPI_Datatype sendtype, int dest, int sendtag,
                  void *recvbuf, int recvcount, MPI_Datatype recvtype, int source, int recvtag,
                  MPI_Comm comm, MPI_Request *request)
                                                                                                              ;
int MPI_Isendrecv_replace(void *buf, int count, MPI_Datatype datatype, int dest, int sendtag,
                          int source, int recvtag, MPI_Comm comm, MPI_Request *request)
                                                                                ;
int MPI_Issend(const void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm,
               MPI_Request *request) ;
int MPI_Mprobe(int source, int tag, MPI_Comm comm, MPI_Message *message, MPI_Status *status)
                    ;
int MPI_Mrecv(void *buf, int count, MPI_Datatype datatype, MPI_Message *message,
              MPI_Status *status) ;
int MPI_Probe(int source, int tag, MPI_Comm comm, MPI_Status *status) ;
int MPI_Recv(void *buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm,
             MPI_Status *status) ;
int MPI_Recv_init(void *buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm,
                  MPI_Request *request) ;
int MPI_Rsend(const void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm)
                                                          ;
int MPI_Rsend_init(const void *buf, int count, MPI_Datatype datatype, int dest, int tag,
                   MPI_Comm comm, MPI_Request *request)
                                                                         ;
int MPI_Send(const void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm)
                                                          ;
int MPI_Send_init(const void *buf, int count, MPI_Datatype datatype, int dest, int tag,
                  MPI_Comm comm, MPI_Request *request)
                                                                        ;
int MPI_Sendrecv(const void *sendbuf, int sendcount, MPI_Datatype sendtype, int dest, int sendtag,
                 void *recvbuf, int recvcount, MPI_Datatype recvtype, int source, int recvtag,
                 MPI_Comm comm, MPI_Status *status)
                                                                                                             ;
int MPI_Sendrecv_replace(void *buf, int count, MPI_Datatype datatype, int dest, int sendtag,
                         int source, int recvtag, MPI_Comm comm, MPI_Status *status)
                                                                               ;
int MPI_Session_attach_buffer(MPI_Session session, void *buffer, int size) ;
int MPI_Session_detach_buffer(MPI_Session session, void *buffer_addr, int *size) ;
int MPI_Session_flush_buffer(MPI_Session session) ;
int MPI_Session_iflush_buffer(MPI_Session session, MPI_Request *request) ;
int MPI_Ssend(const void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm)
                                                          ;
int MPI_Ssend_init(const void *buf, int count, MPI_Datatype datatype, int dest, int tag,
                   MPI_Comm comm, MPI_Request *request)
                                                                         ;
int MPI_Cancel(MPI_Request *request) ;
int MPI_Grequest_complete(MPI_Request request) ;
int MPI_Grequest_start(MPI_Grequest_query_function *query_fn, MPI_Grequest_free_function *free_fn,
                       MPI_Grequest_cancel_function *cancel_fn, void *extra_state,
                       MPI_Request *request) ;
int MPI_Request_free(MPI_Request *request) ;
int MPI_Request_get_status(MPI_Request request, int *flag, MPI_Status *status) ;
int MPI_Request_get_status_all(int count, const MPI_Request array_of_requests[], int *flag,
                               MPI_Status *array_of_statuses) ;
int MPI_Request_get_status_any(int count, const MPI_Request array_of_requests[], int *indx,
                               int *flag, MPI_Status *status) ;
int MPI_Request_get_status_some(int incount, const MPI_Request array_of_requests[], int *outcount,
                                int array_of_indices[], MPI_Status *array_of_statuses)
                                                ;
int MPI_Start(MPI_Request *request) ;
int MPI_Startall(int count, MPI_Request array_of_requests[]) ;
int MPI_Status_get_error(const MPI_Status *status, int *error) ;
int MPI_Status_get_source(const MPI_Status *status, int *source) ;
int MPI_Status_get_tag(const MPI_Status *status, int *tag) ;
int MPI_Status_set_error(MPI_Status *status, int error) ;
int MPI_Status_set_source(MPI_Status *status, int source) ;
int MPI_Status_set_tag(MPI_Status *status, int tag) ;
int MPI_Status_set_cancelled(MPI_Status *status, int flag) ;
int MPI_Test(MPI_Request *request, int *flag, MPI_Status *status) ;
int MPI_Test_cancelled(const MPI_Status *status, int *flag) ;
int MPI_Testall(int count, MPI_Request array_of_requests[], int *flag,
                MPI_Status *array_of_statuses) ;
int MPI_Testany(int count, MPI_Request array_of_requests[], int *indx, int *flag,
                MPI_Status *status) ;
int MPI_Testsome(int incount, MPI_Request array_of_requests[], int *outcount,
                 int array_of_indices[], MPI_Status *array_of_statuses) ;
int MPI_Wait(MPI_Request *request, MPI_Status *status) ;
int MPI_Waitall(int count, MPI_Request array_of_requests[], MPI_Status *array_of_statuses)
                    ;
int MPI_Waitany(int count, MPI_Request array_of_requests[], int *indx, MPI_Status *status)
                    ;
int MPI_Waitsome(int incount, MPI_Request array_of_requests[], int *outcount,
                 int array_of_indices[], MPI_Status *array_of_statuses) ;
int MPIX_Request_is_complete(MPI_Request request) ;
int MPI_Accumulate(const void *origin_addr, int origin_count, MPI_Datatype origin_datatype,
                   int target_rank, MPI_Aint target_disp, int target_count,
                   MPI_Datatype target_datatype, MPI_Op op, MPI_Win win)
                                                                         ;
int MPI_Alloc_mem(MPI_Aint size, MPI_Info info, void *baseptr) ;
int MPI_Compare_and_swap(const void *origin_addr, const void *compare_addr, void *result_addr,
                         MPI_Datatype datatype, int target_rank, MPI_Aint target_disp, MPI_Win win)
                                         ;
int MPI_Fetch_and_op(const void *origin_addr, void *result_addr, MPI_Datatype datatype,
                     int target_rank, MPI_Aint target_disp, MPI_Op op, MPI_Win win)
                                     ;
int MPI_Free_mem(void *base) ;
int MPI_Get(void *origin_addr, int origin_count, MPI_Datatype origin_datatype, int target_rank,
            MPI_Aint target_disp, int target_count, MPI_Datatype target_datatype, MPI_Win win)
                                                                  ;
int MPI_Get_accumulate(const void *origin_addr, int origin_count, MPI_Datatype origin_datatype,
                       void *result_addr, int result_count, MPI_Datatype result_datatype,
                       int target_rank, MPI_Aint target_disp, int target_count,
                       MPI_Datatype target_datatype, MPI_Op op, MPI_Win win)
                                                                                                                   ;
int MPI_Put(const void *origin_addr, int origin_count, MPI_Datatype origin_datatype,
            int target_rank, MPI_Aint target_disp, int target_count, MPI_Datatype target_datatype,
            MPI_Win win) ;
int MPI_Raccumulate(const void *origin_addr, int origin_count, MPI_Datatype origin_datatype,
                    int target_rank, MPI_Aint target_disp, int target_count,
                    MPI_Datatype target_datatype, MPI_Op op, MPI_Win win, MPI_Request *request)
                                                                          ;
int MPI_Rget(void *origin_addr, int origin_count, MPI_Datatype origin_datatype, int target_rank,
             MPI_Aint target_disp, int target_count, MPI_Datatype target_datatype, MPI_Win win,
             MPI_Request *request) ;
int MPI_Rget_accumulate(const void *origin_addr, int origin_count, MPI_Datatype origin_datatype,
                        void *result_addr, int result_count, MPI_Datatype result_datatype,
                        int target_rank, MPI_Aint target_disp, int target_count,
                        MPI_Datatype target_datatype, MPI_Op op, MPI_Win win, MPI_Request *request)
                                                                                                                    ;
int MPI_Rput(const void *origin_addr, int origin_count, MPI_Datatype origin_datatype,
             int target_rank, MPI_Aint target_disp, int target_count, MPI_Datatype target_datatype,
             MPI_Win win, MPI_Request *request)
                                                                   ;
int MPI_Win_allocate(MPI_Aint size, int disp_unit, MPI_Info info, MPI_Comm comm, void *baseptr,
                     MPI_Win *win) ;
int MPI_Win_allocate_shared(MPI_Aint size, int disp_unit, MPI_Info info, MPI_Comm comm,
                            void *baseptr, MPI_Win *win) ;
int MPI_Win_attach(MPI_Win win, void *base, MPI_Aint size) ;
int MPI_Win_complete(MPI_Win win) ;
int MPI_Win_create(void *base, MPI_Aint size, int disp_unit, MPI_Info info, MPI_Comm comm,
                   MPI_Win *win) ;
int MPI_Win_create_dynamic(MPI_Info info, MPI_Comm comm, MPI_Win *win) ;
int MPI_Win_detach(MPI_Win win, const void *base) ;
int MPI_Win_fence(int assert, MPI_Win win) ;
int MPI_Win_flush(int rank, MPI_Win win) ;
int MPI_Win_flush_all(MPI_Win win) ;
int MPI_Win_flush_local(int rank, MPI_Win win) ;
int MPI_Win_flush_local_all(MPI_Win win) ;
int MPI_Win_free(MPI_Win *win) ;
int MPI_Win_get_group(MPI_Win win, MPI_Group *group) ;
int MPI_Win_get_info(MPI_Win win, MPI_Info *info_used) ;
int MPI_Win_get_name(MPI_Win win, char *win_name, int *resultlen) ;
int MPI_Win_lock(int lock_type, int rank, int assert, MPI_Win win) ;
int MPI_Win_lock_all(int assert, MPI_Win win) ;
int MPI_Win_post(MPI_Group group, int assert, MPI_Win win) ;
int MPI_Win_set_info(MPI_Win win, MPI_Info info) ;
int MPI_Win_set_name(MPI_Win win, const char *win_name) ;
int MPI_Win_shared_query(MPI_Win win, int rank, MPI_Aint *size, int *disp_unit, void *baseptr)
                    ;
int MPI_Win_start(MPI_Group group, int assert, MPI_Win win) ;
int MPI_Win_sync(MPI_Win win) ;
int MPI_Win_test(MPI_Win win, int *flag) ;
int MPI_Win_unlock(int rank, MPI_Win win) ;
int MPI_Win_unlock_all(MPI_Win win) ;
int MPI_Win_wait(MPI_Win win) ;
int MPI_Close_port(const char *port_name) ;
int MPI_Comm_accept(const char *port_name, MPI_Info info, int root, MPI_Comm comm,
                    MPI_Comm *newcomm) ;
int MPI_Comm_connect(const char *port_name, MPI_Info info, int root, MPI_Comm comm,
                     MPI_Comm *newcomm) ;
int MPI_Comm_disconnect(MPI_Comm *comm) ;
int MPI_Comm_get_parent(MPI_Comm *parent) ;
int MPI_Comm_join(int fd, MPI_Comm *intercomm) ;
int MPI_Comm_spawn(const char *command, char *argv[], int maxprocs, MPI_Info info, int root,
                   MPI_Comm comm, MPI_Comm *intercomm, int array_of_errcodes[]) ;
int MPI_Comm_spawn_multiple(int count, char *array_of_commands[], char **array_of_argv[],
                            const int array_of_maxprocs[], const MPI_Info array_of_info[], int root,
                            MPI_Comm comm, MPI_Comm *intercomm, int array_of_errcodes[])
                                            ;
int MPI_Lookup_name(const char *service_name, MPI_Info info, char *port_name) ;
int MPI_Open_port(MPI_Info info, char *port_name) ;
int MPI_Publish_name(const char *service_name, MPI_Info info, const char *port_name)
                    ;
int MPI_Unpublish_name(const char *service_name, MPI_Info info, const char *port_name)
                    ;
int MPIX_Stream_create(MPI_Info info, MPIX_Stream *stream) ;
int MPIX_Stream_free(MPIX_Stream *stream) ;
int MPIX_Stream_comm_create(MPI_Comm comm, MPIX_Stream stream, MPI_Comm *newcomm) ;
int MPIX_Stream_comm_create_multiplex(MPI_Comm comm, int count, MPIX_Stream array_of_streams[],
                                      MPI_Comm *newcomm) ;
int MPIX_Comm_get_stream(MPI_Comm comm, int idx, MPIX_Stream *stream) ;
int MPIX_Stream_progress(MPIX_Stream stream) ;
int MPIX_Start_progress_thread(MPIX_Stream stream) ;
int MPIX_Stop_progress_thread(MPIX_Stream stream) ;
int MPIX_Stream_send(const void *buf, int count, MPI_Datatype datatype, int dest, int tag,
                     MPI_Comm comm, int source_stream_index, int dest_stream_index)
                                                                           ;
int MPIX_Stream_isend(const void *buf, int count, MPI_Datatype datatype, int dest, int tag,
                      MPI_Comm comm, int source_stream_index, int dest_stream_index,
                      MPI_Request *request) ;
int MPIX_Stream_recv(void *buf, int count, MPI_Datatype datatype, int source, int tag,
                     MPI_Comm comm, int source_stream_index, int dest_stream_index,
                     MPI_Status *status) ;
int MPIX_Stream_irecv(void *buf, int count, MPI_Datatype datatype, int source, int tag,
                      MPI_Comm comm, int source_stream_index, int dest_stream_index,
                      MPI_Request *request) ;
int MPIX_Send_enqueue(const void *buf, int count, MPI_Datatype datatype, int dest, int tag,
                      MPI_Comm comm) ;
int MPIX_Recv_enqueue(void *buf, int count, MPI_Datatype datatype, int source, int tag,
                      MPI_Comm comm, MPI_Status *status)
                                                                            ;
int MPIX_Isend_enqueue(const void *buf, int count, MPI_Datatype datatype, int dest, int tag,
                       MPI_Comm comm, MPI_Request *request)
                                                                             ;
int MPIX_Irecv_enqueue(void *buf, int count, MPI_Datatype datatype, int source, int tag,
                       MPI_Comm comm, MPI_Request *request)
                                                                             ;
int MPIX_Wait_enqueue(MPI_Request *request, MPI_Status *status) ;
int MPIX_Waitall_enqueue(int count, MPI_Request array_of_requests[], MPI_Status *array_of_statuses)
                    ;
int MPIX_Allreduce_enqueue(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype,
                           MPI_Op op, MPI_Comm comm)
                                                                                 ;
int MPIX_Threadcomm_init(MPI_Comm comm, int num_threads, MPI_Comm *newthreadcomm) ;
int MPIX_Threadcomm_free(MPI_Comm *threadcomm) ;
int MPIX_Threadcomm_start(MPI_Comm threadcomm) ;
int MPIX_Threadcomm_finish(MPI_Comm threadcomm) ;
double MPI_Wtick(void) ;
double MPI_Wtime(void) ;
int MPI_Cart_coords(MPI_Comm comm, int rank, int maxdims, int coords[]) ;
int MPI_Cart_create(MPI_Comm comm_old, int ndims, const int dims[], const int periods[],
                    int reorder, MPI_Comm *comm_cart) ;
int MPI_Cart_get(MPI_Comm comm, int maxdims, int dims[], int periods[], int coords[])
                    ;
int MPI_Cart_map(MPI_Comm comm, int ndims, const int dims[], const int periods[], int *newrank)
                    ;
int MPI_Cart_rank(MPI_Comm comm, const int coords[], int *rank) ;
int MPI_Cart_shift(MPI_Comm comm, int direction, int disp, int *rank_source, int *rank_dest)
                    ;
int MPI_Cart_sub(MPI_Comm comm, const int remain_dims[], MPI_Comm *newcomm) ;
int MPI_Cartdim_get(MPI_Comm comm, int *ndims) ;
int MPI_Dims_create(int nnodes, int ndims, int dims[]) ;
int MPI_Dist_graph_create(MPI_Comm comm_old, int n, const int sources[], const int degrees[],
                          const int destinations[], const int weights[], MPI_Info info, int reorder,
                          MPI_Comm *comm_dist_graph) ;
int MPI_Dist_graph_create_adjacent(MPI_Comm comm_old, int indegree, const int sources[],
                                   const int sourceweights[], int outdegree,
                                   const int destinations[], const int destweights[], MPI_Info info,
                                   int reorder, MPI_Comm *comm_dist_graph) ;
int MPI_Dist_graph_neighbors(MPI_Comm comm, int maxindegree, int sources[], int sourceweights[],
                             int maxoutdegree, int destinations[], int destweights[])
                                             ;
int MPI_Dist_graph_neighbors_count(MPI_Comm comm, int *indegree, int *outdegree, int *weighted)
                    ;
int MPI_Get_hw_resource_info(MPI_Info *hw_info) ;
int MPI_Graph_create(MPI_Comm comm_old, int nnodes, const int indx[], const int edges[],
                     int reorder, MPI_Comm *comm_graph) ;
int MPI_Graph_get(MPI_Comm comm, int maxindex, int maxedges, int indx[], int edges[])
                    ;
int MPI_Graph_map(MPI_Comm comm, int nnodes, const int indx[], const int edges[], int *newrank)
                    ;
int MPI_Graph_neighbors(MPI_Comm comm, int rank, int maxneighbors, int neighbors[])
                    ;
int MPI_Graph_neighbors_count(MPI_Comm comm, int rank, int *nneighbors) ;
int MPI_Graphdims_get(MPI_Comm comm, int *nnodes, int *nedges) ;
int MPI_Topo_test(MPI_Comm comm, int *status) ;
MPI_Fint MPI_File_c2f(MPI_File file) ;
int MPI_File_close(MPI_File *fh) ;
int MPI_File_delete(const char *filename, MPI_Info info) ;
MPI_File MPI_File_f2c(MPI_Fint file) ;
int MPI_File_get_amode(MPI_File fh, int *amode) ;
int MPI_File_get_atomicity(MPI_File fh, int *flag) ;
int MPI_File_get_byte_offset(MPI_File fh, MPI_Offset offset, MPI_Offset *disp) ;
int MPI_File_get_group(MPI_File fh, MPI_Group *group) ;
int MPI_File_get_info(MPI_File fh, MPI_Info *info_used) ;
int MPI_File_get_position(MPI_File fh, MPI_Offset *offset) ;
int MPI_File_get_position_shared(MPI_File fh, MPI_Offset *offset) ;
int MPI_File_get_size(MPI_File fh, MPI_Offset *size) ;
int MPI_File_get_type_extent(MPI_File fh, MPI_Datatype datatype, MPI_Aint *extent)
                    ;
int MPI_File_get_view(MPI_File fh, MPI_Offset *disp, MPI_Datatype *etype, MPI_Datatype *filetype,
                      char *datarep) ;
int MPI_File_iread(MPI_File fh, void *buf, int count, MPI_Datatype datatype, MPI_Request *request)
                                                          ;
int MPI_File_iread_all(MPI_File fh, void *buf, int count, MPI_Datatype datatype,
                       MPI_Request *request)
                                                                             ;
int MPI_File_iread_at(MPI_File fh, MPI_Offset offset, void *buf, int count, MPI_Datatype datatype,
                      MPI_Request *request) ;
int MPI_File_iread_at_all(MPI_File fh, MPI_Offset offset, void *buf, int count,
                          MPI_Datatype datatype, MPI_Request *request)
                                                                                ;
int MPI_File_iread_shared(MPI_File fh, void *buf, int count, MPI_Datatype datatype,
                          MPI_Request *request)
                                                                                ;
int MPI_File_iwrite(MPI_File fh, const void *buf, int count, MPI_Datatype datatype,
                    MPI_Request *request) ;
int MPI_File_iwrite_all(MPI_File fh, const void *buf, int count, MPI_Datatype datatype,
                        MPI_Request *request)
                                                                              ;
int MPI_File_iwrite_at(MPI_File fh, MPI_Offset offset, const void *buf, int count,
                       MPI_Datatype datatype, MPI_Request *request)
                                                                             ;
int MPI_File_iwrite_at_all(MPI_File fh, MPI_Offset offset, const void *buf, int count,
                           MPI_Datatype datatype, MPI_Request *request)
                                                                                 ;
int MPI_File_iwrite_shared(MPI_File fh, const void *buf, int count, MPI_Datatype datatype,
                           MPI_Request *request)
                                                                                 ;
int MPI_File_open(MPI_Comm comm, const char *filename, int amode, MPI_Info info, MPI_File *fh)
                    ;
int MPI_File_preallocate(MPI_File fh, MPI_Offset size) ;
int MPI_File_read(MPI_File fh, void *buf, int count, MPI_Datatype datatype, MPI_Status *status)
                                                          ;
int MPI_File_read_all(MPI_File fh, void *buf, int count, MPI_Datatype datatype, MPI_Status *status)
                                                          ;
int MPI_File_read_all_begin(MPI_File fh, void *buf, int count, MPI_Datatype datatype)
                                                          ;
int MPI_File_read_all_end(MPI_File fh, void *buf, MPI_Status *status) ;
int MPI_File_read_at(MPI_File fh, MPI_Offset offset, void *buf, int count, MPI_Datatype datatype,
                     MPI_Status *status) ;
int MPI_File_read_at_all(MPI_File fh, MPI_Offset offset, void *buf, int count,
                         MPI_Datatype datatype, MPI_Status *status)
                                                                               ;
int MPI_File_read_at_all_begin(MPI_File fh, MPI_Offset offset, void *buf, int count,
                               MPI_Datatype datatype)
                                                                                     ;
int MPI_File_read_at_all_end(MPI_File fh, void *buf, MPI_Status *status) ;
int MPI_File_read_ordered(MPI_File fh, void *buf, int count, MPI_Datatype datatype,
                          MPI_Status *status)
                                                                                ;
int MPI_File_read_ordered_begin(MPI_File fh, void *buf, int count, MPI_Datatype datatype)
                                                          ;
int MPI_File_read_ordered_end(MPI_File fh, void *buf, MPI_Status *status) ;
int MPI_File_read_shared(MPI_File fh, void *buf, int count, MPI_Datatype datatype,
                         MPI_Status *status)
                                                                               ;
int MPI_File_seek(MPI_File fh, MPI_Offset offset, int whence) ;
int MPI_File_seek_shared(MPI_File fh, MPI_Offset offset, int whence) ;
int MPI_File_set_atomicity(MPI_File fh, int flag) ;
int MPI_File_set_info(MPI_File fh, MPI_Info info) ;
int MPI_File_set_size(MPI_File fh, MPI_Offset size) ;
int MPI_File_set_view(MPI_File fh, MPI_Offset disp, MPI_Datatype etype, MPI_Datatype filetype,
                      const char *datarep, MPI_Info info) ;
int MPI_File_sync(MPI_File fh) ;
int MPI_File_write(MPI_File fh, const void *buf, int count, MPI_Datatype datatype,
                   MPI_Status *status) ;
int MPI_File_write_all(MPI_File fh, const void *buf, int count, MPI_Datatype datatype,
                       MPI_Status *status) ;
int MPI_File_write_all_begin(MPI_File fh, const void *buf, int count, MPI_Datatype datatype)
                                                          ;
int MPI_File_write_all_end(MPI_File fh, const void *buf, MPI_Status *status) ;
int MPI_File_write_at(MPI_File fh, MPI_Offset offset, const void *buf, int count,
                      MPI_Datatype datatype, MPI_Status *status)
                                                                            ;
int MPI_File_write_at_all(MPI_File fh, MPI_Offset offset, const void *buf, int count,
                          MPI_Datatype datatype, MPI_Status *status)
                                                                                ;
int MPI_File_write_at_all_begin(MPI_File fh, MPI_Offset offset, const void *buf, int count,
                                MPI_Datatype datatype)
                                                                                      ;
int MPI_File_write_at_all_end(MPI_File fh, const void *buf, MPI_Status *status) ;
int MPI_File_write_ordered(MPI_File fh, const void *buf, int count, MPI_Datatype datatype,
                           MPI_Status *status)
                                                                                 ;
int MPI_File_write_ordered_begin(MPI_File fh, const void *buf, int count, MPI_Datatype datatype)
                                                          ;
int MPI_File_write_ordered_end(MPI_File fh, const void *buf, MPI_Status *status) ;
int MPI_File_write_shared(MPI_File fh, const void *buf, int count, MPI_Datatype datatype,
                          MPI_Status *status)
                                                                                ;
int MPI_Register_datarep(const char *datarep, MPI_Datarep_conversion_function *read_conversion_fn,
                         MPI_Datarep_conversion_function *write_conversion_fn,
                         MPI_Datarep_extent_function *dtype_file_extent_fn, void *extra_state)
                                         ;
int MPI_File_toint(MPI_File file) ;
MPI_File MPI_File_fromint(int file) ;
int MPI_T_category_changed(int *update_number) ;
int MPI_T_category_get_categories(int cat_index, int len, int indices[]) ;
int MPI_T_category_get_cvars(int cat_index, int len, int indices[]) ;
int MPI_T_category_get_events(int cat_index, int len, int indices[]) ;
int MPI_T_category_get_index(const char *name, int *cat_index) ;
int MPI_T_category_get_info(int cat_index, char *name, int *name_len, char *desc, int *desc_len,
                            int *num_cvars, int *num_pvars, int *num_categories) ;
int MPI_T_category_get_num(int *num_cat) ;
int MPI_T_category_get_num_events(int cat_index, int *num_events) ;
int MPI_T_category_get_pvars(int cat_index, int len, int indices[]) ;
int MPI_T_cvar_get_index(const char *name, int *cvar_index) ;
int MPI_T_cvar_get_info(int cvar_index, char *name, int *name_len, int *verbosity,
                        MPI_Datatype *datatype, MPI_T_enum *enumtype, char *desc, int *desc_len,
                        int *bind, int *scope) ;
int MPI_T_cvar_get_num(int *num_cvar) ;
int MPI_T_cvar_handle_alloc(int cvar_index, void *obj_handle, MPI_T_cvar_handle *handle,
                            int *count) ;
int MPI_T_cvar_handle_free(MPI_T_cvar_handle *handle) ;
int MPI_T_cvar_read(MPI_T_cvar_handle handle, void *buf) ;
int MPI_T_cvar_write(MPI_T_cvar_handle handle, const void *buf) ;
int MPI_T_enum_get_info(MPI_T_enum enumtype, int *num, char *name, int *name_len) ;
int MPI_T_enum_get_item(MPI_T_enum enumtype, int indx, int *value, char *name, int *name_len)
                    ;
int MPI_T_event_callback_get_info(MPI_T_event_registration event_registration,
                                  MPI_T_cb_safety cb_safety, MPI_Info *info_used) ;
int MPI_T_event_callback_set_info(MPI_T_event_registration event_registration,
                                  MPI_T_cb_safety cb_safety, MPI_Info info) ;
int MPI_T_event_copy(MPI_T_event_instance event_instance, void *buffer) ;
int MPI_T_event_get_index(const char *name, int *event_index) ;
int MPI_T_event_get_info(int event_index, char *name, int *name_len, int *verbosity,
                         MPI_Datatype array_of_datatypes[], MPI_Aint array_of_displacements[],
                         int *num_elements, MPI_T_enum *enumtype, MPI_Info *info, char *desc,
                         int *desc_len, int *bind) ;
int MPI_T_event_get_num(int *num_events) ;
int MPI_T_event_get_source(MPI_T_event_instance event_instance, int *source_index)
                    ;
int MPI_T_event_get_timestamp(MPI_T_event_instance event_instance, MPI_Count *event_timestamp)
                    ;
int MPI_T_event_handle_alloc(int event_index, void *obj_handle, MPI_Info info,
                             MPI_T_event_registration *event_registration) ;
int MPI_T_event_handle_free(MPI_T_event_registration event_registration, void *user_data,
                            MPI_T_event_free_cb_function free_cb_function) ;
int MPI_T_event_handle_get_info(MPI_T_event_registration event_registration, MPI_Info *info_used)
                    ;
int MPI_T_event_handle_set_info(MPI_T_event_registration event_registration, MPI_Info info)
                    ;
int MPI_T_event_read(MPI_T_event_instance event_instance, int element_index, void *buffer)
                    ;
int MPI_T_event_register_callback(MPI_T_event_registration event_registration,
                                  MPI_T_cb_safety cb_safety, MPI_Info info, void *user_data,
                                  MPI_T_event_cb_function event_cb_function) ;
int MPI_T_event_set_dropped_handler(MPI_T_event_registration event_registration,
                                    MPI_T_event_dropped_cb_function dropped_cb_function)
                                                    ;
int MPI_T_finalize(void) ;
int MPI_T_init_thread(int required, int *provided) ;
int MPI_T_pvar_get_index(const char *name, int var_class, int *pvar_index) ;
int MPI_T_pvar_get_info(int pvar_index, char *name, int *name_len, int *verbosity, int *var_class,
                        MPI_Datatype *datatype, MPI_T_enum *enumtype, char *desc, int *desc_len,
                        int *bind, int *readonly, int *continuous, int *atomic) ;
int MPI_T_pvar_get_num(int *num_pvar) ;
int MPI_T_pvar_handle_alloc(MPI_T_pvar_session session, int pvar_index, void *obj_handle,
                            MPI_T_pvar_handle *handle, int *count) ;
int MPI_T_pvar_handle_free(MPI_T_pvar_session session, MPI_T_pvar_handle *handle) ;
int MPI_T_pvar_read(MPI_T_pvar_session session, MPI_T_pvar_handle handle, void *buf)
                    ;
int MPI_T_pvar_readreset(MPI_T_pvar_session session, MPI_T_pvar_handle handle, void *buf)
                    ;
int MPI_T_pvar_reset(MPI_T_pvar_session session, MPI_T_pvar_handle handle) ;
int MPI_T_pvar_session_create(MPI_T_pvar_session *session) ;
int MPI_T_pvar_session_free(MPI_T_pvar_session *session) ;
int MPI_T_pvar_start(MPI_T_pvar_session session, MPI_T_pvar_handle handle) ;
int MPI_T_pvar_stop(MPI_T_pvar_session session, MPI_T_pvar_handle handle) ;
int MPI_T_pvar_write(MPI_T_pvar_session session, MPI_T_pvar_handle handle, const void *buf)
                    ;
int MPI_T_source_get_info(int source_index, char *name, int *name_len, char *desc, int *desc_len,
                          MPI_T_source_order *ordering, MPI_Count *ticks_per_second,
                          MPI_Count *max_ticks, MPI_Info *info) ;
int MPI_T_source_get_num(int *num_sources) ;
int MPI_T_source_get_timestamp(int source_index, MPI_Count *timestamp) ;
int MPIX_Grequest_start(MPI_Grequest_query_function *query_fn, MPI_Grequest_free_function *free_fn,
                        MPI_Grequest_cancel_function *cancel_fn,
                        MPIX_Grequest_poll_function *poll_fn, MPIX_Grequest_wait_function *wait_fn,
                        void *extra_state, MPI_Request *request) ;
int MPIX_Grequest_class_create(MPI_Grequest_query_function *query_fn,
                               MPI_Grequest_free_function *free_fn,
                               MPI_Grequest_cancel_function *cancel_fn,
                               MPIX_Grequest_poll_function *poll_fn,
                               MPIX_Grequest_wait_function *wait_fn,
                               MPIX_Grequest_class *greq_class) ;
int MPIX_Grequest_class_allocate(MPIX_Grequest_class greq_class, void *extra_state,
                                 MPI_Request *request) ;
int MPI_Allgather_c(const void *sendbuf, MPI_Count sendcount, MPI_Datatype sendtype, void *recvbuf,
                    MPI_Count recvcount, MPI_Datatype recvtype, MPI_Comm comm)
                                                                                                                ;
int MPI_Allgather_init_c(const void *sendbuf, MPI_Count sendcount, MPI_Datatype sendtype,
                         void *recvbuf, MPI_Count recvcount, MPI_Datatype recvtype, MPI_Comm comm,
                         MPI_Info info, MPI_Request *request)
                                                                                                                     ;
int MPI_Allgatherv_c(const void *sendbuf, MPI_Count sendcount, MPI_Datatype sendtype, void *recvbuf,
                     const MPI_Count recvcounts[], const MPI_Aint displs[], MPI_Datatype recvtype,
                     MPI_Comm comm)
                                                                                                                 ;
int MPI_Allgatherv_init_c(const void *sendbuf, MPI_Count sendcount, MPI_Datatype sendtype,
                          void *recvbuf, const MPI_Count recvcounts[], const MPI_Aint displs[],
                          MPI_Datatype recvtype, MPI_Comm comm, MPI_Info info,
                          MPI_Request *request)
                                                                                                                      ;
int MPI_Allreduce_c(const void *sendbuf, void *recvbuf, MPI_Count count, MPI_Datatype datatype,
                    MPI_Op op, MPI_Comm comm)
                                                                                                                ;
int MPI_Allreduce_init_c(const void *sendbuf, void *recvbuf, MPI_Count count, MPI_Datatype datatype,
                         MPI_Op op, MPI_Comm comm, MPI_Info info, MPI_Request *request)
                                                                                                                     ;
int MPI_Alltoall_c(const void *sendbuf, MPI_Count sendcount, MPI_Datatype sendtype, void *recvbuf,
                   MPI_Count recvcount, MPI_Datatype recvtype, MPI_Comm comm)
                                                                                                               ;
int MPI_Alltoall_init_c(const void *sendbuf, MPI_Count sendcount, MPI_Datatype sendtype,
                        void *recvbuf, MPI_Count recvcount, MPI_Datatype recvtype, MPI_Comm comm,
                        MPI_Info info, MPI_Request *request)
                                                                                                                    ;
int MPI_Alltoallv_c(const void *sendbuf, const MPI_Count sendcounts[], const MPI_Aint sdispls[],
                    MPI_Datatype sendtype, void *recvbuf, const MPI_Count recvcounts[],
                    const MPI_Aint rdispls[], MPI_Datatype recvtype, MPI_Comm comm)
                                                                                                                ;
int MPI_Alltoallv_init_c(const void *sendbuf, const MPI_Count sendcounts[],
                         const MPI_Aint sdispls[], MPI_Datatype sendtype, void *recvbuf,
                         const MPI_Count recvcounts[], const MPI_Aint rdispls[],
                         MPI_Datatype recvtype, MPI_Comm comm, MPI_Info info, MPI_Request *request)
                                                                                                                     ;
int MPI_Alltoallw_c(const void *sendbuf, const MPI_Count sendcounts[], const MPI_Aint sdispls[],
                    const MPI_Datatype sendtypes[], void *recvbuf, const MPI_Count recvcounts[],
                    const MPI_Aint rdispls[], const MPI_Datatype recvtypes[], MPI_Comm comm)
                                    ;
int MPI_Alltoallw_init_c(const void *sendbuf, const MPI_Count sendcounts[],
                         const MPI_Aint sdispls[], const MPI_Datatype sendtypes[], void *recvbuf,
                         const MPI_Count recvcounts[], const MPI_Aint rdispls[],
                         const MPI_Datatype recvtypes[], MPI_Comm comm, MPI_Info info,
                         MPI_Request *request) ;
int MPI_Bcast_c(void *buffer, MPI_Count count, MPI_Datatype datatype, int root, MPI_Comm comm)
                                                          ;
int MPI_Bcast_init_c(void *buffer, MPI_Count count, MPI_Datatype datatype, int root, MPI_Comm comm,
                     MPI_Info info, MPI_Request *request)
                                                                           ;
int MPI_Exscan_c(const void *sendbuf, void *recvbuf, MPI_Count count, MPI_Datatype datatype,
                 MPI_Op op, MPI_Comm comm)
                                                                                                             ;
int MPI_Exscan_init_c(const void *sendbuf, void *recvbuf, MPI_Count count, MPI_Datatype datatype,
                      MPI_Op op, MPI_Comm comm, MPI_Info info, MPI_Request *request)
                                                                                                                  ;
int MPI_Gather_c(const void *sendbuf, MPI_Count sendcount, MPI_Datatype sendtype, void *recvbuf,
                 MPI_Count recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm)
                                                                                                             ;
int MPI_Gather_init_c(const void *sendbuf, MPI_Count sendcount, MPI_Datatype sendtype,
                      void *recvbuf, MPI_Count recvcount, MPI_Datatype recvtype, int root,
                      MPI_Comm comm, MPI_Info info, MPI_Request *request)
                                                                                                                  ;
int MPI_Gatherv_c(const void *sendbuf, MPI_Count sendcount, MPI_Datatype sendtype, void *recvbuf,
                  const MPI_Count recvcounts[], const MPI_Aint displs[], MPI_Datatype recvtype,
                  int root, MPI_Comm comm)
                                                                                                              ;
int MPI_Gatherv_init_c(const void *sendbuf, MPI_Count sendcount, MPI_Datatype sendtype,
                       void *recvbuf, const MPI_Count recvcounts[], const MPI_Aint displs[],
                       MPI_Datatype recvtype, int root, MPI_Comm comm, MPI_Info info,
                       MPI_Request *request)
                                                                                                                   ;
int MPI_Iallgather_c(const void *sendbuf, MPI_Count sendcount, MPI_Datatype sendtype, void *recvbuf,
                     MPI_Count recvcount, MPI_Datatype recvtype, MPI_Comm comm,
                     MPI_Request *request)
                                                                                                                 ;
int MPI_Iallgatherv_c(const void *sendbuf, MPI_Count sendcount, MPI_Datatype sendtype,
                      void *recvbuf, const MPI_Count recvcounts[], const MPI_Aint displs[],
                      MPI_Datatype recvtype, MPI_Comm comm, MPI_Request *request)
                                                                                                                  ;
int MPI_Iallreduce_c(const void *sendbuf, void *recvbuf, MPI_Count count, MPI_Datatype datatype,
                     MPI_Op op, MPI_Comm comm, MPI_Request *request)
                                                                                                                 ;
int MPI_Ialltoall_c(const void *sendbuf, MPI_Count sendcount, MPI_Datatype sendtype, void *recvbuf,
                    MPI_Count recvcount, MPI_Datatype recvtype, MPI_Comm comm,
                    MPI_Request *request)
                                                                                                                ;
int MPI_Ialltoallv_c(const void *sendbuf, const MPI_Count sendcounts[], const MPI_Aint sdispls[],
                     MPI_Datatype sendtype, void *recvbuf, const MPI_Count recvcounts[],
                     const MPI_Aint rdispls[], MPI_Datatype recvtype, MPI_Comm comm,
                     MPI_Request *request)
                                                                                                                 ;
int MPI_Ialltoallw_c(const void *sendbuf, const MPI_Count sendcounts[], const MPI_Aint sdispls[],
                     const MPI_Datatype sendtypes[], void *recvbuf, const MPI_Count recvcounts[],
                     const MPI_Aint rdispls[], const MPI_Datatype recvtypes[], MPI_Comm comm,
                     MPI_Request *request) ;
int MPI_Ibcast_c(void *buffer, MPI_Count count, MPI_Datatype datatype, int root, MPI_Comm comm,
                 MPI_Request *request) ;
int MPI_Iexscan_c(const void *sendbuf, void *recvbuf, MPI_Count count, MPI_Datatype datatype,
                  MPI_Op op, MPI_Comm comm, MPI_Request *request)
                                                                                                              ;
int MPI_Igather_c(const void *sendbuf, MPI_Count sendcount, MPI_Datatype sendtype, void *recvbuf,
                  MPI_Count recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm,
                  MPI_Request *request)
                                                                                                              ;
int MPI_Igatherv_c(const void *sendbuf, MPI_Count sendcount, MPI_Datatype sendtype, void *recvbuf,
                   const MPI_Count recvcounts[], const MPI_Aint displs[], MPI_Datatype recvtype,
                   int root, MPI_Comm comm, MPI_Request *request)
                                                                                                               ;
int MPI_Ineighbor_allgather_c(const void *sendbuf, MPI_Count sendcount, MPI_Datatype sendtype,
                              void *recvbuf, MPI_Count recvcount, MPI_Datatype recvtype,
                              MPI_Comm comm, MPI_Request *request)
                                                                                                                          ;
int MPI_Ineighbor_allgatherv_c(const void *sendbuf, MPI_Count sendcount, MPI_Datatype sendtype,
                               void *recvbuf, const MPI_Count recvcounts[], const MPI_Aint displs[],
                               MPI_Datatype recvtype, MPI_Comm comm, MPI_Request *request)
                                                                                                                           ;
int MPI_Ineighbor_alltoall_c(const void *sendbuf, MPI_Count sendcount, MPI_Datatype sendtype,
                             void *recvbuf, MPI_Count recvcount, MPI_Datatype recvtype,
                             MPI_Comm comm, MPI_Request *request)
                                                                                                                         ;
int MPI_Ineighbor_alltoallv_c(const void *sendbuf, const MPI_Count sendcounts[],
                              const MPI_Aint sdispls[], MPI_Datatype sendtype, void *recvbuf,
                              const MPI_Count recvcounts[], const MPI_Aint rdispls[],
                              MPI_Datatype recvtype, MPI_Comm comm, MPI_Request *request)
                                                                                                                          ;
int MPI_Ineighbor_alltoallw_c(const void *sendbuf, const MPI_Count sendcounts[],
                              const MPI_Aint sdispls[], const MPI_Datatype sendtypes[],
                              void *recvbuf, const MPI_Count recvcounts[], const MPI_Aint rdispls[],
                              const MPI_Datatype recvtypes[], MPI_Comm comm, MPI_Request *request)
                                              ;
int MPI_Ireduce_c(const void *sendbuf, void *recvbuf, MPI_Count count, MPI_Datatype datatype,
                  MPI_Op op, int root, MPI_Comm comm, MPI_Request *request)
                                                                                                              ;
int MPI_Ireduce_scatter_c(const void *sendbuf, void *recvbuf, const MPI_Count recvcounts[],
                          MPI_Datatype datatype, MPI_Op op, MPI_Comm comm, MPI_Request *request)
                                                                                                                      ;
int MPI_Ireduce_scatter_block_c(const void *sendbuf, void *recvbuf, MPI_Count recvcount,
                                MPI_Datatype datatype, MPI_Op op, MPI_Comm comm,
                                MPI_Request *request)
                                                                                                                            ;
int MPI_Iscan_c(const void *sendbuf, void *recvbuf, MPI_Count count, MPI_Datatype datatype,
                MPI_Op op, MPI_Comm comm, MPI_Request *request)
                                                                                                            ;
int MPI_Iscatter_c(const void *sendbuf, MPI_Count sendcount, MPI_Datatype sendtype, void *recvbuf,
                   MPI_Count recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm,
                   MPI_Request *request)
                                                                                                               ;
int MPI_Iscatterv_c(const void *sendbuf, const MPI_Count sendcounts[], const MPI_Aint displs[],
                    MPI_Datatype sendtype, void *recvbuf, MPI_Count recvcount,
                    MPI_Datatype recvtype, int root, MPI_Comm comm, MPI_Request *request)
                                                                                                                ;
int MPI_Neighbor_allgather_c(const void *sendbuf, MPI_Count sendcount, MPI_Datatype sendtype,
                             void *recvbuf, MPI_Count recvcount, MPI_Datatype recvtype,
                             MPI_Comm comm)
                                                                                                                         ;
int MPI_Neighbor_allgather_init_c(const void *sendbuf, MPI_Count sendcount, MPI_Datatype sendtype,
                                  void *recvbuf, MPI_Count recvcount, MPI_Datatype recvtype,
                                  MPI_Comm comm, MPI_Info info, MPI_Request *request)
                                                                                                                              ;
int MPI_Neighbor_allgatherv_c(const void *sendbuf, MPI_Count sendcount, MPI_Datatype sendtype,
                              void *recvbuf, const MPI_Count recvcounts[], const MPI_Aint displs[],
                              MPI_Datatype recvtype, MPI_Comm comm)
                                                                                                                          ;
int MPI_Neighbor_allgatherv_init_c(const void *sendbuf, MPI_Count sendcount, MPI_Datatype sendtype,
                                   void *recvbuf, const MPI_Count recvcounts[],
                                   const MPI_Aint displs[], MPI_Datatype recvtype, MPI_Comm comm,
                                   MPI_Info info, MPI_Request *request)
                                                                                                                               ;
int MPI_Neighbor_alltoall_c(const void *sendbuf, MPI_Count sendcount, MPI_Datatype sendtype,
                            void *recvbuf, MPI_Count recvcount, MPI_Datatype recvtype,
                            MPI_Comm comm)
                                                                                                                        ;
int MPI_Neighbor_alltoall_init_c(const void *sendbuf, MPI_Count sendcount, MPI_Datatype sendtype,
                                 void *recvbuf, MPI_Count recvcount, MPI_Datatype recvtype,
                                 MPI_Comm comm, MPI_Info info, MPI_Request *request)
                                                                                                                             ;
int MPI_Neighbor_alltoallv_c(const void *sendbuf, const MPI_Count sendcounts[],
                             const MPI_Aint sdispls[], MPI_Datatype sendtype, void *recvbuf,
                             const MPI_Count recvcounts[], const MPI_Aint rdispls[],
                             MPI_Datatype recvtype, MPI_Comm comm)
                                                                                                                         ;
int MPI_Neighbor_alltoallv_init_c(const void *sendbuf, const MPI_Count sendcounts[],
                                  const MPI_Aint sdispls[], MPI_Datatype sendtype, void *recvbuf,
                                  const MPI_Count recvcounts[], const MPI_Aint rdispls[],
                                  MPI_Datatype recvtype, MPI_Comm comm, MPI_Info info,
                                  MPI_Request *request)
                                                                                                                              ;
int MPI_Neighbor_alltoallw_c(const void *sendbuf, const MPI_Count sendcounts[],
                             const MPI_Aint sdispls[], const MPI_Datatype sendtypes[],
                             void *recvbuf, const MPI_Count recvcounts[], const MPI_Aint rdispls[],
                             const MPI_Datatype recvtypes[], MPI_Comm comm) ;
int MPI_Neighbor_alltoallw_init_c(const void *sendbuf, const MPI_Count sendcounts[],
                                  const MPI_Aint sdispls[], const MPI_Datatype sendtypes[],
                                  void *recvbuf, const MPI_Count recvcounts[],
                                  const MPI_Aint rdispls[], const MPI_Datatype recvtypes[],
                                  MPI_Comm comm, MPI_Info info, MPI_Request *request)
                                                  ;
int MPI_Reduce_c(const void *sendbuf, void *recvbuf, MPI_Count count, MPI_Datatype datatype,
                 MPI_Op op, int root, MPI_Comm comm)
                                                                                                             ;
int MPI_Reduce_init_c(const void *sendbuf, void *recvbuf, MPI_Count count, MPI_Datatype datatype,
                      MPI_Op op, int root, MPI_Comm comm, MPI_Info info, MPI_Request *request)
                                                                                                                  ;
int MPI_Reduce_local_c(const void *inbuf, void *inoutbuf, MPI_Count count, MPI_Datatype datatype,
                       MPI_Op op)
                                                                                                                   ;
int MPI_Reduce_scatter_c(const void *sendbuf, void *recvbuf, const MPI_Count recvcounts[],
                         MPI_Datatype datatype, MPI_Op op, MPI_Comm comm)
                                                                                                                     ;
int MPI_Reduce_scatter_block_c(const void *sendbuf, void *recvbuf, MPI_Count recvcount,
                               MPI_Datatype datatype, MPI_Op op, MPI_Comm comm)
                                                                                                                           ;
int MPI_Reduce_scatter_block_init_c(const void *sendbuf, void *recvbuf, MPI_Count recvcount,
                                    MPI_Datatype datatype, MPI_Op op, MPI_Comm comm, MPI_Info info,
                                    MPI_Request *request)
                                                                                                                                ;
int MPI_Reduce_scatter_init_c(const void *sendbuf, void *recvbuf, const MPI_Count recvcounts[],
                              MPI_Datatype datatype, MPI_Op op, MPI_Comm comm, MPI_Info info,
                              MPI_Request *request)
                                                                                                                          ;
int MPI_Scan_c(const void *sendbuf, void *recvbuf, MPI_Count count, MPI_Datatype datatype,
               MPI_Op op, MPI_Comm comm)
                                                                                                           ;
int MPI_Scan_init_c(const void *sendbuf, void *recvbuf, MPI_Count count, MPI_Datatype datatype,
                    MPI_Op op, MPI_Comm comm, MPI_Info info, MPI_Request *request)
                                                                                                                ;
int MPI_Scatter_c(const void *sendbuf, MPI_Count sendcount, MPI_Datatype sendtype, void *recvbuf,
                  MPI_Count recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm)
                                                                                                              ;
int MPI_Scatter_init_c(const void *sendbuf, MPI_Count sendcount, MPI_Datatype sendtype,
                       void *recvbuf, MPI_Count recvcount, MPI_Datatype recvtype, int root,
                       MPI_Comm comm, MPI_Info info, MPI_Request *request)
                                                                                                                   ;
int MPI_Scatterv_c(const void *sendbuf, const MPI_Count sendcounts[], const MPI_Aint displs[],
                   MPI_Datatype sendtype, void *recvbuf, MPI_Count recvcount, MPI_Datatype recvtype,
                   int root, MPI_Comm comm)
                                                                                                               ;
int MPI_Scatterv_init_c(const void *sendbuf, const MPI_Count sendcounts[], const MPI_Aint displs[],
                        MPI_Datatype sendtype, void *recvbuf, MPI_Count recvcount,
                        MPI_Datatype recvtype, int root, MPI_Comm comm, MPI_Info info,
                        MPI_Request *request)
                                                                                                                    ;
int MPI_Get_count_c(const MPI_Status *status, MPI_Datatype datatype, MPI_Count *count)
                    ;
int MPI_Get_elements_c(const MPI_Status *status, MPI_Datatype datatype, MPI_Count *count)
                    ;
int MPI_Pack_c(const void *inbuf, MPI_Count incount, MPI_Datatype datatype, void *outbuf,
               MPI_Count outsize, MPI_Count *position, MPI_Comm comm) ;
int MPI_Pack_external_c(const char *datarep, const void *inbuf, MPI_Count incount,
                        MPI_Datatype datatype, void *outbuf, MPI_Count outsize,
                        MPI_Count *position) ;
int MPI_Pack_external_size_c(const char *datarep, MPI_Count incount, MPI_Datatype datatype,
                             MPI_Count *size) ;
int MPI_Pack_size_c(MPI_Count incount, MPI_Datatype datatype, MPI_Comm comm, MPI_Count *size)
                    ;
int MPI_Status_set_elements_c(MPI_Status *status, MPI_Datatype datatype, MPI_Count count)
                    ;
int MPI_Type_contiguous_c(MPI_Count count, MPI_Datatype oldtype, MPI_Datatype *newtype)
                    ;
int MPI_Type_create_darray_c(int size, int rank, int ndims, const MPI_Count array_of_gsizes[],
                             const int array_of_distribs[], const int array_of_dargs[],
                             const int array_of_psizes[], int order, MPI_Datatype oldtype,
                             MPI_Datatype *newtype) ;
int MPI_Type_create_hindexed_c(MPI_Count count, const MPI_Count array_of_blocklengths[],
                               const MPI_Count array_of_displacements[], MPI_Datatype oldtype,
                               MPI_Datatype *newtype) ;
int MPI_Type_create_hindexed_block_c(MPI_Count count, MPI_Count blocklength,
                                     const MPI_Count array_of_displacements[], MPI_Datatype oldtype,
                                     MPI_Datatype *newtype) ;
int MPI_Type_create_hvector_c(MPI_Count count, MPI_Count blocklength, MPI_Count stride,
                              MPI_Datatype oldtype, MPI_Datatype *newtype) ;
int MPI_Type_create_indexed_block_c(MPI_Count count, MPI_Count blocklength,
                                    const MPI_Count array_of_displacements[], MPI_Datatype oldtype,
                                    MPI_Datatype *newtype) ;
int MPI_Type_create_resized_c(MPI_Datatype oldtype, MPI_Count lb, MPI_Count extent,
                              MPI_Datatype *newtype) ;
int MPI_Type_create_struct_c(MPI_Count count, const MPI_Count array_of_blocklengths[],
                             const MPI_Count array_of_displacements[],
                             const MPI_Datatype array_of_types[], MPI_Datatype *newtype)
                                             ;
int MPI_Type_create_subarray_c(int ndims, const MPI_Count array_of_sizes[],
                               const MPI_Count array_of_subsizes[],
                               const MPI_Count array_of_starts[], int order, MPI_Datatype oldtype,
                               MPI_Datatype *newtype) ;
int MPI_Type_get_contents_c(MPI_Datatype datatype, MPI_Count max_integers, MPI_Count max_addresses,
                            MPI_Count max_large_counts, MPI_Count max_datatypes,
                            int array_of_integers[], MPI_Aint array_of_addresses[],
                            MPI_Count array_of_large_counts[], MPI_Datatype array_of_datatypes[])
                                            ;
int MPI_Type_get_envelope_c(MPI_Datatype datatype, MPI_Count *num_integers,
                            MPI_Count *num_addresses, MPI_Count *num_large_counts,
                            MPI_Count *num_datatypes, int *combiner) ;
int MPI_Type_get_extent_c(MPI_Datatype datatype, MPI_Count *lb, MPI_Count *extent)
                    ;
int MPI_Type_get_true_extent_c(MPI_Datatype datatype, MPI_Count *true_lb, MPI_Count *true_extent)
                    ;
int MPI_Type_indexed_c(MPI_Count count, const MPI_Count array_of_blocklengths[],
                       const MPI_Count array_of_displacements[], MPI_Datatype oldtype,
                       MPI_Datatype *newtype) ;
int MPI_Type_size_c(MPI_Datatype datatype, MPI_Count *size) ;
int MPI_Type_vector_c(MPI_Count count, MPI_Count blocklength, MPI_Count stride,
                      MPI_Datatype oldtype, MPI_Datatype *newtype) ;
int MPI_Unpack_c(const void *inbuf, MPI_Count insize, MPI_Count *position, void *outbuf,
                 MPI_Count outcount, MPI_Datatype datatype, MPI_Comm comm) ;
int MPI_Unpack_external_c(const char datarep[], const void *inbuf, MPI_Count insize,
                          MPI_Count *position, void *outbuf, MPI_Count outcount,
                          MPI_Datatype datatype) ;
int MPI_Op_create_c(MPI_User_function_c *user_fn, int commute, MPI_Op *op) ;
int MPI_Bsend_c(const void *buf, MPI_Count count, MPI_Datatype datatype, int dest, int tag,
                MPI_Comm comm) ;
int MPI_Bsend_init_c(const void *buf, MPI_Count count, MPI_Datatype datatype, int dest, int tag,
                     MPI_Comm comm, MPI_Request *request)
                                                                           ;
int MPI_Buffer_attach_c(void *buffer, MPI_Count size) ;
int MPI_Buffer_detach_c(void *buffer_addr, MPI_Count *size) ;
int MPI_Comm_attach_buffer_c(MPI_Comm comm, void *buffer, MPI_Count size) ;
int MPI_Comm_detach_buffer_c(MPI_Comm comm, void *buffer_addr, MPI_Count *size) ;
int MPI_Ibsend_c(const void *buf, MPI_Count count, MPI_Datatype datatype, int dest, int tag,
                 MPI_Comm comm, MPI_Request *request)
                                                                       ;
int MPI_Imrecv_c(void *buf, MPI_Count count, MPI_Datatype datatype, MPI_Message *message,
                 MPI_Request *request) ;
int MPI_Irecv_c(void *buf, MPI_Count count, MPI_Datatype datatype, int source, int tag,
                MPI_Comm comm, MPI_Request *request)
                                                                      ;
int MPI_Irsend_c(const void *buf, MPI_Count count, MPI_Datatype datatype, int dest, int tag,
                 MPI_Comm comm, MPI_Request *request)
                                                                       ;
int MPI_Isend_c(const void *buf, MPI_Count count, MPI_Datatype datatype, int dest, int tag,
                MPI_Comm comm, MPI_Request *request)
                                                                      ;
int MPI_Isendrecv_c(const void *sendbuf, MPI_Count sendcount, MPI_Datatype sendtype, int dest,
                    int sendtag, void *recvbuf, MPI_Count recvcount, MPI_Datatype recvtype,
                    int source, int recvtag, MPI_Comm comm, MPI_Request *request)
                                                                                                                ;
int MPI_Isendrecv_replace_c(void *buf, MPI_Count count, MPI_Datatype datatype, int dest,
                            int sendtag, int source, int recvtag, MPI_Comm comm,
                            MPI_Request *request)
                                                                                  ;
int MPI_Issend_c(const void *buf, MPI_Count count, MPI_Datatype datatype, int dest, int tag,
                 MPI_Comm comm, MPI_Request *request)
                                                                       ;
int MPI_Mrecv_c(void *buf, MPI_Count count, MPI_Datatype datatype, MPI_Message *message,
                MPI_Status *status) ;
int MPI_Recv_c(void *buf, MPI_Count count, MPI_Datatype datatype, int source, int tag,
               MPI_Comm comm, MPI_Status *status)
                                                                     ;
int MPI_Recv_init_c(void *buf, MPI_Count count, MPI_Datatype datatype, int source, int tag,
                    MPI_Comm comm, MPI_Request *request)
                                                                          ;
int MPI_Rsend_c(const void *buf, MPI_Count count, MPI_Datatype datatype, int dest, int tag,
                MPI_Comm comm) ;
int MPI_Rsend_init_c(const void *buf, MPI_Count count, MPI_Datatype datatype, int dest, int tag,
                     MPI_Comm comm, MPI_Request *request)
                                                                           ;
int MPI_Send_c(const void *buf, MPI_Count count, MPI_Datatype datatype, int dest, int tag,
               MPI_Comm comm) ;
int MPI_Send_init_c(const void *buf, MPI_Count count, MPI_Datatype datatype, int dest, int tag,
                    MPI_Comm comm, MPI_Request *request)
                                                                          ;
int MPI_Sendrecv_c(const void *sendbuf, MPI_Count sendcount, MPI_Datatype sendtype, int dest,
                   int sendtag, void *recvbuf, MPI_Count recvcount, MPI_Datatype recvtype,
                   int source, int recvtag, MPI_Comm comm, MPI_Status *status)
                                                                                                               ;
int MPI_Sendrecv_replace_c(void *buf, MPI_Count count, MPI_Datatype datatype, int dest, int sendtag,
                           int source, int recvtag, MPI_Comm comm, MPI_Status *status)
                                                                                 ;
int MPI_Session_attach_buffer_c(MPI_Session session, void *buffer, MPI_Count size)
                    ;
int MPI_Session_detach_buffer_c(MPI_Session session, void *buffer_addr, MPI_Count *size)
                    ;
int MPI_Ssend_c(const void *buf, MPI_Count count, MPI_Datatype datatype, int dest, int tag,
                MPI_Comm comm) ;
int MPI_Ssend_init_c(const void *buf, MPI_Count count, MPI_Datatype datatype, int dest, int tag,
                     MPI_Comm comm, MPI_Request *request)
                                                                           ;
int MPI_Accumulate_c(const void *origin_addr, MPI_Count origin_count, MPI_Datatype origin_datatype,
                     int target_rank, MPI_Aint target_disp, MPI_Count target_count,
                     MPI_Datatype target_datatype, MPI_Op op, MPI_Win win)
                                                                           ;
int MPI_Get_c(void *origin_addr, MPI_Count origin_count, MPI_Datatype origin_datatype,
              int target_rank, MPI_Aint target_disp, MPI_Count target_count,
              MPI_Datatype target_datatype, MPI_Win win)
                                                                    ;
int MPI_Get_accumulate_c(const void *origin_addr, MPI_Count origin_count,
                         MPI_Datatype origin_datatype, void *result_addr, MPI_Count result_count,
                         MPI_Datatype result_datatype, int target_rank, MPI_Aint target_disp,
                         MPI_Count target_count, MPI_Datatype target_datatype, MPI_Op op,
                         MPI_Win win)
                                                                                                                     ;
int MPI_Put_c(const void *origin_addr, MPI_Count origin_count, MPI_Datatype origin_datatype,
              int target_rank, MPI_Aint target_disp, MPI_Count target_count,
              MPI_Datatype target_datatype, MPI_Win win)
                                                                    ;
int MPI_Raccumulate_c(const void *origin_addr, MPI_Count origin_count, MPI_Datatype origin_datatype,
                      int target_rank, MPI_Aint target_disp, MPI_Count target_count,
                      MPI_Datatype target_datatype, MPI_Op op, MPI_Win win, MPI_Request *request)
                                                                            ;
int MPI_Rget_c(void *origin_addr, MPI_Count origin_count, MPI_Datatype origin_datatype,
               int target_rank, MPI_Aint target_disp, MPI_Count target_count,
               MPI_Datatype target_datatype, MPI_Win win, MPI_Request *request)
                                                                     ;
int MPI_Rget_accumulate_c(const void *origin_addr, MPI_Count origin_count,
                          MPI_Datatype origin_datatype, void *result_addr, MPI_Count result_count,
                          MPI_Datatype result_datatype, int target_rank, MPI_Aint target_disp,
                          MPI_Count target_count, MPI_Datatype target_datatype, MPI_Op op,
                          MPI_Win win, MPI_Request *request)
                                                                                                                      ;
int MPI_Rput_c(const void *origin_addr, MPI_Count origin_count, MPI_Datatype origin_datatype,
               int target_rank, MPI_Aint target_disp, MPI_Count target_count,
               MPI_Datatype target_datatype, MPI_Win win, MPI_Request *request)
                                                                     ;
int MPI_Win_allocate_c(MPI_Aint size, MPI_Aint disp_unit, MPI_Info info, MPI_Comm comm,
                       void *baseptr, MPI_Win *win) ;
int MPI_Win_allocate_shared_c(MPI_Aint size, MPI_Aint disp_unit, MPI_Info info, MPI_Comm comm,
                              void *baseptr, MPI_Win *win) ;
int MPI_Win_create_c(void *base, MPI_Aint size, MPI_Aint disp_unit, MPI_Info info, MPI_Comm comm,
                     MPI_Win *win) ;
int MPI_Win_shared_query_c(MPI_Win win, int rank, MPI_Aint *size, MPI_Aint *disp_unit,
                           void *baseptr) ;
int MPIX_Stream_send_c(const void *buf, MPI_Count count, MPI_Datatype datatype, int dest, int tag,
                       MPI_Comm comm, int source_stream_index, int dest_stream_index)
                                                                             ;
int MPIX_Stream_isend_c(const void *buf, MPI_Count count, MPI_Datatype datatype, int dest, int tag,
                        MPI_Comm comm, int source_stream_index, int dest_stream_index,
                        MPI_Request *request)
                                                                              ;
int MPIX_Stream_recv_c(void *buf, MPI_Count count, MPI_Datatype datatype, int source, int tag,
                       MPI_Comm comm, int source_stream_index, int dest_stream_index,
                       MPI_Status *status) ;
int MPIX_Stream_irecv_c(void *buf, MPI_Count count, MPI_Datatype datatype, int source, int tag,
                        MPI_Comm comm, int source_stream_index, int dest_stream_index,
                        MPI_Request *request)
                                                                              ;
int MPIX_Send_enqueue_c(const void *buf, MPI_Count count, MPI_Datatype datatype, int dest, int tag,
                        MPI_Comm comm) ;
int MPIX_Recv_enqueue_c(void *buf, MPI_Count count, MPI_Datatype datatype, int source, int tag,
                        MPI_Comm comm, MPI_Status *status)
                                                                              ;
int MPIX_Isend_enqueue_c(const void *buf, MPI_Count count, MPI_Datatype datatype, int dest, int tag,
                         MPI_Comm comm, MPI_Request *request)
                                                                               ;
int MPIX_Irecv_enqueue_c(void *buf, MPI_Count count, MPI_Datatype datatype, int source, int tag,
                         MPI_Comm comm, MPI_Request *request)
                                                                               ;
int MPIX_Allreduce_enqueue_c(const void *sendbuf, void *recvbuf, MPI_Count count,
                             MPI_Datatype datatype, MPI_Op op, MPI_Comm comm)
                                                                                   ;
int MPI_File_get_type_extent_c(MPI_File fh, MPI_Datatype datatype, MPI_Count *extent)
                    ;
int MPI_File_iread_c(MPI_File fh, void *buf, MPI_Count count, MPI_Datatype datatype,
                     MPI_Request *request) ;
int MPI_File_iread_all_c(MPI_File fh, void *buf, MPI_Count count, MPI_Datatype datatype,
                         MPI_Request *request)
                                                                               ;
int MPI_File_iread_at_c(MPI_File fh, MPI_Offset offset, void *buf, MPI_Count count,
                        MPI_Datatype datatype, MPI_Request *request)
                                                                              ;
int MPI_File_iread_at_all_c(MPI_File fh, MPI_Offset offset, void *buf, MPI_Count count,
                            MPI_Datatype datatype, MPI_Request *request)
                                                                                  ;
int MPI_File_iread_shared_c(MPI_File fh, void *buf, MPI_Count count, MPI_Datatype datatype,
                            MPI_Request *request)
                                                                                  ;
int MPI_File_iwrite_c(MPI_File fh, const void *buf, MPI_Count count, MPI_Datatype datatype,
                      MPI_Request *request) ;
int MPI_File_iwrite_all_c(MPI_File fh, const void *buf, MPI_Count count, MPI_Datatype datatype,
                          MPI_Request *request)
                                                                                ;
int MPI_File_iwrite_at_c(MPI_File fh, MPI_Offset offset, const void *buf, MPI_Count count,
                         MPI_Datatype datatype, MPI_Request *request)
                                                                               ;
int MPI_File_iwrite_at_all_c(MPI_File fh, MPI_Offset offset, const void *buf, MPI_Count count,
                             MPI_Datatype datatype, MPI_Request *request)
                                                                                   ;
int MPI_File_iwrite_shared_c(MPI_File fh, const void *buf, MPI_Count count, MPI_Datatype datatype,
                             MPI_Request *request)
                                                                                   ;
int MPI_File_read_c(MPI_File fh, void *buf, MPI_Count count, MPI_Datatype datatype,
                    MPI_Status *status) ;
int MPI_File_read_all_c(MPI_File fh, void *buf, MPI_Count count, MPI_Datatype datatype,
                        MPI_Status *status) ;
int MPI_File_read_all_begin_c(MPI_File fh, void *buf, MPI_Count count, MPI_Datatype datatype)
                                                          ;
int MPI_File_read_at_c(MPI_File fh, MPI_Offset offset, void *buf, MPI_Count count,
                       MPI_Datatype datatype, MPI_Status *status)
                                                                             ;
int MPI_File_read_at_all_c(MPI_File fh, MPI_Offset offset, void *buf, MPI_Count count,
                           MPI_Datatype datatype, MPI_Status *status)
                                                                                 ;
int MPI_File_read_at_all_begin_c(MPI_File fh, MPI_Offset offset, void *buf, MPI_Count count,
                                 MPI_Datatype datatype)
                                                                                       ;
int MPI_File_read_ordered_c(MPI_File fh, void *buf, MPI_Count count, MPI_Datatype datatype,
                            MPI_Status *status)
                                                                                  ;
int MPI_File_read_ordered_begin_c(MPI_File fh, void *buf, MPI_Count count, MPI_Datatype datatype)
                                                          ;
int MPI_File_read_shared_c(MPI_File fh, void *buf, MPI_Count count, MPI_Datatype datatype,
                           MPI_Status *status)
                                                                                 ;
int MPI_File_write_c(MPI_File fh, const void *buf, MPI_Count count, MPI_Datatype datatype,
                     MPI_Status *status) ;
int MPI_File_write_all_c(MPI_File fh, const void *buf, MPI_Count count, MPI_Datatype datatype,
                         MPI_Status *status)
                                                                               ;
int MPI_File_write_all_begin_c(MPI_File fh, const void *buf, MPI_Count count,
                               MPI_Datatype datatype)
                                                                                     ;
int MPI_File_write_at_c(MPI_File fh, MPI_Offset offset, const void *buf, MPI_Count count,
                        MPI_Datatype datatype, MPI_Status *status)
                                                                              ;
int MPI_File_write_at_all_c(MPI_File fh, MPI_Offset offset, const void *buf, MPI_Count count,
                            MPI_Datatype datatype, MPI_Status *status)
                                                                                  ;
int MPI_File_write_at_all_begin_c(MPI_File fh, MPI_Offset offset, const void *buf, MPI_Count count,
                                  MPI_Datatype datatype)
                                                                                        ;
int MPI_File_write_ordered_c(MPI_File fh, const void *buf, MPI_Count count, MPI_Datatype datatype,
                             MPI_Status *status)
                                                                                   ;
int MPI_File_write_ordered_begin_c(MPI_File fh, const void *buf, MPI_Count count,
                                   MPI_Datatype datatype)
                                                                                         ;
int MPI_File_write_shared_c(MPI_File fh, const void *buf, MPI_Count count, MPI_Datatype datatype,
                            MPI_Status *status)
                                                                                  ;
int MPI_Register_datarep_c(const char *datarep,
                           MPI_Datarep_conversion_function_c *read_conversion_fn,
                           MPI_Datarep_conversion_function_c *write_conversion_fn,
                           MPI_Datarep_extent_function *dtype_file_extent_fn, void *extra_state)
                                           ;
int PMPI_Abi_get_fortran_booleans(int logical_size, void *logical_true, void *logical_false,
                                  int *is_set) ;
int PMPI_Abi_get_fortran_info(MPI_Info *info) ;
int PMPI_Abi_get_info(MPI_Info *info) ;
int PMPI_Abi_get_version(int *abi_major, int *abi_minor) ;
int PMPI_Abi_set_fortran_booleans(int logical_size, void *logical_true, void *logical_false)
                    ;
int PMPI_Abi_set_fortran_info(MPI_Info info) ;
int PMPI_Comm_toint(MPI_Comm comm) ;
MPI_Comm PMPI_Comm_fromint(int comm) ;
int PMPI_Errhandler_toint(MPI_Errhandler errhandler) ;
MPI_Errhandler PMPI_Errhandler_fromint(int errhandler) ;
int PMPI_Group_toint(MPI_Group group) ;
MPI_Group PMPI_Group_fromint(int group) ;
int PMPI_Info_toint(MPI_Info info) ;
MPI_Info PMPI_Info_fromint(int info) ;
int PMPI_Message_toint(MPI_Message message) ;
MPI_Message PMPI_Message_fromint(int message) ;
int PMPI_Op_toint(MPI_Op op) ;
MPI_Op PMPI_Op_fromint(int op) ;
int PMPI_Request_toint(MPI_Request request) ;
MPI_Request PMPI_Request_fromint(int request) ;
int PMPI_Session_toint(MPI_Session session) ;
MPI_Session PMPI_Session_fromint(int session) ;
int PMPI_Type_toint(MPI_Datatype datatype) ;
MPI_Datatype PMPI_Type_fromint(int datatype) ;
int PMPI_Win_toint(MPI_Win win) ;
MPI_Win PMPI_Win_fromint(int win) ;
int PMPIX_Async_start(MPIX_Async_poll_function *poll_fn, void *extra_state, MPIX_Stream stream)
                    ;
void * PMPIX_Async_get_state(MPIX_Async_thing async_thing) ;
int PMPIX_Async_spawn(MPIX_Async_thing async_thing, MPIX_Async_poll_function *poll_fn,
                      void *extra_state, MPIX_Stream stream) ;
int PMPI_Comm_create_keyval(MPI_Comm_copy_attr_function *comm_copy_attr_fn,
                            MPI_Comm_delete_attr_function *comm_delete_attr_fn, int *comm_keyval,
                            void *extra_state) ;
int PMPI_Keyval_create(MPI_Copy_function *copy_fn, MPI_Delete_function *delete_fn, int *keyval,
                       void *extra_state) ;
int PMPI_Comm_delete_attr(MPI_Comm comm, int comm_keyval) ;
int PMPI_Attr_delete(MPI_Comm comm, int keyval) ;
int PMPI_Comm_free_keyval(int *comm_keyval) ;
int PMPI_Keyval_free(int *keyval) ;
int PMPI_Comm_get_attr(MPI_Comm comm, int comm_keyval, void *attribute_val, int *flag)
                    ;
int PMPI_Attr_get(MPI_Comm comm, int keyval, void *attribute_val, int *flag) ;
int PMPI_Comm_set_attr(MPI_Comm comm, int comm_keyval, void *attribute_val) ;
int PMPI_Attr_put(MPI_Comm comm, int keyval, void *attribute_val) ;
int PMPI_Type_create_keyval(MPI_Type_copy_attr_function *type_copy_attr_fn,
                            MPI_Type_delete_attr_function *type_delete_attr_fn, int *type_keyval,
                            void *extra_state) ;
int PMPI_Type_delete_attr(MPI_Datatype datatype, int type_keyval) ;
int PMPI_Type_free_keyval(int *type_keyval) ;
int PMPI_Type_get_attr(MPI_Datatype datatype, int type_keyval, void *attribute_val, int *flag)
                    ;
int PMPI_Type_set_attr(MPI_Datatype datatype, int type_keyval, void *attribute_val)
                    ;
int PMPI_Win_create_keyval(MPI_Win_copy_attr_function *win_copy_attr_fn,
                           MPI_Win_delete_attr_function *win_delete_attr_fn, int *win_keyval,
                           void *extra_state) ;
int PMPI_Win_delete_attr(MPI_Win win, int win_keyval) ;
int PMPI_Win_free_keyval(int *win_keyval) ;
int PMPI_Win_get_attr(MPI_Win win, int win_keyval, void *attribute_val, int *flag)
                    ;
int PMPI_Win_set_attr(MPI_Win win, int win_keyval, void *attribute_val) ;
int PMPIX_Op_create_x(MPIX_User_function_x *user_fn_x, MPIX_Destructor_function *destructor_fn,
                      int commute, void *extra_state, MPI_Op *op) ;
int PMPIX_Comm_create_errhandler_x(MPIX_Comm_errhandler_function_x *comm_errhandler_fn_x,
                                   MPIX_Destructor_function *destructor_fn, void *extra_state,
                                   MPI_Errhandler *errhandler) ;
int PMPIX_Win_create_errhandler_x(MPIX_Win_errhandler_function_x *comm_errhandler_fn_x,
                                  MPIX_Destructor_function *destructor_fn, void *extra_state,
                                  MPI_Errhandler *errhandler) ;
int PMPIX_File_create_errhandler_x(MPIX_File_errhandler_function_x *comm_errhandler_fn_x,
                                   MPIX_Destructor_function *destructor_fn, void *extra_state,
                                   MPI_Errhandler *errhandler) ;
int PMPIX_Session_create_errhandler_x(MPIX_Session_errhandler_function_x *comm_errhandler_fn_x,
                                      MPIX_Destructor_function *destructor_fn, void *extra_state,
                                      MPI_Errhandler *errhandler) ;
int PMPI_Allgather(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
                   int recvcount, MPI_Datatype recvtype, MPI_Comm comm)
                                                                                                               ;
int PMPI_Allgather_c(const void *sendbuf, MPI_Count sendcount, MPI_Datatype sendtype, void *recvbuf,
                     MPI_Count recvcount, MPI_Datatype recvtype, MPI_Comm comm)
                                                                                                                 ;
int PMPI_Allgather_init(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
                        int recvcount, MPI_Datatype recvtype, MPI_Comm comm, MPI_Info info,
                        MPI_Request *request)
                                                                                                                    ;
int PMPI_Allgather_init_c(const void *sendbuf, MPI_Count sendcount, MPI_Datatype sendtype,
                          void *recvbuf, MPI_Count recvcount, MPI_Datatype recvtype, MPI_Comm comm,
                          MPI_Info info, MPI_Request *request)
                                                                                                                      ;
int PMPI_Allgatherv(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
                    const int recvcounts[], const int displs[], MPI_Datatype recvtype,
                    MPI_Comm comm)
                                                                                                                ;
int PMPI_Allgatherv_c(const void *sendbuf, MPI_Count sendcount, MPI_Datatype sendtype,
                      void *recvbuf, const MPI_Count recvcounts[], const MPI_Aint displs[],
                      MPI_Datatype recvtype, MPI_Comm comm)
                                                                                                                  ;
int PMPI_Allgatherv_init(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
                         const int recvcounts[], const int displs[], MPI_Datatype recvtype,
                         MPI_Comm comm, MPI_Info info, MPI_Request *request)
                                                                                                                     ;
int PMPI_Allgatherv_init_c(const void *sendbuf, MPI_Count sendcount, MPI_Datatype sendtype,
                           void *recvbuf, const MPI_Count recvcounts[], const MPI_Aint displs[],
                           MPI_Datatype recvtype, MPI_Comm comm, MPI_Info info,
                           MPI_Request *request)
                                                                                                                       ;
int PMPI_Allreduce(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op,
                   MPI_Comm comm)
                                                                                                               ;
int PMPI_Allreduce_c(const void *sendbuf, void *recvbuf, MPI_Count count, MPI_Datatype datatype,
                     MPI_Op op, MPI_Comm comm)
                                                                                                                 ;
int PMPI_Allreduce_init(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype,
                        MPI_Op op, MPI_Comm comm, MPI_Info info, MPI_Request *request)
                                                                                                                    ;
int PMPI_Allreduce_init_c(const void *sendbuf, void *recvbuf, MPI_Count count,
                          MPI_Datatype datatype, MPI_Op op, MPI_Comm comm, MPI_Info info,
                          MPI_Request *request)
                                                                                                                      ;
int PMPI_Alltoall(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
                  int recvcount, MPI_Datatype recvtype, MPI_Comm comm)
                                                                                                              ;
int PMPI_Alltoall_c(const void *sendbuf, MPI_Count sendcount, MPI_Datatype sendtype, void *recvbuf,
                    MPI_Count recvcount, MPI_Datatype recvtype, MPI_Comm comm)
                                                                                                                ;
int PMPI_Alltoall_init(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
                       int recvcount, MPI_Datatype recvtype, MPI_Comm comm, MPI_Info info,
                       MPI_Request *request)
                                                                                                                   ;
int PMPI_Alltoall_init_c(const void *sendbuf, MPI_Count sendcount, MPI_Datatype sendtype,
                         void *recvbuf, MPI_Count recvcount, MPI_Datatype recvtype, MPI_Comm comm,
                         MPI_Info info, MPI_Request *request)
                                                                                                                     ;
int PMPI_Alltoallv(const void *sendbuf, const int sendcounts[], const int sdispls[],
                   MPI_Datatype sendtype, void *recvbuf, const int recvcounts[],
                   const int rdispls[], MPI_Datatype recvtype, MPI_Comm comm)
                                                                                                               ;
int PMPI_Alltoallv_c(const void *sendbuf, const MPI_Count sendcounts[], const MPI_Aint sdispls[],
                     MPI_Datatype sendtype, void *recvbuf, const MPI_Count recvcounts[],
                     const MPI_Aint rdispls[], MPI_Datatype recvtype, MPI_Comm comm)
                                                                                                                 ;
int PMPI_Alltoallv_init(const void *sendbuf, const int sendcounts[], const int sdispls[],
                        MPI_Datatype sendtype, void *recvbuf, const int recvcounts[],
                        const int rdispls[], MPI_Datatype recvtype, MPI_Comm comm, MPI_Info info,
                        MPI_Request *request)
                                                                                                                    ;
int PMPI_Alltoallv_init_c(const void *sendbuf, const MPI_Count sendcounts[],
                          const MPI_Aint sdispls[], MPI_Datatype sendtype, void *recvbuf,
                          const MPI_Count recvcounts[], const MPI_Aint rdispls[],
                          MPI_Datatype recvtype, MPI_Comm comm, MPI_Info info,
                          MPI_Request *request)
                                                                                                                      ;
int PMPI_Alltoallw(const void *sendbuf, const int sendcounts[], const int sdispls[],
                   const MPI_Datatype sendtypes[], void *recvbuf, const int recvcounts[],
                   const int rdispls[], const MPI_Datatype recvtypes[], MPI_Comm comm)
                                   ;
int PMPI_Alltoallw_c(const void *sendbuf, const MPI_Count sendcounts[], const MPI_Aint sdispls[],
                     const MPI_Datatype sendtypes[], void *recvbuf, const MPI_Count recvcounts[],
                     const MPI_Aint rdispls[], const MPI_Datatype recvtypes[], MPI_Comm comm)
                                     ;
int PMPI_Alltoallw_init(const void *sendbuf, const int sendcounts[], const int sdispls[],
                        const MPI_Datatype sendtypes[], void *recvbuf, const int recvcounts[],
                        const int rdispls[], const MPI_Datatype recvtypes[], MPI_Comm comm,
                        MPI_Info info, MPI_Request *request) ;
int PMPI_Alltoallw_init_c(const void *sendbuf, const MPI_Count sendcounts[],
                          const MPI_Aint sdispls[], const MPI_Datatype sendtypes[], void *recvbuf,
                          const MPI_Count recvcounts[], const MPI_Aint rdispls[],
                          const MPI_Datatype recvtypes[], MPI_Comm comm, MPI_Info info,
                          MPI_Request *request) ;
int PMPI_Barrier(MPI_Comm comm) ;
int PMPI_Barrier_init(MPI_Comm comm, MPI_Info info, MPI_Request *request) ;
int PMPI_Bcast(void *buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm)
                                                          ;
int PMPI_Bcast_c(void *buffer, MPI_Count count, MPI_Datatype datatype, int root, MPI_Comm comm)
                                                          ;
int PMPI_Bcast_init(void *buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm,
                    MPI_Info info, MPI_Request *request)
                                                                          ;
int PMPI_Bcast_init_c(void *buffer, MPI_Count count, MPI_Datatype datatype, int root, MPI_Comm comm,
                      MPI_Info info, MPI_Request *request)
                                                                            ;
int PMPI_Exscan(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op,
                MPI_Comm comm)
                                                                                                            ;
int PMPI_Exscan_c(const void *sendbuf, void *recvbuf, MPI_Count count, MPI_Datatype datatype,
                  MPI_Op op, MPI_Comm comm)
                                                                                                              ;
int PMPI_Exscan_init(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype,
                     MPI_Op op, MPI_Comm comm, MPI_Info info, MPI_Request *request)
                                                                                                                 ;
int PMPI_Exscan_init_c(const void *sendbuf, void *recvbuf, MPI_Count count, MPI_Datatype datatype,
                       MPI_Op op, MPI_Comm comm, MPI_Info info, MPI_Request *request)
                                                                                                                   ;
int PMPI_Gather(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
                int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm)
                                                                                                            ;
int PMPI_Gather_c(const void *sendbuf, MPI_Count sendcount, MPI_Datatype sendtype, void *recvbuf,
                  MPI_Count recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm)
                                                                                                              ;
int PMPI_Gather_init(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
                     int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm, MPI_Info info,
                     MPI_Request *request)
                                                                                                                 ;
int PMPI_Gather_init_c(const void *sendbuf, MPI_Count sendcount, MPI_Datatype sendtype,
                       void *recvbuf, MPI_Count recvcount, MPI_Datatype recvtype, int root,
                       MPI_Comm comm, MPI_Info info, MPI_Request *request)
                                                                                                                   ;
int PMPI_Gatherv(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
                 const int recvcounts[], const int displs[], MPI_Datatype recvtype, int root,
                 MPI_Comm comm)
                                                                                                             ;
int PMPI_Gatherv_c(const void *sendbuf, MPI_Count sendcount, MPI_Datatype sendtype, void *recvbuf,
                   const MPI_Count recvcounts[], const MPI_Aint displs[], MPI_Datatype recvtype,
                   int root, MPI_Comm comm)
                                                                                                               ;
int PMPI_Gatherv_init(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
                      const int recvcounts[], const int displs[], MPI_Datatype recvtype, int root,
                      MPI_Comm comm, MPI_Info info, MPI_Request *request)
                                                                                                                  ;
int PMPI_Gatherv_init_c(const void *sendbuf, MPI_Count sendcount, MPI_Datatype sendtype,
                        void *recvbuf, const MPI_Count recvcounts[], const MPI_Aint displs[],
                        MPI_Datatype recvtype, int root, MPI_Comm comm, MPI_Info info,
                        MPI_Request *request)
                                                                                                                    ;
int PMPI_Iallgather(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
                    int recvcount, MPI_Datatype recvtype, MPI_Comm comm, MPI_Request *request)
                                                                                                                ;
int PMPI_Iallgather_c(const void *sendbuf, MPI_Count sendcount, MPI_Datatype sendtype,
                      void *recvbuf, MPI_Count recvcount, MPI_Datatype recvtype, MPI_Comm comm,
                      MPI_Request *request)
                                                                                                                  ;
int PMPI_Iallgatherv(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
                     const int recvcounts[], const int displs[], MPI_Datatype recvtype,
                     MPI_Comm comm, MPI_Request *request)
                                                                                                                 ;
int PMPI_Iallgatherv_c(const void *sendbuf, MPI_Count sendcount, MPI_Datatype sendtype,
                       void *recvbuf, const MPI_Count recvcounts[], const MPI_Aint displs[],
                       MPI_Datatype recvtype, MPI_Comm comm, MPI_Request *request)
                                                                                                                   ;
int PMPI_Iallreduce(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op,
                    MPI_Comm comm, MPI_Request *request)
                                                                                                                ;
int PMPI_Iallreduce_c(const void *sendbuf, void *recvbuf, MPI_Count count, MPI_Datatype datatype,
                      MPI_Op op, MPI_Comm comm, MPI_Request *request)
                                                                                                                  ;
int PMPI_Ialltoall(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
                   int recvcount, MPI_Datatype recvtype, MPI_Comm comm, MPI_Request *request)
                                                                                                               ;
int PMPI_Ialltoall_c(const void *sendbuf, MPI_Count sendcount, MPI_Datatype sendtype, void *recvbuf,
                     MPI_Count recvcount, MPI_Datatype recvtype, MPI_Comm comm,
                     MPI_Request *request)
                                                                                                                 ;
int PMPI_Ialltoallv(const void *sendbuf, const int sendcounts[], const int sdispls[],
                    MPI_Datatype sendtype, void *recvbuf, const int recvcounts[],
                    const int rdispls[], MPI_Datatype recvtype, MPI_Comm comm,
                    MPI_Request *request)
                                                                                                                ;
int PMPI_Ialltoallv_c(const void *sendbuf, const MPI_Count sendcounts[], const MPI_Aint sdispls[],
                      MPI_Datatype sendtype, void *recvbuf, const MPI_Count recvcounts[],
                      const MPI_Aint rdispls[], MPI_Datatype recvtype, MPI_Comm comm,
                      MPI_Request *request)
                                                                                                                  ;
int PMPI_Ialltoallw(const void *sendbuf, const int sendcounts[], const int sdispls[],
                    const MPI_Datatype sendtypes[], void *recvbuf, const int recvcounts[],
                    const int rdispls[], const MPI_Datatype recvtypes[], MPI_Comm comm,
                    MPI_Request *request) ;
int PMPI_Ialltoallw_c(const void *sendbuf, const MPI_Count sendcounts[], const MPI_Aint sdispls[],
                      const MPI_Datatype sendtypes[], void *recvbuf, const MPI_Count recvcounts[],
                      const MPI_Aint rdispls[], const MPI_Datatype recvtypes[], MPI_Comm comm,
                      MPI_Request *request) ;
int PMPI_Ibarrier(MPI_Comm comm, MPI_Request *request) ;
int PMPI_Ibcast(void *buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm,
                MPI_Request *request) ;
int PMPI_Ibcast_c(void *buffer, MPI_Count count, MPI_Datatype datatype, int root, MPI_Comm comm,
                  MPI_Request *request) ;
int PMPI_Iexscan(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op,
                 MPI_Comm comm, MPI_Request *request)
                                                                                                             ;
int PMPI_Iexscan_c(const void *sendbuf, void *recvbuf, MPI_Count count, MPI_Datatype datatype,
                   MPI_Op op, MPI_Comm comm, MPI_Request *request)
                                                                                                               ;
int PMPI_Igather(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
                 int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm,
                 MPI_Request *request)
                                                                                                             ;
int PMPI_Igather_c(const void *sendbuf, MPI_Count sendcount, MPI_Datatype sendtype, void *recvbuf,
                   MPI_Count recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm,
                   MPI_Request *request)
                                                                                                               ;
int PMPI_Igatherv(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
                  const int recvcounts[], const int displs[], MPI_Datatype recvtype, int root,
                  MPI_Comm comm, MPI_Request *request)
                                                                                                              ;
int PMPI_Igatherv_c(const void *sendbuf, MPI_Count sendcount, MPI_Datatype sendtype, void *recvbuf,
                    const MPI_Count recvcounts[], const MPI_Aint displs[], MPI_Datatype recvtype,
                    int root, MPI_Comm comm, MPI_Request *request)
                                                                                                                ;
int PMPI_Ineighbor_allgather(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
                             void *recvbuf, int recvcount, MPI_Datatype recvtype, MPI_Comm comm,
                             MPI_Request *request)
                                                                                                                         ;
int PMPI_Ineighbor_allgather_c(const void *sendbuf, MPI_Count sendcount, MPI_Datatype sendtype,
                               void *recvbuf, MPI_Count recvcount, MPI_Datatype recvtype,
                               MPI_Comm comm, MPI_Request *request)
                                                                                                                           ;
int PMPI_Ineighbor_allgatherv(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
                              void *recvbuf, const int recvcounts[], const int displs[],
                              MPI_Datatype recvtype, MPI_Comm comm, MPI_Request *request)
                                                                                                                          ;
int PMPI_Ineighbor_allgatherv_c(const void *sendbuf, MPI_Count sendcount, MPI_Datatype sendtype,
                                void *recvbuf, const MPI_Count recvcounts[],
                                const MPI_Aint displs[], MPI_Datatype recvtype, MPI_Comm comm,
                                MPI_Request *request)
                                                                                                                            ;
int PMPI_Ineighbor_alltoall(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
                            void *recvbuf, int recvcount, MPI_Datatype recvtype, MPI_Comm comm,
                            MPI_Request *request)
                                                                                                                        ;
int PMPI_Ineighbor_alltoall_c(const void *sendbuf, MPI_Count sendcount, MPI_Datatype sendtype,
                              void *recvbuf, MPI_Count recvcount, MPI_Datatype recvtype,
                              MPI_Comm comm, MPI_Request *request)
                                                                                                                          ;
int PMPI_Ineighbor_alltoallv(const void *sendbuf, const int sendcounts[], const int sdispls[],
                             MPI_Datatype sendtype, void *recvbuf, const int recvcounts[],
                             const int rdispls[], MPI_Datatype recvtype, MPI_Comm comm,
                             MPI_Request *request)
                                                                                                                         ;
int PMPI_Ineighbor_alltoallv_c(const void *sendbuf, const MPI_Count sendcounts[],
                               const MPI_Aint sdispls[], MPI_Datatype sendtype, void *recvbuf,
                               const MPI_Count recvcounts[], const MPI_Aint rdispls[],
                               MPI_Datatype recvtype, MPI_Comm comm, MPI_Request *request)
                                                                                                                           ;
int PMPI_Ineighbor_alltoallw(const void *sendbuf, const int sendcounts[], const MPI_Aint sdispls[],
                             const MPI_Datatype sendtypes[], void *recvbuf, const int recvcounts[],
                             const MPI_Aint rdispls[], const MPI_Datatype recvtypes[],
                             MPI_Comm comm, MPI_Request *request) ;
int PMPI_Ineighbor_alltoallw_c(const void *sendbuf, const MPI_Count sendcounts[],
                               const MPI_Aint sdispls[], const MPI_Datatype sendtypes[],
                               void *recvbuf, const MPI_Count recvcounts[],
                               const MPI_Aint rdispls[], const MPI_Datatype recvtypes[],
                               MPI_Comm comm, MPI_Request *request) ;
int PMPI_Ireduce(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op,
                 int root, MPI_Comm comm, MPI_Request *request)
                                                                                                             ;
int PMPI_Ireduce_c(const void *sendbuf, void *recvbuf, MPI_Count count, MPI_Datatype datatype,
                   MPI_Op op, int root, MPI_Comm comm, MPI_Request *request)
                                                                                                               ;
int PMPI_Ireduce_scatter(const void *sendbuf, void *recvbuf, const int recvcounts[],
                         MPI_Datatype datatype, MPI_Op op, MPI_Comm comm, MPI_Request *request)
                                                                                                                     ;
int PMPI_Ireduce_scatter_c(const void *sendbuf, void *recvbuf, const MPI_Count recvcounts[],
                           MPI_Datatype datatype, MPI_Op op, MPI_Comm comm, MPI_Request *request)
                                                                                                                       ;
int PMPI_Ireduce_scatter_block(const void *sendbuf, void *recvbuf, int recvcount,
                               MPI_Datatype datatype, MPI_Op op, MPI_Comm comm,
                               MPI_Request *request)
                                                                                                                           ;
int PMPI_Ireduce_scatter_block_c(const void *sendbuf, void *recvbuf, MPI_Count recvcount,
                                 MPI_Datatype datatype, MPI_Op op, MPI_Comm comm,
                                 MPI_Request *request)
                                                                                                                             ;
int PMPI_Iscan(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op,
               MPI_Comm comm, MPI_Request *request)
                                                                                                           ;
int PMPI_Iscan_c(const void *sendbuf, void *recvbuf, MPI_Count count, MPI_Datatype datatype,
                 MPI_Op op, MPI_Comm comm, MPI_Request *request)
                                                                                                             ;
int PMPI_Iscatter(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
                  int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm,
                  MPI_Request *request)
                                                                                                              ;
int PMPI_Iscatter_c(const void *sendbuf, MPI_Count sendcount, MPI_Datatype sendtype, void *recvbuf,
                    MPI_Count recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm,
                    MPI_Request *request)
                                                                                                                ;
int PMPI_Iscatterv(const void *sendbuf, const int sendcounts[], const int displs[],
                   MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype,
                   int root, MPI_Comm comm, MPI_Request *request)
                                                                                                               ;
int PMPI_Iscatterv_c(const void *sendbuf, const MPI_Count sendcounts[], const MPI_Aint displs[],
                     MPI_Datatype sendtype, void *recvbuf, MPI_Count recvcount,
                     MPI_Datatype recvtype, int root, MPI_Comm comm, MPI_Request *request)
                                                                                                                 ;
int PMPI_Neighbor_allgather(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
                            void *recvbuf, int recvcount, MPI_Datatype recvtype, MPI_Comm comm)
                                                                                                                        ;
int PMPI_Neighbor_allgather_c(const void *sendbuf, MPI_Count sendcount, MPI_Datatype sendtype,
                              void *recvbuf, MPI_Count recvcount, MPI_Datatype recvtype,
                              MPI_Comm comm)
                                                                                                                          ;
int PMPI_Neighbor_allgather_init(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
                                 void *recvbuf, int recvcount, MPI_Datatype recvtype, MPI_Comm comm,
                                 MPI_Info info, MPI_Request *request)
                                                                                                                             ;
int PMPI_Neighbor_allgather_init_c(const void *sendbuf, MPI_Count sendcount, MPI_Datatype sendtype,
                                   void *recvbuf, MPI_Count recvcount, MPI_Datatype recvtype,
                                   MPI_Comm comm, MPI_Info info, MPI_Request *request)
                                                                                                                               ;
int PMPI_Neighbor_allgatherv(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
                             void *recvbuf, const int recvcounts[], const int displs[],
                             MPI_Datatype recvtype, MPI_Comm comm)
                                                                                                                         ;
int PMPI_Neighbor_allgatherv_c(const void *sendbuf, MPI_Count sendcount, MPI_Datatype sendtype,
                               void *recvbuf, const MPI_Count recvcounts[], const MPI_Aint displs[],
                               MPI_Datatype recvtype, MPI_Comm comm)
                                                                                                                           ;
int PMPI_Neighbor_allgatherv_init(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
                                  void *recvbuf, const int recvcounts[], const int displs[],
                                  MPI_Datatype recvtype, MPI_Comm comm, MPI_Info info,
                                  MPI_Request *request)
                                                                                                                              ;
int PMPI_Neighbor_allgatherv_init_c(const void *sendbuf, MPI_Count sendcount, MPI_Datatype sendtype,
                                    void *recvbuf, const MPI_Count recvcounts[],
                                    const MPI_Aint displs[], MPI_Datatype recvtype, MPI_Comm comm,
                                    MPI_Info info, MPI_Request *request)
                                                                                                                                ;
int PMPI_Neighbor_alltoall(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
                           int recvcount, MPI_Datatype recvtype, MPI_Comm comm)
                                                                                                                       ;
int PMPI_Neighbor_alltoall_c(const void *sendbuf, MPI_Count sendcount, MPI_Datatype sendtype,
                             void *recvbuf, MPI_Count recvcount, MPI_Datatype recvtype,
                             MPI_Comm comm)
                                                                                                                         ;
int PMPI_Neighbor_alltoall_init(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
                                void *recvbuf, int recvcount, MPI_Datatype recvtype, MPI_Comm comm,
                                MPI_Info info, MPI_Request *request)
                                                                                                                            ;
int PMPI_Neighbor_alltoall_init_c(const void *sendbuf, MPI_Count sendcount, MPI_Datatype sendtype,
                                  void *recvbuf, MPI_Count recvcount, MPI_Datatype recvtype,
                                  MPI_Comm comm, MPI_Info info, MPI_Request *request)
                                                                                                                              ;
int PMPI_Neighbor_alltoallv(const void *sendbuf, const int sendcounts[], const int sdispls[],
                            MPI_Datatype sendtype, void *recvbuf, const int recvcounts[],
                            const int rdispls[], MPI_Datatype recvtype, MPI_Comm comm)
                                                                                                                        ;
int PMPI_Neighbor_alltoallv_c(const void *sendbuf, const MPI_Count sendcounts[],
                              const MPI_Aint sdispls[], MPI_Datatype sendtype, void *recvbuf,
                              const MPI_Count recvcounts[], const MPI_Aint rdispls[],
                              MPI_Datatype recvtype, MPI_Comm comm)
                                                                                                                          ;
int PMPI_Neighbor_alltoallv_init(const void *sendbuf, const int sendcounts[], const int sdispls[],
                                 MPI_Datatype sendtype, void *recvbuf, const int recvcounts[],
                                 const int rdispls[], MPI_Datatype recvtype, MPI_Comm comm,
                                 MPI_Info info, MPI_Request *request)
                                                                                                                             ;
int PMPI_Neighbor_alltoallv_init_c(const void *sendbuf, const MPI_Count sendcounts[],
                                   const MPI_Aint sdispls[], MPI_Datatype sendtype, void *recvbuf,
                                   const MPI_Count recvcounts[], const MPI_Aint rdispls[],
                                   MPI_Datatype recvtype, MPI_Comm comm, MPI_Info info,
                                   MPI_Request *request)
                                                                                                                               ;
int PMPI_Neighbor_alltoallw(const void *sendbuf, const int sendcounts[], const MPI_Aint sdispls[],
                            const MPI_Datatype sendtypes[], void *recvbuf, const int recvcounts[],
                            const MPI_Aint rdispls[], const MPI_Datatype recvtypes[],
                            MPI_Comm comm) ;
int PMPI_Neighbor_alltoallw_c(const void *sendbuf, const MPI_Count sendcounts[],
                              const MPI_Aint sdispls[], const MPI_Datatype sendtypes[],
                              void *recvbuf, const MPI_Count recvcounts[], const MPI_Aint rdispls[],
                              const MPI_Datatype recvtypes[], MPI_Comm comm) ;
int PMPI_Neighbor_alltoallw_init(const void *sendbuf, const int sendcounts[],
                                 const MPI_Aint sdispls[], const MPI_Datatype sendtypes[],
                                 void *recvbuf, const int recvcounts[], const MPI_Aint rdispls[],
                                 const MPI_Datatype recvtypes[], MPI_Comm comm, MPI_Info info,
                                 MPI_Request *request) ;
int PMPI_Neighbor_alltoallw_init_c(const void *sendbuf, const MPI_Count sendcounts[],
                                   const MPI_Aint sdispls[], const MPI_Datatype sendtypes[],
                                   void *recvbuf, const MPI_Count recvcounts[],
                                   const MPI_Aint rdispls[], const MPI_Datatype recvtypes[],
                                   MPI_Comm comm, MPI_Info info, MPI_Request *request)
                                                   ;
int PMPI_Reduce(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op,
                int root, MPI_Comm comm)
                                                                                                            ;
int PMPI_Reduce_c(const void *sendbuf, void *recvbuf, MPI_Count count, MPI_Datatype datatype,
                  MPI_Op op, int root, MPI_Comm comm)
                                                                                                              ;
int PMPI_Reduce_init(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype,
                     MPI_Op op, int root, MPI_Comm comm, MPI_Info info, MPI_Request *request)
                                                                                                                 ;
int PMPI_Reduce_init_c(const void *sendbuf, void *recvbuf, MPI_Count count, MPI_Datatype datatype,
                       MPI_Op op, int root, MPI_Comm comm, MPI_Info info, MPI_Request *request)
                                                                                                                   ;
int PMPI_Reduce_local(const void *inbuf, void *inoutbuf, int count, MPI_Datatype datatype,
                      MPI_Op op)
                                                                                                                  ;
int PMPI_Reduce_local_c(const void *inbuf, void *inoutbuf, MPI_Count count, MPI_Datatype datatype,
                        MPI_Op op)
                                                                                                                    ;
int PMPI_Reduce_scatter(const void *sendbuf, void *recvbuf, const int recvcounts[],
                        MPI_Datatype datatype, MPI_Op op, MPI_Comm comm)
                                                                                                                    ;
int PMPI_Reduce_scatter_c(const void *sendbuf, void *recvbuf, const MPI_Count recvcounts[],
                          MPI_Datatype datatype, MPI_Op op, MPI_Comm comm)
                                                                                                                      ;
int PMPI_Reduce_scatter_block(const void *sendbuf, void *recvbuf, int recvcount,
                              MPI_Datatype datatype, MPI_Op op, MPI_Comm comm)
                                                                                                                          ;
int PMPI_Reduce_scatter_block_c(const void *sendbuf, void *recvbuf, MPI_Count recvcount,
                                MPI_Datatype datatype, MPI_Op op, MPI_Comm comm)
                                                                                                                            ;
int PMPI_Reduce_scatter_block_init(const void *sendbuf, void *recvbuf, int recvcount,
                                   MPI_Datatype datatype, MPI_Op op, MPI_Comm comm, MPI_Info info,
                                   MPI_Request *request)
                                                                                                                               ;
int PMPI_Reduce_scatter_block_init_c(const void *sendbuf, void *recvbuf, MPI_Count recvcount,
                                     MPI_Datatype datatype, MPI_Op op, MPI_Comm comm, MPI_Info info,
                                     MPI_Request *request)
                                                                                                                                 ;
int PMPI_Reduce_scatter_init(const void *sendbuf, void *recvbuf, const int recvcounts[],
                             MPI_Datatype datatype, MPI_Op op, MPI_Comm comm, MPI_Info info,
                             MPI_Request *request)
                                                                                                                         ;
int PMPI_Reduce_scatter_init_c(const void *sendbuf, void *recvbuf, const MPI_Count recvcounts[],
                               MPI_Datatype datatype, MPI_Op op, MPI_Comm comm, MPI_Info info,
                               MPI_Request *request)
                                                                                                                           ;
int PMPI_Scan(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op,
              MPI_Comm comm)
                                                                                                          ;
int PMPI_Scan_c(const void *sendbuf, void *recvbuf, MPI_Count count, MPI_Datatype datatype,
                MPI_Op op, MPI_Comm comm)
                                                                                                            ;
int PMPI_Scan_init(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op,
                   MPI_Comm comm, MPI_Info info, MPI_Request *request)
                                                                                                               ;
int PMPI_Scan_init_c(const void *sendbuf, void *recvbuf, MPI_Count count, MPI_Datatype datatype,
                     MPI_Op op, MPI_Comm comm, MPI_Info info, MPI_Request *request)
                                                                                                                 ;
int PMPI_Scatter(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
                 int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm)
                                                                                                             ;
int PMPI_Scatter_c(const void *sendbuf, MPI_Count sendcount, MPI_Datatype sendtype, void *recvbuf,
                   MPI_Count recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm)
                                                                                                               ;
int PMPI_Scatter_init(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
                      int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm, MPI_Info info,
                      MPI_Request *request)
                                                                                                                  ;
int PMPI_Scatter_init_c(const void *sendbuf, MPI_Count sendcount, MPI_Datatype sendtype,
                        void *recvbuf, MPI_Count recvcount, MPI_Datatype recvtype, int root,
                        MPI_Comm comm, MPI_Info info, MPI_Request *request)
                                                                                                                    ;
int PMPI_Scatterv(const void *sendbuf, const int sendcounts[], const int displs[],
                  MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype,
                  int root, MPI_Comm comm)
                                                                                                              ;
int PMPI_Scatterv_c(const void *sendbuf, const MPI_Count sendcounts[], const MPI_Aint displs[],
                    MPI_Datatype sendtype, void *recvbuf, MPI_Count recvcount,
                    MPI_Datatype recvtype, int root, MPI_Comm comm)
                                                                                                                ;
int PMPI_Scatterv_init(const void *sendbuf, const int sendcounts[], const int displs[],
                       MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype,
                       int root, MPI_Comm comm, MPI_Info info, MPI_Request *request)
                                                                                                                   ;
int PMPI_Scatterv_init_c(const void *sendbuf, const MPI_Count sendcounts[], const MPI_Aint displs[],
                         MPI_Datatype sendtype, void *recvbuf, MPI_Count recvcount,
                         MPI_Datatype recvtype, int root, MPI_Comm comm, MPI_Info info,
                         MPI_Request *request)
                                                                                                                     ;
int PMPI_Comm_compare(MPI_Comm comm1, MPI_Comm comm2, int *result) ;
int PMPI_Comm_create(MPI_Comm comm, MPI_Group group, MPI_Comm *newcomm) ;
int PMPI_Comm_create_group(MPI_Comm comm, MPI_Group group, int tag, MPI_Comm *newcomm)
                    ;
int PMPI_Comm_dup(MPI_Comm comm, MPI_Comm *newcomm) ;
int PMPI_Comm_dup_with_info(MPI_Comm comm, MPI_Info info, MPI_Comm *newcomm) ;
int PMPI_Comm_free(MPI_Comm *comm) ;
int PMPI_Comm_get_info(MPI_Comm comm, MPI_Info *info_used) ;
int PMPI_Comm_get_name(MPI_Comm comm, char *comm_name, int *resultlen) ;
int PMPI_Comm_group(MPI_Comm comm, MPI_Group *group) ;
int PMPI_Comm_idup(MPI_Comm comm, MPI_Comm *newcomm, MPI_Request *request) ;
int PMPI_Comm_idup_with_info(MPI_Comm comm, MPI_Info info, MPI_Comm *newcomm, MPI_Request *request)
                    ;
int PMPI_Comm_rank(MPI_Comm comm, int *rank) ;
int PMPI_Comm_remote_group(MPI_Comm comm, MPI_Group *group) ;
int PMPI_Comm_remote_size(MPI_Comm comm, int *size) ;
int PMPI_Comm_set_info(MPI_Comm comm, MPI_Info info) ;
int PMPI_Comm_set_name(MPI_Comm comm, const char *comm_name) ;
int PMPI_Comm_size(MPI_Comm comm, int *size) ;
int PMPI_Comm_split(MPI_Comm comm, int color, int key, MPI_Comm *newcomm) ;
int PMPI_Comm_split_type(MPI_Comm comm, int split_type, int key, MPI_Info info, MPI_Comm *newcomm)
                    ;
int PMPI_Comm_test_inter(MPI_Comm comm, int *flag) ;
int PMPI_Intercomm_create(MPI_Comm local_comm, int local_leader, MPI_Comm peer_comm,
                          int remote_leader, int tag, MPI_Comm *newintercomm) ;
int PMPI_Intercomm_create_from_groups(MPI_Group local_group, int local_leader,
                                      MPI_Group remote_group, int remote_leader,
                                      const char *stringtag, MPI_Info info,
                                      MPI_Errhandler errhandler, MPI_Comm *newintercomm)
                                                      ;
int PMPI_Intercomm_merge(MPI_Comm intercomm, int high, MPI_Comm *newintracomm) ;
int PMPIX_Comm_test_threadcomm(MPI_Comm comm, int *flag) ;
int PMPIX_Comm_revoke(MPI_Comm comm) ;
int PMPIX_Comm_shrink(MPI_Comm comm, MPI_Comm *newcomm) ;
int PMPIX_Comm_failure_ack(MPI_Comm comm) ;
int PMPIX_Comm_failure_get_acked(MPI_Comm comm, MPI_Group *failedgrp) ;
int PMPIX_Comm_agree(MPI_Comm comm, int *flag) ;
int PMPIX_Comm_get_failed(MPI_Comm comm, MPI_Group *failedgrp) ;
int PMPI_Get_address(const void *location, MPI_Aint *address) ;
int PMPI_Address(void *location, MPI_Aint *address) ;
int PMPI_Get_count(const MPI_Status *status, MPI_Datatype datatype, int *count) ;
int PMPI_Get_count_c(const MPI_Status *status, MPI_Datatype datatype, MPI_Count *count)
                    ;
int PMPI_Get_elements(const MPI_Status *status, MPI_Datatype datatype, int *count)
                    ;
int PMPI_Get_elements_c(const MPI_Status *status, MPI_Datatype datatype, MPI_Count *count)
                    ;
int PMPI_Get_elements_x(const MPI_Status *status, MPI_Datatype datatype, MPI_Count *count)
                    ;
int PMPI_Pack(const void *inbuf, int incount, MPI_Datatype datatype, void *outbuf, int outsize,
              int *position, MPI_Comm comm) ;
int PMPI_Pack_c(const void *inbuf, MPI_Count incount, MPI_Datatype datatype, void *outbuf,
                MPI_Count outsize, MPI_Count *position, MPI_Comm comm) ;
int PMPI_Pack_external(const char *datarep, const void *inbuf, int incount, MPI_Datatype datatype,
                       void *outbuf, MPI_Aint outsize, MPI_Aint *position) ;
int PMPI_Pack_external_c(const char *datarep, const void *inbuf, MPI_Count incount,
                         MPI_Datatype datatype, void *outbuf, MPI_Count outsize,
                         MPI_Count *position) ;
int PMPI_Pack_external_size(const char *datarep, int incount, MPI_Datatype datatype,
                            MPI_Aint *size) ;
int PMPI_Pack_external_size_c(const char *datarep, MPI_Count incount, MPI_Datatype datatype,
                              MPI_Count *size) ;
int PMPI_Pack_size(int incount, MPI_Datatype datatype, MPI_Comm comm, int *size) ;
int PMPI_Pack_size_c(MPI_Count incount, MPI_Datatype datatype, MPI_Comm comm, MPI_Count *size)
                    ;
int PMPI_Status_set_elements(MPI_Status *status, MPI_Datatype datatype, int count)
                    ;
int PMPI_Status_set_elements_c(MPI_Status *status, MPI_Datatype datatype, MPI_Count count)
                    ;
int PMPI_Status_set_elements_x(MPI_Status *status, MPI_Datatype datatype, MPI_Count count)
                    ;
int PMPI_Type_commit(MPI_Datatype *datatype) ;
int PMPI_Type_contiguous(int count, MPI_Datatype oldtype, MPI_Datatype *newtype) ;
int PMPI_Type_contiguous_c(MPI_Count count, MPI_Datatype oldtype, MPI_Datatype *newtype)
                    ;
int PMPI_Type_create_darray(int size, int rank, int ndims, const int array_of_gsizes[],
                            const int array_of_distribs[], const int array_of_dargs[],
                            const int array_of_psizes[], int order, MPI_Datatype oldtype,
                            MPI_Datatype *newtype) ;
int PMPI_Type_create_darray_c(int size, int rank, int ndims, const MPI_Count array_of_gsizes[],
                              const int array_of_distribs[], const int array_of_dargs[],
                              const int array_of_psizes[], int order, MPI_Datatype oldtype,
                              MPI_Datatype *newtype) ;
int PMPI_Type_create_f90_complex(int p, int r, MPI_Datatype *newtype) ;
int PMPI_Type_create_f90_integer(int r, MPI_Datatype *newtype) ;
int PMPI_Type_create_f90_real(int p, int r, MPI_Datatype *newtype) ;
int PMPI_Type_create_hindexed(int count, const int array_of_blocklengths[],
                              const MPI_Aint array_of_displacements[], MPI_Datatype oldtype,
                              MPI_Datatype *newtype) ;
int PMPI_Type_create_hindexed_c(MPI_Count count, const MPI_Count array_of_blocklengths[],
                                const MPI_Count array_of_displacements[], MPI_Datatype oldtype,
                                MPI_Datatype *newtype) ;
int PMPI_Type_hindexed(int count, int array_of_blocklengths[], MPI_Aint array_of_displacements[],
                       MPI_Datatype oldtype, MPI_Datatype *newtype) ;
int PMPI_Type_create_hindexed_block(int count, int blocklength,
                                    const MPI_Aint array_of_displacements[], MPI_Datatype oldtype,
                                    MPI_Datatype *newtype) ;
int PMPI_Type_create_hindexed_block_c(MPI_Count count, MPI_Count blocklength,
                                      const MPI_Count array_of_displacements[],
                                      MPI_Datatype oldtype, MPI_Datatype *newtype)
                                                      ;
int PMPI_Type_create_hvector(int count, int blocklength, MPI_Aint stride, MPI_Datatype oldtype,
                             MPI_Datatype *newtype) ;
int PMPI_Type_create_hvector_c(MPI_Count count, MPI_Count blocklength, MPI_Count stride,
                               MPI_Datatype oldtype, MPI_Datatype *newtype) ;
int PMPI_Type_hvector(int count, int blocklength, MPI_Aint stride, MPI_Datatype oldtype,
                      MPI_Datatype *newtype) ;
int PMPI_Type_create_indexed_block(int count, int blocklength, const int array_of_displacements[],
                                   MPI_Datatype oldtype, MPI_Datatype *newtype) ;
int PMPI_Type_create_indexed_block_c(MPI_Count count, MPI_Count blocklength,
                                     const MPI_Count array_of_displacements[], MPI_Datatype oldtype,
                                     MPI_Datatype *newtype) ;
int PMPI_Type_create_resized(MPI_Datatype oldtype, MPI_Aint lb, MPI_Aint extent,
                             MPI_Datatype *newtype) ;
int PMPI_Type_create_resized_c(MPI_Datatype oldtype, MPI_Count lb, MPI_Count extent,
                               MPI_Datatype *newtype) ;
int PMPI_Type_create_struct(int count, const int array_of_blocklengths[],
                            const MPI_Aint array_of_displacements[],
                            const MPI_Datatype array_of_types[], MPI_Datatype *newtype)
                                            ;
int PMPI_Type_create_struct_c(MPI_Count count, const MPI_Count array_of_blocklengths[],
                              const MPI_Count array_of_displacements[],
                              const MPI_Datatype array_of_types[], MPI_Datatype *newtype)
                                              ;
int PMPI_Type_struct(int count, int array_of_blocklengths[], MPI_Aint array_of_displacements[],
                     MPI_Datatype array_of_types[], MPI_Datatype *newtype) ;
int PMPI_Type_create_subarray(int ndims, const int array_of_sizes[], const int array_of_subsizes[],
                              const int array_of_starts[], int order, MPI_Datatype oldtype,
                              MPI_Datatype *newtype) ;
int PMPI_Type_create_subarray_c(int ndims, const MPI_Count array_of_sizes[],
                                const MPI_Count array_of_subsizes[],
                                const MPI_Count array_of_starts[], int order, MPI_Datatype oldtype,
                                MPI_Datatype *newtype) ;
int PMPI_Type_dup(MPI_Datatype oldtype, MPI_Datatype *newtype) ;
int PMPI_Type_free(MPI_Datatype *datatype) ;
int PMPI_Type_get_contents(MPI_Datatype datatype, int max_integers, int max_addresses,
                           int max_datatypes, int array_of_integers[],
                           MPI_Aint array_of_addresses[], MPI_Datatype array_of_datatypes[])
                                           ;
int PMPI_Type_get_contents_c(MPI_Datatype datatype, MPI_Count max_integers, MPI_Count max_addresses,
                             MPI_Count max_large_counts, MPI_Count max_datatypes,
                             int array_of_integers[], MPI_Aint array_of_addresses[],
                             MPI_Count array_of_large_counts[], MPI_Datatype array_of_datatypes[])
                                             ;
int PMPI_Type_get_envelope(MPI_Datatype datatype, int *num_integers, int *num_addresses,
                           int *num_datatypes, int *combiner) ;
int PMPI_Type_get_envelope_c(MPI_Datatype datatype, MPI_Count *num_integers,
                             MPI_Count *num_addresses, MPI_Count *num_large_counts,
                             MPI_Count *num_datatypes, int *combiner) ;
int PMPI_Type_get_extent(MPI_Datatype datatype, MPI_Aint *lb, MPI_Aint *extent) ;
int PMPI_Type_get_extent_c(MPI_Datatype datatype, MPI_Count *lb, MPI_Count *extent)
                    ;
int PMPI_Type_get_extent_x(MPI_Datatype datatype, MPI_Count *lb, MPI_Count *extent)
                    ;
int PMPI_Type_get_name(MPI_Datatype datatype, char *type_name, int *resultlen) ;
int PMPI_Type_get_true_extent(MPI_Datatype datatype, MPI_Aint *true_lb, MPI_Aint *true_extent)
                    ;
int PMPI_Type_get_true_extent_c(MPI_Datatype datatype, MPI_Count *true_lb, MPI_Count *true_extent)
                    ;
int PMPI_Type_get_true_extent_x(MPI_Datatype datatype, MPI_Count *true_lb, MPI_Count *true_extent)
                    ;
int PMPI_Type_get_value_index(MPI_Datatype value_type, MPI_Datatype index_type,
                              MPI_Datatype *pair_type) ;
int PMPI_Type_indexed(int count, const int array_of_blocklengths[],
                      const int array_of_displacements[], MPI_Datatype oldtype,
                      MPI_Datatype *newtype) ;
int PMPI_Type_indexed_c(MPI_Count count, const MPI_Count array_of_blocklengths[],
                        const MPI_Count array_of_displacements[], MPI_Datatype oldtype,
                        MPI_Datatype *newtype) ;
int PMPI_Type_match_size(int typeclass, int size, MPI_Datatype *datatype) ;
int PMPI_Type_set_name(MPI_Datatype datatype, const char *type_name) ;
int PMPI_Type_size(MPI_Datatype datatype, int *size) ;
int PMPI_Type_size_c(MPI_Datatype datatype, MPI_Count *size) ;
int PMPI_Type_size_x(MPI_Datatype datatype, MPI_Count *size) ;
int PMPI_Type_vector(int count, int blocklength, int stride, MPI_Datatype oldtype,
                     MPI_Datatype *newtype) ;
int PMPI_Type_vector_c(MPI_Count count, MPI_Count blocklength, MPI_Count stride,
                       MPI_Datatype oldtype, MPI_Datatype *newtype) ;
int PMPI_Unpack(const void *inbuf, int insize, int *position, void *outbuf, int outcount,
                MPI_Datatype datatype, MPI_Comm comm) ;
int PMPI_Unpack_c(const void *inbuf, MPI_Count insize, MPI_Count *position, void *outbuf,
                  MPI_Count outcount, MPI_Datatype datatype, MPI_Comm comm) ;
int PMPI_Unpack_external(const char datarep[], const void *inbuf, MPI_Aint insize,
                         MPI_Aint *position, void *outbuf, int outcount, MPI_Datatype datatype)
                                         ;
int PMPI_Unpack_external_c(const char datarep[], const void *inbuf, MPI_Count insize,
                           MPI_Count *position, void *outbuf, MPI_Count outcount,
                           MPI_Datatype datatype) ;
int PMPI_Type_extent(MPI_Datatype datatype, MPI_Aint *extent) ;
int PMPI_Type_lb(MPI_Datatype datatype, MPI_Aint *displacement) ;
int PMPI_Type_ub(MPI_Datatype datatype, MPI_Aint *displacement) ;
int PMPIX_Type_iov_len(MPI_Datatype datatype, MPI_Count max_iov_bytes, MPI_Count *iov_len,
                       MPI_Count *actual_iov_bytes) ;
int PMPIX_Type_iov(MPI_Datatype datatype, MPI_Count iov_offset, MPIX_Iov *iov,
                   MPI_Count max_iov_len, MPI_Count *actual_iov_len) ;
int PMPI_Add_error_class(int *errorclass) ;
int PMPI_Add_error_code(int errorclass, int *errorcode) ;
int PMPI_Add_error_string(int errorcode, const char *string) ;
int PMPI_Comm_call_errhandler(MPI_Comm comm, int errorcode) ;
int PMPI_Comm_create_errhandler(MPI_Comm_errhandler_function *comm_errhandler_fn,
                                MPI_Errhandler *errhandler) ;
int PMPI_Errhandler_create(MPI_Comm_errhandler_function *comm_errhandler_fn,
                           MPI_Errhandler *errhandler) ;
int PMPI_Comm_get_errhandler(MPI_Comm comm, MPI_Errhandler *errhandler) ;
int PMPI_Errhandler_get(MPI_Comm comm, MPI_Errhandler *errhandler) ;
int PMPI_Comm_set_errhandler(MPI_Comm comm, MPI_Errhandler errhandler) ;
int PMPI_Errhandler_set(MPI_Comm comm, MPI_Errhandler errhandler) ;
int PMPI_Errhandler_free(MPI_Errhandler *errhandler) ;
int PMPI_Error_class(int errorcode, int *errorclass) ;
int PMPI_Error_string(int errorcode, char *string, int *resultlen) ;
int PMPI_File_call_errhandler(MPI_File fh, int errorcode) ;
int PMPI_File_create_errhandler(MPI_File_errhandler_function *file_errhandler_fn,
                                MPI_Errhandler *errhandler) ;
int PMPI_File_get_errhandler(MPI_File file, MPI_Errhandler *errhandler) ;
int PMPI_File_set_errhandler(MPI_File file, MPI_Errhandler errhandler) ;
int PMPI_Remove_error_class(int errorclass) ;
int PMPI_Remove_error_code(int errorcode) ;
int PMPI_Remove_error_string(int errorcode) ;
int PMPI_Session_call_errhandler(MPI_Session session, int errorcode) ;
int PMPI_Session_create_errhandler(MPI_Session_errhandler_function *session_errhandler_fn,
                                   MPI_Errhandler *errhandler) ;
int PMPI_Session_get_errhandler(MPI_Session session, MPI_Errhandler *errhandler) ;
int PMPI_Session_set_errhandler(MPI_Session session, MPI_Errhandler errhandler) ;
int PMPI_Win_call_errhandler(MPI_Win win, int errorcode) ;
int PMPI_Win_create_errhandler(MPI_Win_errhandler_function *win_errhandler_fn,
                               MPI_Errhandler *errhandler) ;
int PMPI_Win_get_errhandler(MPI_Win win, MPI_Errhandler *errhandler) ;
int PMPI_Win_set_errhandler(MPI_Win win, MPI_Errhandler errhandler) ;
MPI_Fint PMPI_Comm_c2f(MPI_Comm comm) ;
MPI_Comm PMPI_Comm_f2c(MPI_Fint comm) ;
MPI_Fint PMPI_Errhandler_c2f(MPI_Errhandler errhandler) ;
MPI_Errhandler PMPI_Errhandler_f2c(MPI_Fint errhandler) ;
MPI_Fint PMPI_Group_c2f(MPI_Group group) ;
MPI_Group PMPI_Group_f2c(MPI_Fint group) ;
MPI_Fint PMPI_Info_c2f(MPI_Info info) ;
MPI_Info PMPI_Info_f2c(MPI_Fint info) ;
MPI_Fint PMPI_Message_c2f(MPI_Message message) ;
MPI_Message PMPI_Message_f2c(MPI_Fint message) ;
MPI_Fint PMPI_Op_c2f(MPI_Op op) ;
MPI_Op PMPI_Op_f2c(MPI_Fint op) ;
MPI_Fint PMPI_Request_c2f(MPI_Request request) ;
MPI_Request PMPI_Request_f2c(MPI_Fint request) ;
MPI_Fint PMPI_Session_c2f(MPI_Session session) ;
MPI_Session PMPI_Session_f2c(MPI_Fint session) ;
int PMPI_Status_c2f(const MPI_Status *c_status, MPI_Fint *f_status) ;
int PMPI_Status_c2f08(const MPI_Status *c_status, MPI_F08_status *f08_status) ;
int PMPI_Status_f082c(const MPI_F08_status *f08_status, MPI_Status *c_status) ;
int PMPI_Status_f082f(const MPI_F08_status *f08_status, MPI_Fint *f_status) ;
int PMPI_Status_f2c(const MPI_Fint *f_status, MPI_Status *c_status) ;
int PMPI_Status_f2f08(const MPI_Fint *f_status, MPI_F08_status *f08_status) ;
MPI_Fint PMPI_Type_c2f(MPI_Datatype datatype) ;
MPI_Datatype PMPI_Type_f2c(MPI_Fint datatype) ;
MPI_Fint PMPI_Win_c2f(MPI_Win win) ;
MPI_Win PMPI_Win_f2c(MPI_Fint win) ;
int PMPI_Group_compare(MPI_Group group1, MPI_Group group2, int *result) ;
int PMPI_Group_difference(MPI_Group group1, MPI_Group group2, MPI_Group *newgroup)
                    ;
int PMPI_Group_excl(MPI_Group group, int n, const int ranks[], MPI_Group *newgroup)
                    ;
int PMPI_Group_free(MPI_Group *group) ;
int PMPI_Group_incl(MPI_Group group, int n, const int ranks[], MPI_Group *newgroup)
                    ;
int PMPI_Group_intersection(MPI_Group group1, MPI_Group group2, MPI_Group *newgroup)
                    ;
int PMPI_Group_range_excl(MPI_Group group, int n, int ranges[][3], MPI_Group *newgroup)
                    ;
int PMPI_Group_range_incl(MPI_Group group, int n, int ranges[][3], MPI_Group *newgroup)
                    ;
int PMPI_Group_rank(MPI_Group group, int *rank) ;
int PMPI_Group_size(MPI_Group group, int *size) ;
int PMPI_Group_translate_ranks(MPI_Group group1, int n, const int ranks1[], MPI_Group group2,
                               int ranks2[]) ;
int PMPI_Group_union(MPI_Group group1, MPI_Group group2, MPI_Group *newgroup) ;
int PMPI_Info_create(MPI_Info *info) ;
int PMPI_Info_create_env(int argc, char *argv[], MPI_Info *info) ;
int PMPI_Info_delete(MPI_Info info, const char *key) ;
int PMPI_Info_dup(MPI_Info info, MPI_Info *newinfo) ;
int PMPI_Info_free(MPI_Info *info) ;
int PMPI_Info_get(MPI_Info info, const char *key, int valuelen, char *value, int *flag)
                    ;
int PMPI_Info_get_nkeys(MPI_Info info, int *nkeys) ;
int PMPI_Info_get_nthkey(MPI_Info info, int n, char *key) ;
int PMPI_Info_get_string(MPI_Info info, const char *key, int *buflen, char *value, int *flag)
                    ;
int PMPI_Info_get_valuelen(MPI_Info info, const char *key, int *valuelen, int *flag)
                    ;
int PMPI_Info_set(MPI_Info info, const char *key, const char *value) ;
int PMPIX_Info_set_hex(MPI_Info info, const char *key, const void *value, int value_size)
                    ;
int PMPI_Abort(MPI_Comm comm, int errorcode) ;
int PMPI_Comm_create_from_group(MPI_Group group, const char *stringtag, MPI_Info info,
                                MPI_Errhandler errhandler, MPI_Comm *newcomm) ;
int PMPI_Finalize(void) ;
int PMPI_Finalized(int *flag) ;
int PMPI_Group_from_session_pset(MPI_Session session, const char *pset_name, MPI_Group *newgroup)
                    ;
int PMPI_Init(int *argc, char ***argv) ;
int PMPI_Init_thread(int *argc, char ***argv, int required, int *provided) ;
int PMPI_Initialized(int *flag) ;
int PMPI_Is_thread_main(int *flag) ;
int PMPI_Query_thread(int *provided) ;
int PMPI_Session_finalize(MPI_Session *session) ;
int PMPI_Session_get_info(MPI_Session session, MPI_Info *info_used) ;
int PMPI_Session_get_nth_pset(MPI_Session session, MPI_Info info, int n, int *pset_len,
                              char *pset_name) ;
int PMPI_Session_get_num_psets(MPI_Session session, MPI_Info info, int *npset_names)
                    ;
int PMPI_Session_get_pset_info(MPI_Session session, const char *pset_name, MPI_Info *info)
                    ;
int PMPI_Session_init(MPI_Info info, MPI_Errhandler errhandler, MPI_Session *session)
                    ;
MPI_Aint PMPI_Aint_add(MPI_Aint base, MPI_Aint disp) ;
MPI_Aint PMPI_Aint_diff(MPI_Aint addr1, MPI_Aint addr2) ;
int PMPI_Get_library_version(char *version, int *resultlen) ;
int PMPI_Get_processor_name(char *name, int *resultlen) ;
int PMPI_Get_version(int *version, int *subversion) ;
int PMPI_Pcontrol(const int level, ...) ;
int PMPIX_GPU_query_support(int gpu_type, int *is_supported) ;
int PMPIX_Query_cuda_support(void) ;
int PMPIX_Query_ze_support(void) ;
int PMPIX_Query_hip_support(void) ;
int PMPI_T_category_changed(int *update_number) ;
int PMPI_T_category_get_categories(int cat_index, int len, int indices[]) ;
int PMPI_T_category_get_cvars(int cat_index, int len, int indices[]) ;
int PMPI_T_category_get_events(int cat_index, int len, int indices[]) ;
int PMPI_T_category_get_index(const char *name, int *cat_index) ;
int PMPI_T_category_get_info(int cat_index, char *name, int *name_len, char *desc, int *desc_len,
                             int *num_cvars, int *num_pvars, int *num_categories) ;
int PMPI_T_category_get_num(int *num_cat) ;
int PMPI_T_category_get_num_events(int cat_index, int *num_events) ;
int PMPI_T_category_get_pvars(int cat_index, int len, int indices[]) ;
int PMPI_T_cvar_get_index(const char *name, int *cvar_index) ;
int PMPI_T_cvar_get_info(int cvar_index, char *name, int *name_len, int *verbosity,
                         MPI_Datatype *datatype, MPI_T_enum *enumtype, char *desc, int *desc_len,
                         int *bind, int *scope) ;
int PMPI_T_cvar_get_num(int *num_cvar) ;
int PMPI_T_cvar_handle_alloc(int cvar_index, void *obj_handle, MPI_T_cvar_handle *handle,
                             int *count) ;
int PMPI_T_cvar_handle_free(MPI_T_cvar_handle *handle) ;
int PMPI_T_cvar_read(MPI_T_cvar_handle handle, void *buf) ;
int PMPI_T_cvar_write(MPI_T_cvar_handle handle, const void *buf) ;
int PMPI_T_enum_get_info(MPI_T_enum enumtype, int *num, char *name, int *name_len)
                    ;
int PMPI_T_enum_get_item(MPI_T_enum enumtype, int indx, int *value, char *name, int *name_len)
                    ;
int PMPI_T_event_callback_get_info(MPI_T_event_registration event_registration,
                                   MPI_T_cb_safety cb_safety, MPI_Info *info_used)
                                                   ;
int PMPI_T_event_callback_set_info(MPI_T_event_registration event_registration,
                                   MPI_T_cb_safety cb_safety, MPI_Info info) ;
int PMPI_T_event_copy(MPI_T_event_instance event_instance, void *buffer) ;
int PMPI_T_event_get_index(const char *name, int *event_index) ;
int PMPI_T_event_get_info(int event_index, char *name, int *name_len, int *verbosity,
                          MPI_Datatype array_of_datatypes[], MPI_Aint array_of_displacements[],
                          int *num_elements, MPI_T_enum *enumtype, MPI_Info *info, char *desc,
                          int *desc_len, int *bind) ;
int PMPI_T_event_get_num(int *num_events) ;
int PMPI_T_event_get_source(MPI_T_event_instance event_instance, int *source_index)
                    ;
int PMPI_T_event_get_timestamp(MPI_T_event_instance event_instance, MPI_Count *event_timestamp)
                    ;
int PMPI_T_event_handle_alloc(int event_index, void *obj_handle, MPI_Info info,
                              MPI_T_event_registration *event_registration) ;
int PMPI_T_event_handle_free(MPI_T_event_registration event_registration, void *user_data,
                             MPI_T_event_free_cb_function free_cb_function) ;
int PMPI_T_event_handle_get_info(MPI_T_event_registration event_registration, MPI_Info *info_used)
                    ;
int PMPI_T_event_handle_set_info(MPI_T_event_registration event_registration, MPI_Info info)
                    ;
int PMPI_T_event_read(MPI_T_event_instance event_instance, int element_index, void *buffer)
                    ;
int PMPI_T_event_register_callback(MPI_T_event_registration event_registration,
                                   MPI_T_cb_safety cb_safety, MPI_Info info, void *user_data,
                                   MPI_T_event_cb_function event_cb_function) ;
int PMPI_T_event_set_dropped_handler(MPI_T_event_registration event_registration,
                                     MPI_T_event_dropped_cb_function dropped_cb_function)
                                                     ;
int PMPI_T_finalize(void) ;
int PMPI_T_init_thread(int required, int *provided) ;
int PMPI_T_pvar_get_index(const char *name, int var_class, int *pvar_index) ;
int PMPI_T_pvar_get_info(int pvar_index, char *name, int *name_len, int *verbosity, int *var_class,
                         MPI_Datatype *datatype, MPI_T_enum *enumtype, char *desc, int *desc_len,
                         int *bind, int *readonly, int *continuous, int *atomic) ;
int PMPI_T_pvar_get_num(int *num_pvar) ;
int PMPI_T_pvar_handle_alloc(MPI_T_pvar_session session, int pvar_index, void *obj_handle,
                             MPI_T_pvar_handle *handle, int *count) ;
int PMPI_T_pvar_handle_free(MPI_T_pvar_session session, MPI_T_pvar_handle *handle)
                    ;
int PMPI_T_pvar_read(MPI_T_pvar_session session, MPI_T_pvar_handle handle, void *buf)
                    ;
int PMPI_T_pvar_readreset(MPI_T_pvar_session session, MPI_T_pvar_handle handle, void *buf)
                    ;
int PMPI_T_pvar_reset(MPI_T_pvar_session session, MPI_T_pvar_handle handle) ;
int PMPI_T_pvar_session_create(MPI_T_pvar_session *session) ;
int PMPI_T_pvar_session_free(MPI_T_pvar_session *session) ;
int PMPI_T_pvar_start(MPI_T_pvar_session session, MPI_T_pvar_handle handle) ;
int PMPI_T_pvar_stop(MPI_T_pvar_session session, MPI_T_pvar_handle handle) ;
int PMPI_T_pvar_write(MPI_T_pvar_session session, MPI_T_pvar_handle handle, const void *buf)
                    ;
int PMPI_T_source_get_info(int source_index, char *name, int *name_len, char *desc, int *desc_len,
                           MPI_T_source_order *ordering, MPI_Count *ticks_per_second,
                           MPI_Count *max_ticks, MPI_Info *info) ;
int PMPI_T_source_get_num(int *num_sources) ;
int PMPI_T_source_get_timestamp(int source_index, MPI_Count *timestamp) ;
int PMPI_Op_commutative(MPI_Op op, int *commute) ;
int PMPI_Op_create(MPI_User_function *user_fn, int commute, MPI_Op *op) ;
int PMPI_Op_create_c(MPI_User_function_c *user_fn, int commute, MPI_Op *op) ;
int PMPI_Op_free(MPI_Op *op) ;
int PMPI_Parrived(MPI_Request request, int partition, int *flag) ;
int PMPI_Pready(int partition, MPI_Request request) ;
int PMPI_Pready_list(int length, const int array_of_partitions[], MPI_Request request)
                    ;
int PMPI_Pready_range(int partition_low, int partition_high, MPI_Request request) ;
int PMPI_Precv_init(void *buf, int partitions, MPI_Count count, MPI_Datatype datatype, int dest,
                    int tag, MPI_Comm comm, MPI_Info info, MPI_Request *request)
                                                                          ;
int PMPI_Psend_init(const void *buf, int partitions, MPI_Count count, MPI_Datatype datatype,
                    int dest, int tag, MPI_Comm comm, MPI_Info info, MPI_Request *request)
                                                                          ;
int PMPI_Bsend(const void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm)
                                                          ;
int PMPI_Bsend_c(const void *buf, MPI_Count count, MPI_Datatype datatype, int dest, int tag,
                 MPI_Comm comm) ;
int PMPI_Bsend_init(const void *buf, int count, MPI_Datatype datatype, int dest, int tag,
                    MPI_Comm comm, MPI_Request *request)
                                                                          ;
int PMPI_Bsend_init_c(const void *buf, MPI_Count count, MPI_Datatype datatype, int dest, int tag,
                      MPI_Comm comm, MPI_Request *request)
                                                                            ;
int PMPI_Buffer_attach(void *buffer, int size) ;
int PMPI_Buffer_attach_c(void *buffer, MPI_Count size) ;
int PMPI_Buffer_detach(void *buffer_addr, int *size) ;
int PMPI_Buffer_detach_c(void *buffer_addr, MPI_Count *size) ;
int PMPI_Buffer_flush(void) ;
int PMPI_Buffer_iflush(MPI_Request *request) ;
int PMPI_Comm_attach_buffer(MPI_Comm comm, void *buffer, int size) ;
int PMPI_Comm_attach_buffer_c(MPI_Comm comm, void *buffer, MPI_Count size) ;
int PMPI_Comm_detach_buffer(MPI_Comm comm, void *buffer_addr, int *size) ;
int PMPI_Comm_detach_buffer_c(MPI_Comm comm, void *buffer_addr, MPI_Count *size) ;
int PMPI_Comm_flush_buffer(MPI_Comm comm) ;
int PMPI_Comm_iflush_buffer(MPI_Comm comm, MPI_Request *request) ;
int PMPI_Ibsend(const void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm,
                MPI_Request *request) ;
int PMPI_Ibsend_c(const void *buf, MPI_Count count, MPI_Datatype datatype, int dest, int tag,
                  MPI_Comm comm, MPI_Request *request)
                                                                        ;
int PMPI_Improbe(int source, int tag, MPI_Comm comm, int *flag, MPI_Message *message,
                 MPI_Status *status) ;
int PMPI_Imrecv(void *buf, int count, MPI_Datatype datatype, MPI_Message *message,
                MPI_Request *request) ;
int PMPI_Imrecv_c(void *buf, MPI_Count count, MPI_Datatype datatype, MPI_Message *message,
                  MPI_Request *request) ;
int PMPI_Iprobe(int source, int tag, MPI_Comm comm, int *flag, MPI_Status *status)
                    ;
int PMPI_Irecv(void *buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm,
               MPI_Request *request) ;
int PMPI_Irecv_c(void *buf, MPI_Count count, MPI_Datatype datatype, int source, int tag,
                 MPI_Comm comm, MPI_Request *request)
                                                                       ;
int PMPI_Irsend(const void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm,
                MPI_Request *request) ;
int PMPI_Irsend_c(const void *buf, MPI_Count count, MPI_Datatype datatype, int dest, int tag,
                  MPI_Comm comm, MPI_Request *request)
                                                                        ;
int PMPI_Isend(const void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm,
               MPI_Request *request) ;
int PMPI_Isend_c(const void *buf, MPI_Count count, MPI_Datatype datatype, int dest, int tag,
                 MPI_Comm comm, MPI_Request *request)
                                                                       ;
int PMPI_Isendrecv(const void *sendbuf, int sendcount, MPI_Datatype sendtype, int dest, int sendtag,
                   void *recvbuf, int recvcount, MPI_Datatype recvtype, int source, int recvtag,
                   MPI_Comm comm, MPI_Request *request)
                                                                                                               ;
int PMPI_Isendrecv_c(const void *sendbuf, MPI_Count sendcount, MPI_Datatype sendtype, int dest,
                     int sendtag, void *recvbuf, MPI_Count recvcount, MPI_Datatype recvtype,
                     int source, int recvtag, MPI_Comm comm, MPI_Request *request)
                                                                                                                 ;
int PMPI_Isendrecv_replace(void *buf, int count, MPI_Datatype datatype, int dest, int sendtag,
                           int source, int recvtag, MPI_Comm comm, MPI_Request *request)
                                                                                 ;
int PMPI_Isendrecv_replace_c(void *buf, MPI_Count count, MPI_Datatype datatype, int dest,
                             int sendtag, int source, int recvtag, MPI_Comm comm,
                             MPI_Request *request)
                                                                                   ;
int PMPI_Issend(const void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm,
                MPI_Request *request) ;
int PMPI_Issend_c(const void *buf, MPI_Count count, MPI_Datatype datatype, int dest, int tag,
                  MPI_Comm comm, MPI_Request *request)
                                                                        ;
int PMPI_Mprobe(int source, int tag, MPI_Comm comm, MPI_Message *message, MPI_Status *status)
                    ;
int PMPI_Mrecv(void *buf, int count, MPI_Datatype datatype, MPI_Message *message,
               MPI_Status *status) ;
int PMPI_Mrecv_c(void *buf, MPI_Count count, MPI_Datatype datatype, MPI_Message *message,
                 MPI_Status *status) ;
int PMPI_Probe(int source, int tag, MPI_Comm comm, MPI_Status *status) ;
int PMPI_Recv(void *buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm,
              MPI_Status *status) ;
int PMPI_Recv_c(void *buf, MPI_Count count, MPI_Datatype datatype, int source, int tag,
                MPI_Comm comm, MPI_Status *status)
                                                                      ;
int PMPI_Recv_init(void *buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm,
                   MPI_Request *request) ;
int PMPI_Recv_init_c(void *buf, MPI_Count count, MPI_Datatype datatype, int source, int tag,
                     MPI_Comm comm, MPI_Request *request)
                                                                           ;
int PMPI_Rsend(const void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm)
                                                          ;
int PMPI_Rsend_c(const void *buf, MPI_Count count, MPI_Datatype datatype, int dest, int tag,
                 MPI_Comm comm) ;
int PMPI_Rsend_init(const void *buf, int count, MPI_Datatype datatype, int dest, int tag,
                    MPI_Comm comm, MPI_Request *request)
                                                                          ;
int PMPI_Rsend_init_c(const void *buf, MPI_Count count, MPI_Datatype datatype, int dest, int tag,
                      MPI_Comm comm, MPI_Request *request)
                                                                            ;
int PMPI_Send(const void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm)
                                                          ;
int PMPI_Send_c(const void *buf, MPI_Count count, MPI_Datatype datatype, int dest, int tag,
                MPI_Comm comm) ;
int PMPI_Send_init(const void *buf, int count, MPI_Datatype datatype, int dest, int tag,
                   MPI_Comm comm, MPI_Request *request)
                                                                         ;
int PMPI_Send_init_c(const void *buf, MPI_Count count, MPI_Datatype datatype, int dest, int tag,
                     MPI_Comm comm, MPI_Request *request)
                                                                           ;
int PMPI_Sendrecv(const void *sendbuf, int sendcount, MPI_Datatype sendtype, int dest, int sendtag,
                  void *recvbuf, int recvcount, MPI_Datatype recvtype, int source, int recvtag,
                  MPI_Comm comm, MPI_Status *status)
                                                                                                              ;
int PMPI_Sendrecv_c(const void *sendbuf, MPI_Count sendcount, MPI_Datatype sendtype, int dest,
                    int sendtag, void *recvbuf, MPI_Count recvcount, MPI_Datatype recvtype,
                    int source, int recvtag, MPI_Comm comm, MPI_Status *status)
                                                                                                                ;
int PMPI_Sendrecv_replace(void *buf, int count, MPI_Datatype datatype, int dest, int sendtag,
                          int source, int recvtag, MPI_Comm comm, MPI_Status *status)
                                                                                ;
int PMPI_Sendrecv_replace_c(void *buf, MPI_Count count, MPI_Datatype datatype, int dest,
                            int sendtag, int source, int recvtag, MPI_Comm comm,
                            MPI_Status *status)
                                                                                  ;
int PMPI_Session_attach_buffer(MPI_Session session, void *buffer, int size) ;
int PMPI_Session_attach_buffer_c(MPI_Session session, void *buffer, MPI_Count size)
                    ;
int PMPI_Session_detach_buffer(MPI_Session session, void *buffer_addr, int *size) ;
int PMPI_Session_detach_buffer_c(MPI_Session session, void *buffer_addr, MPI_Count *size)
                    ;
int PMPI_Session_flush_buffer(MPI_Session session) ;
int PMPI_Session_iflush_buffer(MPI_Session session, MPI_Request *request) ;
int PMPI_Ssend(const void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm)
                                                          ;
int PMPI_Ssend_c(const void *buf, MPI_Count count, MPI_Datatype datatype, int dest, int tag,
                 MPI_Comm comm) ;
int PMPI_Ssend_init(const void *buf, int count, MPI_Datatype datatype, int dest, int tag,
                    MPI_Comm comm, MPI_Request *request)
                                                                          ;
int PMPI_Ssend_init_c(const void *buf, MPI_Count count, MPI_Datatype datatype, int dest, int tag,
                      MPI_Comm comm, MPI_Request *request)
                                                                            ;
int PMPI_Cancel(MPI_Request *request) ;
int PMPI_Grequest_complete(MPI_Request request) ;
int PMPI_Grequest_start(MPI_Grequest_query_function *query_fn, MPI_Grequest_free_function *free_fn,
                        MPI_Grequest_cancel_function *cancel_fn, void *extra_state,
                        MPI_Request *request) ;
int PMPI_Request_free(MPI_Request *request) ;
int PMPI_Request_get_status(MPI_Request request, int *flag, MPI_Status *status) ;
int PMPI_Request_get_status_all(int count, const MPI_Request array_of_requests[], int *flag,
                                MPI_Status *array_of_statuses) ;
int PMPI_Request_get_status_any(int count, const MPI_Request array_of_requests[], int *indx,
                                int *flag, MPI_Status *status) ;
int PMPI_Request_get_status_some(int incount, const MPI_Request array_of_requests[], int *outcount,
                                 int array_of_indices[], MPI_Status *array_of_statuses)
                                                 ;
int PMPI_Start(MPI_Request *request) ;
int PMPI_Startall(int count, MPI_Request array_of_requests[]) ;
int PMPI_Status_get_error(const MPI_Status *status, int *error) ;
int PMPI_Status_get_source(const MPI_Status *status, int *source) ;
int PMPI_Status_get_tag(const MPI_Status *status, int *tag) ;
int PMPI_Status_set_error(MPI_Status *status, int error) ;
int PMPI_Status_set_source(MPI_Status *status, int source) ;
int PMPI_Status_set_tag(MPI_Status *status, int tag) ;
int PMPI_Status_set_cancelled(MPI_Status *status, int flag) ;
int PMPI_Test(MPI_Request *request, int *flag, MPI_Status *status) ;
int PMPI_Test_cancelled(const MPI_Status *status, int *flag) ;
int PMPI_Testall(int count, MPI_Request array_of_requests[], int *flag,
                 MPI_Status *array_of_statuses) ;
int PMPI_Testany(int count, MPI_Request array_of_requests[], int *indx, int *flag,
                 MPI_Status *status) ;
int PMPI_Testsome(int incount, MPI_Request array_of_requests[], int *outcount,
                  int array_of_indices[], MPI_Status *array_of_statuses) ;
int PMPI_Wait(MPI_Request *request, MPI_Status *status) ;
int PMPI_Waitall(int count, MPI_Request array_of_requests[], MPI_Status *array_of_statuses)
                    ;
int PMPI_Waitany(int count, MPI_Request array_of_requests[], int *indx, MPI_Status *status)
                    ;
int PMPI_Waitsome(int incount, MPI_Request array_of_requests[], int *outcount,
                  int array_of_indices[], MPI_Status *array_of_statuses) ;
int PMPIX_Grequest_start(MPI_Grequest_query_function *query_fn, MPI_Grequest_free_function *free_fn,
                         MPI_Grequest_cancel_function *cancel_fn,
                         MPIX_Grequest_poll_function *poll_fn, MPIX_Grequest_wait_function *wait_fn,
                         void *extra_state, MPI_Request *request) ;
int PMPIX_Grequest_class_create(MPI_Grequest_query_function *query_fn,
                                MPI_Grequest_free_function *free_fn,
                                MPI_Grequest_cancel_function *cancel_fn,
                                MPIX_Grequest_poll_function *poll_fn,
                                MPIX_Grequest_wait_function *wait_fn,
                                MPIX_Grequest_class *greq_class) ;
int PMPIX_Grequest_class_allocate(MPIX_Grequest_class greq_class, void *extra_state,
                                  MPI_Request *request) ;
int PMPIX_Request_is_complete(MPI_Request request) ;
int PMPI_Accumulate(const void *origin_addr, int origin_count, MPI_Datatype origin_datatype,
                    int target_rank, MPI_Aint target_disp, int target_count,
                    MPI_Datatype target_datatype, MPI_Op op, MPI_Win win)
                                                                          ;
int PMPI_Accumulate_c(const void *origin_addr, MPI_Count origin_count, MPI_Datatype origin_datatype,
                      int target_rank, MPI_Aint target_disp, MPI_Count target_count,
                      MPI_Datatype target_datatype, MPI_Op op, MPI_Win win)
                                                                            ;
int PMPI_Alloc_mem(MPI_Aint size, MPI_Info info, void *baseptr) ;
int PMPI_Compare_and_swap(const void *origin_addr, const void *compare_addr, void *result_addr,
                          MPI_Datatype datatype, int target_rank, MPI_Aint target_disp,
                          MPI_Win win) ;
int PMPI_Fetch_and_op(const void *origin_addr, void *result_addr, MPI_Datatype datatype,
                      int target_rank, MPI_Aint target_disp, MPI_Op op, MPI_Win win)
                                      ;
int PMPI_Free_mem(void *base) ;
int PMPI_Get(void *origin_addr, int origin_count, MPI_Datatype origin_datatype, int target_rank,
             MPI_Aint target_disp, int target_count, MPI_Datatype target_datatype, MPI_Win win)
                                                                   ;
int PMPI_Get_c(void *origin_addr, MPI_Count origin_count, MPI_Datatype origin_datatype,
               int target_rank, MPI_Aint target_disp, MPI_Count target_count,
               MPI_Datatype target_datatype, MPI_Win win)
                                                                     ;
int PMPI_Get_accumulate(const void *origin_addr, int origin_count, MPI_Datatype origin_datatype,
                        void *result_addr, int result_count, MPI_Datatype result_datatype,
                        int target_rank, MPI_Aint target_disp, int target_count,
                        MPI_Datatype target_datatype, MPI_Op op, MPI_Win win)
                                                                                                                    ;
int PMPI_Get_accumulate_c(const void *origin_addr, MPI_Count origin_count,
                          MPI_Datatype origin_datatype, void *result_addr, MPI_Count result_count,
                          MPI_Datatype result_datatype, int target_rank, MPI_Aint target_disp,
                          MPI_Count target_count, MPI_Datatype target_datatype, MPI_Op op,
                          MPI_Win win)
                                                                                                                      ;
int PMPI_Put(const void *origin_addr, int origin_count, MPI_Datatype origin_datatype,
             int target_rank, MPI_Aint target_disp, int target_count, MPI_Datatype target_datatype,
             MPI_Win win) ;
int PMPI_Put_c(const void *origin_addr, MPI_Count origin_count, MPI_Datatype origin_datatype,
               int target_rank, MPI_Aint target_disp, MPI_Count target_count,
               MPI_Datatype target_datatype, MPI_Win win)
                                                                     ;
int PMPI_Raccumulate(const void *origin_addr, int origin_count, MPI_Datatype origin_datatype,
                     int target_rank, MPI_Aint target_disp, int target_count,
                     MPI_Datatype target_datatype, MPI_Op op, MPI_Win win, MPI_Request *request)
                                                                           ;
int PMPI_Raccumulate_c(const void *origin_addr, MPI_Count origin_count,
                       MPI_Datatype origin_datatype, int target_rank, MPI_Aint target_disp,
                       MPI_Count target_count, MPI_Datatype target_datatype, MPI_Op op, MPI_Win win,
                       MPI_Request *request)
                                                                             ;
int PMPI_Rget(void *origin_addr, int origin_count, MPI_Datatype origin_datatype, int target_rank,
              MPI_Aint target_disp, int target_count, MPI_Datatype target_datatype, MPI_Win win,
              MPI_Request *request) ;
int PMPI_Rget_c(void *origin_addr, MPI_Count origin_count, MPI_Datatype origin_datatype,
                int target_rank, MPI_Aint target_disp, MPI_Count target_count,
                MPI_Datatype target_datatype, MPI_Win win, MPI_Request *request)
                                                                      ;
int PMPI_Rget_accumulate(const void *origin_addr, int origin_count, MPI_Datatype origin_datatype,
                         void *result_addr, int result_count, MPI_Datatype result_datatype,
                         int target_rank, MPI_Aint target_disp, int target_count,
                         MPI_Datatype target_datatype, MPI_Op op, MPI_Win win,
                         MPI_Request *request)
                                                                                                                     ;
int PMPI_Rget_accumulate_c(const void *origin_addr, MPI_Count origin_count,
                           MPI_Datatype origin_datatype, void *result_addr, MPI_Count result_count,
                           MPI_Datatype result_datatype, int target_rank, MPI_Aint target_disp,
                           MPI_Count target_count, MPI_Datatype target_datatype, MPI_Op op,
                           MPI_Win win, MPI_Request *request)
                                                                                                                       ;
int PMPI_Rput(const void *origin_addr, int origin_count, MPI_Datatype origin_datatype,
              int target_rank, MPI_Aint target_disp, int target_count, MPI_Datatype target_datatype,
              MPI_Win win, MPI_Request *request)
                                                                    ;
int PMPI_Rput_c(const void *origin_addr, MPI_Count origin_count, MPI_Datatype origin_datatype,
                int target_rank, MPI_Aint target_disp, MPI_Count target_count,
                MPI_Datatype target_datatype, MPI_Win win, MPI_Request *request)
                                                                      ;
int PMPI_Win_allocate(MPI_Aint size, int disp_unit, MPI_Info info, MPI_Comm comm, void *baseptr,
                      MPI_Win *win) ;
int PMPI_Win_allocate_c(MPI_Aint size, MPI_Aint disp_unit, MPI_Info info, MPI_Comm comm,
                        void *baseptr, MPI_Win *win) ;
int PMPI_Win_allocate_shared(MPI_Aint size, int disp_unit, MPI_Info info, MPI_Comm comm,
                             void *baseptr, MPI_Win *win) ;
int PMPI_Win_allocate_shared_c(MPI_Aint size, MPI_Aint disp_unit, MPI_Info info, MPI_Comm comm,
                               void *baseptr, MPI_Win *win) ;
int PMPI_Win_attach(MPI_Win win, void *base, MPI_Aint size) ;
int PMPI_Win_complete(MPI_Win win) ;
int PMPI_Win_create(void *base, MPI_Aint size, int disp_unit, MPI_Info info, MPI_Comm comm,
                    MPI_Win *win) ;
int PMPI_Win_create_c(void *base, MPI_Aint size, MPI_Aint disp_unit, MPI_Info info, MPI_Comm comm,
                      MPI_Win *win) ;
int PMPI_Win_create_dynamic(MPI_Info info, MPI_Comm comm, MPI_Win *win) ;
int PMPI_Win_detach(MPI_Win win, const void *base) ;
int PMPI_Win_fence(int assert, MPI_Win win) ;
int PMPI_Win_flush(int rank, MPI_Win win) ;
int PMPI_Win_flush_all(MPI_Win win) ;
int PMPI_Win_flush_local(int rank, MPI_Win win) ;
int PMPI_Win_flush_local_all(MPI_Win win) ;
int PMPI_Win_free(MPI_Win *win) ;
int PMPI_Win_get_group(MPI_Win win, MPI_Group *group) ;
int PMPI_Win_get_info(MPI_Win win, MPI_Info *info_used) ;
int PMPI_Win_get_name(MPI_Win win, char *win_name, int *resultlen) ;
int PMPI_Win_lock(int lock_type, int rank, int assert, MPI_Win win) ;
int PMPI_Win_lock_all(int assert, MPI_Win win) ;
int PMPI_Win_post(MPI_Group group, int assert, MPI_Win win) ;
int PMPI_Win_set_info(MPI_Win win, MPI_Info info) ;
int PMPI_Win_set_name(MPI_Win win, const char *win_name) ;
int PMPI_Win_shared_query(MPI_Win win, int rank, MPI_Aint *size, int *disp_unit, void *baseptr)
                    ;
int PMPI_Win_shared_query_c(MPI_Win win, int rank, MPI_Aint *size, MPI_Aint *disp_unit,
                            void *baseptr) ;
int PMPI_Win_start(MPI_Group group, int assert, MPI_Win win) ;
int PMPI_Win_sync(MPI_Win win) ;
int PMPI_Win_test(MPI_Win win, int *flag) ;
int PMPI_Win_unlock(int rank, MPI_Win win) ;
int PMPI_Win_unlock_all(MPI_Win win) ;
int PMPI_Win_wait(MPI_Win win) ;
int PMPI_Close_port(const char *port_name) ;
int PMPI_Comm_accept(const char *port_name, MPI_Info info, int root, MPI_Comm comm,
                     MPI_Comm *newcomm) ;
int PMPI_Comm_connect(const char *port_name, MPI_Info info, int root, MPI_Comm comm,
                      MPI_Comm *newcomm) ;
int PMPI_Comm_disconnect(MPI_Comm *comm) ;
int PMPI_Comm_get_parent(MPI_Comm *parent) ;
int PMPI_Comm_join(int fd, MPI_Comm *intercomm) ;
int PMPI_Comm_spawn(const char *command, char *argv[], int maxprocs, MPI_Info info, int root,
                    MPI_Comm comm, MPI_Comm *intercomm, int array_of_errcodes[]) ;
int PMPI_Comm_spawn_multiple(int count, char *array_of_commands[], char **array_of_argv[],
                             const int array_of_maxprocs[], const MPI_Info array_of_info[],
                             int root, MPI_Comm comm, MPI_Comm *intercomm, int array_of_errcodes[])
                                             ;
int PMPI_Lookup_name(const char *service_name, MPI_Info info, char *port_name) ;
int PMPI_Open_port(MPI_Info info, char *port_name) ;
int PMPI_Publish_name(const char *service_name, MPI_Info info, const char *port_name)
                    ;
int PMPI_Unpublish_name(const char *service_name, MPI_Info info, const char *port_name)
                    ;
int PMPIX_Stream_create(MPI_Info info, MPIX_Stream *stream) ;
int PMPIX_Stream_free(MPIX_Stream *stream) ;
int PMPIX_Stream_comm_create(MPI_Comm comm, MPIX_Stream stream, MPI_Comm *newcomm)
                    ;
int PMPIX_Stream_comm_create_multiplex(MPI_Comm comm, int count, MPIX_Stream array_of_streams[],
                                       MPI_Comm *newcomm) ;
int PMPIX_Comm_get_stream(MPI_Comm comm, int idx, MPIX_Stream *stream) ;
int PMPIX_Stream_progress(MPIX_Stream stream) ;
int PMPIX_Start_progress_thread(MPIX_Stream stream) ;
int PMPIX_Stop_progress_thread(MPIX_Stream stream) ;
int PMPIX_Stream_send(const void *buf, int count, MPI_Datatype datatype, int dest, int tag,
                      MPI_Comm comm, int source_stream_index, int dest_stream_index)
                                                                            ;
int PMPIX_Stream_send_c(const void *buf, MPI_Count count, MPI_Datatype datatype, int dest, int tag,
                        MPI_Comm comm, int source_stream_index, int dest_stream_index)
                                                                              ;
int PMPIX_Stream_isend(const void *buf, int count, MPI_Datatype datatype, int dest, int tag,
                       MPI_Comm comm, int source_stream_index, int dest_stream_index,
                       MPI_Request *request)
                                                                             ;
int PMPIX_Stream_isend_c(const void *buf, MPI_Count count, MPI_Datatype datatype, int dest, int tag,
                         MPI_Comm comm, int source_stream_index, int dest_stream_index,
                         MPI_Request *request)
                                                                               ;
int PMPIX_Stream_recv(void *buf, int count, MPI_Datatype datatype, int source, int tag,
                      MPI_Comm comm, int source_stream_index, int dest_stream_index,
                      MPI_Status *status) ;
int PMPIX_Stream_recv_c(void *buf, MPI_Count count, MPI_Datatype datatype, int source, int tag,
                        MPI_Comm comm, int source_stream_index, int dest_stream_index,
                        MPI_Status *status) ;
int PMPIX_Stream_irecv(void *buf, int count, MPI_Datatype datatype, int source, int tag,
                       MPI_Comm comm, int source_stream_index, int dest_stream_index,
                       MPI_Request *request)
                                                                             ;
int PMPIX_Stream_irecv_c(void *buf, MPI_Count count, MPI_Datatype datatype, int source, int tag,
                         MPI_Comm comm, int source_stream_index, int dest_stream_index,
                         MPI_Request *request)
                                                                               ;
int PMPIX_Send_enqueue(const void *buf, int count, MPI_Datatype datatype, int dest, int tag,
                       MPI_Comm comm) ;
int PMPIX_Send_enqueue_c(const void *buf, MPI_Count count, MPI_Datatype datatype, int dest, int tag,
                         MPI_Comm comm) ;
int PMPIX_Recv_enqueue(void *buf, int count, MPI_Datatype datatype, int source, int tag,
                       MPI_Comm comm, MPI_Status *status)
                                                                             ;
int PMPIX_Recv_enqueue_c(void *buf, MPI_Count count, MPI_Datatype datatype, int source, int tag,
                         MPI_Comm comm, MPI_Status *status)
                                                                               ;
int PMPIX_Isend_enqueue(const void *buf, int count, MPI_Datatype datatype, int dest, int tag,
                        MPI_Comm comm, MPI_Request *request)
                                                                              ;
int PMPIX_Isend_enqueue_c(const void *buf, MPI_Count count, MPI_Datatype datatype, int dest,
                          int tag, MPI_Comm comm, MPI_Request *request)
                                                                                ;
int PMPIX_Irecv_enqueue(void *buf, int count, MPI_Datatype datatype, int source, int tag,
                        MPI_Comm comm, MPI_Request *request)
                                                                              ;
int PMPIX_Irecv_enqueue_c(void *buf, MPI_Count count, MPI_Datatype datatype, int source, int tag,
                          MPI_Comm comm, MPI_Request *request)
                                                                                ;
int PMPIX_Wait_enqueue(MPI_Request *request, MPI_Status *status) ;
int PMPIX_Waitall_enqueue(int count, MPI_Request array_of_requests[],
                          MPI_Status *array_of_statuses) ;
int PMPIX_Allreduce_enqueue(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype,
                            MPI_Op op, MPI_Comm comm)
                                                                                  ;
int PMPIX_Allreduce_enqueue_c(const void *sendbuf, void *recvbuf, MPI_Count count,
                              MPI_Datatype datatype, MPI_Op op, MPI_Comm comm)
                                                                                    ;
int PMPIX_Threadcomm_init(MPI_Comm comm, int num_threads, MPI_Comm *newthreadcomm)
                    ;
int PMPIX_Threadcomm_free(MPI_Comm *threadcomm) ;
int PMPIX_Threadcomm_start(MPI_Comm threadcomm) ;
int PMPIX_Threadcomm_finish(MPI_Comm threadcomm) ;
double PMPI_Wtick(void) ;
double PMPI_Wtime(void) ;
int PMPI_Cart_coords(MPI_Comm comm, int rank, int maxdims, int coords[]) ;
int PMPI_Cart_create(MPI_Comm comm_old, int ndims, const int dims[], const int periods[],
                     int reorder, MPI_Comm *comm_cart) ;
int PMPI_Cart_get(MPI_Comm comm, int maxdims, int dims[], int periods[], int coords[])
                    ;
int PMPI_Cart_map(MPI_Comm comm, int ndims, const int dims[], const int periods[], int *newrank)
                    ;
int PMPI_Cart_rank(MPI_Comm comm, const int coords[], int *rank) ;
int PMPI_Cart_shift(MPI_Comm comm, int direction, int disp, int *rank_source, int *rank_dest)
                    ;
int PMPI_Cart_sub(MPI_Comm comm, const int remain_dims[], MPI_Comm *newcomm) ;
int PMPI_Cartdim_get(MPI_Comm comm, int *ndims) ;
int PMPI_Dims_create(int nnodes, int ndims, int dims[]) ;
int PMPI_Dist_graph_create(MPI_Comm comm_old, int n, const int sources[], const int degrees[],
                           const int destinations[], const int weights[], MPI_Info info,
                           int reorder, MPI_Comm *comm_dist_graph) ;
int PMPI_Dist_graph_create_adjacent(MPI_Comm comm_old, int indegree, const int sources[],
                                    const int sourceweights[], int outdegree,
                                    const int destinations[], const int destweights[],
                                    MPI_Info info, int reorder, MPI_Comm *comm_dist_graph)
                                                    ;
int PMPI_Dist_graph_neighbors(MPI_Comm comm, int maxindegree, int sources[], int sourceweights[],
                              int maxoutdegree, int destinations[], int destweights[])
                                              ;
int PMPI_Dist_graph_neighbors_count(MPI_Comm comm, int *indegree, int *outdegree, int *weighted)
                    ;
int PMPI_Get_hw_resource_info(MPI_Info *hw_info) ;
int PMPI_Graph_create(MPI_Comm comm_old, int nnodes, const int indx[], const int edges[],
                      int reorder, MPI_Comm *comm_graph) ;
int PMPI_Graph_get(MPI_Comm comm, int maxindex, int maxedges, int indx[], int edges[])
                    ;
int PMPI_Graph_map(MPI_Comm comm, int nnodes, const int indx[], const int edges[], int *newrank)
                    ;
int PMPI_Graph_neighbors(MPI_Comm comm, int rank, int maxneighbors, int neighbors[])
                    ;
int PMPI_Graph_neighbors_count(MPI_Comm comm, int rank, int *nneighbors) ;
int PMPI_Graphdims_get(MPI_Comm comm, int *nnodes, int *nedges) ;
int PMPI_Topo_test(MPI_Comm comm, int *status) ;
MPI_Fint PMPI_File_c2f(MPI_File file) ;
int PMPI_File_close(MPI_File *fh) ;
int PMPI_File_delete(const char *filename, MPI_Info info) ;
MPI_File PMPI_File_f2c(MPI_Fint file) ;
int PMPI_File_get_amode(MPI_File fh, int *amode) ;
int PMPI_File_get_atomicity(MPI_File fh, int *flag) ;
int PMPI_File_get_byte_offset(MPI_File fh, MPI_Offset offset, MPI_Offset *disp) ;
int PMPI_File_get_group(MPI_File fh, MPI_Group *group) ;
int PMPI_File_get_info(MPI_File fh, MPI_Info *info_used) ;
int PMPI_File_get_position(MPI_File fh, MPI_Offset *offset) ;
int PMPI_File_get_position_shared(MPI_File fh, MPI_Offset *offset) ;
int PMPI_File_get_size(MPI_File fh, MPI_Offset *size) ;
int PMPI_File_get_type_extent(MPI_File fh, MPI_Datatype datatype, MPI_Aint *extent)
                    ;
int PMPI_File_get_type_extent_c(MPI_File fh, MPI_Datatype datatype, MPI_Count *extent)
                    ;
int PMPI_File_get_view(MPI_File fh, MPI_Offset *disp, MPI_Datatype *etype, MPI_Datatype *filetype,
                       char *datarep) ;
int PMPI_File_iread(MPI_File fh, void *buf, int count, MPI_Datatype datatype, MPI_Request *request)
                                                          ;
int PMPI_File_iread_c(MPI_File fh, void *buf, MPI_Count count, MPI_Datatype datatype,
                      MPI_Request *request) ;
int PMPI_File_iread_all(MPI_File fh, void *buf, int count, MPI_Datatype datatype,
                        MPI_Request *request)
                                                                              ;
int PMPI_File_iread_all_c(MPI_File fh, void *buf, MPI_Count count, MPI_Datatype datatype,
                          MPI_Request *request)
                                                                                ;
int PMPI_File_iread_at(MPI_File fh, MPI_Offset offset, void *buf, int count, MPI_Datatype datatype,
                       MPI_Request *request)
                                                                             ;
int PMPI_File_iread_at_c(MPI_File fh, MPI_Offset offset, void *buf, MPI_Count count,
                         MPI_Datatype datatype, MPI_Request *request)
                                                                               ;
int PMPI_File_iread_at_all(MPI_File fh, MPI_Offset offset, void *buf, int count,
                           MPI_Datatype datatype, MPI_Request *request)
                                                                                 ;
int PMPI_File_iread_at_all_c(MPI_File fh, MPI_Offset offset, void *buf, MPI_Count count,
                             MPI_Datatype datatype, MPI_Request *request)
                                                                                   ;
int PMPI_File_iread_shared(MPI_File fh, void *buf, int count, MPI_Datatype datatype,
                           MPI_Request *request)
                                                                                 ;
int PMPI_File_iread_shared_c(MPI_File fh, void *buf, MPI_Count count, MPI_Datatype datatype,
                             MPI_Request *request)
                                                                                   ;
int PMPI_File_iwrite(MPI_File fh, const void *buf, int count, MPI_Datatype datatype,
                     MPI_Request *request) ;
int PMPI_File_iwrite_c(MPI_File fh, const void *buf, MPI_Count count, MPI_Datatype datatype,
                       MPI_Request *request)
                                                                             ;
int PMPI_File_iwrite_all(MPI_File fh, const void *buf, int count, MPI_Datatype datatype,
                         MPI_Request *request)
                                                                               ;
int PMPI_File_iwrite_all_c(MPI_File fh, const void *buf, MPI_Count count, MPI_Datatype datatype,
                           MPI_Request *request)
                                                                                 ;
int PMPI_File_iwrite_at(MPI_File fh, MPI_Offset offset, const void *buf, int count,
                        MPI_Datatype datatype, MPI_Request *request)
                                                                              ;
int PMPI_File_iwrite_at_c(MPI_File fh, MPI_Offset offset, const void *buf, MPI_Count count,
                          MPI_Datatype datatype, MPI_Request *request)
                                                                                ;
int PMPI_File_iwrite_at_all(MPI_File fh, MPI_Offset offset, const void *buf, int count,
                            MPI_Datatype datatype, MPI_Request *request)
                                                                                  ;
int PMPI_File_iwrite_at_all_c(MPI_File fh, MPI_Offset offset, const void *buf, MPI_Count count,
                              MPI_Datatype datatype, MPI_Request *request)
                                                                                    ;
int PMPI_File_iwrite_shared(MPI_File fh, const void *buf, int count, MPI_Datatype datatype,
                            MPI_Request *request)
                                                                                  ;
int PMPI_File_iwrite_shared_c(MPI_File fh, const void *buf, MPI_Count count, MPI_Datatype datatype,
                              MPI_Request *request)
                                                                                    ;
int PMPI_File_open(MPI_Comm comm, const char *filename, int amode, MPI_Info info, MPI_File *fh)
                    ;
int PMPI_File_preallocate(MPI_File fh, MPI_Offset size) ;
int PMPI_File_read(MPI_File fh, void *buf, int count, MPI_Datatype datatype, MPI_Status *status)
                                                          ;
int PMPI_File_read_c(MPI_File fh, void *buf, MPI_Count count, MPI_Datatype datatype,
                     MPI_Status *status) ;
int PMPI_File_read_all(MPI_File fh, void *buf, int count, MPI_Datatype datatype,
                       MPI_Status *status) ;
int PMPI_File_read_all_c(MPI_File fh, void *buf, MPI_Count count, MPI_Datatype datatype,
                         MPI_Status *status)
                                                                               ;
int PMPI_File_read_all_begin(MPI_File fh, void *buf, int count, MPI_Datatype datatype)
                                                          ;
int PMPI_File_read_all_begin_c(MPI_File fh, void *buf, MPI_Count count, MPI_Datatype datatype)
                                                          ;
int PMPI_File_read_all_end(MPI_File fh, void *buf, MPI_Status *status) ;
int PMPI_File_read_at(MPI_File fh, MPI_Offset offset, void *buf, int count, MPI_Datatype datatype,
                      MPI_Status *status) ;
int PMPI_File_read_at_c(MPI_File fh, MPI_Offset offset, void *buf, MPI_Count count,
                        MPI_Datatype datatype, MPI_Status *status)
                                                                              ;
int PMPI_File_read_at_all(MPI_File fh, MPI_Offset offset, void *buf, int count,
                          MPI_Datatype datatype, MPI_Status *status)
                                                                                ;
int PMPI_File_read_at_all_c(MPI_File fh, MPI_Offset offset, void *buf, MPI_Count count,
                            MPI_Datatype datatype, MPI_Status *status)
                                                                                  ;
int PMPI_File_read_at_all_begin(MPI_File fh, MPI_Offset offset, void *buf, int count,
                                MPI_Datatype datatype)
                                                                                      ;
int PMPI_File_read_at_all_begin_c(MPI_File fh, MPI_Offset offset, void *buf, MPI_Count count,
                                  MPI_Datatype datatype)
                                                                                        ;
int PMPI_File_read_at_all_end(MPI_File fh, void *buf, MPI_Status *status) ;
int PMPI_File_read_ordered(MPI_File fh, void *buf, int count, MPI_Datatype datatype,
                           MPI_Status *status)
                                                                                 ;
int PMPI_File_read_ordered_c(MPI_File fh, void *buf, MPI_Count count, MPI_Datatype datatype,
                             MPI_Status *status)
                                                                                   ;
int PMPI_File_read_ordered_begin(MPI_File fh, void *buf, int count, MPI_Datatype datatype)
                                                          ;
int PMPI_File_read_ordered_begin_c(MPI_File fh, void *buf, MPI_Count count, MPI_Datatype datatype)
                                                          ;
int PMPI_File_read_ordered_end(MPI_File fh, void *buf, MPI_Status *status) ;
int PMPI_File_read_shared(MPI_File fh, void *buf, int count, MPI_Datatype datatype,
                          MPI_Status *status)
                                                                                ;
int PMPI_File_read_shared_c(MPI_File fh, void *buf, MPI_Count count, MPI_Datatype datatype,
                            MPI_Status *status)
                                                                                  ;
int PMPI_File_seek(MPI_File fh, MPI_Offset offset, int whence) ;
int PMPI_File_seek_shared(MPI_File fh, MPI_Offset offset, int whence) ;
int PMPI_File_set_atomicity(MPI_File fh, int flag) ;
int PMPI_File_set_info(MPI_File fh, MPI_Info info) ;
int PMPI_File_set_size(MPI_File fh, MPI_Offset size) ;
int PMPI_File_set_view(MPI_File fh, MPI_Offset disp, MPI_Datatype etype, MPI_Datatype filetype,
                       const char *datarep, MPI_Info info) ;
int PMPI_File_sync(MPI_File fh) ;
int PMPI_File_write(MPI_File fh, const void *buf, int count, MPI_Datatype datatype,
                    MPI_Status *status) ;
int PMPI_File_write_c(MPI_File fh, const void *buf, MPI_Count count, MPI_Datatype datatype,
                      MPI_Status *status) ;
int PMPI_File_write_all(MPI_File fh, const void *buf, int count, MPI_Datatype datatype,
                        MPI_Status *status) ;
int PMPI_File_write_all_c(MPI_File fh, const void *buf, MPI_Count count, MPI_Datatype datatype,
                          MPI_Status *status)
                                                                                ;
int PMPI_File_write_all_begin(MPI_File fh, const void *buf, int count, MPI_Datatype datatype)
                                                          ;
int PMPI_File_write_all_begin_c(MPI_File fh, const void *buf, MPI_Count count,
                                MPI_Datatype datatype)
                                                                                      ;
int PMPI_File_write_all_end(MPI_File fh, const void *buf, MPI_Status *status) ;
int PMPI_File_write_at(MPI_File fh, MPI_Offset offset, const void *buf, int count,
                       MPI_Datatype datatype, MPI_Status *status)
                                                                             ;
int PMPI_File_write_at_c(MPI_File fh, MPI_Offset offset, const void *buf, MPI_Count count,
                         MPI_Datatype datatype, MPI_Status *status)
                                                                               ;
int PMPI_File_write_at_all(MPI_File fh, MPI_Offset offset, const void *buf, int count,
                           MPI_Datatype datatype, MPI_Status *status)
                                                                                 ;
int PMPI_File_write_at_all_c(MPI_File fh, MPI_Offset offset, const void *buf, MPI_Count count,
                             MPI_Datatype datatype, MPI_Status *status)
                                                                                   ;
int PMPI_File_write_at_all_begin(MPI_File fh, MPI_Offset offset, const void *buf, int count,
                                 MPI_Datatype datatype)
                                                                                       ;
int PMPI_File_write_at_all_begin_c(MPI_File fh, MPI_Offset offset, const void *buf, MPI_Count count,
                                   MPI_Datatype datatype)
                                                                                         ;
int PMPI_File_write_at_all_end(MPI_File fh, const void *buf, MPI_Status *status) ;
int PMPI_File_write_ordered(MPI_File fh, const void *buf, int count, MPI_Datatype datatype,
                            MPI_Status *status)
                                                                                  ;
int PMPI_File_write_ordered_c(MPI_File fh, const void *buf, MPI_Count count, MPI_Datatype datatype,
                              MPI_Status *status)
                                                                                    ;
int PMPI_File_write_ordered_begin(MPI_File fh, const void *buf, int count, MPI_Datatype datatype)
                                                          ;
int PMPI_File_write_ordered_begin_c(MPI_File fh, const void *buf, MPI_Count count,
                                    MPI_Datatype datatype)
                                                                                          ;
int PMPI_File_write_ordered_end(MPI_File fh, const void *buf, MPI_Status *status) ;
int PMPI_File_write_shared(MPI_File fh, const void *buf, int count, MPI_Datatype datatype,
                           MPI_Status *status)
                                                                                 ;
int PMPI_File_write_shared_c(MPI_File fh, const void *buf, MPI_Count count, MPI_Datatype datatype,
                             MPI_Status *status)
                                                                                   ;
int PMPI_Register_datarep(const char *datarep, MPI_Datarep_conversion_function *read_conversion_fn,
                          MPI_Datarep_conversion_function *write_conversion_fn,
                          MPI_Datarep_extent_function *dtype_file_extent_fn, void *extra_state)
                                          ;
int PMPI_Register_datarep_c(const char *datarep,
                            MPI_Datarep_conversion_function_c *read_conversion_fn,
                            MPI_Datarep_conversion_function_c *write_conversion_fn,
                            MPI_Datarep_extent_function *dtype_file_extent_fn, void *extra_state)
                                            ;
int PMPI_File_toint(MPI_File file) ;
MPI_File PMPI_File_fromint(int file) ;
enum QMPI_Functions_enum {
    MPI_ABI_GET_FORTRAN_BOOLEANS_T,
    MPI_ABI_GET_FORTRAN_INFO_T,
    MPI_ABI_GET_INFO_T,
    MPI_ABI_GET_VERSION_T,
    MPI_ABI_SET_FORTRAN_BOOLEANS_T,
    MPI_ABI_SET_FORTRAN_INFO_T,
    MPI_COMM_TOINT_T,
    MPI_COMM_FROMINT_T,
    MPI_ERRHANDLER_TOINT_T,
    MPI_ERRHANDLER_FROMINT_T,
    MPI_GROUP_TOINT_T,
    MPI_GROUP_FROMINT_T,
    MPI_INFO_TOINT_T,
    MPI_INFO_FROMINT_T,
    MPI_MESSAGE_TOINT_T,
    MPI_MESSAGE_FROMINT_T,
    MPI_OP_TOINT_T,
    MPI_OP_FROMINT_T,
    MPI_REQUEST_TOINT_T,
    MPI_REQUEST_FROMINT_T,
    MPI_SESSION_TOINT_T,
    MPI_SESSION_FROMINT_T,
    MPI_TYPE_TOINT_T,
    MPI_TYPE_FROMINT_T,
    MPI_WIN_TOINT_T,
    MPI_WIN_FROMINT_T,
    MPIX_ASYNC_START_T,
    MPIX_ASYNC_GET_STATE_T,
    MPIX_ASYNC_SPAWN_T,
    MPI_COMM_CREATE_KEYVAL_T,
    MPI_KEYVAL_CREATE_T,
    MPI_COMM_DELETE_ATTR_T,
    MPI_ATTR_DELETE_T,
    MPI_COMM_FREE_KEYVAL_T,
    MPI_KEYVAL_FREE_T,
    MPI_COMM_GET_ATTR_T,
    MPI_ATTR_GET_T,
    MPI_COMM_SET_ATTR_T,
    MPI_ATTR_PUT_T,
    MPI_TYPE_CREATE_KEYVAL_T,
    MPI_TYPE_DELETE_ATTR_T,
    MPI_TYPE_FREE_KEYVAL_T,
    MPI_TYPE_GET_ATTR_T,
    MPI_TYPE_SET_ATTR_T,
    MPI_WIN_CREATE_KEYVAL_T,
    MPI_WIN_DELETE_ATTR_T,
    MPI_WIN_FREE_KEYVAL_T,
    MPI_WIN_GET_ATTR_T,
    MPI_WIN_SET_ATTR_T,
    MPIX_OP_CREATE_X_T,
    MPIX_COMM_CREATE_ERRHANDLER_X_T,
    MPIX_WIN_CREATE_ERRHANDLER_X_T,
    MPIX_FILE_CREATE_ERRHANDLER_X_T,
    MPIX_SESSION_CREATE_ERRHANDLER_X_T,
    MPI_ALLGATHER_T,
    MPI_ALLGATHER_C_T,
    MPI_ALLGATHER_INIT_T,
    MPI_ALLGATHER_INIT_C_T,
    MPI_ALLGATHERV_T,
    MPI_ALLGATHERV_C_T,
    MPI_ALLGATHERV_INIT_T,
    MPI_ALLGATHERV_INIT_C_T,
    MPI_ALLREDUCE_T,
    MPI_ALLREDUCE_C_T,
    MPI_ALLREDUCE_INIT_T,
    MPI_ALLREDUCE_INIT_C_T,
    MPI_ALLTOALL_T,
    MPI_ALLTOALL_C_T,
    MPI_ALLTOALL_INIT_T,
    MPI_ALLTOALL_INIT_C_T,
    MPI_ALLTOALLV_T,
    MPI_ALLTOALLV_C_T,
    MPI_ALLTOALLV_INIT_T,
    MPI_ALLTOALLV_INIT_C_T,
    MPI_ALLTOALLW_T,
    MPI_ALLTOALLW_C_T,
    MPI_ALLTOALLW_INIT_T,
    MPI_ALLTOALLW_INIT_C_T,
    MPI_BARRIER_T,
    MPI_BARRIER_INIT_T,
    MPI_BCAST_T,
    MPI_BCAST_C_T,
    MPI_BCAST_INIT_T,
    MPI_BCAST_INIT_C_T,
    MPI_EXSCAN_T,
    MPI_EXSCAN_C_T,
    MPI_EXSCAN_INIT_T,
    MPI_EXSCAN_INIT_C_T,
    MPI_GATHER_T,
    MPI_GATHER_C_T,
    MPI_GATHER_INIT_T,
    MPI_GATHER_INIT_C_T,
    MPI_GATHERV_T,
    MPI_GATHERV_C_T,
    MPI_GATHERV_INIT_T,
    MPI_GATHERV_INIT_C_T,
    MPI_IALLGATHER_T,
    MPI_IALLGATHER_C_T,
    MPI_IALLGATHERV_T,
    MPI_IALLGATHERV_C_T,
    MPI_IALLREDUCE_T,
    MPI_IALLREDUCE_C_T,
    MPI_IALLTOALL_T,
    MPI_IALLTOALL_C_T,
    MPI_IALLTOALLV_T,
    MPI_IALLTOALLV_C_T,
    MPI_IALLTOALLW_T,
    MPI_IALLTOALLW_C_T,
    MPI_IBARRIER_T,
    MPI_IBCAST_T,
    MPI_IBCAST_C_T,
    MPI_IEXSCAN_T,
    MPI_IEXSCAN_C_T,
    MPI_IGATHER_T,
    MPI_IGATHER_C_T,
    MPI_IGATHERV_T,
    MPI_IGATHERV_C_T,
    MPI_INEIGHBOR_ALLGATHER_T,
    MPI_INEIGHBOR_ALLGATHER_C_T,
    MPI_INEIGHBOR_ALLGATHERV_T,
    MPI_INEIGHBOR_ALLGATHERV_C_T,
    MPI_INEIGHBOR_ALLTOALL_T,
    MPI_INEIGHBOR_ALLTOALL_C_T,
    MPI_INEIGHBOR_ALLTOALLV_T,
    MPI_INEIGHBOR_ALLTOALLV_C_T,
    MPI_INEIGHBOR_ALLTOALLW_T,
    MPI_INEIGHBOR_ALLTOALLW_C_T,
    MPI_IREDUCE_T,
    MPI_IREDUCE_C_T,
    MPI_IREDUCE_SCATTER_T,
    MPI_IREDUCE_SCATTER_C_T,
    MPI_IREDUCE_SCATTER_BLOCK_T,
    MPI_IREDUCE_SCATTER_BLOCK_C_T,
    MPI_ISCAN_T,
    MPI_ISCAN_C_T,
    MPI_ISCATTER_T,
    MPI_ISCATTER_C_T,
    MPI_ISCATTERV_T,
    MPI_ISCATTERV_C_T,
    MPI_NEIGHBOR_ALLGATHER_T,
    MPI_NEIGHBOR_ALLGATHER_C_T,
    MPI_NEIGHBOR_ALLGATHER_INIT_T,
    MPI_NEIGHBOR_ALLGATHER_INIT_C_T,
    MPI_NEIGHBOR_ALLGATHERV_T,
    MPI_NEIGHBOR_ALLGATHERV_C_T,
    MPI_NEIGHBOR_ALLGATHERV_INIT_T,
    MPI_NEIGHBOR_ALLGATHERV_INIT_C_T,
    MPI_NEIGHBOR_ALLTOALL_T,
    MPI_NEIGHBOR_ALLTOALL_C_T,
    MPI_NEIGHBOR_ALLTOALL_INIT_T,
    MPI_NEIGHBOR_ALLTOALL_INIT_C_T,
    MPI_NEIGHBOR_ALLTOALLV_T,
    MPI_NEIGHBOR_ALLTOALLV_C_T,
    MPI_NEIGHBOR_ALLTOALLV_INIT_T,
    MPI_NEIGHBOR_ALLTOALLV_INIT_C_T,
    MPI_NEIGHBOR_ALLTOALLW_T,
    MPI_NEIGHBOR_ALLTOALLW_C_T,
    MPI_NEIGHBOR_ALLTOALLW_INIT_T,
    MPI_NEIGHBOR_ALLTOALLW_INIT_C_T,
    MPI_REDUCE_T,
    MPI_REDUCE_C_T,
    MPI_REDUCE_INIT_T,
    MPI_REDUCE_INIT_C_T,
    MPI_REDUCE_LOCAL_T,
    MPI_REDUCE_LOCAL_C_T,
    MPI_REDUCE_SCATTER_T,
    MPI_REDUCE_SCATTER_C_T,
    MPI_REDUCE_SCATTER_BLOCK_T,
    MPI_REDUCE_SCATTER_BLOCK_C_T,
    MPI_REDUCE_SCATTER_BLOCK_INIT_T,
    MPI_REDUCE_SCATTER_BLOCK_INIT_C_T,
    MPI_REDUCE_SCATTER_INIT_T,
    MPI_REDUCE_SCATTER_INIT_C_T,
    MPI_SCAN_T,
    MPI_SCAN_C_T,
    MPI_SCAN_INIT_T,
    MPI_SCAN_INIT_C_T,
    MPI_SCATTER_T,
    MPI_SCATTER_C_T,
    MPI_SCATTER_INIT_T,
    MPI_SCATTER_INIT_C_T,
    MPI_SCATTERV_T,
    MPI_SCATTERV_C_T,
    MPI_SCATTERV_INIT_T,
    MPI_SCATTERV_INIT_C_T,
    MPI_COMM_COMPARE_T,
    MPI_COMM_CREATE_T,
    MPI_COMM_CREATE_GROUP_T,
    MPI_COMM_DUP_T,
    MPI_COMM_DUP_WITH_INFO_T,
    MPI_COMM_FREE_T,
    MPI_COMM_GET_INFO_T,
    MPI_COMM_GET_NAME_T,
    MPI_COMM_GROUP_T,
    MPI_COMM_IDUP_T,
    MPI_COMM_IDUP_WITH_INFO_T,
    MPI_COMM_RANK_T,
    MPI_COMM_REMOTE_GROUP_T,
    MPI_COMM_REMOTE_SIZE_T,
    MPI_COMM_SET_INFO_T,
    MPI_COMM_SET_NAME_T,
    MPI_COMM_SIZE_T,
    MPI_COMM_SPLIT_T,
    MPI_COMM_SPLIT_TYPE_T,
    MPI_COMM_TEST_INTER_T,
    MPI_INTERCOMM_CREATE_T,
    MPI_INTERCOMM_CREATE_FROM_GROUPS_T,
    MPI_INTERCOMM_MERGE_T,
    MPIX_COMM_TEST_THREADCOMM_T,
    MPIX_COMM_REVOKE_T,
    MPIX_COMM_SHRINK_T,
    MPIX_COMM_FAILURE_ACK_T,
    MPIX_COMM_FAILURE_GET_ACKED_T,
    MPIX_COMM_AGREE_T,
    MPIX_COMM_GET_FAILED_T,
    MPI_GET_ADDRESS_T,
    MPI_ADDRESS_T,
    MPI_GET_COUNT_T,
    MPI_GET_COUNT_C_T,
    MPI_GET_ELEMENTS_T,
    MPI_GET_ELEMENTS_C_T,
    MPI_GET_ELEMENTS_X_T,
    MPI_PACK_T,
    MPI_PACK_C_T,
    MPI_PACK_EXTERNAL_T,
    MPI_PACK_EXTERNAL_C_T,
    MPI_PACK_EXTERNAL_SIZE_T,
    MPI_PACK_EXTERNAL_SIZE_C_T,
    MPI_PACK_SIZE_T,
    MPI_PACK_SIZE_C_T,
    MPI_STATUS_SET_ELEMENTS_T,
    MPI_STATUS_SET_ELEMENTS_C_T,
    MPI_STATUS_SET_ELEMENTS_X_T,
    MPI_TYPE_COMMIT_T,
    MPI_TYPE_CONTIGUOUS_T,
    MPI_TYPE_CONTIGUOUS_C_T,
    MPI_TYPE_CREATE_DARRAY_T,
    MPI_TYPE_CREATE_DARRAY_C_T,
    MPI_TYPE_CREATE_F90_COMPLEX_T,
    MPI_TYPE_CREATE_F90_INTEGER_T,
    MPI_TYPE_CREATE_F90_REAL_T,
    MPI_TYPE_CREATE_HINDEXED_T,
    MPI_TYPE_CREATE_HINDEXED_C_T,
    MPI_TYPE_HINDEXED_T,
    MPI_TYPE_CREATE_HINDEXED_BLOCK_T,
    MPI_TYPE_CREATE_HINDEXED_BLOCK_C_T,
    MPI_TYPE_CREATE_HVECTOR_T,
    MPI_TYPE_CREATE_HVECTOR_C_T,
    MPI_TYPE_HVECTOR_T,
    MPI_TYPE_CREATE_INDEXED_BLOCK_T,
    MPI_TYPE_CREATE_INDEXED_BLOCK_C_T,
    MPI_TYPE_CREATE_RESIZED_T,
    MPI_TYPE_CREATE_RESIZED_C_T,
    MPI_TYPE_CREATE_STRUCT_T,
    MPI_TYPE_CREATE_STRUCT_C_T,
    MPI_TYPE_STRUCT_T,
    MPI_TYPE_CREATE_SUBARRAY_T,
    MPI_TYPE_CREATE_SUBARRAY_C_T,
    MPI_TYPE_DUP_T,
    MPI_TYPE_FREE_T,
    MPI_TYPE_GET_CONTENTS_T,
    MPI_TYPE_GET_CONTENTS_C_T,
    MPI_TYPE_GET_ENVELOPE_T,
    MPI_TYPE_GET_ENVELOPE_C_T,
    MPI_TYPE_GET_EXTENT_T,
    MPI_TYPE_GET_EXTENT_C_T,
    MPI_TYPE_GET_EXTENT_X_T,
    MPI_TYPE_GET_NAME_T,
    MPI_TYPE_GET_TRUE_EXTENT_T,
    MPI_TYPE_GET_TRUE_EXTENT_C_T,
    MPI_TYPE_GET_TRUE_EXTENT_X_T,
    MPI_TYPE_GET_VALUE_INDEX_T,
    MPI_TYPE_INDEXED_T,
    MPI_TYPE_INDEXED_C_T,
    MPI_TYPE_MATCH_SIZE_T,
    MPI_TYPE_SET_NAME_T,
    MPI_TYPE_SIZE_T,
    MPI_TYPE_SIZE_C_T,
    MPI_TYPE_SIZE_X_T,
    MPI_TYPE_VECTOR_T,
    MPI_TYPE_VECTOR_C_T,
    MPI_UNPACK_T,
    MPI_UNPACK_C_T,
    MPI_UNPACK_EXTERNAL_T,
    MPI_UNPACK_EXTERNAL_C_T,
    MPI_TYPE_EXTENT_T,
    MPI_TYPE_LB_T,
    MPI_TYPE_UB_T,
    MPIX_TYPE_IOV_LEN_T,
    MPIX_TYPE_IOV_T,
    MPI_ADD_ERROR_CLASS_T,
    MPI_ADD_ERROR_CODE_T,
    MPI_ADD_ERROR_STRING_T,
    MPI_COMM_CALL_ERRHANDLER_T,
    MPI_COMM_CREATE_ERRHANDLER_T,
    MPI_ERRHANDLER_CREATE_T,
    MPI_COMM_GET_ERRHANDLER_T,
    MPI_ERRHANDLER_GET_T,
    MPI_COMM_SET_ERRHANDLER_T,
    MPI_ERRHANDLER_SET_T,
    MPI_ERRHANDLER_FREE_T,
    MPI_ERROR_CLASS_T,
    MPI_ERROR_STRING_T,
    MPI_FILE_CALL_ERRHANDLER_T,
    MPI_FILE_CREATE_ERRHANDLER_T,
    MPI_FILE_GET_ERRHANDLER_T,
    MPI_FILE_SET_ERRHANDLER_T,
    MPI_REMOVE_ERROR_CLASS_T,
    MPI_REMOVE_ERROR_CODE_T,
    MPI_REMOVE_ERROR_STRING_T,
    MPI_SESSION_CALL_ERRHANDLER_T,
    MPI_SESSION_CREATE_ERRHANDLER_T,
    MPI_SESSION_GET_ERRHANDLER_T,
    MPI_SESSION_SET_ERRHANDLER_T,
    MPI_WIN_CALL_ERRHANDLER_T,
    MPI_WIN_CREATE_ERRHANDLER_T,
    MPI_WIN_GET_ERRHANDLER_T,
    MPI_WIN_SET_ERRHANDLER_T,
    MPI_COMM_C2F_T,
    MPI_COMM_F2C_T,
    MPI_ERRHANDLER_C2F_T,
    MPI_ERRHANDLER_F2C_T,
    MPI_GROUP_C2F_T,
    MPI_GROUP_F2C_T,
    MPI_INFO_C2F_T,
    MPI_INFO_F2C_T,
    MPI_MESSAGE_C2F_T,
    MPI_MESSAGE_F2C_T,
    MPI_OP_C2F_T,
    MPI_OP_F2C_T,
    MPI_REQUEST_C2F_T,
    MPI_REQUEST_F2C_T,
    MPI_SESSION_C2F_T,
    MPI_SESSION_F2C_T,
    MPI_STATUS_C2F_T,
    MPI_STATUS_C2F08_T,
    MPI_STATUS_F082C_T,
    MPI_STATUS_F082F_T,
    MPI_STATUS_F2C_T,
    MPI_STATUS_F2F08_T,
    MPI_TYPE_C2F_T,
    MPI_TYPE_F2C_T,
    MPI_WIN_C2F_T,
    MPI_WIN_F2C_T,
    MPI_GROUP_COMPARE_T,
    MPI_GROUP_DIFFERENCE_T,
    MPI_GROUP_EXCL_T,
    MPI_GROUP_FREE_T,
    MPI_GROUP_INCL_T,
    MPI_GROUP_INTERSECTION_T,
    MPI_GROUP_RANGE_EXCL_T,
    MPI_GROUP_RANGE_INCL_T,
    MPI_GROUP_RANK_T,
    MPI_GROUP_SIZE_T,
    MPI_GROUP_TRANSLATE_RANKS_T,
    MPI_GROUP_UNION_T,
    MPI_INFO_CREATE_T,
    MPI_INFO_CREATE_ENV_T,
    MPI_INFO_DELETE_T,
    MPI_INFO_DUP_T,
    MPI_INFO_FREE_T,
    MPI_INFO_GET_T,
    MPI_INFO_GET_NKEYS_T,
    MPI_INFO_GET_NTHKEY_T,
    MPI_INFO_GET_STRING_T,
    MPI_INFO_GET_VALUELEN_T,
    MPI_INFO_SET_T,
    MPIX_INFO_SET_HEX_T,
    MPI_ABORT_T,
    MPI_COMM_CREATE_FROM_GROUP_T,
    MPI_FINALIZE_T,
    MPI_FINALIZED_T,
    MPI_GROUP_FROM_SESSION_PSET_T,
    MPI_INIT_T,
    MPI_INIT_THREAD_T,
    MPI_INITIALIZED_T,
    MPI_IS_THREAD_MAIN_T,
    MPI_QUERY_THREAD_T,
    MPI_SESSION_FINALIZE_T,
    MPI_SESSION_GET_INFO_T,
    MPI_SESSION_GET_NTH_PSET_T,
    MPI_SESSION_GET_NUM_PSETS_T,
    MPI_SESSION_GET_PSET_INFO_T,
    MPI_SESSION_INIT_T,
    MPI_AINT_ADD_T,
    MPI_AINT_DIFF_T,
    MPI_GET_LIBRARY_VERSION_T,
    MPI_GET_PROCESSOR_NAME_T,
    MPI_GET_VERSION_T,
    MPI_PCONTROL_T,
    MPIX_GPU_QUERY_SUPPORT_T,
    MPIX_QUERY_CUDA_SUPPORT_T,
    MPIX_QUERY_ZE_SUPPORT_T,
    MPIX_QUERY_HIP_SUPPORT_T,
    MPI_T_CATEGORY_CHANGED_T,
    MPI_T_CATEGORY_GET_CATEGORIES_T,
    MPI_T_CATEGORY_GET_CVARS_T,
    MPI_T_CATEGORY_GET_EVENTS_T,
    MPI_T_CATEGORY_GET_INDEX_T,
    MPI_T_CATEGORY_GET_INFO_T,
    MPI_T_CATEGORY_GET_NUM_T,
    MPI_T_CATEGORY_GET_NUM_EVENTS_T,
    MPI_T_CATEGORY_GET_PVARS_T,
    MPI_T_CVAR_GET_INDEX_T,
    MPI_T_CVAR_GET_INFO_T,
    MPI_T_CVAR_GET_NUM_T,
    MPI_T_CVAR_HANDLE_ALLOC_T,
    MPI_T_CVAR_HANDLE_FREE_T,
    MPI_T_CVAR_READ_T,
    MPI_T_CVAR_WRITE_T,
    MPI_T_ENUM_GET_INFO_T,
    MPI_T_ENUM_GET_ITEM_T,
    MPI_T_EVENT_CALLBACK_GET_INFO_T,
    MPI_T_EVENT_CALLBACK_SET_INFO_T,
    MPI_T_EVENT_COPY_T,
    MPI_T_EVENT_GET_INDEX_T,
    MPI_T_EVENT_GET_INFO_T,
    MPI_T_EVENT_GET_NUM_T,
    MPI_T_EVENT_GET_SOURCE_T,
    MPI_T_EVENT_GET_TIMESTAMP_T,
    MPI_T_EVENT_HANDLE_ALLOC_T,
    MPI_T_EVENT_HANDLE_FREE_T,
    MPI_T_EVENT_HANDLE_GET_INFO_T,
    MPI_T_EVENT_HANDLE_SET_INFO_T,
    MPI_T_EVENT_READ_T,
    MPI_T_EVENT_REGISTER_CALLBACK_T,
    MPI_T_EVENT_SET_DROPPED_HANDLER_T,
    MPI_T_FINALIZE_T,
    MPI_T_INIT_THREAD_T,
    MPI_T_PVAR_GET_INDEX_T,
    MPI_T_PVAR_GET_INFO_T,
    MPI_T_PVAR_GET_NUM_T,
    MPI_T_PVAR_HANDLE_ALLOC_T,
    MPI_T_PVAR_HANDLE_FREE_T,
    MPI_T_PVAR_READ_T,
    MPI_T_PVAR_READRESET_T,
    MPI_T_PVAR_RESET_T,
    MPI_T_PVAR_SESSION_CREATE_T,
    MPI_T_PVAR_SESSION_FREE_T,
    MPI_T_PVAR_START_T,
    MPI_T_PVAR_STOP_T,
    MPI_T_PVAR_WRITE_T,
    MPI_T_SOURCE_GET_INFO_T,
    MPI_T_SOURCE_GET_NUM_T,
    MPI_T_SOURCE_GET_TIMESTAMP_T,
    MPI_OP_COMMUTATIVE_T,
    MPI_OP_CREATE_T,
    MPI_OP_CREATE_C_T,
    MPI_OP_FREE_T,
    MPI_PARRIVED_T,
    MPI_PREADY_T,
    MPI_PREADY_LIST_T,
    MPI_PREADY_RANGE_T,
    MPI_PRECV_INIT_T,
    MPI_PSEND_INIT_T,
    MPI_BSEND_T,
    MPI_BSEND_C_T,
    MPI_BSEND_INIT_T,
    MPI_BSEND_INIT_C_T,
    MPI_BUFFER_ATTACH_T,
    MPI_BUFFER_ATTACH_C_T,
    MPI_BUFFER_DETACH_T,
    MPI_BUFFER_DETACH_C_T,
    MPI_BUFFER_FLUSH_T,
    MPI_BUFFER_IFLUSH_T,
    MPI_COMM_ATTACH_BUFFER_T,
    MPI_COMM_ATTACH_BUFFER_C_T,
    MPI_COMM_DETACH_BUFFER_T,
    MPI_COMM_DETACH_BUFFER_C_T,
    MPI_COMM_FLUSH_BUFFER_T,
    MPI_COMM_IFLUSH_BUFFER_T,
    MPI_IBSEND_T,
    MPI_IBSEND_C_T,
    MPI_IMPROBE_T,
    MPI_IMRECV_T,
    MPI_IMRECV_C_T,
    MPI_IPROBE_T,
    MPI_IRECV_T,
    MPI_IRECV_C_T,
    MPI_IRSEND_T,
    MPI_IRSEND_C_T,
    MPI_ISEND_T,
    MPI_ISEND_C_T,
    MPI_ISENDRECV_T,
    MPI_ISENDRECV_C_T,
    MPI_ISENDRECV_REPLACE_T,
    MPI_ISENDRECV_REPLACE_C_T,
    MPI_ISSEND_T,
    MPI_ISSEND_C_T,
    MPI_MPROBE_T,
    MPI_MRECV_T,
    MPI_MRECV_C_T,
    MPI_PROBE_T,
    MPI_RECV_T,
    MPI_RECV_C_T,
    MPI_RECV_INIT_T,
    MPI_RECV_INIT_C_T,
    MPI_RSEND_T,
    MPI_RSEND_C_T,
    MPI_RSEND_INIT_T,
    MPI_RSEND_INIT_C_T,
    MPI_SEND_T,
    MPI_SEND_C_T,
    MPI_SEND_INIT_T,
    MPI_SEND_INIT_C_T,
    MPI_SENDRECV_T,
    MPI_SENDRECV_C_T,
    MPI_SENDRECV_REPLACE_T,
    MPI_SENDRECV_REPLACE_C_T,
    MPI_SESSION_ATTACH_BUFFER_T,
    MPI_SESSION_ATTACH_BUFFER_C_T,
    MPI_SESSION_DETACH_BUFFER_T,
    MPI_SESSION_DETACH_BUFFER_C_T,
    MPI_SESSION_FLUSH_BUFFER_T,
    MPI_SESSION_IFLUSH_BUFFER_T,
    MPI_SSEND_T,
    MPI_SSEND_C_T,
    MPI_SSEND_INIT_T,
    MPI_SSEND_INIT_C_T,
    MPI_CANCEL_T,
    MPI_GREQUEST_COMPLETE_T,
    MPI_GREQUEST_START_T,
    MPI_REQUEST_FREE_T,
    MPI_REQUEST_GET_STATUS_T,
    MPI_REQUEST_GET_STATUS_ALL_T,
    MPI_REQUEST_GET_STATUS_ANY_T,
    MPI_REQUEST_GET_STATUS_SOME_T,
    MPI_START_T,
    MPI_STARTALL_T,
    MPI_STATUS_GET_ERROR_T,
    MPI_STATUS_GET_SOURCE_T,
    MPI_STATUS_GET_TAG_T,
    MPI_STATUS_SET_ERROR_T,
    MPI_STATUS_SET_SOURCE_T,
    MPI_STATUS_SET_TAG_T,
    MPI_STATUS_SET_CANCELLED_T,
    MPI_TEST_T,
    MPI_TEST_CANCELLED_T,
    MPI_TESTALL_T,
    MPI_TESTANY_T,
    MPI_TESTSOME_T,
    MPI_WAIT_T,
    MPI_WAITALL_T,
    MPI_WAITANY_T,
    MPI_WAITSOME_T,
    MPIX_GREQUEST_START_T,
    MPIX_GREQUEST_CLASS_CREATE_T,
    MPIX_GREQUEST_CLASS_ALLOCATE_T,
    MPIX_REQUEST_IS_COMPLETE_T,
    MPI_ACCUMULATE_T,
    MPI_ACCUMULATE_C_T,
    MPI_ALLOC_MEM_T,
    MPI_COMPARE_AND_SWAP_T,
    MPI_FETCH_AND_OP_T,
    MPI_FREE_MEM_T,
    MPI_GET_T,
    MPI_GET_C_T,
    MPI_GET_ACCUMULATE_T,
    MPI_GET_ACCUMULATE_C_T,
    MPI_PUT_T,
    MPI_PUT_C_T,
    MPI_RACCUMULATE_T,
    MPI_RACCUMULATE_C_T,
    MPI_RGET_T,
    MPI_RGET_C_T,
    MPI_RGET_ACCUMULATE_T,
    MPI_RGET_ACCUMULATE_C_T,
    MPI_RPUT_T,
    MPI_RPUT_C_T,
    MPI_WIN_ALLOCATE_T,
    MPI_WIN_ALLOCATE_C_T,
    MPI_WIN_ALLOCATE_SHARED_T,
    MPI_WIN_ALLOCATE_SHARED_C_T,
    MPI_WIN_ATTACH_T,
    MPI_WIN_COMPLETE_T,
    MPI_WIN_CREATE_T,
    MPI_WIN_CREATE_C_T,
    MPI_WIN_CREATE_DYNAMIC_T,
    MPI_WIN_DETACH_T,
    MPI_WIN_FENCE_T,
    MPI_WIN_FLUSH_T,
    MPI_WIN_FLUSH_ALL_T,
    MPI_WIN_FLUSH_LOCAL_T,
    MPI_WIN_FLUSH_LOCAL_ALL_T,
    MPI_WIN_FREE_T,
    MPI_WIN_GET_GROUP_T,
    MPI_WIN_GET_INFO_T,
    MPI_WIN_GET_NAME_T,
    MPI_WIN_LOCK_T,
    MPI_WIN_LOCK_ALL_T,
    MPI_WIN_POST_T,
    MPI_WIN_SET_INFO_T,
    MPI_WIN_SET_NAME_T,
    MPI_WIN_SHARED_QUERY_T,
    MPI_WIN_SHARED_QUERY_C_T,
    MPI_WIN_START_T,
    MPI_WIN_SYNC_T,
    MPI_WIN_TEST_T,
    MPI_WIN_UNLOCK_T,
    MPI_WIN_UNLOCK_ALL_T,
    MPI_WIN_WAIT_T,
    MPI_CLOSE_PORT_T,
    MPI_COMM_ACCEPT_T,
    MPI_COMM_CONNECT_T,
    MPI_COMM_DISCONNECT_T,
    MPI_COMM_GET_PARENT_T,
    MPI_COMM_JOIN_T,
    MPI_COMM_SPAWN_T,
    MPI_COMM_SPAWN_MULTIPLE_T,
    MPI_LOOKUP_NAME_T,
    MPI_OPEN_PORT_T,
    MPI_PUBLISH_NAME_T,
    MPI_UNPUBLISH_NAME_T,
    MPIX_STREAM_CREATE_T,
    MPIX_STREAM_FREE_T,
    MPIX_STREAM_COMM_CREATE_T,
    MPIX_STREAM_COMM_CREATE_MULTIPLEX_T,
    MPIX_COMM_GET_STREAM_T,
    MPIX_STREAM_PROGRESS_T,
    MPIX_START_PROGRESS_THREAD_T,
    MPIX_STOP_PROGRESS_THREAD_T,
    MPIX_STREAM_SEND_T,
    MPIX_STREAM_SEND_C_T,
    MPIX_STREAM_ISEND_T,
    MPIX_STREAM_ISEND_C_T,
    MPIX_STREAM_RECV_T,
    MPIX_STREAM_RECV_C_T,
    MPIX_STREAM_IRECV_T,
    MPIX_STREAM_IRECV_C_T,
    MPIX_SEND_ENQUEUE_T,
    MPIX_SEND_ENQUEUE_C_T,
    MPIX_RECV_ENQUEUE_T,
    MPIX_RECV_ENQUEUE_C_T,
    MPIX_ISEND_ENQUEUE_T,
    MPIX_ISEND_ENQUEUE_C_T,
    MPIX_IRECV_ENQUEUE_T,
    MPIX_IRECV_ENQUEUE_C_T,
    MPIX_WAIT_ENQUEUE_T,
    MPIX_WAITALL_ENQUEUE_T,
    MPIX_ALLREDUCE_ENQUEUE_T,
    MPIX_ALLREDUCE_ENQUEUE_C_T,
    MPIX_THREADCOMM_INIT_T,
    MPIX_THREADCOMM_FREE_T,
    MPIX_THREADCOMM_START_T,
    MPIX_THREADCOMM_FINISH_T,
    MPI_WTICK_T,
    MPI_WTIME_T,
    MPI_CART_COORDS_T,
    MPI_CART_CREATE_T,
    MPI_CART_GET_T,
    MPI_CART_MAP_T,
    MPI_CART_RANK_T,
    MPI_CART_SHIFT_T,
    MPI_CART_SUB_T,
    MPI_CARTDIM_GET_T,
    MPI_DIMS_CREATE_T,
    MPI_DIST_GRAPH_CREATE_T,
    MPI_DIST_GRAPH_CREATE_ADJACENT_T,
    MPI_DIST_GRAPH_NEIGHBORS_T,
    MPI_DIST_GRAPH_NEIGHBORS_COUNT_T,
    MPI_GET_HW_RESOURCE_INFO_T,
    MPI_GRAPH_CREATE_T,
    MPI_GRAPH_GET_T,
    MPI_GRAPH_MAP_T,
    MPI_GRAPH_NEIGHBORS_T,
    MPI_GRAPH_NEIGHBORS_COUNT_T,
    MPI_GRAPHDIMS_GET_T,
    MPI_TOPO_TEST_T,
    MPI_FILE_C2F_T,
    MPI_FILE_CLOSE_T,
    MPI_FILE_DELETE_T,
    MPI_FILE_F2C_T,
    MPI_FILE_GET_AMODE_T,
    MPI_FILE_GET_ATOMICITY_T,
    MPI_FILE_GET_BYTE_OFFSET_T,
    MPI_FILE_GET_GROUP_T,
    MPI_FILE_GET_INFO_T,
    MPI_FILE_GET_POSITION_T,
    MPI_FILE_GET_POSITION_SHARED_T,
    MPI_FILE_GET_SIZE_T,
    MPI_FILE_GET_TYPE_EXTENT_T,
    MPI_FILE_GET_TYPE_EXTENT_C_T,
    MPI_FILE_GET_VIEW_T,
    MPI_FILE_IREAD_T,
    MPI_FILE_IREAD_C_T,
    MPI_FILE_IREAD_ALL_T,
    MPI_FILE_IREAD_ALL_C_T,
    MPI_FILE_IREAD_AT_T,
    MPI_FILE_IREAD_AT_C_T,
    MPI_FILE_IREAD_AT_ALL_T,
    MPI_FILE_IREAD_AT_ALL_C_T,
    MPI_FILE_IREAD_SHARED_T,
    MPI_FILE_IREAD_SHARED_C_T,
    MPI_FILE_IWRITE_T,
    MPI_FILE_IWRITE_C_T,
    MPI_FILE_IWRITE_ALL_T,
    MPI_FILE_IWRITE_ALL_C_T,
    MPI_FILE_IWRITE_AT_T,
    MPI_FILE_IWRITE_AT_C_T,
    MPI_FILE_IWRITE_AT_ALL_T,
    MPI_FILE_IWRITE_AT_ALL_C_T,
    MPI_FILE_IWRITE_SHARED_T,
    MPI_FILE_IWRITE_SHARED_C_T,
    MPI_FILE_OPEN_T,
    MPI_FILE_PREALLOCATE_T,
    MPI_FILE_READ_T,
    MPI_FILE_READ_C_T,
    MPI_FILE_READ_ALL_T,
    MPI_FILE_READ_ALL_C_T,
    MPI_FILE_READ_ALL_BEGIN_T,
    MPI_FILE_READ_ALL_BEGIN_C_T,
    MPI_FILE_READ_ALL_END_T,
    MPI_FILE_READ_AT_T,
    MPI_FILE_READ_AT_C_T,
    MPI_FILE_READ_AT_ALL_T,
    MPI_FILE_READ_AT_ALL_C_T,
    MPI_FILE_READ_AT_ALL_BEGIN_T,
    MPI_FILE_READ_AT_ALL_BEGIN_C_T,
    MPI_FILE_READ_AT_ALL_END_T,
    MPI_FILE_READ_ORDERED_T,
    MPI_FILE_READ_ORDERED_C_T,
    MPI_FILE_READ_ORDERED_BEGIN_T,
    MPI_FILE_READ_ORDERED_BEGIN_C_T,
    MPI_FILE_READ_ORDERED_END_T,
    MPI_FILE_READ_SHARED_T,
    MPI_FILE_READ_SHARED_C_T,
    MPI_FILE_SEEK_T,
    MPI_FILE_SEEK_SHARED_T,
    MPI_FILE_SET_ATOMICITY_T,
    MPI_FILE_SET_INFO_T,
    MPI_FILE_SET_SIZE_T,
    MPI_FILE_SET_VIEW_T,
    MPI_FILE_SYNC_T,
    MPI_FILE_WRITE_T,
    MPI_FILE_WRITE_C_T,
    MPI_FILE_WRITE_ALL_T,
    MPI_FILE_WRITE_ALL_C_T,
    MPI_FILE_WRITE_ALL_BEGIN_T,
    MPI_FILE_WRITE_ALL_BEGIN_C_T,
    MPI_FILE_WRITE_ALL_END_T,
    MPI_FILE_WRITE_AT_T,
    MPI_FILE_WRITE_AT_C_T,
    MPI_FILE_WRITE_AT_ALL_T,
    MPI_FILE_WRITE_AT_ALL_C_T,
    MPI_FILE_WRITE_AT_ALL_BEGIN_T,
    MPI_FILE_WRITE_AT_ALL_BEGIN_C_T,
    MPI_FILE_WRITE_AT_ALL_END_T,
    MPI_FILE_WRITE_ORDERED_T,
    MPI_FILE_WRITE_ORDERED_C_T,
    MPI_FILE_WRITE_ORDERED_BEGIN_T,
    MPI_FILE_WRITE_ORDERED_BEGIN_C_T,
    MPI_FILE_WRITE_ORDERED_END_T,
    MPI_FILE_WRITE_SHARED_T,
    MPI_FILE_WRITE_SHARED_C_T,
    MPI_REGISTER_DATAREP_T,
    MPI_REGISTER_DATAREP_C_T,
    MPI_FILE_TOINT_T,
    MPI_FILE_FROMINT_T,
    MPI_LAST_FUNC_T
};
int QMPI_Abi_get_fortran_booleans(QMPI_Context context, int tool_id, int logical_size,
                                  void *logical_true, void *logical_false, int *is_set)
                                                  ;
int QMPI_Abi_get_fortran_info(QMPI_Context context, int tool_id, MPI_Info *info) ;
int QMPI_Abi_get_info(QMPI_Context context, int tool_id, MPI_Info *info) ;
int QMPI_Abi_get_version(QMPI_Context context, int tool_id, int *abi_major, int *abi_minor)
                    ;
int QMPI_Abi_set_fortran_booleans(QMPI_Context context, int tool_id, int logical_size,
                                  void *logical_true, void *logical_false) ;
int QMPI_Abi_set_fortran_info(QMPI_Context context, int tool_id, MPI_Info info) ;
int QMPI_Comm_toint(QMPI_Context context, int tool_id, MPI_Comm comm) ;
MPI_Comm QMPI_Comm_fromint(QMPI_Context context, int tool_id, int comm) ;
int QMPI_Errhandler_toint(QMPI_Context context, int tool_id, MPI_Errhandler errhandler)
                    ;
MPI_Errhandler QMPI_Errhandler_fromint(QMPI_Context context, int tool_id, int errhandler)
                    ;
int QMPI_Group_toint(QMPI_Context context, int tool_id, MPI_Group group) ;
MPI_Group QMPI_Group_fromint(QMPI_Context context, int tool_id, int group) ;
int QMPI_Info_toint(QMPI_Context context, int tool_id, MPI_Info info) ;
MPI_Info QMPI_Info_fromint(QMPI_Context context, int tool_id, int info) ;
int QMPI_Message_toint(QMPI_Context context, int tool_id, MPI_Message message) ;
MPI_Message QMPI_Message_fromint(QMPI_Context context, int tool_id, int message) ;
int QMPI_Op_toint(QMPI_Context context, int tool_id, MPI_Op op) ;
MPI_Op QMPI_Op_fromint(QMPI_Context context, int tool_id, int op) ;
int QMPI_Request_toint(QMPI_Context context, int tool_id, MPI_Request request) ;
MPI_Request QMPI_Request_fromint(QMPI_Context context, int tool_id, int request) ;
int QMPI_Session_toint(QMPI_Context context, int tool_id, MPI_Session session) ;
MPI_Session QMPI_Session_fromint(QMPI_Context context, int tool_id, int session) ;
int QMPI_Type_toint(QMPI_Context context, int tool_id, MPI_Datatype datatype) ;
MPI_Datatype QMPI_Type_fromint(QMPI_Context context, int tool_id, int datatype) ;
int QMPI_Win_toint(QMPI_Context context, int tool_id, MPI_Win win) ;
MPI_Win QMPI_Win_fromint(QMPI_Context context, int tool_id, int win) ;
int QMPIX_Async_start(QMPI_Context context, int tool_id, MPIX_Async_poll_function *poll_fn,
                      void *extra_state, MPIX_Stream stream) ;
void * QMPIX_Async_get_state(QMPI_Context context, int tool_id, MPIX_Async_thing async_thing)
                    ;
int QMPIX_Async_spawn(QMPI_Context context, int tool_id, MPIX_Async_thing async_thing,
                      MPIX_Async_poll_function *poll_fn, void *extra_state, MPIX_Stream stream)
                                      ;
int QMPI_Comm_create_keyval(QMPI_Context context, int tool_id,
                            MPI_Comm_copy_attr_function *comm_copy_attr_fn,
                            MPI_Comm_delete_attr_function *comm_delete_attr_fn, int *comm_keyval,
                            void *extra_state) ;
int QMPI_Keyval_create(QMPI_Context context, int tool_id, MPI_Copy_function *copy_fn,
                       MPI_Delete_function *delete_fn, int *keyval, void *extra_state)
                                       ;
int QMPI_Comm_delete_attr(QMPI_Context context, int tool_id, MPI_Comm comm, int comm_keyval)
                    ;
int QMPI_Attr_delete(QMPI_Context context, int tool_id, MPI_Comm comm, int keyval)
                    ;
int QMPI_Comm_free_keyval(QMPI_Context context, int tool_id, int *comm_keyval) ;
int QMPI_Keyval_free(QMPI_Context context, int tool_id, int *keyval) ;
int QMPI_Comm_get_attr(QMPI_Context context, int tool_id, MPI_Comm comm, int comm_keyval,
                       void *attribute_val, int *flag) ;
int QMPI_Attr_get(QMPI_Context context, int tool_id, MPI_Comm comm, int keyval, void *attribute_val,
                  int *flag) ;
int QMPI_Comm_set_attr(QMPI_Context context, int tool_id, MPI_Comm comm, int comm_keyval,
                       void *attribute_val) ;
int QMPI_Attr_put(QMPI_Context context, int tool_id, MPI_Comm comm, int keyval,
                  void *attribute_val) ;
int QMPI_Type_create_keyval(QMPI_Context context, int tool_id,
                            MPI_Type_copy_attr_function *type_copy_attr_fn,
                            MPI_Type_delete_attr_function *type_delete_attr_fn, int *type_keyval,
                            void *extra_state) ;
int QMPI_Type_delete_attr(QMPI_Context context, int tool_id, MPI_Datatype datatype,
                          int type_keyval) ;
int QMPI_Type_free_keyval(QMPI_Context context, int tool_id, int *type_keyval) ;
int QMPI_Type_get_attr(QMPI_Context context, int tool_id, MPI_Datatype datatype, int type_keyval,
                       void *attribute_val, int *flag) ;
int QMPI_Type_set_attr(QMPI_Context context, int tool_id, MPI_Datatype datatype, int type_keyval,
                       void *attribute_val) ;
int QMPI_Win_create_keyval(QMPI_Context context, int tool_id,
                           MPI_Win_copy_attr_function *win_copy_attr_fn,
                           MPI_Win_delete_attr_function *win_delete_attr_fn, int *win_keyval,
                           void *extra_state) ;
int QMPI_Win_delete_attr(QMPI_Context context, int tool_id, MPI_Win win, int win_keyval)
                    ;
int QMPI_Win_free_keyval(QMPI_Context context, int tool_id, int *win_keyval) ;
int QMPI_Win_get_attr(QMPI_Context context, int tool_id, MPI_Win win, int win_keyval,
                      void *attribute_val, int *flag) ;
int QMPI_Win_set_attr(QMPI_Context context, int tool_id, MPI_Win win, int win_keyval,
                      void *attribute_val) ;
int QMPIX_Op_create_x(QMPI_Context context, int tool_id, MPIX_User_function_x *user_fn_x,
                      MPIX_Destructor_function *destructor_fn, int commute, void *extra_state,
                      MPI_Op *op) ;
int QMPIX_Comm_create_errhandler_x(QMPI_Context context, int tool_id,
                                   MPIX_Comm_errhandler_function_x *comm_errhandler_fn_x,
                                   MPIX_Destructor_function *destructor_fn, void *extra_state,
                                   MPI_Errhandler *errhandler) ;
int QMPIX_Win_create_errhandler_x(QMPI_Context context, int tool_id,
                                  MPIX_Win_errhandler_function_x *comm_errhandler_fn_x,
                                  MPIX_Destructor_function *destructor_fn, void *extra_state,
                                  MPI_Errhandler *errhandler) ;
int QMPIX_File_create_errhandler_x(QMPI_Context context, int tool_id,
                                   MPIX_File_errhandler_function_x *comm_errhandler_fn_x,
                                   MPIX_Destructor_function *destructor_fn, void *extra_state,
                                   MPI_Errhandler *errhandler) ;
int QMPIX_Session_create_errhandler_x(QMPI_Context context, int tool_id,
                                      MPIX_Session_errhandler_function_x *comm_errhandler_fn_x,
                                      MPIX_Destructor_function *destructor_fn, void *extra_state,
                                      MPI_Errhandler *errhandler) ;
int QMPI_Allgather(QMPI_Context context, int tool_id, const void *sendbuf, int sendcount,
                   MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype,
                   MPI_Comm comm)
                                                                                                               ;
int QMPI_Allgather_c(QMPI_Context context, int tool_id, const void *sendbuf, MPI_Count sendcount,
                     MPI_Datatype sendtype, void *recvbuf, MPI_Count recvcount,
                     MPI_Datatype recvtype, MPI_Comm comm)
                                                                                                                 ;
int QMPI_Allgather_init(QMPI_Context context, int tool_id, const void *sendbuf, int sendcount,
                        MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype,
                        MPI_Comm comm, MPI_Info info, MPI_Request *request)
                                                                                                                    ;
int QMPI_Allgather_init_c(QMPI_Context context, int tool_id, const void *sendbuf,
                          MPI_Count sendcount, MPI_Datatype sendtype, void *recvbuf,
                          MPI_Count recvcount, MPI_Datatype recvtype, MPI_Comm comm, MPI_Info info,
                          MPI_Request *request)
                                                                                                                      ;
int QMPI_Allgatherv(QMPI_Context context, int tool_id, const void *sendbuf, int sendcount,
                    MPI_Datatype sendtype, void *recvbuf, const int recvcounts[],
                    const int displs[], MPI_Datatype recvtype, MPI_Comm comm)
                                                                                                                ;
int QMPI_Allgatherv_c(QMPI_Context context, int tool_id, const void *sendbuf, MPI_Count sendcount,
                      MPI_Datatype sendtype, void *recvbuf, const MPI_Count recvcounts[],
                      const MPI_Aint displs[], MPI_Datatype recvtype, MPI_Comm comm)
                                                                                                                  ;
int QMPI_Allgatherv_init(QMPI_Context context, int tool_id, const void *sendbuf, int sendcount,
                         MPI_Datatype sendtype, void *recvbuf, const int recvcounts[],
                         const int displs[], MPI_Datatype recvtype, MPI_Comm comm, MPI_Info info,
                         MPI_Request *request)
                                                                                                                     ;
int QMPI_Allgatherv_init_c(QMPI_Context context, int tool_id, const void *sendbuf,
                           MPI_Count sendcount, MPI_Datatype sendtype, void *recvbuf,
                           const MPI_Count recvcounts[], const MPI_Aint displs[],
                           MPI_Datatype recvtype, MPI_Comm comm, MPI_Info info,
                           MPI_Request *request)
                                                                                                                       ;
int QMPI_Allreduce(QMPI_Context context, int tool_id, const void *sendbuf, void *recvbuf, int count,
                   MPI_Datatype datatype, MPI_Op op, MPI_Comm comm)
                                                                                                               ;
int QMPI_Allreduce_c(QMPI_Context context, int tool_id, const void *sendbuf, void *recvbuf,
                     MPI_Count count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm)
                                                                                                                 ;
int QMPI_Allreduce_init(QMPI_Context context, int tool_id, const void *sendbuf, void *recvbuf,
                        int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm, MPI_Info info,
                        MPI_Request *request)
                                                                                                                    ;
int QMPI_Allreduce_init_c(QMPI_Context context, int tool_id, const void *sendbuf, void *recvbuf,
                          MPI_Count count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm,
                          MPI_Info info, MPI_Request *request)
                                                                                                                      ;
int QMPI_Alltoall(QMPI_Context context, int tool_id, const void *sendbuf, int sendcount,
                  MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype,
                  MPI_Comm comm)
                                                                                                              ;
int QMPI_Alltoall_c(QMPI_Context context, int tool_id, const void *sendbuf, MPI_Count sendcount,
                    MPI_Datatype sendtype, void *recvbuf, MPI_Count recvcount,
                    MPI_Datatype recvtype, MPI_Comm comm)
                                                                                                                ;
int QMPI_Alltoall_init(QMPI_Context context, int tool_id, const void *sendbuf, int sendcount,
                       MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype,
                       MPI_Comm comm, MPI_Info info, MPI_Request *request)
                                                                                                                   ;
int QMPI_Alltoall_init_c(QMPI_Context context, int tool_id, const void *sendbuf,
                         MPI_Count sendcount, MPI_Datatype sendtype, void *recvbuf,
                         MPI_Count recvcount, MPI_Datatype recvtype, MPI_Comm comm, MPI_Info info,
                         MPI_Request *request)
                                                                                                                     ;
int QMPI_Alltoallv(QMPI_Context context, int tool_id, const void *sendbuf, const int sendcounts[],
                   const int sdispls[], MPI_Datatype sendtype, void *recvbuf,
                   const int recvcounts[], const int rdispls[], MPI_Datatype recvtype,
                   MPI_Comm comm)
                                                                                                                ;
int QMPI_Alltoallv_c(QMPI_Context context, int tool_id, const void *sendbuf,
                     const MPI_Count sendcounts[], const MPI_Aint sdispls[], MPI_Datatype sendtype,
                     void *recvbuf, const MPI_Count recvcounts[], const MPI_Aint rdispls[],
                     MPI_Datatype recvtype, MPI_Comm comm)
                                                                                                                  ;
int QMPI_Alltoallv_init(QMPI_Context context, int tool_id, const void *sendbuf,
                        const int sendcounts[], const int sdispls[], MPI_Datatype sendtype,
                        void *recvbuf, const int recvcounts[], const int rdispls[],
                        MPI_Datatype recvtype, MPI_Comm comm, MPI_Info info, MPI_Request *request)
                                                                                                                     ;
int QMPI_Alltoallv_init_c(QMPI_Context context, int tool_id, const void *sendbuf,
                          const MPI_Count sendcounts[], const MPI_Aint sdispls[],
                          MPI_Datatype sendtype, void *recvbuf, const MPI_Count recvcounts[],
                          const MPI_Aint rdispls[], MPI_Datatype recvtype, MPI_Comm comm,
                          MPI_Info info, MPI_Request *request)
                                                                                                                       ;
int QMPI_Alltoallw(QMPI_Context context, int tool_id, const void *sendbuf, const int sendcounts[],
                   const int sdispls[], const MPI_Datatype sendtypes[], void *recvbuf,
                   const int recvcounts[], const int rdispls[], const MPI_Datatype recvtypes[],
                   MPI_Comm comm) ;
int QMPI_Alltoallw_c(QMPI_Context context, int tool_id, const void *sendbuf,
                     const MPI_Count sendcounts[], const MPI_Aint sdispls[],
                     const MPI_Datatype sendtypes[], void *recvbuf, const MPI_Count recvcounts[],
                     const MPI_Aint rdispls[], const MPI_Datatype recvtypes[], MPI_Comm comm)
                                     ;
int QMPI_Alltoallw_init(QMPI_Context context, int tool_id, const void *sendbuf,
                        const int sendcounts[], const int sdispls[], const MPI_Datatype sendtypes[],
                        void *recvbuf, const int recvcounts[], const int rdispls[],
                        const MPI_Datatype recvtypes[], MPI_Comm comm, MPI_Info info,
                        MPI_Request *request) ;
int QMPI_Alltoallw_init_c(QMPI_Context context, int tool_id, const void *sendbuf,
                          const MPI_Count sendcounts[], const MPI_Aint sdispls[],
                          const MPI_Datatype sendtypes[], void *recvbuf,
                          const MPI_Count recvcounts[], const MPI_Aint rdispls[],
                          const MPI_Datatype recvtypes[], MPI_Comm comm, MPI_Info info,
                          MPI_Request *request) ;
int QMPI_Barrier(QMPI_Context context, int tool_id, MPI_Comm comm) ;
int QMPI_Barrier_init(QMPI_Context context, int tool_id, MPI_Comm comm, MPI_Info info,
                      MPI_Request *request) ;
int QMPI_Bcast(QMPI_Context context, int tool_id, void *buffer, int count, MPI_Datatype datatype,
               int root, MPI_Comm comm) ;
int QMPI_Bcast_c(QMPI_Context context, int tool_id, void *buffer, MPI_Count count,
                 MPI_Datatype datatype, int root, MPI_Comm comm)
                                                                       ;
int QMPI_Bcast_init(QMPI_Context context, int tool_id, void *buffer, int count,
                    MPI_Datatype datatype, int root, MPI_Comm comm, MPI_Info info,
                    MPI_Request *request) ;
int QMPI_Bcast_init_c(QMPI_Context context, int tool_id, void *buffer, MPI_Count count,
                      MPI_Datatype datatype, int root, MPI_Comm comm, MPI_Info info,
                      MPI_Request *request) ;
int QMPI_Exscan(QMPI_Context context, int tool_id, const void *sendbuf, void *recvbuf, int count,
                MPI_Datatype datatype, MPI_Op op, MPI_Comm comm)
                                                                                                            ;
int QMPI_Exscan_c(QMPI_Context context, int tool_id, const void *sendbuf, void *recvbuf,
                  MPI_Count count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm)
                                                                                                              ;
int QMPI_Exscan_init(QMPI_Context context, int tool_id, const void *sendbuf, void *recvbuf,
                     int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm, MPI_Info info,
                     MPI_Request *request)
                                                                                                                 ;
int QMPI_Exscan_init_c(QMPI_Context context, int tool_id, const void *sendbuf, void *recvbuf,
                       MPI_Count count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm,
                       MPI_Info info, MPI_Request *request)
                                                                                                                   ;
int QMPI_Gather(QMPI_Context context, int tool_id, const void *sendbuf, int sendcount,
                MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype,
                int root, MPI_Comm comm)
                                                                                                            ;
int QMPI_Gather_c(QMPI_Context context, int tool_id, const void *sendbuf, MPI_Count sendcount,
                  MPI_Datatype sendtype, void *recvbuf, MPI_Count recvcount, MPI_Datatype recvtype,
                  int root, MPI_Comm comm)
                                                                                                              ;
int QMPI_Gather_init(QMPI_Context context, int tool_id, const void *sendbuf, int sendcount,
                     MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype,
                     int root, MPI_Comm comm, MPI_Info info, MPI_Request *request)
                                                                                                                 ;
int QMPI_Gather_init_c(QMPI_Context context, int tool_id, const void *sendbuf, MPI_Count sendcount,
                       MPI_Datatype sendtype, void *recvbuf, MPI_Count recvcount,
                       MPI_Datatype recvtype, int root, MPI_Comm comm, MPI_Info info,
                       MPI_Request *request)
                                                                                                                   ;
int QMPI_Gatherv(QMPI_Context context, int tool_id, const void *sendbuf, int sendcount,
                 MPI_Datatype sendtype, void *recvbuf, const int recvcounts[], const int displs[],
                 MPI_Datatype recvtype, int root, MPI_Comm comm)
                                                                                                             ;
int QMPI_Gatherv_c(QMPI_Context context, int tool_id, const void *sendbuf, MPI_Count sendcount,
                   MPI_Datatype sendtype, void *recvbuf, const MPI_Count recvcounts[],
                   const MPI_Aint displs[], MPI_Datatype recvtype, int root, MPI_Comm comm)
                                                                                                               ;
int QMPI_Gatherv_init(QMPI_Context context, int tool_id, const void *sendbuf, int sendcount,
                      MPI_Datatype sendtype, void *recvbuf, const int recvcounts[],
                      const int displs[], MPI_Datatype recvtype, int root, MPI_Comm comm,
                      MPI_Info info, MPI_Request *request)
                                                                                                                  ;
int QMPI_Gatherv_init_c(QMPI_Context context, int tool_id, const void *sendbuf, MPI_Count sendcount,
                        MPI_Datatype sendtype, void *recvbuf, const MPI_Count recvcounts[],
                        const MPI_Aint displs[], MPI_Datatype recvtype, int root, MPI_Comm comm,
                        MPI_Info info, MPI_Request *request)
                                                                                                                    ;
int QMPI_Iallgather(QMPI_Context context, int tool_id, const void *sendbuf, int sendcount,
                    MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype,
                    MPI_Comm comm, MPI_Request *request)
                                                                                                                ;
int QMPI_Iallgather_c(QMPI_Context context, int tool_id, const void *sendbuf, MPI_Count sendcount,
                      MPI_Datatype sendtype, void *recvbuf, MPI_Count recvcount,
                      MPI_Datatype recvtype, MPI_Comm comm, MPI_Request *request)
                                                                                                                  ;
int QMPI_Iallgatherv(QMPI_Context context, int tool_id, const void *sendbuf, int sendcount,
                     MPI_Datatype sendtype, void *recvbuf, const int recvcounts[],
                     const int displs[], MPI_Datatype recvtype, MPI_Comm comm,
                     MPI_Request *request)
                                                                                                                 ;
int QMPI_Iallgatherv_c(QMPI_Context context, int tool_id, const void *sendbuf, MPI_Count sendcount,
                       MPI_Datatype sendtype, void *recvbuf, const MPI_Count recvcounts[],
                       const MPI_Aint displs[], MPI_Datatype recvtype, MPI_Comm comm,
                       MPI_Request *request)
                                                                                                                   ;
int QMPI_Iallreduce(QMPI_Context context, int tool_id, const void *sendbuf, void *recvbuf,
                    int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm,
                    MPI_Request *request)
                                                                                                                ;
int QMPI_Iallreduce_c(QMPI_Context context, int tool_id, const void *sendbuf, void *recvbuf,
                      MPI_Count count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm,
                      MPI_Request *request)
                                                                                                                  ;
int QMPI_Ialltoall(QMPI_Context context, int tool_id, const void *sendbuf, int sendcount,
                   MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype,
                   MPI_Comm comm, MPI_Request *request)
                                                                                                               ;
int QMPI_Ialltoall_c(QMPI_Context context, int tool_id, const void *sendbuf, MPI_Count sendcount,
                     MPI_Datatype sendtype, void *recvbuf, MPI_Count recvcount,
                     MPI_Datatype recvtype, MPI_Comm comm, MPI_Request *request)
                                                                                                                 ;
int QMPI_Ialltoallv(QMPI_Context context, int tool_id, const void *sendbuf, const int sendcounts[],
                    const int sdispls[], MPI_Datatype sendtype, void *recvbuf,
                    const int recvcounts[], const int rdispls[], MPI_Datatype recvtype,
                    MPI_Comm comm, MPI_Request *request)
                                                                                                                 ;
int QMPI_Ialltoallv_c(QMPI_Context context, int tool_id, const void *sendbuf,
                      const MPI_Count sendcounts[], const MPI_Aint sdispls[], MPI_Datatype sendtype,
                      void *recvbuf, const MPI_Count recvcounts[], const MPI_Aint rdispls[],
                      MPI_Datatype recvtype, MPI_Comm comm, MPI_Request *request)
                                                                                                                   ;
int QMPI_Ialltoallw(QMPI_Context context, int tool_id, const void *sendbuf, const int sendcounts[],
                    const int sdispls[], const MPI_Datatype sendtypes[], void *recvbuf,
                    const int recvcounts[], const int rdispls[], const MPI_Datatype recvtypes[],
                    MPI_Comm comm, MPI_Request *request) ;
int QMPI_Ialltoallw_c(QMPI_Context context, int tool_id, const void *sendbuf,
                      const MPI_Count sendcounts[], const MPI_Aint sdispls[],
                      const MPI_Datatype sendtypes[], void *recvbuf, const MPI_Count recvcounts[],
                      const MPI_Aint rdispls[], const MPI_Datatype recvtypes[], MPI_Comm comm,
                      MPI_Request *request) ;
int QMPI_Ibarrier(QMPI_Context context, int tool_id, MPI_Comm comm, MPI_Request *request)
                    ;
int QMPI_Ibcast(QMPI_Context context, int tool_id, void *buffer, int count, MPI_Datatype datatype,
                int root, MPI_Comm comm, MPI_Request *request)
                                                                      ;
int QMPI_Ibcast_c(QMPI_Context context, int tool_id, void *buffer, MPI_Count count,
                  MPI_Datatype datatype, int root, MPI_Comm comm, MPI_Request *request)
                                                                        ;
int QMPI_Iexscan(QMPI_Context context, int tool_id, const void *sendbuf, void *recvbuf, int count,
                 MPI_Datatype datatype, MPI_Op op, MPI_Comm comm, MPI_Request *request)
                                                                                                             ;
int QMPI_Iexscan_c(QMPI_Context context, int tool_id, const void *sendbuf, void *recvbuf,
                   MPI_Count count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm,
                   MPI_Request *request)
                                                                                                               ;
int QMPI_Igather(QMPI_Context context, int tool_id, const void *sendbuf, int sendcount,
                 MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype,
                 int root, MPI_Comm comm, MPI_Request *request)
                                                                                                             ;
int QMPI_Igather_c(QMPI_Context context, int tool_id, const void *sendbuf, MPI_Count sendcount,
                   MPI_Datatype sendtype, void *recvbuf, MPI_Count recvcount, MPI_Datatype recvtype,
                   int root, MPI_Comm comm, MPI_Request *request)
                                                                                                               ;
int QMPI_Igatherv(QMPI_Context context, int tool_id, const void *sendbuf, int sendcount,
                  MPI_Datatype sendtype, void *recvbuf, const int recvcounts[], const int displs[],
                  MPI_Datatype recvtype, int root, MPI_Comm comm, MPI_Request *request)
                                                                                                              ;
int QMPI_Igatherv_c(QMPI_Context context, int tool_id, const void *sendbuf, MPI_Count sendcount,
                    MPI_Datatype sendtype, void *recvbuf, const MPI_Count recvcounts[],
                    const MPI_Aint displs[], MPI_Datatype recvtype, int root, MPI_Comm comm,
                    MPI_Request *request)
                                                                                                                ;
int QMPI_Ineighbor_allgather(QMPI_Context context, int tool_id, const void *sendbuf, int sendcount,
                             MPI_Datatype sendtype, void *recvbuf, int recvcount,
                             MPI_Datatype recvtype, MPI_Comm comm, MPI_Request *request)
                                                                                                                         ;
int QMPI_Ineighbor_allgather_c(QMPI_Context context, int tool_id, const void *sendbuf,
                               MPI_Count sendcount, MPI_Datatype sendtype, void *recvbuf,
                               MPI_Count recvcount, MPI_Datatype recvtype, MPI_Comm comm,
                               MPI_Request *request)
                                                                                                                           ;
int QMPI_Ineighbor_allgatherv(QMPI_Context context, int tool_id, const void *sendbuf, int sendcount,
                              MPI_Datatype sendtype, void *recvbuf, const int recvcounts[],
                              const int displs[], MPI_Datatype recvtype, MPI_Comm comm,
                              MPI_Request *request)
                                                                                                                          ;
int QMPI_Ineighbor_allgatherv_c(QMPI_Context context, int tool_id, const void *sendbuf,
                                MPI_Count sendcount, MPI_Datatype sendtype, void *recvbuf,
                                const MPI_Count recvcounts[], const MPI_Aint displs[],
                                MPI_Datatype recvtype, MPI_Comm comm, MPI_Request *request)
                                                                                                                            ;
int QMPI_Ineighbor_alltoall(QMPI_Context context, int tool_id, const void *sendbuf, int sendcount,
                            MPI_Datatype sendtype, void *recvbuf, int recvcount,
                            MPI_Datatype recvtype, MPI_Comm comm, MPI_Request *request)
                                                                                                                        ;
int QMPI_Ineighbor_alltoall_c(QMPI_Context context, int tool_id, const void *sendbuf,
                              MPI_Count sendcount, MPI_Datatype sendtype, void *recvbuf,
                              MPI_Count recvcount, MPI_Datatype recvtype, MPI_Comm comm,
                              MPI_Request *request)
                                                                                                                          ;
int QMPI_Ineighbor_alltoallv(QMPI_Context context, int tool_id, const void *sendbuf,
                             const int sendcounts[], const int sdispls[], MPI_Datatype sendtype,
                             void *recvbuf, const int recvcounts[], const int rdispls[],
                             MPI_Datatype recvtype, MPI_Comm comm, MPI_Request *request)
                                                                                                                          ;
int QMPI_Ineighbor_alltoallv_c(QMPI_Context context, int tool_id, const void *sendbuf,
                               const MPI_Count sendcounts[], const MPI_Aint sdispls[],
                               MPI_Datatype sendtype, void *recvbuf, const MPI_Count recvcounts[],
                               const MPI_Aint rdispls[], MPI_Datatype recvtype, MPI_Comm comm,
                               MPI_Request *request)
                                                                                                                            ;
int QMPI_Ineighbor_alltoallw(QMPI_Context context, int tool_id, const void *sendbuf,
                             const int sendcounts[], const MPI_Aint sdispls[],
                             const MPI_Datatype sendtypes[], void *recvbuf, const int recvcounts[],
                             const MPI_Aint rdispls[], const MPI_Datatype recvtypes[],
                             MPI_Comm comm, MPI_Request *request) ;
int QMPI_Ineighbor_alltoallw_c(QMPI_Context context, int tool_id, const void *sendbuf,
                               const MPI_Count sendcounts[], const MPI_Aint sdispls[],
                               const MPI_Datatype sendtypes[], void *recvbuf,
                               const MPI_Count recvcounts[], const MPI_Aint rdispls[],
                               const MPI_Datatype recvtypes[], MPI_Comm comm, MPI_Request *request)
                                               ;
int QMPI_Ireduce(QMPI_Context context, int tool_id, const void *sendbuf, void *recvbuf, int count,
                 MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm, MPI_Request *request)
                                                                                                             ;
int QMPI_Ireduce_c(QMPI_Context context, int tool_id, const void *sendbuf, void *recvbuf,
                   MPI_Count count, MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm,
                   MPI_Request *request)
                                                                                                               ;
int QMPI_Ireduce_scatter(QMPI_Context context, int tool_id, const void *sendbuf, void *recvbuf,
                         const int recvcounts[], MPI_Datatype datatype, MPI_Op op, MPI_Comm comm,
                         MPI_Request *request)
                                                                                                                     ;
int QMPI_Ireduce_scatter_c(QMPI_Context context, int tool_id, const void *sendbuf, void *recvbuf,
                           const MPI_Count recvcounts[], MPI_Datatype datatype, MPI_Op op,
                           MPI_Comm comm, MPI_Request *request)
                                                                                                                       ;
int QMPI_Ireduce_scatter_block(QMPI_Context context, int tool_id, const void *sendbuf,
                               void *recvbuf, int recvcount, MPI_Datatype datatype, MPI_Op op,
                               MPI_Comm comm, MPI_Request *request)
                                                                                                                           ;
int QMPI_Ireduce_scatter_block_c(QMPI_Context context, int tool_id, const void *sendbuf,
                                 void *recvbuf, MPI_Count recvcount, MPI_Datatype datatype,
                                 MPI_Op op, MPI_Comm comm, MPI_Request *request)
                                                                                                                             ;
int QMPI_Iscan(QMPI_Context context, int tool_id, const void *sendbuf, void *recvbuf, int count,
               MPI_Datatype datatype, MPI_Op op, MPI_Comm comm, MPI_Request *request)
                                                                                                           ;
int QMPI_Iscan_c(QMPI_Context context, int tool_id, const void *sendbuf, void *recvbuf,
                 MPI_Count count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm,
                 MPI_Request *request)
                                                                                                             ;
int QMPI_Iscatter(QMPI_Context context, int tool_id, const void *sendbuf, int sendcount,
                  MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype,
                  int root, MPI_Comm comm, MPI_Request *request)
                                                                                                              ;
int QMPI_Iscatter_c(QMPI_Context context, int tool_id, const void *sendbuf, MPI_Count sendcount,
                    MPI_Datatype sendtype, void *recvbuf, MPI_Count recvcount,
                    MPI_Datatype recvtype, int root, MPI_Comm comm, MPI_Request *request)
                                                                                                                ;
int QMPI_Iscatterv(QMPI_Context context, int tool_id, const void *sendbuf, const int sendcounts[],
                   const int displs[], MPI_Datatype sendtype, void *recvbuf, int recvcount,
                   MPI_Datatype recvtype, int root, MPI_Comm comm, MPI_Request *request)
                                                                                                               ;
int QMPI_Iscatterv_c(QMPI_Context context, int tool_id, const void *sendbuf,
                     const MPI_Count sendcounts[], const MPI_Aint displs[], MPI_Datatype sendtype,
                     void *recvbuf, MPI_Count recvcount, MPI_Datatype recvtype, int root,
                     MPI_Comm comm, MPI_Request *request)
                                                                                                                 ;
int QMPI_Neighbor_allgather(QMPI_Context context, int tool_id, const void *sendbuf, int sendcount,
                            MPI_Datatype sendtype, void *recvbuf, int recvcount,
                            MPI_Datatype recvtype, MPI_Comm comm)
                                                                                                                        ;
int QMPI_Neighbor_allgather_c(QMPI_Context context, int tool_id, const void *sendbuf,
                              MPI_Count sendcount, MPI_Datatype sendtype, void *recvbuf,
                              MPI_Count recvcount, MPI_Datatype recvtype, MPI_Comm comm)
                                                                                                                          ;
int QMPI_Neighbor_allgather_init(QMPI_Context context, int tool_id, const void *sendbuf,
                                 int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount,
                                 MPI_Datatype recvtype, MPI_Comm comm, MPI_Info info,
                                 MPI_Request *request)
                                                                                                                             ;
int QMPI_Neighbor_allgather_init_c(QMPI_Context context, int tool_id, const void *sendbuf,
                                   MPI_Count sendcount, MPI_Datatype sendtype, void *recvbuf,
                                   MPI_Count recvcount, MPI_Datatype recvtype, MPI_Comm comm,
                                   MPI_Info info, MPI_Request *request)
                                                                                                                               ;
int QMPI_Neighbor_allgatherv(QMPI_Context context, int tool_id, const void *sendbuf, int sendcount,
                             MPI_Datatype sendtype, void *recvbuf, const int recvcounts[],
                             const int displs[], MPI_Datatype recvtype, MPI_Comm comm)
                                                                                                                         ;
int QMPI_Neighbor_allgatherv_c(QMPI_Context context, int tool_id, const void *sendbuf,
                               MPI_Count sendcount, MPI_Datatype sendtype, void *recvbuf,
                               const MPI_Count recvcounts[], const MPI_Aint displs[],
                               MPI_Datatype recvtype, MPI_Comm comm)
                                                                                                                           ;
int QMPI_Neighbor_allgatherv_init(QMPI_Context context, int tool_id, const void *sendbuf,
                                  int sendcount, MPI_Datatype sendtype, void *recvbuf,
                                  const int recvcounts[], const int displs[], MPI_Datatype recvtype,
                                  MPI_Comm comm, MPI_Info info, MPI_Request *request)
                                                                                                                              ;
int QMPI_Neighbor_allgatherv_init_c(QMPI_Context context, int tool_id, const void *sendbuf,
                                    MPI_Count sendcount, MPI_Datatype sendtype, void *recvbuf,
                                    const MPI_Count recvcounts[], const MPI_Aint displs[],
                                    MPI_Datatype recvtype, MPI_Comm comm, MPI_Info info,
                                    MPI_Request *request)
                                                                                                                                ;
int QMPI_Neighbor_alltoall(QMPI_Context context, int tool_id, const void *sendbuf, int sendcount,
                           MPI_Datatype sendtype, void *recvbuf, int recvcount,
                           MPI_Datatype recvtype, MPI_Comm comm)
                                                                                                                       ;
int QMPI_Neighbor_alltoall_c(QMPI_Context context, int tool_id, const void *sendbuf,
                             MPI_Count sendcount, MPI_Datatype sendtype, void *recvbuf,
                             MPI_Count recvcount, MPI_Datatype recvtype, MPI_Comm comm)
                                                                                                                         ;
int QMPI_Neighbor_alltoall_init(QMPI_Context context, int tool_id, const void *sendbuf,
                                int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount,
                                MPI_Datatype recvtype, MPI_Comm comm, MPI_Info info,
                                MPI_Request *request)
                                                                                                                            ;
int QMPI_Neighbor_alltoall_init_c(QMPI_Context context, int tool_id, const void *sendbuf,
                                  MPI_Count sendcount, MPI_Datatype sendtype, void *recvbuf,
                                  MPI_Count recvcount, MPI_Datatype recvtype, MPI_Comm comm,
                                  MPI_Info info, MPI_Request *request)
                                                                                                                              ;
int QMPI_Neighbor_alltoallv(QMPI_Context context, int tool_id, const void *sendbuf,
                            const int sendcounts[], const int sdispls[], MPI_Datatype sendtype,
                            void *recvbuf, const int recvcounts[], const int rdispls[],
                            MPI_Datatype recvtype, MPI_Comm comm)
                                                                                                                         ;
int QMPI_Neighbor_alltoallv_c(QMPI_Context context, int tool_id, const void *sendbuf,
                              const MPI_Count sendcounts[], const MPI_Aint sdispls[],
                              MPI_Datatype sendtype, void *recvbuf, const MPI_Count recvcounts[],
                              const MPI_Aint rdispls[], MPI_Datatype recvtype, MPI_Comm comm)
                                                                                                                           ;
int QMPI_Neighbor_alltoallv_init(QMPI_Context context, int tool_id, const void *sendbuf,
                                 const int sendcounts[], const int sdispls[], MPI_Datatype sendtype,
                                 void *recvbuf, const int recvcounts[], const int rdispls[],
                                 MPI_Datatype recvtype, MPI_Comm comm, MPI_Info info,
                                 MPI_Request *request)
                                                                                                                              ;
int QMPI_Neighbor_alltoallv_init_c(QMPI_Context context, int tool_id, const void *sendbuf,
                                   const MPI_Count sendcounts[], const MPI_Aint sdispls[],
                                   MPI_Datatype sendtype, void *recvbuf,
                                   const MPI_Count recvcounts[], const MPI_Aint rdispls[],
                                   MPI_Datatype recvtype, MPI_Comm comm, MPI_Info info,
                                   MPI_Request *request)
                                                                                                                                ;
int QMPI_Neighbor_alltoallw(QMPI_Context context, int tool_id, const void *sendbuf,
                            const int sendcounts[], const MPI_Aint sdispls[],
                            const MPI_Datatype sendtypes[], void *recvbuf, const int recvcounts[],
                            const MPI_Aint rdispls[], const MPI_Datatype recvtypes[],
                            MPI_Comm comm) ;
int QMPI_Neighbor_alltoallw_c(QMPI_Context context, int tool_id, const void *sendbuf,
                              const MPI_Count sendcounts[], const MPI_Aint sdispls[],
                              const MPI_Datatype sendtypes[], void *recvbuf,
                              const MPI_Count recvcounts[], const MPI_Aint rdispls[],
                              const MPI_Datatype recvtypes[], MPI_Comm comm) ;
int QMPI_Neighbor_alltoallw_init(QMPI_Context context, int tool_id, const void *sendbuf,
                                 const int sendcounts[], const MPI_Aint sdispls[],
                                 const MPI_Datatype sendtypes[], void *recvbuf,
                                 const int recvcounts[], const MPI_Aint rdispls[],
                                 const MPI_Datatype recvtypes[], MPI_Comm comm, MPI_Info info,
                                 MPI_Request *request) ;
int QMPI_Neighbor_alltoallw_init_c(QMPI_Context context, int tool_id, const void *sendbuf,
                                   const MPI_Count sendcounts[], const MPI_Aint sdispls[],
                                   const MPI_Datatype sendtypes[], void *recvbuf,
                                   const MPI_Count recvcounts[], const MPI_Aint rdispls[],
                                   const MPI_Datatype recvtypes[], MPI_Comm comm, MPI_Info info,
                                   MPI_Request *request) ;
int QMPI_Reduce(QMPI_Context context, int tool_id, const void *sendbuf, void *recvbuf, int count,
                MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm)
                                                                                                            ;
int QMPI_Reduce_c(QMPI_Context context, int tool_id, const void *sendbuf, void *recvbuf,
                  MPI_Count count, MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm)
                                                                                                              ;
int QMPI_Reduce_init(QMPI_Context context, int tool_id, const void *sendbuf, void *recvbuf,
                     int count, MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm,
                     MPI_Info info, MPI_Request *request)
                                                                                                                 ;
int QMPI_Reduce_init_c(QMPI_Context context, int tool_id, const void *sendbuf, void *recvbuf,
                       MPI_Count count, MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm,
                       MPI_Info info, MPI_Request *request)
                                                                                                                   ;
int QMPI_Reduce_local(QMPI_Context context, int tool_id, const void *inbuf, void *inoutbuf,
                      int count, MPI_Datatype datatype, MPI_Op op)
                                                                                                                  ;
int QMPI_Reduce_local_c(QMPI_Context context, int tool_id, const void *inbuf, void *inoutbuf,
                        MPI_Count count, MPI_Datatype datatype, MPI_Op op)
                                                                                                                    ;
int QMPI_Reduce_scatter(QMPI_Context context, int tool_id, const void *sendbuf, void *recvbuf,
                        const int recvcounts[], MPI_Datatype datatype, MPI_Op op, MPI_Comm comm)
                                                                                                                    ;
int QMPI_Reduce_scatter_c(QMPI_Context context, int tool_id, const void *sendbuf, void *recvbuf,
                          const MPI_Count recvcounts[], MPI_Datatype datatype, MPI_Op op,
                          MPI_Comm comm)
                                                                                                                      ;
int QMPI_Reduce_scatter_block(QMPI_Context context, int tool_id, const void *sendbuf, void *recvbuf,
                              int recvcount, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm)
                                                                                                                          ;
int QMPI_Reduce_scatter_block_c(QMPI_Context context, int tool_id, const void *sendbuf,
                                void *recvbuf, MPI_Count recvcount, MPI_Datatype datatype,
                                MPI_Op op, MPI_Comm comm)
                                                                                                                            ;
int QMPI_Reduce_scatter_block_init(QMPI_Context context, int tool_id, const void *sendbuf,
                                   void *recvbuf, int recvcount, MPI_Datatype datatype, MPI_Op op,
                                   MPI_Comm comm, MPI_Info info, MPI_Request *request)
                                                                                                                               ;
int QMPI_Reduce_scatter_block_init_c(QMPI_Context context, int tool_id, const void *sendbuf,
                                     void *recvbuf, MPI_Count recvcount, MPI_Datatype datatype,
                                     MPI_Op op, MPI_Comm comm, MPI_Info info, MPI_Request *request)
                                                                                                                                 ;
int QMPI_Reduce_scatter_init(QMPI_Context context, int tool_id, const void *sendbuf, void *recvbuf,
                             const int recvcounts[], MPI_Datatype datatype, MPI_Op op,
                             MPI_Comm comm, MPI_Info info, MPI_Request *request)
                                                                                                                         ;
int QMPI_Reduce_scatter_init_c(QMPI_Context context, int tool_id, const void *sendbuf,
                               void *recvbuf, const MPI_Count recvcounts[], MPI_Datatype datatype,
                               MPI_Op op, MPI_Comm comm, MPI_Info info, MPI_Request *request)
                                                                                                                           ;
int QMPI_Scan(QMPI_Context context, int tool_id, const void *sendbuf, void *recvbuf, int count,
              MPI_Datatype datatype, MPI_Op op, MPI_Comm comm)
                                                                                                          ;
int QMPI_Scan_c(QMPI_Context context, int tool_id, const void *sendbuf, void *recvbuf,
                MPI_Count count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm)
                                                                                                            ;
int QMPI_Scan_init(QMPI_Context context, int tool_id, const void *sendbuf, void *recvbuf, int count,
                   MPI_Datatype datatype, MPI_Op op, MPI_Comm comm, MPI_Info info,
                   MPI_Request *request)
                                                                                                               ;
int QMPI_Scan_init_c(QMPI_Context context, int tool_id, const void *sendbuf, void *recvbuf,
                     MPI_Count count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm,
                     MPI_Info info, MPI_Request *request)
                                                                                                                 ;
int QMPI_Scatter(QMPI_Context context, int tool_id, const void *sendbuf, int sendcount,
                 MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype,
                 int root, MPI_Comm comm)
                                                                                                             ;
int QMPI_Scatter_c(QMPI_Context context, int tool_id, const void *sendbuf, MPI_Count sendcount,
                   MPI_Datatype sendtype, void *recvbuf, MPI_Count recvcount, MPI_Datatype recvtype,
                   int root, MPI_Comm comm)
                                                                                                               ;
int QMPI_Scatter_init(QMPI_Context context, int tool_id, const void *sendbuf, int sendcount,
                      MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype,
                      int root, MPI_Comm comm, MPI_Info info, MPI_Request *request)
                                                                                                                  ;
int QMPI_Scatter_init_c(QMPI_Context context, int tool_id, const void *sendbuf, MPI_Count sendcount,
                        MPI_Datatype sendtype, void *recvbuf, MPI_Count recvcount,
                        MPI_Datatype recvtype, int root, MPI_Comm comm, MPI_Info info,
                        MPI_Request *request)
                                                                                                                    ;
int QMPI_Scatterv(QMPI_Context context, int tool_id, const void *sendbuf, const int sendcounts[],
                  const int displs[], MPI_Datatype sendtype, void *recvbuf, int recvcount,
                  MPI_Datatype recvtype, int root, MPI_Comm comm)
                                                                                                              ;
int QMPI_Scatterv_c(QMPI_Context context, int tool_id, const void *sendbuf,
                    const MPI_Count sendcounts[], const MPI_Aint displs[], MPI_Datatype sendtype,
                    void *recvbuf, MPI_Count recvcount, MPI_Datatype recvtype, int root,
                    MPI_Comm comm)
                                                                                                                ;
int QMPI_Scatterv_init(QMPI_Context context, int tool_id, const void *sendbuf,
                       const int sendcounts[], const int displs[], MPI_Datatype sendtype,
                       void *recvbuf, int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm,
                       MPI_Info info, MPI_Request *request)
                                                                                                                   ;
int QMPI_Scatterv_init_c(QMPI_Context context, int tool_id, const void *sendbuf,
                         const MPI_Count sendcounts[], const MPI_Aint displs[],
                         MPI_Datatype sendtype, void *recvbuf, MPI_Count recvcount,
                         MPI_Datatype recvtype, int root, MPI_Comm comm, MPI_Info info,
                         MPI_Request *request)
                                                                                                                     ;
int QMPI_Comm_compare(QMPI_Context context, int tool_id, MPI_Comm comm1, MPI_Comm comm2,
                      int *result) ;
int QMPI_Comm_create(QMPI_Context context, int tool_id, MPI_Comm comm, MPI_Group group,
                     MPI_Comm *newcomm) ;
int QMPI_Comm_create_group(QMPI_Context context, int tool_id, MPI_Comm comm, MPI_Group group,
                           int tag, MPI_Comm *newcomm) ;
int QMPI_Comm_dup(QMPI_Context context, int tool_id, MPI_Comm comm, MPI_Comm *newcomm)
                    ;
int QMPI_Comm_dup_with_info(QMPI_Context context, int tool_id, MPI_Comm comm, MPI_Info info,
                            MPI_Comm *newcomm) ;
int QMPI_Comm_free(QMPI_Context context, int tool_id, MPI_Comm *comm) ;
int QMPI_Comm_get_info(QMPI_Context context, int tool_id, MPI_Comm comm, MPI_Info *info_used)
                    ;
int QMPI_Comm_get_name(QMPI_Context context, int tool_id, MPI_Comm comm, char *comm_name,
                       int *resultlen) ;
int QMPI_Comm_group(QMPI_Context context, int tool_id, MPI_Comm comm, MPI_Group *group)
                    ;
int QMPI_Comm_idup(QMPI_Context context, int tool_id, MPI_Comm comm, MPI_Comm *newcomm,
                   MPI_Request *request) ;
int QMPI_Comm_idup_with_info(QMPI_Context context, int tool_id, MPI_Comm comm, MPI_Info info,
                             MPI_Comm *newcomm, MPI_Request *request) ;
int QMPI_Comm_rank(QMPI_Context context, int tool_id, MPI_Comm comm, int *rank) ;
int QMPI_Comm_remote_group(QMPI_Context context, int tool_id, MPI_Comm comm, MPI_Group *group)
                    ;
int QMPI_Comm_remote_size(QMPI_Context context, int tool_id, MPI_Comm comm, int *size)
                    ;
int QMPI_Comm_set_info(QMPI_Context context, int tool_id, MPI_Comm comm, MPI_Info info)
                    ;
int QMPI_Comm_set_name(QMPI_Context context, int tool_id, MPI_Comm comm, const char *comm_name)
                    ;
int QMPI_Comm_size(QMPI_Context context, int tool_id, MPI_Comm comm, int *size) ;
int QMPI_Comm_split(QMPI_Context context, int tool_id, MPI_Comm comm, int color, int key,
                    MPI_Comm *newcomm) ;
int QMPI_Comm_split_type(QMPI_Context context, int tool_id, MPI_Comm comm, int split_type, int key,
                         MPI_Info info, MPI_Comm *newcomm) ;
int QMPI_Comm_test_inter(QMPI_Context context, int tool_id, MPI_Comm comm, int *flag)
                    ;
int QMPI_Intercomm_create(QMPI_Context context, int tool_id, MPI_Comm local_comm, int local_leader,
                          MPI_Comm peer_comm, int remote_leader, int tag, MPI_Comm *newintercomm)
                                          ;
int QMPI_Intercomm_create_from_groups(QMPI_Context context, int tool_id, MPI_Group local_group,
                                      int local_leader, MPI_Group remote_group, int remote_leader,
                                      const char *stringtag, MPI_Info info,
                                      MPI_Errhandler errhandler, MPI_Comm *newintercomm)
                                                      ;
int QMPI_Intercomm_merge(QMPI_Context context, int tool_id, MPI_Comm intercomm, int high,
                         MPI_Comm *newintracomm) ;
int QMPIX_Comm_test_threadcomm(QMPI_Context context, int tool_id, MPI_Comm comm, int *flag)
                    ;
int QMPIX_Comm_revoke(QMPI_Context context, int tool_id, MPI_Comm comm) ;
int QMPIX_Comm_shrink(QMPI_Context context, int tool_id, MPI_Comm comm, MPI_Comm *newcomm)
                    ;
int QMPIX_Comm_failure_ack(QMPI_Context context, int tool_id, MPI_Comm comm) ;
int QMPIX_Comm_failure_get_acked(QMPI_Context context, int tool_id, MPI_Comm comm,
                                 MPI_Group *failedgrp) ;
int QMPIX_Comm_agree(QMPI_Context context, int tool_id, MPI_Comm comm, int *flag) ;
int QMPIX_Comm_get_failed(QMPI_Context context, int tool_id, MPI_Comm comm, MPI_Group *failedgrp)
                    ;
int QMPI_Get_address(QMPI_Context context, int tool_id, const void *location, MPI_Aint *address)
                    ;
int QMPI_Address(QMPI_Context context, int tool_id, void *location, MPI_Aint *address)
                    ;
int QMPI_Get_count(QMPI_Context context, int tool_id, const MPI_Status *status,
                   MPI_Datatype datatype, int *count) ;
int QMPI_Get_count_c(QMPI_Context context, int tool_id, const MPI_Status *status,
                     MPI_Datatype datatype, MPI_Count *count) ;
int QMPI_Get_elements(QMPI_Context context, int tool_id, const MPI_Status *status,
                      MPI_Datatype datatype, int *count) ;
int QMPI_Get_elements_c(QMPI_Context context, int tool_id, const MPI_Status *status,
                        MPI_Datatype datatype, MPI_Count *count) ;
int QMPI_Get_elements_x(QMPI_Context context, int tool_id, const MPI_Status *status,
                        MPI_Datatype datatype, MPI_Count *count) ;
int QMPI_Pack(QMPI_Context context, int tool_id, const void *inbuf, int incount,
              MPI_Datatype datatype, void *outbuf, int outsize, int *position, MPI_Comm comm)
                              ;
int QMPI_Pack_c(QMPI_Context context, int tool_id, const void *inbuf, MPI_Count incount,
                MPI_Datatype datatype, void *outbuf, MPI_Count outsize, MPI_Count *position,
                MPI_Comm comm) ;
int QMPI_Pack_external(QMPI_Context context, int tool_id, const char *datarep, const void *inbuf,
                       int incount, MPI_Datatype datatype, void *outbuf, MPI_Aint outsize,
                       MPI_Aint *position) ;
int QMPI_Pack_external_c(QMPI_Context context, int tool_id, const char *datarep, const void *inbuf,
                         MPI_Count incount, MPI_Datatype datatype, void *outbuf, MPI_Count outsize,
                         MPI_Count *position) ;
int QMPI_Pack_external_size(QMPI_Context context, int tool_id, const char *datarep, int incount,
                            MPI_Datatype datatype, MPI_Aint *size) ;
int QMPI_Pack_external_size_c(QMPI_Context context, int tool_id, const char *datarep,
                              MPI_Count incount, MPI_Datatype datatype, MPI_Count *size)
                                              ;
int QMPI_Pack_size(QMPI_Context context, int tool_id, int incount, MPI_Datatype datatype,
                   MPI_Comm comm, int *size) ;
int QMPI_Pack_size_c(QMPI_Context context, int tool_id, MPI_Count incount, MPI_Datatype datatype,
                     MPI_Comm comm, MPI_Count *size) ;
int QMPI_Status_set_elements(QMPI_Context context, int tool_id, MPI_Status *status,
                             MPI_Datatype datatype, int count) ;
int QMPI_Status_set_elements_c(QMPI_Context context, int tool_id, MPI_Status *status,
                               MPI_Datatype datatype, MPI_Count count) ;
int QMPI_Status_set_elements_x(QMPI_Context context, int tool_id, MPI_Status *status,
                               MPI_Datatype datatype, MPI_Count count) ;
int QMPI_Type_commit(QMPI_Context context, int tool_id, MPI_Datatype *datatype) ;
int QMPI_Type_contiguous(QMPI_Context context, int tool_id, int count, MPI_Datatype oldtype,
                         MPI_Datatype *newtype) ;
int QMPI_Type_contiguous_c(QMPI_Context context, int tool_id, MPI_Count count, MPI_Datatype oldtype,
                           MPI_Datatype *newtype) ;
int QMPI_Type_create_darray(QMPI_Context context, int tool_id, int size, int rank, int ndims,
                            const int array_of_gsizes[], const int array_of_distribs[],
                            const int array_of_dargs[], const int array_of_psizes[], int order,
                            MPI_Datatype oldtype, MPI_Datatype *newtype) ;
int QMPI_Type_create_darray_c(QMPI_Context context, int tool_id, int size, int rank, int ndims,
                              const MPI_Count array_of_gsizes[], const int array_of_distribs[],
                              const int array_of_dargs[], const int array_of_psizes[], int order,
                              MPI_Datatype oldtype, MPI_Datatype *newtype) ;
int QMPI_Type_create_f90_complex(QMPI_Context context, int tool_id, int p, int r,
                                 MPI_Datatype *newtype) ;
int QMPI_Type_create_f90_integer(QMPI_Context context, int tool_id, int r, MPI_Datatype *newtype)
                    ;
int QMPI_Type_create_f90_real(QMPI_Context context, int tool_id, int p, int r,
                              MPI_Datatype *newtype) ;
int QMPI_Type_create_hindexed(QMPI_Context context, int tool_id, int count,
                              const int array_of_blocklengths[],
                              const MPI_Aint array_of_displacements[], MPI_Datatype oldtype,
                              MPI_Datatype *newtype) ;
int QMPI_Type_create_hindexed_c(QMPI_Context context, int tool_id, MPI_Count count,
                                const MPI_Count array_of_blocklengths[],
                                const MPI_Count array_of_displacements[], MPI_Datatype oldtype,
                                MPI_Datatype *newtype) ;
int QMPI_Type_hindexed(QMPI_Context context, int tool_id, int count, int array_of_blocklengths[],
                       MPI_Aint array_of_displacements[], MPI_Datatype oldtype,
                       MPI_Datatype *newtype) ;
int QMPI_Type_create_hindexed_block(QMPI_Context context, int tool_id, int count, int blocklength,
                                    const MPI_Aint array_of_displacements[], MPI_Datatype oldtype,
                                    MPI_Datatype *newtype) ;
int QMPI_Type_create_hindexed_block_c(QMPI_Context context, int tool_id, MPI_Count count,
                                      MPI_Count blocklength,
                                      const MPI_Count array_of_displacements[],
                                      MPI_Datatype oldtype, MPI_Datatype *newtype)
                                                      ;
int QMPI_Type_create_hvector(QMPI_Context context, int tool_id, int count, int blocklength,
                             MPI_Aint stride, MPI_Datatype oldtype, MPI_Datatype *newtype)
                                             ;
int QMPI_Type_create_hvector_c(QMPI_Context context, int tool_id, MPI_Count count,
                               MPI_Count blocklength, MPI_Count stride, MPI_Datatype oldtype,
                               MPI_Datatype *newtype) ;
int QMPI_Type_hvector(QMPI_Context context, int tool_id, int count, int blocklength,
                      MPI_Aint stride, MPI_Datatype oldtype, MPI_Datatype *newtype)
                                      ;
int QMPI_Type_create_indexed_block(QMPI_Context context, int tool_id, int count, int blocklength,
                                   const int array_of_displacements[], MPI_Datatype oldtype,
                                   MPI_Datatype *newtype) ;
int QMPI_Type_create_indexed_block_c(QMPI_Context context, int tool_id, MPI_Count count,
                                     MPI_Count blocklength,
                                     const MPI_Count array_of_displacements[], MPI_Datatype oldtype,
                                     MPI_Datatype *newtype) ;
int QMPI_Type_create_resized(QMPI_Context context, int tool_id, MPI_Datatype oldtype, MPI_Aint lb,
                             MPI_Aint extent, MPI_Datatype *newtype) ;
int QMPI_Type_create_resized_c(QMPI_Context context, int tool_id, MPI_Datatype oldtype,
                               MPI_Count lb, MPI_Count extent, MPI_Datatype *newtype)
                                               ;
int QMPI_Type_create_struct(QMPI_Context context, int tool_id, int count,
                            const int array_of_blocklengths[],
                            const MPI_Aint array_of_displacements[],
                            const MPI_Datatype array_of_types[], MPI_Datatype *newtype)
                                            ;
int QMPI_Type_create_struct_c(QMPI_Context context, int tool_id, MPI_Count count,
                              const MPI_Count array_of_blocklengths[],
                              const MPI_Count array_of_displacements[],
                              const MPI_Datatype array_of_types[], MPI_Datatype *newtype)
                                              ;
int QMPI_Type_struct(QMPI_Context context, int tool_id, int count, int array_of_blocklengths[],
                     MPI_Aint array_of_displacements[], MPI_Datatype array_of_types[],
                     MPI_Datatype *newtype) ;
int QMPI_Type_create_subarray(QMPI_Context context, int tool_id, int ndims,
                              const int array_of_sizes[], const int array_of_subsizes[],
                              const int array_of_starts[], int order, MPI_Datatype oldtype,
                              MPI_Datatype *newtype) ;
int QMPI_Type_create_subarray_c(QMPI_Context context, int tool_id, int ndims,
                                const MPI_Count array_of_sizes[],
                                const MPI_Count array_of_subsizes[],
                                const MPI_Count array_of_starts[], int order, MPI_Datatype oldtype,
                                MPI_Datatype *newtype) ;
int QMPI_Type_dup(QMPI_Context context, int tool_id, MPI_Datatype oldtype, MPI_Datatype *newtype)
                    ;
int QMPI_Type_free(QMPI_Context context, int tool_id, MPI_Datatype *datatype) ;
int QMPI_Type_get_contents(QMPI_Context context, int tool_id, MPI_Datatype datatype,
                           int max_integers, int max_addresses, int max_datatypes,
                           int array_of_integers[], MPI_Aint array_of_addresses[],
                           MPI_Datatype array_of_datatypes[]) ;
int QMPI_Type_get_contents_c(QMPI_Context context, int tool_id, MPI_Datatype datatype,
                             MPI_Count max_integers, MPI_Count max_addresses,
                             MPI_Count max_large_counts, MPI_Count max_datatypes,
                             int array_of_integers[], MPI_Aint array_of_addresses[],
                             MPI_Count array_of_large_counts[], MPI_Datatype array_of_datatypes[])
                                             ;
int QMPI_Type_get_envelope(QMPI_Context context, int tool_id, MPI_Datatype datatype,
                           int *num_integers, int *num_addresses, int *num_datatypes,
                           int *combiner) ;
int QMPI_Type_get_envelope_c(QMPI_Context context, int tool_id, MPI_Datatype datatype,
                             MPI_Count *num_integers, MPI_Count *num_addresses,
                             MPI_Count *num_large_counts, MPI_Count *num_datatypes, int *combiner)
                                             ;
int QMPI_Type_get_extent(QMPI_Context context, int tool_id, MPI_Datatype datatype, MPI_Aint *lb,
                         MPI_Aint *extent) ;
int QMPI_Type_get_extent_c(QMPI_Context context, int tool_id, MPI_Datatype datatype, MPI_Count *lb,
                           MPI_Count *extent) ;
int QMPI_Type_get_extent_x(QMPI_Context context, int tool_id, MPI_Datatype datatype, MPI_Count *lb,
                           MPI_Count *extent) ;
int QMPI_Type_get_name(QMPI_Context context, int tool_id, MPI_Datatype datatype, char *type_name,
                       int *resultlen) ;
int QMPI_Type_get_true_extent(QMPI_Context context, int tool_id, MPI_Datatype datatype,
                              MPI_Aint *true_lb, MPI_Aint *true_extent) ;
int QMPI_Type_get_true_extent_c(QMPI_Context context, int tool_id, MPI_Datatype datatype,
                                MPI_Count *true_lb, MPI_Count *true_extent) ;
int QMPI_Type_get_true_extent_x(QMPI_Context context, int tool_id, MPI_Datatype datatype,
                                MPI_Count *true_lb, MPI_Count *true_extent) ;
int QMPI_Type_get_value_index(QMPI_Context context, int tool_id, MPI_Datatype value_type,
                              MPI_Datatype index_type, MPI_Datatype *pair_type) ;
int QMPI_Type_indexed(QMPI_Context context, int tool_id, int count,
                      const int array_of_blocklengths[], const int array_of_displacements[],
                      MPI_Datatype oldtype, MPI_Datatype *newtype) ;
int QMPI_Type_indexed_c(QMPI_Context context, int tool_id, MPI_Count count,
                        const MPI_Count array_of_blocklengths[],
                        const MPI_Count array_of_displacements[], MPI_Datatype oldtype,
                        MPI_Datatype *newtype) ;
int QMPI_Type_match_size(QMPI_Context context, int tool_id, int typeclass, int size,
                         MPI_Datatype *datatype) ;
int QMPI_Type_set_name(QMPI_Context context, int tool_id, MPI_Datatype datatype,
                       const char *type_name) ;
int QMPI_Type_size(QMPI_Context context, int tool_id, MPI_Datatype datatype, int *size)
                    ;
int QMPI_Type_size_c(QMPI_Context context, int tool_id, MPI_Datatype datatype, MPI_Count *size)
                    ;
int QMPI_Type_size_x(QMPI_Context context, int tool_id, MPI_Datatype datatype, MPI_Count *size)
                    ;
int QMPI_Type_vector(QMPI_Context context, int tool_id, int count, int blocklength, int stride,
                     MPI_Datatype oldtype, MPI_Datatype *newtype) ;
int QMPI_Type_vector_c(QMPI_Context context, int tool_id, MPI_Count count, MPI_Count blocklength,
                       MPI_Count stride, MPI_Datatype oldtype, MPI_Datatype *newtype)
                                       ;
int QMPI_Unpack(QMPI_Context context, int tool_id, const void *inbuf, int insize, int *position,
                void *outbuf, int outcount, MPI_Datatype datatype, MPI_Comm comm) ;
int QMPI_Unpack_c(QMPI_Context context, int tool_id, const void *inbuf, MPI_Count insize,
                  MPI_Count *position, void *outbuf, MPI_Count outcount, MPI_Datatype datatype,
                  MPI_Comm comm) ;
int QMPI_Unpack_external(QMPI_Context context, int tool_id, const char datarep[], const void *inbuf,
                         MPI_Aint insize, MPI_Aint *position, void *outbuf, int outcount,
                         MPI_Datatype datatype) ;
int QMPI_Unpack_external_c(QMPI_Context context, int tool_id, const char datarep[],
                           const void *inbuf, MPI_Count insize, MPI_Count *position, void *outbuf,
                           MPI_Count outcount, MPI_Datatype datatype) ;
int QMPI_Type_extent(QMPI_Context context, int tool_id, MPI_Datatype datatype, MPI_Aint *extent)
                    ;
int QMPI_Type_lb(QMPI_Context context, int tool_id, MPI_Datatype datatype, MPI_Aint *displacement)
                    ;
int QMPI_Type_ub(QMPI_Context context, int tool_id, MPI_Datatype datatype, MPI_Aint *displacement)
                    ;
int QMPIX_Type_iov_len(QMPI_Context context, int tool_id, MPI_Datatype datatype,
                       MPI_Count max_iov_bytes, MPI_Count *iov_len, MPI_Count *actual_iov_bytes)
                                       ;
int QMPIX_Type_iov(QMPI_Context context, int tool_id, MPI_Datatype datatype, MPI_Count iov_offset,
                   MPIX_Iov *iov, MPI_Count max_iov_len, MPI_Count *actual_iov_len)
                                   ;
int QMPI_Add_error_class(QMPI_Context context, int tool_id, int *errorclass) ;
int QMPI_Add_error_code(QMPI_Context context, int tool_id, int errorclass, int *errorcode)
                    ;
int QMPI_Add_error_string(QMPI_Context context, int tool_id, int errorcode, const char *string)
                    ;
int QMPI_Comm_call_errhandler(QMPI_Context context, int tool_id, MPI_Comm comm, int errorcode)
                    ;
int QMPI_Comm_create_errhandler(QMPI_Context context, int tool_id,
                                MPI_Comm_errhandler_function *comm_errhandler_fn,
                                MPI_Errhandler *errhandler) ;
int QMPI_Errhandler_create(QMPI_Context context, int tool_id,
                           MPI_Comm_errhandler_function *comm_errhandler_fn,
                           MPI_Errhandler *errhandler) ;
int QMPI_Comm_get_errhandler(QMPI_Context context, int tool_id, MPI_Comm comm,
                             MPI_Errhandler *errhandler) ;
int QMPI_Errhandler_get(QMPI_Context context, int tool_id, MPI_Comm comm,
                        MPI_Errhandler *errhandler) ;
int QMPI_Comm_set_errhandler(QMPI_Context context, int tool_id, MPI_Comm comm,
                             MPI_Errhandler errhandler) ;
int QMPI_Errhandler_set(QMPI_Context context, int tool_id, MPI_Comm comm,
                        MPI_Errhandler errhandler) ;
int QMPI_Errhandler_free(QMPI_Context context, int tool_id, MPI_Errhandler *errhandler)
                    ;
int QMPI_Error_class(QMPI_Context context, int tool_id, int errorcode, int *errorclass)
                    ;
int QMPI_Error_string(QMPI_Context context, int tool_id, int errorcode, char *string,
                      int *resultlen) ;
int QMPI_File_call_errhandler(QMPI_Context context, int tool_id, MPI_File fh, int errorcode)
                    ;
int QMPI_File_create_errhandler(QMPI_Context context, int tool_id,
                                MPI_File_errhandler_function *file_errhandler_fn,
                                MPI_Errhandler *errhandler) ;
int QMPI_File_get_errhandler(QMPI_Context context, int tool_id, MPI_File file,
                             MPI_Errhandler *errhandler) ;
int QMPI_File_set_errhandler(QMPI_Context context, int tool_id, MPI_File file,
                             MPI_Errhandler errhandler) ;
int QMPI_Remove_error_class(QMPI_Context context, int tool_id, int errorclass) ;
int QMPI_Remove_error_code(QMPI_Context context, int tool_id, int errorcode) ;
int QMPI_Remove_error_string(QMPI_Context context, int tool_id, int errorcode) ;
int QMPI_Session_call_errhandler(QMPI_Context context, int tool_id, MPI_Session session,
                                 int errorcode) ;
int QMPI_Session_create_errhandler(QMPI_Context context, int tool_id,
                                   MPI_Session_errhandler_function *session_errhandler_fn,
                                   MPI_Errhandler *errhandler) ;
int QMPI_Session_get_errhandler(QMPI_Context context, int tool_id, MPI_Session session,
                                MPI_Errhandler *errhandler) ;
int QMPI_Session_set_errhandler(QMPI_Context context, int tool_id, MPI_Session session,
                                MPI_Errhandler errhandler) ;
int QMPI_Win_call_errhandler(QMPI_Context context, int tool_id, MPI_Win win, int errorcode)
                    ;
int QMPI_Win_create_errhandler(QMPI_Context context, int tool_id,
                               MPI_Win_errhandler_function *win_errhandler_fn,
                               MPI_Errhandler *errhandler) ;
int QMPI_Win_get_errhandler(QMPI_Context context, int tool_id, MPI_Win win,
                            MPI_Errhandler *errhandler) ;
int QMPI_Win_set_errhandler(QMPI_Context context, int tool_id, MPI_Win win,
                            MPI_Errhandler errhandler) ;
MPI_Fint QMPI_Comm_c2f(QMPI_Context context, int tool_id, MPI_Comm comm) ;
MPI_Comm QMPI_Comm_f2c(QMPI_Context context, int tool_id, MPI_Fint comm) ;
MPI_Fint QMPI_Errhandler_c2f(QMPI_Context context, int tool_id, MPI_Errhandler errhandler)
                    ;
MPI_Errhandler QMPI_Errhandler_f2c(QMPI_Context context, int tool_id, MPI_Fint errhandler)
                    ;
MPI_Fint QMPI_Group_c2f(QMPI_Context context, int tool_id, MPI_Group group) ;
MPI_Group QMPI_Group_f2c(QMPI_Context context, int tool_id, MPI_Fint group) ;
MPI_Fint QMPI_Info_c2f(QMPI_Context context, int tool_id, MPI_Info info) ;
MPI_Info QMPI_Info_f2c(QMPI_Context context, int tool_id, MPI_Fint info) ;
MPI_Fint QMPI_Message_c2f(QMPI_Context context, int tool_id, MPI_Message message) ;
MPI_Message QMPI_Message_f2c(QMPI_Context context, int tool_id, MPI_Fint message) ;
MPI_Fint QMPI_Op_c2f(QMPI_Context context, int tool_id, MPI_Op op) ;
MPI_Op QMPI_Op_f2c(QMPI_Context context, int tool_id, MPI_Fint op) ;
MPI_Fint QMPI_Request_c2f(QMPI_Context context, int tool_id, MPI_Request request) ;
MPI_Request QMPI_Request_f2c(QMPI_Context context, int tool_id, MPI_Fint request) ;
MPI_Fint QMPI_Session_c2f(QMPI_Context context, int tool_id, MPI_Session session) ;
MPI_Session QMPI_Session_f2c(QMPI_Context context, int tool_id, MPI_Fint session) ;
int QMPI_Status_c2f(QMPI_Context context, int tool_id, const MPI_Status *c_status,
                    MPI_Fint *f_status) ;
int QMPI_Status_c2f08(QMPI_Context context, int tool_id, const MPI_Status *c_status,
                      MPI_F08_status *f08_status) ;
int QMPI_Status_f082c(QMPI_Context context, int tool_id, const MPI_F08_status *f08_status,
                      MPI_Status *c_status) ;
int QMPI_Status_f082f(QMPI_Context context, int tool_id, const MPI_F08_status *f08_status,
                      MPI_Fint *f_status) ;
int QMPI_Status_f2c(QMPI_Context context, int tool_id, const MPI_Fint *f_status,
                    MPI_Status *c_status) ;
int QMPI_Status_f2f08(QMPI_Context context, int tool_id, const MPI_Fint *f_status,
                      MPI_F08_status *f08_status) ;
MPI_Fint QMPI_Type_c2f(QMPI_Context context, int tool_id, MPI_Datatype datatype) ;
MPI_Datatype QMPI_Type_f2c(QMPI_Context context, int tool_id, MPI_Fint datatype) ;
MPI_Fint QMPI_Win_c2f(QMPI_Context context, int tool_id, MPI_Win win) ;
MPI_Win QMPI_Win_f2c(QMPI_Context context, int tool_id, MPI_Fint win) ;
int QMPI_Group_compare(QMPI_Context context, int tool_id, MPI_Group group1, MPI_Group group2,
                       int *result) ;
int QMPI_Group_difference(QMPI_Context context, int tool_id, MPI_Group group1, MPI_Group group2,
                          MPI_Group *newgroup) ;
int QMPI_Group_excl(QMPI_Context context, int tool_id, MPI_Group group, int n, const int ranks[],
                    MPI_Group *newgroup) ;
int QMPI_Group_free(QMPI_Context context, int tool_id, MPI_Group *group) ;
int QMPI_Group_incl(QMPI_Context context, int tool_id, MPI_Group group, int n, const int ranks[],
                    MPI_Group *newgroup) ;
int QMPI_Group_intersection(QMPI_Context context, int tool_id, MPI_Group group1, MPI_Group group2,
                            MPI_Group *newgroup) ;
int QMPI_Group_range_excl(QMPI_Context context, int tool_id, MPI_Group group, int n,
                          int ranges[][3], MPI_Group *newgroup) ;
int QMPI_Group_range_incl(QMPI_Context context, int tool_id, MPI_Group group, int n,
                          int ranges[][3], MPI_Group *newgroup) ;
int QMPI_Group_rank(QMPI_Context context, int tool_id, MPI_Group group, int *rank)
                    ;
int QMPI_Group_size(QMPI_Context context, int tool_id, MPI_Group group, int *size)
                    ;
int QMPI_Group_translate_ranks(QMPI_Context context, int tool_id, MPI_Group group1, int n,
                               const int ranks1[], MPI_Group group2, int ranks2[])
                                               ;
int QMPI_Group_union(QMPI_Context context, int tool_id, MPI_Group group1, MPI_Group group2,
                     MPI_Group *newgroup) ;
int QMPI_Info_create(QMPI_Context context, int tool_id, MPI_Info *info) ;
int QMPI_Info_create_env(QMPI_Context context, int tool_id, int argc, char *argv[], MPI_Info *info)
                    ;
int QMPI_Info_delete(QMPI_Context context, int tool_id, MPI_Info info, const char *key)
                    ;
int QMPI_Info_dup(QMPI_Context context, int tool_id, MPI_Info info, MPI_Info *newinfo)
                    ;
int QMPI_Info_free(QMPI_Context context, int tool_id, MPI_Info *info) ;
int QMPI_Info_get(QMPI_Context context, int tool_id, MPI_Info info, const char *key, int valuelen,
                  char *value, int *flag) ;
int QMPI_Info_get_nkeys(QMPI_Context context, int tool_id, MPI_Info info, int *nkeys)
                    ;
int QMPI_Info_get_nthkey(QMPI_Context context, int tool_id, MPI_Info info, int n, char *key)
                    ;
int QMPI_Info_get_string(QMPI_Context context, int tool_id, MPI_Info info, const char *key,
                         int *buflen, char *value, int *flag) ;
int QMPI_Info_get_valuelen(QMPI_Context context, int tool_id, MPI_Info info, const char *key,
                           int *valuelen, int *flag) ;
int QMPI_Info_set(QMPI_Context context, int tool_id, MPI_Info info, const char *key,
                  const char *value) ;
int QMPIX_Info_set_hex(QMPI_Context context, int tool_id, MPI_Info info, const char *key,
                       const void *value, int value_size) ;
int QMPI_Abort(QMPI_Context context, int tool_id, MPI_Comm comm, int errorcode) ;
int QMPI_Comm_create_from_group(QMPI_Context context, int tool_id, MPI_Group group,
                                const char *stringtag, MPI_Info info, MPI_Errhandler errhandler,
                                MPI_Comm *newcomm) ;
int QMPI_Finalize(QMPI_Context context, int tool_id) ;
int QMPI_Finalized(QMPI_Context context, int tool_id, int *flag) ;
int QMPI_Group_from_session_pset(QMPI_Context context, int tool_id, MPI_Session session,
                                 const char *pset_name, MPI_Group *newgroup) ;
int QMPI_Init(QMPI_Context context, int tool_id, int *argc, char ***argv) ;
int QMPI_Init_thread(QMPI_Context context, int tool_id, int *argc, char ***argv, int required,
                     int *provided) ;
int QMPI_Initialized(QMPI_Context context, int tool_id, int *flag) ;
int QMPI_Is_thread_main(QMPI_Context context, int tool_id, int *flag) ;
int QMPI_Query_thread(QMPI_Context context, int tool_id, int *provided) ;
int QMPI_Session_finalize(QMPI_Context context, int tool_id, MPI_Session *session)
                    ;
int QMPI_Session_get_info(QMPI_Context context, int tool_id, MPI_Session session,
                          MPI_Info *info_used) ;
int QMPI_Session_get_nth_pset(QMPI_Context context, int tool_id, MPI_Session session, MPI_Info info,
                              int n, int *pset_len, char *pset_name) ;
int QMPI_Session_get_num_psets(QMPI_Context context, int tool_id, MPI_Session session,
                               MPI_Info info, int *npset_names) ;
int QMPI_Session_get_pset_info(QMPI_Context context, int tool_id, MPI_Session session,
                               const char *pset_name, MPI_Info *info) ;
int QMPI_Session_init(QMPI_Context context, int tool_id, MPI_Info info, MPI_Errhandler errhandler,
                      MPI_Session *session) ;
MPI_Aint QMPI_Aint_add(QMPI_Context context, int tool_id, MPI_Aint base, MPI_Aint disp)
                    ;
MPI_Aint QMPI_Aint_diff(QMPI_Context context, int tool_id, MPI_Aint addr1, MPI_Aint addr2)
                    ;
int QMPI_Get_library_version(QMPI_Context context, int tool_id, char *version, int *resultlen)
                    ;
int QMPI_Get_processor_name(QMPI_Context context, int tool_id, char *name, int *resultlen)
                    ;
int QMPI_Get_version(QMPI_Context context, int tool_id, int *version, int *subversion)
                    ;
int QMPI_Pcontrol(QMPI_Context context, int tool_id, const int level, ...) ;
int QMPIX_GPU_query_support(QMPI_Context context, int tool_id, int gpu_type, int *is_supported)
                    ;
int QMPIX_Query_cuda_support(QMPI_Context context, int tool_id) ;
int QMPIX_Query_ze_support(QMPI_Context context, int tool_id) ;
int QMPIX_Query_hip_support(QMPI_Context context, int tool_id) ;
int QMPI_T_category_changed(QMPI_Context context, int tool_id, int *update_number)
                    ;
int QMPI_T_category_get_categories(QMPI_Context context, int tool_id, int cat_index, int len,
                                   int indices[]) ;
int QMPI_T_category_get_cvars(QMPI_Context context, int tool_id, int cat_index, int len,
                              int indices[]) ;
int QMPI_T_category_get_events(QMPI_Context context, int tool_id, int cat_index, int len,
                               int indices[]) ;
int QMPI_T_category_get_index(QMPI_Context context, int tool_id, const char *name, int *cat_index)
                    ;
int QMPI_T_category_get_info(QMPI_Context context, int tool_id, int cat_index, char *name,
                             int *name_len, char *desc, int *desc_len, int *num_cvars,
                             int *num_pvars, int *num_categories) ;
int QMPI_T_category_get_num(QMPI_Context context, int tool_id, int *num_cat) ;
int QMPI_T_category_get_num_events(QMPI_Context context, int tool_id, int cat_index,
                                   int *num_events) ;
int QMPI_T_category_get_pvars(QMPI_Context context, int tool_id, int cat_index, int len,
                              int indices[]) ;
int QMPI_T_cvar_get_index(QMPI_Context context, int tool_id, const char *name, int *cvar_index)
                    ;
int QMPI_T_cvar_get_info(QMPI_Context context, int tool_id, int cvar_index, char *name,
                         int *name_len, int *verbosity, MPI_Datatype *datatype,
                         MPI_T_enum *enumtype, char *desc, int *desc_len, int *bind, int *scope)
                                         ;
int QMPI_T_cvar_get_num(QMPI_Context context, int tool_id, int *num_cvar) ;
int QMPI_T_cvar_handle_alloc(QMPI_Context context, int tool_id, int cvar_index, void *obj_handle,
                             MPI_T_cvar_handle *handle, int *count) ;
int QMPI_T_cvar_handle_free(QMPI_Context context, int tool_id, MPI_T_cvar_handle *handle)
                    ;
int QMPI_T_cvar_read(QMPI_Context context, int tool_id, MPI_T_cvar_handle handle, void *buf)
                    ;
int QMPI_T_cvar_write(QMPI_Context context, int tool_id, MPI_T_cvar_handle handle, const void *buf)
                    ;
int QMPI_T_enum_get_info(QMPI_Context context, int tool_id, MPI_T_enum enumtype, int *num,
                         char *name, int *name_len) ;
int QMPI_T_enum_get_item(QMPI_Context context, int tool_id, MPI_T_enum enumtype, int indx,
                         int *value, char *name, int *name_len) ;
int QMPI_T_event_callback_get_info(QMPI_Context context, int tool_id,
                                   MPI_T_event_registration event_registration,
                                   MPI_T_cb_safety cb_safety, MPI_Info *info_used)
                                                   ;
int QMPI_T_event_callback_set_info(QMPI_Context context, int tool_id,
                                   MPI_T_event_registration event_registration,
                                   MPI_T_cb_safety cb_safety, MPI_Info info) ;
int QMPI_T_event_copy(QMPI_Context context, int tool_id, MPI_T_event_instance event_instance,
                      void *buffer) ;
int QMPI_T_event_get_index(QMPI_Context context, int tool_id, const char *name, int *event_index)
                    ;
int QMPI_T_event_get_info(QMPI_Context context, int tool_id, int event_index, char *name,
                          int *name_len, int *verbosity, MPI_Datatype array_of_datatypes[],
                          MPI_Aint array_of_displacements[], int *num_elements,
                          MPI_T_enum *enumtype, MPI_Info *info, char *desc, int *desc_len,
                          int *bind) ;
int QMPI_T_event_get_num(QMPI_Context context, int tool_id, int *num_events) ;
int QMPI_T_event_get_source(QMPI_Context context, int tool_id, MPI_T_event_instance event_instance,
                            int *source_index) ;
int QMPI_T_event_get_timestamp(QMPI_Context context, int tool_id,
                               MPI_T_event_instance event_instance, MPI_Count *event_timestamp)
                                               ;
int QMPI_T_event_handle_alloc(QMPI_Context context, int tool_id, int event_index, void *obj_handle,
                              MPI_Info info, MPI_T_event_registration *event_registration)
                                              ;
int QMPI_T_event_handle_free(QMPI_Context context, int tool_id,
                             MPI_T_event_registration event_registration, void *user_data,
                             MPI_T_event_free_cb_function free_cb_function) ;
int QMPI_T_event_handle_get_info(QMPI_Context context, int tool_id,
                                 MPI_T_event_registration event_registration, MPI_Info *info_used)
                                                 ;
int QMPI_T_event_handle_set_info(QMPI_Context context, int tool_id,
                                 MPI_T_event_registration event_registration, MPI_Info info)
                                                 ;
int QMPI_T_event_read(QMPI_Context context, int tool_id, MPI_T_event_instance event_instance,
                      int element_index, void *buffer) ;
int QMPI_T_event_register_callback(QMPI_Context context, int tool_id,
                                   MPI_T_event_registration event_registration,
                                   MPI_T_cb_safety cb_safety, MPI_Info info, void *user_data,
                                   MPI_T_event_cb_function event_cb_function) ;
int QMPI_T_event_set_dropped_handler(QMPI_Context context, int tool_id,
                                     MPI_T_event_registration event_registration,
                                     MPI_T_event_dropped_cb_function dropped_cb_function)
                                                     ;
int QMPI_T_finalize(QMPI_Context context, int tool_id) ;
int QMPI_T_init_thread(QMPI_Context context, int tool_id, int required, int *provided)
                    ;
int QMPI_T_pvar_get_index(QMPI_Context context, int tool_id, const char *name, int var_class,
                          int *pvar_index) ;
int QMPI_T_pvar_get_info(QMPI_Context context, int tool_id, int pvar_index, char *name,
                         int *name_len, int *verbosity, int *var_class, MPI_Datatype *datatype,
                         MPI_T_enum *enumtype, char *desc, int *desc_len, int *bind, int *readonly,
                         int *continuous, int *atomic) ;
int QMPI_T_pvar_get_num(QMPI_Context context, int tool_id, int *num_pvar) ;
int QMPI_T_pvar_handle_alloc(QMPI_Context context, int tool_id, MPI_T_pvar_session session,
                             int pvar_index, void *obj_handle, MPI_T_pvar_handle *handle,
                             int *count) ;
int QMPI_T_pvar_handle_free(QMPI_Context context, int tool_id, MPI_T_pvar_session session,
                            MPI_T_pvar_handle *handle) ;
int QMPI_T_pvar_read(QMPI_Context context, int tool_id, MPI_T_pvar_session session,
                     MPI_T_pvar_handle handle, void *buf) ;
int QMPI_T_pvar_readreset(QMPI_Context context, int tool_id, MPI_T_pvar_session session,
                          MPI_T_pvar_handle handle, void *buf) ;
int QMPI_T_pvar_reset(QMPI_Context context, int tool_id, MPI_T_pvar_session session,
                      MPI_T_pvar_handle handle) ;
int QMPI_T_pvar_session_create(QMPI_Context context, int tool_id, MPI_T_pvar_session *session)
                    ;
int QMPI_T_pvar_session_free(QMPI_Context context, int tool_id, MPI_T_pvar_session *session)
                    ;
int QMPI_T_pvar_start(QMPI_Context context, int tool_id, MPI_T_pvar_session session,
                      MPI_T_pvar_handle handle) ;
int QMPI_T_pvar_stop(QMPI_Context context, int tool_id, MPI_T_pvar_session session,
                     MPI_T_pvar_handle handle) ;
int QMPI_T_pvar_write(QMPI_Context context, int tool_id, MPI_T_pvar_session session,
                      MPI_T_pvar_handle handle, const void *buf) ;
int QMPI_T_source_get_info(QMPI_Context context, int tool_id, int source_index, char *name,
                           int *name_len, char *desc, int *desc_len, MPI_T_source_order *ordering,
                           MPI_Count *ticks_per_second, MPI_Count *max_ticks, MPI_Info *info)
                                           ;
int QMPI_T_source_get_num(QMPI_Context context, int tool_id, int *num_sources) ;
int QMPI_T_source_get_timestamp(QMPI_Context context, int tool_id, int source_index,
                                MPI_Count *timestamp) ;
int QMPI_Op_commutative(QMPI_Context context, int tool_id, MPI_Op op, int *commute)
                    ;
int QMPI_Op_create(QMPI_Context context, int tool_id, MPI_User_function *user_fn, int commute,
                   MPI_Op *op) ;
int QMPI_Op_create_c(QMPI_Context context, int tool_id, MPI_User_function_c *user_fn, int commute,
                     MPI_Op *op) ;
int QMPI_Op_free(QMPI_Context context, int tool_id, MPI_Op *op) ;
int QMPI_Parrived(QMPI_Context context, int tool_id, MPI_Request request, int partition, int *flag)
                    ;
int QMPI_Pready(QMPI_Context context, int tool_id, int partition, MPI_Request request)
                    ;
int QMPI_Pready_list(QMPI_Context context, int tool_id, int length, const int array_of_partitions[],
                     MPI_Request request) ;
int QMPI_Pready_range(QMPI_Context context, int tool_id, int partition_low, int partition_high,
                      MPI_Request request) ;
int QMPI_Precv_init(QMPI_Context context, int tool_id, void *buf, int partitions, MPI_Count count,
                    MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Info info,
                    MPI_Request *request) ;
int QMPI_Psend_init(QMPI_Context context, int tool_id, const void *buf, int partitions,
                    MPI_Count count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm,
                    MPI_Info info, MPI_Request *request)
                                                                          ;
int QMPI_Bsend(QMPI_Context context, int tool_id, const void *buf, int count, MPI_Datatype datatype,
               int dest, int tag, MPI_Comm comm)
                                                                     ;
int QMPI_Bsend_c(QMPI_Context context, int tool_id, const void *buf, MPI_Count count,
                 MPI_Datatype datatype, int dest, int tag, MPI_Comm comm)
                                                                       ;
int QMPI_Bsend_init(QMPI_Context context, int tool_id, const void *buf, int count,
                    MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Request *request)
                                                                          ;
int QMPI_Bsend_init_c(QMPI_Context context, int tool_id, const void *buf, MPI_Count count,
                      MPI_Datatype datatype, int dest, int tag, MPI_Comm comm,
                      MPI_Request *request) ;
int QMPI_Buffer_attach(QMPI_Context context, int tool_id, void *buffer, int size) ;
int QMPI_Buffer_attach_c(QMPI_Context context, int tool_id, void *buffer, MPI_Count size)
                    ;
int QMPI_Buffer_detach(QMPI_Context context, int tool_id, void *buffer_addr, int *size)
                    ;
int QMPI_Buffer_detach_c(QMPI_Context context, int tool_id, void *buffer_addr, MPI_Count *size)
                    ;
int QMPI_Buffer_flush(QMPI_Context context, int tool_id) ;
int QMPI_Buffer_iflush(QMPI_Context context, int tool_id, MPI_Request *request) ;
int QMPI_Comm_attach_buffer(QMPI_Context context, int tool_id, MPI_Comm comm, void *buffer,
                            int size) ;
int QMPI_Comm_attach_buffer_c(QMPI_Context context, int tool_id, MPI_Comm comm, void *buffer,
                              MPI_Count size) ;
int QMPI_Comm_detach_buffer(QMPI_Context context, int tool_id, MPI_Comm comm, void *buffer_addr,
                            int *size) ;
int QMPI_Comm_detach_buffer_c(QMPI_Context context, int tool_id, MPI_Comm comm, void *buffer_addr,
                              MPI_Count *size) ;
int QMPI_Comm_flush_buffer(QMPI_Context context, int tool_id, MPI_Comm comm) ;
int QMPI_Comm_iflush_buffer(QMPI_Context context, int tool_id, MPI_Comm comm, MPI_Request *request)
                    ;
int QMPI_Ibsend(QMPI_Context context, int tool_id, const void *buf, int count,
                MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Request *request)
                                                                      ;
int QMPI_Ibsend_c(QMPI_Context context, int tool_id, const void *buf, MPI_Count count,
                  MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Request *request)
                                                                        ;
int QMPI_Improbe(QMPI_Context context, int tool_id, int source, int tag, MPI_Comm comm, int *flag,
                 MPI_Message *message, MPI_Status *status) ;
int QMPI_Imrecv(QMPI_Context context, int tool_id, void *buf, int count, MPI_Datatype datatype,
                MPI_Message *message, MPI_Request *request)
                                                                      ;
int QMPI_Imrecv_c(QMPI_Context context, int tool_id, void *buf, MPI_Count count,
                  MPI_Datatype datatype, MPI_Message *message, MPI_Request *request)
                                                                        ;
int QMPI_Iprobe(QMPI_Context context, int tool_id, int source, int tag, MPI_Comm comm, int *flag,
                MPI_Status *status) ;
int QMPI_Irecv(QMPI_Context context, int tool_id, void *buf, int count, MPI_Datatype datatype,
               int source, int tag, MPI_Comm comm, MPI_Request *request)
                                                                     ;
int QMPI_Irecv_c(QMPI_Context context, int tool_id, void *buf, MPI_Count count,
                 MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Request *request)
                                                                       ;
int QMPI_Irsend(QMPI_Context context, int tool_id, const void *buf, int count,
                MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Request *request)
                                                                      ;
int QMPI_Irsend_c(QMPI_Context context, int tool_id, const void *buf, MPI_Count count,
                  MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Request *request)
                                                                        ;
int QMPI_Isend(QMPI_Context context, int tool_id, const void *buf, int count, MPI_Datatype datatype,
               int dest, int tag, MPI_Comm comm, MPI_Request *request)
                                                                     ;
int QMPI_Isend_c(QMPI_Context context, int tool_id, const void *buf, MPI_Count count,
                 MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Request *request)
                                                                       ;
int QMPI_Isendrecv(QMPI_Context context, int tool_id, const void *sendbuf, int sendcount,
                   MPI_Datatype sendtype, int dest, int sendtag, void *recvbuf, int recvcount,
                   MPI_Datatype recvtype, int source, int recvtag, MPI_Comm comm,
                   MPI_Request *request)
                                                                                                                ;
int QMPI_Isendrecv_c(QMPI_Context context, int tool_id, const void *sendbuf, MPI_Count sendcount,
                     MPI_Datatype sendtype, int dest, int sendtag, void *recvbuf,
                     MPI_Count recvcount, MPI_Datatype recvtype, int source, int recvtag,
                     MPI_Comm comm, MPI_Request *request)
                                                                                                                  ;
int QMPI_Isendrecv_replace(QMPI_Context context, int tool_id, void *buf, int count,
                           MPI_Datatype datatype, int dest, int sendtag, int source, int recvtag,
                           MPI_Comm comm, MPI_Request *request)
                                                                                 ;
int QMPI_Isendrecv_replace_c(QMPI_Context context, int tool_id, void *buf, MPI_Count count,
                             MPI_Datatype datatype, int dest, int sendtag, int source, int recvtag,
                             MPI_Comm comm, MPI_Request *request)
                                                                                   ;
int QMPI_Issend(QMPI_Context context, int tool_id, const void *buf, int count,
                MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Request *request)
                                                                      ;
int QMPI_Issend_c(QMPI_Context context, int tool_id, const void *buf, MPI_Count count,
                  MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Request *request)
                                                                        ;
int QMPI_Mprobe(QMPI_Context context, int tool_id, int source, int tag, MPI_Comm comm,
                MPI_Message *message, MPI_Status *status) ;
int QMPI_Mrecv(QMPI_Context context, int tool_id, void *buf, int count, MPI_Datatype datatype,
               MPI_Message *message, MPI_Status *status)
                                                                     ;
int QMPI_Mrecv_c(QMPI_Context context, int tool_id, void *buf, MPI_Count count,
                 MPI_Datatype datatype, MPI_Message *message, MPI_Status *status)
                                                                       ;
int QMPI_Probe(QMPI_Context context, int tool_id, int source, int tag, MPI_Comm comm,
               MPI_Status *status) ;
int QMPI_Recv(QMPI_Context context, int tool_id, void *buf, int count, MPI_Datatype datatype,
              int source, int tag, MPI_Comm comm, MPI_Status *status)
                                                                    ;
int QMPI_Recv_c(QMPI_Context context, int tool_id, void *buf, MPI_Count count,
                MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Status *status)
                                                                      ;
int QMPI_Recv_init(QMPI_Context context, int tool_id, void *buf, int count, MPI_Datatype datatype,
                   int source, int tag, MPI_Comm comm, MPI_Request *request)
                                                                         ;
int QMPI_Recv_init_c(QMPI_Context context, int tool_id, void *buf, MPI_Count count,
                     MPI_Datatype datatype, int source, int tag, MPI_Comm comm,
                     MPI_Request *request) ;
int QMPI_Rsend(QMPI_Context context, int tool_id, const void *buf, int count, MPI_Datatype datatype,
               int dest, int tag, MPI_Comm comm)
                                                                     ;
int QMPI_Rsend_c(QMPI_Context context, int tool_id, const void *buf, MPI_Count count,
                 MPI_Datatype datatype, int dest, int tag, MPI_Comm comm)
                                                                       ;
int QMPI_Rsend_init(QMPI_Context context, int tool_id, const void *buf, int count,
                    MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Request *request)
                                                                          ;
int QMPI_Rsend_init_c(QMPI_Context context, int tool_id, const void *buf, MPI_Count count,
                      MPI_Datatype datatype, int dest, int tag, MPI_Comm comm,
                      MPI_Request *request) ;
int QMPI_Send(QMPI_Context context, int tool_id, const void *buf, int count, MPI_Datatype datatype,
              int dest, int tag, MPI_Comm comm)
                                                                    ;
int QMPI_Send_c(QMPI_Context context, int tool_id, const void *buf, MPI_Count count,
                MPI_Datatype datatype, int dest, int tag, MPI_Comm comm)
                                                                      ;
int QMPI_Send_init(QMPI_Context context, int tool_id, const void *buf, int count,
                   MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Request *request)
                                                                         ;
int QMPI_Send_init_c(QMPI_Context context, int tool_id, const void *buf, MPI_Count count,
                     MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Request *request)
                                                                           ;
int QMPI_Sendrecv(QMPI_Context context, int tool_id, const void *sendbuf, int sendcount,
                  MPI_Datatype sendtype, int dest, int sendtag, void *recvbuf, int recvcount,
                  MPI_Datatype recvtype, int source, int recvtag, MPI_Comm comm,
                  MPI_Status *status)
                                                                                                               ;
int QMPI_Sendrecv_c(QMPI_Context context, int tool_id, const void *sendbuf, MPI_Count sendcount,
                    MPI_Datatype sendtype, int dest, int sendtag, void *recvbuf,
                    MPI_Count recvcount, MPI_Datatype recvtype, int source, int recvtag,
                    MPI_Comm comm, MPI_Status *status)
                                                                                                                 ;
int QMPI_Sendrecv_replace(QMPI_Context context, int tool_id, void *buf, int count,
                          MPI_Datatype datatype, int dest, int sendtag, int source, int recvtag,
                          MPI_Comm comm, MPI_Status *status)
                                                                                ;
int QMPI_Sendrecv_replace_c(QMPI_Context context, int tool_id, void *buf, MPI_Count count,
                            MPI_Datatype datatype, int dest, int sendtag, int source, int recvtag,
                            MPI_Comm comm, MPI_Status *status)
                                                                                  ;
int QMPI_Session_attach_buffer(QMPI_Context context, int tool_id, MPI_Session session, void *buffer,
                               int size) ;
int QMPI_Session_attach_buffer_c(QMPI_Context context, int tool_id, MPI_Session session,
                                 void *buffer, MPI_Count size) ;
int QMPI_Session_detach_buffer(QMPI_Context context, int tool_id, MPI_Session session,
                               void *buffer_addr, int *size) ;
int QMPI_Session_detach_buffer_c(QMPI_Context context, int tool_id, MPI_Session session,
                                 void *buffer_addr, MPI_Count *size) ;
int QMPI_Session_flush_buffer(QMPI_Context context, int tool_id, MPI_Session session)
                    ;
int QMPI_Session_iflush_buffer(QMPI_Context context, int tool_id, MPI_Session session,
                               MPI_Request *request) ;
int QMPI_Ssend(QMPI_Context context, int tool_id, const void *buf, int count, MPI_Datatype datatype,
               int dest, int tag, MPI_Comm comm)
                                                                     ;
int QMPI_Ssend_c(QMPI_Context context, int tool_id, const void *buf, MPI_Count count,
                 MPI_Datatype datatype, int dest, int tag, MPI_Comm comm)
                                                                       ;
int QMPI_Ssend_init(QMPI_Context context, int tool_id, const void *buf, int count,
                    MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Request *request)
                                                                          ;
int QMPI_Ssend_init_c(QMPI_Context context, int tool_id, const void *buf, MPI_Count count,
                      MPI_Datatype datatype, int dest, int tag, MPI_Comm comm,
                      MPI_Request *request) ;
int QMPI_Cancel(QMPI_Context context, int tool_id, MPI_Request *request) ;
int QMPI_Grequest_complete(QMPI_Context context, int tool_id, MPI_Request request)
                    ;
int QMPI_Grequest_start(QMPI_Context context, int tool_id, MPI_Grequest_query_function *query_fn,
                        MPI_Grequest_free_function *free_fn,
                        MPI_Grequest_cancel_function *cancel_fn, void *extra_state,
                        MPI_Request *request) ;
int QMPI_Request_free(QMPI_Context context, int tool_id, MPI_Request *request) ;
int QMPI_Request_get_status(QMPI_Context context, int tool_id, MPI_Request request, int *flag,
                            MPI_Status *status) ;
int QMPI_Request_get_status_all(QMPI_Context context, int tool_id, int count,
                                const MPI_Request array_of_requests[], int *flag,
                                MPI_Status *array_of_statuses) ;
int QMPI_Request_get_status_any(QMPI_Context context, int tool_id, int count,
                                const MPI_Request array_of_requests[], int *indx, int *flag,
                                MPI_Status *status) ;
int QMPI_Request_get_status_some(QMPI_Context context, int tool_id, int incount,
                                 const MPI_Request array_of_requests[], int *outcount,
                                 int array_of_indices[], MPI_Status *array_of_statuses)
                                                 ;
int QMPI_Start(QMPI_Context context, int tool_id, MPI_Request *request) ;
int QMPI_Startall(QMPI_Context context, int tool_id, int count, MPI_Request array_of_requests[])
                    ;
int QMPI_Status_get_error(QMPI_Context context, int tool_id, const MPI_Status *status, int *error)
                    ;
int QMPI_Status_get_source(QMPI_Context context, int tool_id, const MPI_Status *status,
                           int *source) ;
int QMPI_Status_get_tag(QMPI_Context context, int tool_id, const MPI_Status *status, int *tag)
                    ;
int QMPI_Status_set_error(QMPI_Context context, int tool_id, MPI_Status *status, int error)
                    ;
int QMPI_Status_set_source(QMPI_Context context, int tool_id, MPI_Status *status, int source)
                    ;
int QMPI_Status_set_tag(QMPI_Context context, int tool_id, MPI_Status *status, int tag)
                    ;
int QMPI_Status_set_cancelled(QMPI_Context context, int tool_id, MPI_Status *status, int flag)
                    ;
int QMPI_Test(QMPI_Context context, int tool_id, MPI_Request *request, int *flag,
              MPI_Status *status) ;
int QMPI_Test_cancelled(QMPI_Context context, int tool_id, const MPI_Status *status, int *flag)
                    ;
int QMPI_Testall(QMPI_Context context, int tool_id, int count, MPI_Request array_of_requests[],
                 int *flag, MPI_Status *array_of_statuses) ;
int QMPI_Testany(QMPI_Context context, int tool_id, int count, MPI_Request array_of_requests[],
                 int *indx, int *flag, MPI_Status *status) ;
int QMPI_Testsome(QMPI_Context context, int tool_id, int incount, MPI_Request array_of_requests[],
                  int *outcount, int array_of_indices[], MPI_Status *array_of_statuses)
                                  ;
int QMPI_Wait(QMPI_Context context, int tool_id, MPI_Request *request, MPI_Status *status)
                    ;
int QMPI_Waitall(QMPI_Context context, int tool_id, int count, MPI_Request array_of_requests[],
                 MPI_Status *array_of_statuses) ;
int QMPI_Waitany(QMPI_Context context, int tool_id, int count, MPI_Request array_of_requests[],
                 int *indx, MPI_Status *status) ;
int QMPI_Waitsome(QMPI_Context context, int tool_id, int incount, MPI_Request array_of_requests[],
                  int *outcount, int array_of_indices[], MPI_Status *array_of_statuses)
                                  ;
int QMPIX_Grequest_start(QMPI_Context context, int tool_id, MPI_Grequest_query_function *query_fn,
                         MPI_Grequest_free_function *free_fn,
                         MPI_Grequest_cancel_function *cancel_fn,
                         MPIX_Grequest_poll_function *poll_fn, MPIX_Grequest_wait_function *wait_fn,
                         void *extra_state, MPI_Request *request) ;
int QMPIX_Grequest_class_create(QMPI_Context context, int tool_id,
                                MPI_Grequest_query_function *query_fn,
                                MPI_Grequest_free_function *free_fn,
                                MPI_Grequest_cancel_function *cancel_fn,
                                MPIX_Grequest_poll_function *poll_fn,
                                MPIX_Grequest_wait_function *wait_fn,
                                MPIX_Grequest_class *greq_class) ;
int QMPIX_Grequest_class_allocate(QMPI_Context context, int tool_id, MPIX_Grequest_class greq_class,
                                  void *extra_state, MPI_Request *request) ;
int QMPIX_Request_is_complete(QMPI_Context context, int tool_id, MPI_Request request)
                    ;
int QMPI_Accumulate(QMPI_Context context, int tool_id, const void *origin_addr, int origin_count,
                    MPI_Datatype origin_datatype, int target_rank, MPI_Aint target_disp,
                    int target_count, MPI_Datatype target_datatype, MPI_Op op, MPI_Win win)
                                                                          ;
int QMPI_Accumulate_c(QMPI_Context context, int tool_id, const void *origin_addr,
                      MPI_Count origin_count, MPI_Datatype origin_datatype, int target_rank,
                      MPI_Aint target_disp, MPI_Count target_count, MPI_Datatype target_datatype,
                      MPI_Op op, MPI_Win win)
                                                                            ;
int QMPI_Alloc_mem(QMPI_Context context, int tool_id, MPI_Aint size, MPI_Info info, void *baseptr)
                    ;
int QMPI_Compare_and_swap(QMPI_Context context, int tool_id, const void *origin_addr,
                          const void *compare_addr, void *result_addr, MPI_Datatype datatype,
                          int target_rank, MPI_Aint target_disp, MPI_Win win) ;
int QMPI_Fetch_and_op(QMPI_Context context, int tool_id, const void *origin_addr, void *result_addr,
                      MPI_Datatype datatype, int target_rank, MPI_Aint target_disp, MPI_Op op,
                      MPI_Win win) ;
int QMPI_Free_mem(QMPI_Context context, int tool_id, void *base) ;
int QMPI_Get(QMPI_Context context, int tool_id, void *origin_addr, int origin_count,
             MPI_Datatype origin_datatype, int target_rank, MPI_Aint target_disp, int target_count,
             MPI_Datatype target_datatype, MPI_Win win)
                                                                   ;
int QMPI_Get_c(QMPI_Context context, int tool_id, void *origin_addr, MPI_Count origin_count,
               MPI_Datatype origin_datatype, int target_rank, MPI_Aint target_disp,
               MPI_Count target_count, MPI_Datatype target_datatype, MPI_Win win)
                                                                     ;
int QMPI_Get_accumulate(QMPI_Context context, int tool_id, const void *origin_addr,
                        int origin_count, MPI_Datatype origin_datatype, void *result_addr,
                        int result_count, MPI_Datatype result_datatype, int target_rank,
                        MPI_Aint target_disp, int target_count, MPI_Datatype target_datatype,
                        MPI_Op op, MPI_Win win)
                                                                                                                    ;
int QMPI_Get_accumulate_c(QMPI_Context context, int tool_id, const void *origin_addr,
                          MPI_Count origin_count, MPI_Datatype origin_datatype, void *result_addr,
                          MPI_Count result_count, MPI_Datatype result_datatype, int target_rank,
                          MPI_Aint target_disp, MPI_Count target_count,
                          MPI_Datatype target_datatype, MPI_Op op, MPI_Win win)
                                                                                                                      ;
int QMPI_Put(QMPI_Context context, int tool_id, const void *origin_addr, int origin_count,
             MPI_Datatype origin_datatype, int target_rank, MPI_Aint target_disp, int target_count,
             MPI_Datatype target_datatype, MPI_Win win)
                                                                   ;
int QMPI_Put_c(QMPI_Context context, int tool_id, const void *origin_addr, MPI_Count origin_count,
               MPI_Datatype origin_datatype, int target_rank, MPI_Aint target_disp,
               MPI_Count target_count, MPI_Datatype target_datatype, MPI_Win win)
                                                                     ;
int QMPI_Raccumulate(QMPI_Context context, int tool_id, const void *origin_addr, int origin_count,
                     MPI_Datatype origin_datatype, int target_rank, MPI_Aint target_disp,
                     int target_count, MPI_Datatype target_datatype, MPI_Op op, MPI_Win win,
                     MPI_Request *request) ;
int QMPI_Raccumulate_c(QMPI_Context context, int tool_id, const void *origin_addr,
                       MPI_Count origin_count, MPI_Datatype origin_datatype, int target_rank,
                       MPI_Aint target_disp, MPI_Count target_count, MPI_Datatype target_datatype,
                       MPI_Op op, MPI_Win win, MPI_Request *request)
                                                                             ;
int QMPI_Rget(QMPI_Context context, int tool_id, void *origin_addr, int origin_count,
              MPI_Datatype origin_datatype, int target_rank, MPI_Aint target_disp, int target_count,
              MPI_Datatype target_datatype, MPI_Win win, MPI_Request *request)
                                                                    ;
int QMPI_Rget_c(QMPI_Context context, int tool_id, void *origin_addr, MPI_Count origin_count,
                MPI_Datatype origin_datatype, int target_rank, MPI_Aint target_disp,
                MPI_Count target_count, MPI_Datatype target_datatype, MPI_Win win,
                MPI_Request *request) ;
int QMPI_Rget_accumulate(QMPI_Context context, int tool_id, const void *origin_addr,
                         int origin_count, MPI_Datatype origin_datatype, void *result_addr,
                         int result_count, MPI_Datatype result_datatype, int target_rank,
                         MPI_Aint target_disp, int target_count, MPI_Datatype target_datatype,
                         MPI_Op op, MPI_Win win, MPI_Request *request)
                                                                                                                     ;
int QMPI_Rget_accumulate_c(QMPI_Context context, int tool_id, const void *origin_addr,
                           MPI_Count origin_count, MPI_Datatype origin_datatype, void *result_addr,
                           MPI_Count result_count, MPI_Datatype result_datatype, int target_rank,
                           MPI_Aint target_disp, MPI_Count target_count,
                           MPI_Datatype target_datatype, MPI_Op op, MPI_Win win,
                           MPI_Request *request)
                                                                                                                       ;
int QMPI_Rput(QMPI_Context context, int tool_id, const void *origin_addr, int origin_count,
              MPI_Datatype origin_datatype, int target_rank, MPI_Aint target_disp, int target_count,
              MPI_Datatype target_datatype, MPI_Win win, MPI_Request *request)
                                                                    ;
int QMPI_Rput_c(QMPI_Context context, int tool_id, const void *origin_addr, MPI_Count origin_count,
                MPI_Datatype origin_datatype, int target_rank, MPI_Aint target_disp,
                MPI_Count target_count, MPI_Datatype target_datatype, MPI_Win win,
                MPI_Request *request) ;
int QMPI_Win_allocate(QMPI_Context context, int tool_id, MPI_Aint size, int disp_unit,
                      MPI_Info info, MPI_Comm comm, void *baseptr, MPI_Win *win) ;
int QMPI_Win_allocate_c(QMPI_Context context, int tool_id, MPI_Aint size, MPI_Aint disp_unit,
                        MPI_Info info, MPI_Comm comm, void *baseptr, MPI_Win *win)
                                        ;
int QMPI_Win_allocate_shared(QMPI_Context context, int tool_id, MPI_Aint size, int disp_unit,
                             MPI_Info info, MPI_Comm comm, void *baseptr, MPI_Win *win)
                                             ;
int QMPI_Win_allocate_shared_c(QMPI_Context context, int tool_id, MPI_Aint size, MPI_Aint disp_unit,
                               MPI_Info info, MPI_Comm comm, void *baseptr, MPI_Win *win)
                                               ;
int QMPI_Win_attach(QMPI_Context context, int tool_id, MPI_Win win, void *base, MPI_Aint size)
                    ;
int QMPI_Win_complete(QMPI_Context context, int tool_id, MPI_Win win) ;
int QMPI_Win_create(QMPI_Context context, int tool_id, void *base, MPI_Aint size, int disp_unit,
                    MPI_Info info, MPI_Comm comm, MPI_Win *win) ;
int QMPI_Win_create_c(QMPI_Context context, int tool_id, void *base, MPI_Aint size,
                      MPI_Aint disp_unit, MPI_Info info, MPI_Comm comm, MPI_Win *win)
                                      ;
int QMPI_Win_create_dynamic(QMPI_Context context, int tool_id, MPI_Info info, MPI_Comm comm,
                            MPI_Win *win) ;
int QMPI_Win_detach(QMPI_Context context, int tool_id, MPI_Win win, const void *base)
                    ;
int QMPI_Win_fence(QMPI_Context context, int tool_id, int assert, MPI_Win win) ;
int QMPI_Win_flush(QMPI_Context context, int tool_id, int rank, MPI_Win win) ;
int QMPI_Win_flush_all(QMPI_Context context, int tool_id, MPI_Win win) ;
int QMPI_Win_flush_local(QMPI_Context context, int tool_id, int rank, MPI_Win win)
                    ;
int QMPI_Win_flush_local_all(QMPI_Context context, int tool_id, MPI_Win win) ;
int QMPI_Win_free(QMPI_Context context, int tool_id, MPI_Win *win) ;
int QMPI_Win_get_group(QMPI_Context context, int tool_id, MPI_Win win, MPI_Group *group)
                    ;
int QMPI_Win_get_info(QMPI_Context context, int tool_id, MPI_Win win, MPI_Info *info_used)
                    ;
int QMPI_Win_get_name(QMPI_Context context, int tool_id, MPI_Win win, char *win_name,
                      int *resultlen) ;
int QMPI_Win_lock(QMPI_Context context, int tool_id, int lock_type, int rank, int assert,
                  MPI_Win win) ;
int QMPI_Win_lock_all(QMPI_Context context, int tool_id, int assert, MPI_Win win) ;
int QMPI_Win_post(QMPI_Context context, int tool_id, MPI_Group group, int assert, MPI_Win win)
                    ;
int QMPI_Win_set_info(QMPI_Context context, int tool_id, MPI_Win win, MPI_Info info)
                    ;
int QMPI_Win_set_name(QMPI_Context context, int tool_id, MPI_Win win, const char *win_name)
                    ;
int QMPI_Win_shared_query(QMPI_Context context, int tool_id, MPI_Win win, int rank, MPI_Aint *size,
                          int *disp_unit, void *baseptr) ;
int QMPI_Win_shared_query_c(QMPI_Context context, int tool_id, MPI_Win win, int rank,
                            MPI_Aint *size, MPI_Aint *disp_unit, void *baseptr) ;
int QMPI_Win_start(QMPI_Context context, int tool_id, MPI_Group group, int assert, MPI_Win win)
                    ;
int QMPI_Win_sync(QMPI_Context context, int tool_id, MPI_Win win) ;
int QMPI_Win_test(QMPI_Context context, int tool_id, MPI_Win win, int *flag) ;
int QMPI_Win_unlock(QMPI_Context context, int tool_id, int rank, MPI_Win win) ;
int QMPI_Win_unlock_all(QMPI_Context context, int tool_id, MPI_Win win) ;
int QMPI_Win_wait(QMPI_Context context, int tool_id, MPI_Win win) ;
int QMPI_Close_port(QMPI_Context context, int tool_id, const char *port_name) ;
int QMPI_Comm_accept(QMPI_Context context, int tool_id, const char *port_name, MPI_Info info,
                     int root, MPI_Comm comm, MPI_Comm *newcomm) ;
int QMPI_Comm_connect(QMPI_Context context, int tool_id, const char *port_name, MPI_Info info,
                      int root, MPI_Comm comm, MPI_Comm *newcomm) ;
int QMPI_Comm_disconnect(QMPI_Context context, int tool_id, MPI_Comm *comm) ;
int QMPI_Comm_get_parent(QMPI_Context context, int tool_id, MPI_Comm *parent) ;
int QMPI_Comm_join(QMPI_Context context, int tool_id, int fd, MPI_Comm *intercomm)
                    ;
int QMPI_Comm_spawn(QMPI_Context context, int tool_id, const char *command, char *argv[],
                    int maxprocs, MPI_Info info, int root, MPI_Comm comm, MPI_Comm *intercomm,
                    int array_of_errcodes[]) ;
int QMPI_Comm_spawn_multiple(QMPI_Context context, int tool_id, int count,
                             char *array_of_commands[], char **array_of_argv[],
                             const int array_of_maxprocs[], const MPI_Info array_of_info[],
                             int root, MPI_Comm comm, MPI_Comm *intercomm, int array_of_errcodes[])
                                             ;
int QMPI_Lookup_name(QMPI_Context context, int tool_id, const char *service_name, MPI_Info info,
                     char *port_name) ;
int QMPI_Open_port(QMPI_Context context, int tool_id, MPI_Info info, char *port_name)
                    ;
int QMPI_Publish_name(QMPI_Context context, int tool_id, const char *service_name, MPI_Info info,
                      const char *port_name) ;
int QMPI_Unpublish_name(QMPI_Context context, int tool_id, const char *service_name, MPI_Info info,
                        const char *port_name) ;
int QMPIX_Stream_create(QMPI_Context context, int tool_id, MPI_Info info, MPIX_Stream *stream)
                    ;
int QMPIX_Stream_free(QMPI_Context context, int tool_id, MPIX_Stream *stream) ;
int QMPIX_Stream_comm_create(QMPI_Context context, int tool_id, MPI_Comm comm, MPIX_Stream stream,
                             MPI_Comm *newcomm) ;
int QMPIX_Stream_comm_create_multiplex(QMPI_Context context, int tool_id, MPI_Comm comm, int count,
                                       MPIX_Stream array_of_streams[], MPI_Comm *newcomm)
                                                       ;
int QMPIX_Comm_get_stream(QMPI_Context context, int tool_id, MPI_Comm comm, int idx,
                          MPIX_Stream *stream) ;
int QMPIX_Stream_progress(QMPI_Context context, int tool_id, MPIX_Stream stream) ;
int QMPIX_Start_progress_thread(QMPI_Context context, int tool_id, MPIX_Stream stream)
                    ;
int QMPIX_Stop_progress_thread(QMPI_Context context, int tool_id, MPIX_Stream stream)
                    ;
int QMPIX_Stream_send(QMPI_Context context, int tool_id, const void *buf, int count,
                      MPI_Datatype datatype, int dest, int tag, MPI_Comm comm,
                      int source_stream_index, int dest_stream_index)
                                                                            ;
int QMPIX_Stream_send_c(QMPI_Context context, int tool_id, const void *buf, MPI_Count count,
                        MPI_Datatype datatype, int dest, int tag, MPI_Comm comm,
                        int source_stream_index, int dest_stream_index)
                                                                              ;
int QMPIX_Stream_isend(QMPI_Context context, int tool_id, const void *buf, int count,
                       MPI_Datatype datatype, int dest, int tag, MPI_Comm comm,
                       int source_stream_index, int dest_stream_index, MPI_Request *request)
                                                                             ;
int QMPIX_Stream_isend_c(QMPI_Context context, int tool_id, const void *buf, MPI_Count count,
                         MPI_Datatype datatype, int dest, int tag, MPI_Comm comm,
                         int source_stream_index, int dest_stream_index, MPI_Request *request)
                                                                               ;
int QMPIX_Stream_recv(QMPI_Context context, int tool_id, void *buf, int count,
                      MPI_Datatype datatype, int source, int tag, MPI_Comm comm,
                      int source_stream_index, int dest_stream_index, MPI_Status *status)
                                                                            ;
int QMPIX_Stream_recv_c(QMPI_Context context, int tool_id, void *buf, MPI_Count count,
                        MPI_Datatype datatype, int source, int tag, MPI_Comm comm,
                        int source_stream_index, int dest_stream_index, MPI_Status *status)
                                                                              ;
int QMPIX_Stream_irecv(QMPI_Context context, int tool_id, void *buf, int count,
                       MPI_Datatype datatype, int source, int tag, MPI_Comm comm,
                       int source_stream_index, int dest_stream_index, MPI_Request *request)
                                                                             ;
int QMPIX_Stream_irecv_c(QMPI_Context context, int tool_id, void *buf, MPI_Count count,
                         MPI_Datatype datatype, int source, int tag, MPI_Comm comm,
                         int source_stream_index, int dest_stream_index, MPI_Request *request)
                                                                               ;
int QMPIX_Send_enqueue(QMPI_Context context, int tool_id, const void *buf, int count,
                       MPI_Datatype datatype, int dest, int tag, MPI_Comm comm)
                                                                             ;
int QMPIX_Send_enqueue_c(QMPI_Context context, int tool_id, const void *buf, MPI_Count count,
                         MPI_Datatype datatype, int dest, int tag, MPI_Comm comm)
                                                                               ;
int QMPIX_Recv_enqueue(QMPI_Context context, int tool_id, void *buf, int count,
                       MPI_Datatype datatype, int source, int tag, MPI_Comm comm,
                       MPI_Status *status) ;
int QMPIX_Recv_enqueue_c(QMPI_Context context, int tool_id, void *buf, MPI_Count count,
                         MPI_Datatype datatype, int source, int tag, MPI_Comm comm,
                         MPI_Status *status)
                                                                               ;
int QMPIX_Isend_enqueue(QMPI_Context context, int tool_id, const void *buf, int count,
                        MPI_Datatype datatype, int dest, int tag, MPI_Comm comm,
                        MPI_Request *request)
                                                                              ;
int QMPIX_Isend_enqueue_c(QMPI_Context context, int tool_id, const void *buf, MPI_Count count,
                          MPI_Datatype datatype, int dest, int tag, MPI_Comm comm,
                          MPI_Request *request)
                                                                                ;
int QMPIX_Irecv_enqueue(QMPI_Context context, int tool_id, void *buf, int count,
                        MPI_Datatype datatype, int source, int tag, MPI_Comm comm,
                        MPI_Request *request)
                                                                              ;
int QMPIX_Irecv_enqueue_c(QMPI_Context context, int tool_id, void *buf, MPI_Count count,
                          MPI_Datatype datatype, int source, int tag, MPI_Comm comm,
                          MPI_Request *request)
                                                                                ;
int QMPIX_Wait_enqueue(QMPI_Context context, int tool_id, MPI_Request *request, MPI_Status *status)
                    ;
int QMPIX_Waitall_enqueue(QMPI_Context context, int tool_id, int count,
                          MPI_Request array_of_requests[], MPI_Status *array_of_statuses)
                                          ;
int QMPIX_Allreduce_enqueue(QMPI_Context context, int tool_id, const void *sendbuf, void *recvbuf,
                            int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm)
                                                                                  ;
int QMPIX_Allreduce_enqueue_c(QMPI_Context context, int tool_id, const void *sendbuf, void *recvbuf,
                              MPI_Count count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm)
                                                                                    ;
int QMPIX_Threadcomm_init(QMPI_Context context, int tool_id, MPI_Comm comm, int num_threads,
                          MPI_Comm *newthreadcomm) ;
int QMPIX_Threadcomm_free(QMPI_Context context, int tool_id, MPI_Comm *threadcomm)
                    ;
int QMPIX_Threadcomm_start(QMPI_Context context, int tool_id, MPI_Comm threadcomm)
                    ;
int QMPIX_Threadcomm_finish(QMPI_Context context, int tool_id, MPI_Comm threadcomm)
                    ;
double QMPI_Wtick(QMPI_Context context, int tool_id) ;
double QMPI_Wtime(QMPI_Context context, int tool_id) ;
int QMPI_Cart_coords(QMPI_Context context, int tool_id, MPI_Comm comm, int rank, int maxdims,
                     int coords[]) ;
int QMPI_Cart_create(QMPI_Context context, int tool_id, MPI_Comm comm_old, int ndims,
                     const int dims[], const int periods[], int reorder, MPI_Comm *comm_cart)
                                     ;
int QMPI_Cart_get(QMPI_Context context, int tool_id, MPI_Comm comm, int maxdims, int dims[],
                  int periods[], int coords[]) ;
int QMPI_Cart_map(QMPI_Context context, int tool_id, MPI_Comm comm, int ndims, const int dims[],
                  const int periods[], int *newrank) ;
int QMPI_Cart_rank(QMPI_Context context, int tool_id, MPI_Comm comm, const int coords[], int *rank)
                    ;
int QMPI_Cart_shift(QMPI_Context context, int tool_id, MPI_Comm comm, int direction, int disp,
                    int *rank_source, int *rank_dest) ;
int QMPI_Cart_sub(QMPI_Context context, int tool_id, MPI_Comm comm, const int remain_dims[],
                  MPI_Comm *newcomm) ;
int QMPI_Cartdim_get(QMPI_Context context, int tool_id, MPI_Comm comm, int *ndims)
                    ;
int QMPI_Dims_create(QMPI_Context context, int tool_id, int nnodes, int ndims, int dims[])
                    ;
int QMPI_Dist_graph_create(QMPI_Context context, int tool_id, MPI_Comm comm_old, int n,
                           const int sources[], const int degrees[], const int destinations[],
                           const int weights[], MPI_Info info, int reorder,
                           MPI_Comm *comm_dist_graph) ;
int QMPI_Dist_graph_create_adjacent(QMPI_Context context, int tool_id, MPI_Comm comm_old,
                                    int indegree, const int sources[], const int sourceweights[],
                                    int outdegree, const int destinations[],
                                    const int destweights[], MPI_Info info, int reorder,
                                    MPI_Comm *comm_dist_graph) ;
int QMPI_Dist_graph_neighbors(QMPI_Context context, int tool_id, MPI_Comm comm, int maxindegree,
                              int sources[], int sourceweights[], int maxoutdegree,
                              int destinations[], int destweights[]) ;
int QMPI_Dist_graph_neighbors_count(QMPI_Context context, int tool_id, MPI_Comm comm, int *indegree,
                                    int *outdegree, int *weighted) ;
int QMPI_Get_hw_resource_info(QMPI_Context context, int tool_id, MPI_Info *hw_info)
                    ;
int QMPI_Graph_create(QMPI_Context context, int tool_id, MPI_Comm comm_old, int nnodes,
                      const int indx[], const int edges[], int reorder, MPI_Comm *comm_graph)
                                      ;
int QMPI_Graph_get(QMPI_Context context, int tool_id, MPI_Comm comm, int maxindex, int maxedges,
                   int indx[], int edges[]) ;
int QMPI_Graph_map(QMPI_Context context, int tool_id, MPI_Comm comm, int nnodes, const int indx[],
                   const int edges[], int *newrank) ;
int QMPI_Graph_neighbors(QMPI_Context context, int tool_id, MPI_Comm comm, int rank,
                         int maxneighbors, int neighbors[]) ;
int QMPI_Graph_neighbors_count(QMPI_Context context, int tool_id, MPI_Comm comm, int rank,
                               int *nneighbors) ;
int QMPI_Graphdims_get(QMPI_Context context, int tool_id, MPI_Comm comm, int *nnodes, int *nedges)
                    ;
int QMPI_Topo_test(QMPI_Context context, int tool_id, MPI_Comm comm, int *status) ;
MPI_Fint QMPI_File_c2f(QMPI_Context context, int tool_id, MPI_File file) ;
int QMPI_File_close(QMPI_Context context, int tool_id, MPI_File *fh) ;
int QMPI_File_delete(QMPI_Context context, int tool_id, const char *filename, MPI_Info info)
                    ;
MPI_File QMPI_File_f2c(QMPI_Context context, int tool_id, MPI_Fint file) ;
int QMPI_File_get_amode(QMPI_Context context, int tool_id, MPI_File fh, int *amode)
                    ;
int QMPI_File_get_atomicity(QMPI_Context context, int tool_id, MPI_File fh, int *flag)
                    ;
int QMPI_File_get_byte_offset(QMPI_Context context, int tool_id, MPI_File fh, MPI_Offset offset,
                              MPI_Offset *disp) ;
int QMPI_File_get_group(QMPI_Context context, int tool_id, MPI_File fh, MPI_Group *group)
                    ;
int QMPI_File_get_info(QMPI_Context context, int tool_id, MPI_File fh, MPI_Info *info_used)
                    ;
int QMPI_File_get_position(QMPI_Context context, int tool_id, MPI_File fh, MPI_Offset *offset)
                    ;
int QMPI_File_get_position_shared(QMPI_Context context, int tool_id, MPI_File fh,
                                  MPI_Offset *offset) ;
int QMPI_File_get_size(QMPI_Context context, int tool_id, MPI_File fh, MPI_Offset *size)
                    ;
int QMPI_File_get_type_extent(QMPI_Context context, int tool_id, MPI_File fh, MPI_Datatype datatype,
                              MPI_Aint *extent) ;
int QMPI_File_get_type_extent_c(QMPI_Context context, int tool_id, MPI_File fh,
                                MPI_Datatype datatype, MPI_Count *extent) ;
int QMPI_File_get_view(QMPI_Context context, int tool_id, MPI_File fh, MPI_Offset *disp,
                       MPI_Datatype *etype, MPI_Datatype *filetype, char *datarep)
                                       ;
int QMPI_File_iread(QMPI_Context context, int tool_id, MPI_File fh, void *buf, int count,
                    MPI_Datatype datatype, MPI_Request *request)
                                                                          ;
int QMPI_File_iread_c(QMPI_Context context, int tool_id, MPI_File fh, void *buf, MPI_Count count,
                      MPI_Datatype datatype, MPI_Request *request)
                                                                            ;
int QMPI_File_iread_all(QMPI_Context context, int tool_id, MPI_File fh, void *buf, int count,
                        MPI_Datatype datatype, MPI_Request *request)
                                                                              ;
int QMPI_File_iread_all_c(QMPI_Context context, int tool_id, MPI_File fh, void *buf,
                          MPI_Count count, MPI_Datatype datatype, MPI_Request *request)
                                                                                ;
int QMPI_File_iread_at(QMPI_Context context, int tool_id, MPI_File fh, MPI_Offset offset, void *buf,
                       int count, MPI_Datatype datatype, MPI_Request *request)
                                                                             ;
int QMPI_File_iread_at_c(QMPI_Context context, int tool_id, MPI_File fh, MPI_Offset offset,
                         void *buf, MPI_Count count, MPI_Datatype datatype, MPI_Request *request)
                                                                               ;
int QMPI_File_iread_at_all(QMPI_Context context, int tool_id, MPI_File fh, MPI_Offset offset,
                           void *buf, int count, MPI_Datatype datatype, MPI_Request *request)
                                                                                 ;
int QMPI_File_iread_at_all_c(QMPI_Context context, int tool_id, MPI_File fh, MPI_Offset offset,
                             void *buf, MPI_Count count, MPI_Datatype datatype,
                             MPI_Request *request)
                                                                                   ;
int QMPI_File_iread_shared(QMPI_Context context, int tool_id, MPI_File fh, void *buf, int count,
                           MPI_Datatype datatype, MPI_Request *request)
                                                                                 ;
int QMPI_File_iread_shared_c(QMPI_Context context, int tool_id, MPI_File fh, void *buf,
                             MPI_Count count, MPI_Datatype datatype, MPI_Request *request)
                                                                                   ;
int QMPI_File_iwrite(QMPI_Context context, int tool_id, MPI_File fh, const void *buf, int count,
                     MPI_Datatype datatype, MPI_Request *request)
                                                                           ;
int QMPI_File_iwrite_c(QMPI_Context context, int tool_id, MPI_File fh, const void *buf,
                       MPI_Count count, MPI_Datatype datatype, MPI_Request *request)
                                                                             ;
int QMPI_File_iwrite_all(QMPI_Context context, int tool_id, MPI_File fh, const void *buf, int count,
                         MPI_Datatype datatype, MPI_Request *request)
                                                                               ;
int QMPI_File_iwrite_all_c(QMPI_Context context, int tool_id, MPI_File fh, const void *buf,
                           MPI_Count count, MPI_Datatype datatype, MPI_Request *request)
                                                                                 ;
int QMPI_File_iwrite_at(QMPI_Context context, int tool_id, MPI_File fh, MPI_Offset offset,
                        const void *buf, int count, MPI_Datatype datatype, MPI_Request *request)
                                                                              ;
int QMPI_File_iwrite_at_c(QMPI_Context context, int tool_id, MPI_File fh, MPI_Offset offset,
                          const void *buf, MPI_Count count, MPI_Datatype datatype,
                          MPI_Request *request)
                                                                                ;
int QMPI_File_iwrite_at_all(QMPI_Context context, int tool_id, MPI_File fh, MPI_Offset offset,
                            const void *buf, int count, MPI_Datatype datatype,
                            MPI_Request *request)
                                                                                  ;
int QMPI_File_iwrite_at_all_c(QMPI_Context context, int tool_id, MPI_File fh, MPI_Offset offset,
                              const void *buf, MPI_Count count, MPI_Datatype datatype,
                              MPI_Request *request)
                                                                                    ;
int QMPI_File_iwrite_shared(QMPI_Context context, int tool_id, MPI_File fh, const void *buf,
                            int count, MPI_Datatype datatype, MPI_Request *request)
                                                                                  ;
int QMPI_File_iwrite_shared_c(QMPI_Context context, int tool_id, MPI_File fh, const void *buf,
                              MPI_Count count, MPI_Datatype datatype, MPI_Request *request)
                                                                                    ;
int QMPI_File_open(QMPI_Context context, int tool_id, MPI_Comm comm, const char *filename,
                   int amode, MPI_Info info, MPI_File *fh) ;
int QMPI_File_preallocate(QMPI_Context context, int tool_id, MPI_File fh, MPI_Offset size)
                    ;
int QMPI_File_read(QMPI_Context context, int tool_id, MPI_File fh, void *buf, int count,
                   MPI_Datatype datatype, MPI_Status *status)
                                                                         ;
int QMPI_File_read_c(QMPI_Context context, int tool_id, MPI_File fh, void *buf, MPI_Count count,
                     MPI_Datatype datatype, MPI_Status *status)
                                                                           ;
int QMPI_File_read_all(QMPI_Context context, int tool_id, MPI_File fh, void *buf, int count,
                       MPI_Datatype datatype, MPI_Status *status)
                                                                             ;
int QMPI_File_read_all_c(QMPI_Context context, int tool_id, MPI_File fh, void *buf, MPI_Count count,
                         MPI_Datatype datatype, MPI_Status *status)
                                                                               ;
int QMPI_File_read_all_begin(QMPI_Context context, int tool_id, MPI_File fh, void *buf, int count,
                             MPI_Datatype datatype)
                                                                                   ;
int QMPI_File_read_all_begin_c(QMPI_Context context, int tool_id, MPI_File fh, void *buf,
                               MPI_Count count, MPI_Datatype datatype)
                                                                                     ;
int QMPI_File_read_all_end(QMPI_Context context, int tool_id, MPI_File fh, void *buf,
                           MPI_Status *status) ;
int QMPI_File_read_at(QMPI_Context context, int tool_id, MPI_File fh, MPI_Offset offset, void *buf,
                      int count, MPI_Datatype datatype, MPI_Status *status)
                                                                            ;
int QMPI_File_read_at_c(QMPI_Context context, int tool_id, MPI_File fh, MPI_Offset offset,
                        void *buf, MPI_Count count, MPI_Datatype datatype, MPI_Status *status)
                                                                              ;
int QMPI_File_read_at_all(QMPI_Context context, int tool_id, MPI_File fh, MPI_Offset offset,
                          void *buf, int count, MPI_Datatype datatype, MPI_Status *status)
                                                                                ;
int QMPI_File_read_at_all_c(QMPI_Context context, int tool_id, MPI_File fh, MPI_Offset offset,
                            void *buf, MPI_Count count, MPI_Datatype datatype, MPI_Status *status)
                                                                                  ;
int QMPI_File_read_at_all_begin(QMPI_Context context, int tool_id, MPI_File fh, MPI_Offset offset,
                                void *buf, int count, MPI_Datatype datatype)
                                                                                      ;
int QMPI_File_read_at_all_begin_c(QMPI_Context context, int tool_id, MPI_File fh, MPI_Offset offset,
                                  void *buf, MPI_Count count, MPI_Datatype datatype)
                                                                                        ;
int QMPI_File_read_at_all_end(QMPI_Context context, int tool_id, MPI_File fh, void *buf,
                              MPI_Status *status) ;
int QMPI_File_read_ordered(QMPI_Context context, int tool_id, MPI_File fh, void *buf, int count,
                           MPI_Datatype datatype, MPI_Status *status)
                                                                                 ;
int QMPI_File_read_ordered_c(QMPI_Context context, int tool_id, MPI_File fh, void *buf,
                             MPI_Count count, MPI_Datatype datatype, MPI_Status *status)
                                                                                   ;
int QMPI_File_read_ordered_begin(QMPI_Context context, int tool_id, MPI_File fh, void *buf,
                                 int count, MPI_Datatype datatype)
                                                                                       ;
int QMPI_File_read_ordered_begin_c(QMPI_Context context, int tool_id, MPI_File fh, void *buf,
                                   MPI_Count count, MPI_Datatype datatype)
                                                                                         ;
int QMPI_File_read_ordered_end(QMPI_Context context, int tool_id, MPI_File fh, void *buf,
                               MPI_Status *status) ;
int QMPI_File_read_shared(QMPI_Context context, int tool_id, MPI_File fh, void *buf, int count,
                          MPI_Datatype datatype, MPI_Status *status)
                                                                                ;
int QMPI_File_read_shared_c(QMPI_Context context, int tool_id, MPI_File fh, void *buf,
                            MPI_Count count, MPI_Datatype datatype, MPI_Status *status)
                                                                                  ;
int QMPI_File_seek(QMPI_Context context, int tool_id, MPI_File fh, MPI_Offset offset, int whence)
                    ;
int QMPI_File_seek_shared(QMPI_Context context, int tool_id, MPI_File fh, MPI_Offset offset,
                          int whence) ;
int QMPI_File_set_atomicity(QMPI_Context context, int tool_id, MPI_File fh, int flag)
                    ;
int QMPI_File_set_info(QMPI_Context context, int tool_id, MPI_File fh, MPI_Info info)
                    ;
int QMPI_File_set_size(QMPI_Context context, int tool_id, MPI_File fh, MPI_Offset size)
                    ;
int QMPI_File_set_view(QMPI_Context context, int tool_id, MPI_File fh, MPI_Offset disp,
                       MPI_Datatype etype, MPI_Datatype filetype, const char *datarep,
                       MPI_Info info) ;
int QMPI_File_sync(QMPI_Context context, int tool_id, MPI_File fh) ;
int QMPI_File_write(QMPI_Context context, int tool_id, MPI_File fh, const void *buf, int count,
                    MPI_Datatype datatype, MPI_Status *status)
                                                                          ;
int QMPI_File_write_c(QMPI_Context context, int tool_id, MPI_File fh, const void *buf,
                      MPI_Count count, MPI_Datatype datatype, MPI_Status *status)
                                                                            ;
int QMPI_File_write_all(QMPI_Context context, int tool_id, MPI_File fh, const void *buf, int count,
                        MPI_Datatype datatype, MPI_Status *status)
                                                                              ;
int QMPI_File_write_all_c(QMPI_Context context, int tool_id, MPI_File fh, const void *buf,
                          MPI_Count count, MPI_Datatype datatype, MPI_Status *status)
                                                                                ;
int QMPI_File_write_all_begin(QMPI_Context context, int tool_id, MPI_File fh, const void *buf,
                              int count, MPI_Datatype datatype)
                                                                                    ;
int QMPI_File_write_all_begin_c(QMPI_Context context, int tool_id, MPI_File fh, const void *buf,
                                MPI_Count count, MPI_Datatype datatype)
                                                                                      ;
int QMPI_File_write_all_end(QMPI_Context context, int tool_id, MPI_File fh, const void *buf,
                            MPI_Status *status) ;
int QMPI_File_write_at(QMPI_Context context, int tool_id, MPI_File fh, MPI_Offset offset,
                       const void *buf, int count, MPI_Datatype datatype, MPI_Status *status)
                                                                             ;
int QMPI_File_write_at_c(QMPI_Context context, int tool_id, MPI_File fh, MPI_Offset offset,
                         const void *buf, MPI_Count count, MPI_Datatype datatype,
                         MPI_Status *status)
                                                                               ;
int QMPI_File_write_at_all(QMPI_Context context, int tool_id, MPI_File fh, MPI_Offset offset,
                           const void *buf, int count, MPI_Datatype datatype, MPI_Status *status)
                                                                                 ;
int QMPI_File_write_at_all_c(QMPI_Context context, int tool_id, MPI_File fh, MPI_Offset offset,
                             const void *buf, MPI_Count count, MPI_Datatype datatype,
                             MPI_Status *status)
                                                                                   ;
int QMPI_File_write_at_all_begin(QMPI_Context context, int tool_id, MPI_File fh, MPI_Offset offset,
                                 const void *buf, int count, MPI_Datatype datatype)
                                                                                       ;
int QMPI_File_write_at_all_begin_c(QMPI_Context context, int tool_id, MPI_File fh,
                                   MPI_Offset offset, const void *buf, MPI_Count count,
                                   MPI_Datatype datatype)
                                                                                         ;
int QMPI_File_write_at_all_end(QMPI_Context context, int tool_id, MPI_File fh, const void *buf,
                               MPI_Status *status) ;
int QMPI_File_write_ordered(QMPI_Context context, int tool_id, MPI_File fh, const void *buf,
                            int count, MPI_Datatype datatype, MPI_Status *status)
                                                                                  ;
int QMPI_File_write_ordered_c(QMPI_Context context, int tool_id, MPI_File fh, const void *buf,
                              MPI_Count count, MPI_Datatype datatype, MPI_Status *status)
                                                                                    ;
int QMPI_File_write_ordered_begin(QMPI_Context context, int tool_id, MPI_File fh, const void *buf,
                                  int count, MPI_Datatype datatype)
                                                                                        ;
int QMPI_File_write_ordered_begin_c(QMPI_Context context, int tool_id, MPI_File fh, const void *buf,
                                    MPI_Count count, MPI_Datatype datatype)
                                                                                          ;
int QMPI_File_write_ordered_end(QMPI_Context context, int tool_id, MPI_File fh, const void *buf,
                                MPI_Status *status) ;
int QMPI_File_write_shared(QMPI_Context context, int tool_id, MPI_File fh, const void *buf,
                           int count, MPI_Datatype datatype, MPI_Status *status)
                                                                                 ;
int QMPI_File_write_shared_c(QMPI_Context context, int tool_id, MPI_File fh, const void *buf,
                             MPI_Count count, MPI_Datatype datatype, MPI_Status *status)
                                                                                   ;
int QMPI_Register_datarep(QMPI_Context context, int tool_id, const char *datarep,
                          MPI_Datarep_conversion_function *read_conversion_fn,
                          MPI_Datarep_conversion_function *write_conversion_fn,
                          MPI_Datarep_extent_function *dtype_file_extent_fn, void *extra_state)
                                          ;
int QMPI_Register_datarep_c(QMPI_Context context, int tool_id, const char *datarep,
                            MPI_Datarep_conversion_function_c *read_conversion_fn,
                            MPI_Datarep_conversion_function_c *write_conversion_fn,
                            MPI_Datarep_extent_function *dtype_file_extent_fn, void *extra_state)
                                            ;
int QMPI_File_toint(QMPI_Context context, int tool_id, MPI_File file) ;
MPI_File QMPI_File_fromint(QMPI_Context context, int tool_id, int file) ;
typedef int (QMPI_Abi_get_fortran_booleans_t) (QMPI_Context context, int tool_id, int logical_size,
             void *logical_true, void *logical_false, int *is_set);
typedef int (QMPI_Abi_get_fortran_info_t) (QMPI_Context context, int tool_id, MPI_Info *info);
typedef int (QMPI_Abi_get_info_t) (QMPI_Context context, int tool_id, MPI_Info *info);
typedef int (QMPI_Abi_get_version_t) (QMPI_Context context, int tool_id, int *abi_major,
             int *abi_minor);
typedef int (QMPI_Abi_set_fortran_booleans_t) (QMPI_Context context, int tool_id, int logical_size,
             void *logical_true, void *logical_false);
typedef int (QMPI_Abi_set_fortran_info_t) (QMPI_Context context, int tool_id, MPI_Info info);
typedef int (QMPI_Comm_toint_t) (QMPI_Context context, int tool_id, MPI_Comm comm);
typedef MPI_Comm (QMPI_Comm_fromint_t) (QMPI_Context context, int tool_id, int comm);
typedef int (QMPI_Errhandler_toint_t) (QMPI_Context context, int tool_id,
             MPI_Errhandler errhandler);
typedef MPI_Errhandler (QMPI_Errhandler_fromint_t) (QMPI_Context context, int tool_id,
                        int errhandler);
typedef int (QMPI_Group_toint_t) (QMPI_Context context, int tool_id, MPI_Group group);
typedef MPI_Group (QMPI_Group_fromint_t) (QMPI_Context context, int tool_id, int group);
typedef int (QMPI_Info_toint_t) (QMPI_Context context, int tool_id, MPI_Info info);
typedef MPI_Info (QMPI_Info_fromint_t) (QMPI_Context context, int tool_id, int info);
typedef int (QMPI_Message_toint_t) (QMPI_Context context, int tool_id, MPI_Message message);
typedef MPI_Message (QMPI_Message_fromint_t) (QMPI_Context context, int tool_id, int message);
typedef int (QMPI_Op_toint_t) (QMPI_Context context, int tool_id, MPI_Op op);
typedef MPI_Op (QMPI_Op_fromint_t) (QMPI_Context context, int tool_id, int op);
typedef int (QMPI_Request_toint_t) (QMPI_Context context, int tool_id, MPI_Request request);
typedef MPI_Request (QMPI_Request_fromint_t) (QMPI_Context context, int tool_id, int request);
typedef int (QMPI_Session_toint_t) (QMPI_Context context, int tool_id, MPI_Session session);
typedef MPI_Session (QMPI_Session_fromint_t) (QMPI_Context context, int tool_id, int session);
typedef int (QMPI_Type_toint_t) (QMPI_Context context, int tool_id, MPI_Datatype datatype);
typedef MPI_Datatype (QMPI_Type_fromint_t) (QMPI_Context context, int tool_id, int datatype);
typedef int (QMPI_Win_toint_t) (QMPI_Context context, int tool_id, MPI_Win win);
typedef MPI_Win (QMPI_Win_fromint_t) (QMPI_Context context, int tool_id, int win);
typedef int (QMPIX_Async_start_t) (QMPI_Context context, int tool_id,
             MPIX_Async_poll_function *poll_fn, void *extra_state, MPIX_Stream stream);
typedef void * (QMPIX_Async_get_state_t) (QMPI_Context context, int tool_id,
                MPIX_Async_thing async_thing);
typedef int (QMPIX_Async_spawn_t) (QMPI_Context context, int tool_id, MPIX_Async_thing async_thing,
             MPIX_Async_poll_function *poll_fn, void *extra_state, MPIX_Stream stream);
typedef int (QMPI_Comm_create_keyval_t) (QMPI_Context context, int tool_id,
             MPI_Comm_copy_attr_function *comm_copy_attr_fn,
             MPI_Comm_delete_attr_function *comm_delete_attr_fn, int *comm_keyval,
             void *extra_state);
typedef int (QMPI_Keyval_create_t) (QMPI_Context context, int tool_id, MPI_Copy_function *copy_fn,
             MPI_Delete_function *delete_fn, int *keyval, void *extra_state);
typedef int (QMPI_Comm_delete_attr_t) (QMPI_Context context, int tool_id, MPI_Comm comm,
             int comm_keyval);
typedef int (QMPI_Attr_delete_t) (QMPI_Context context, int tool_id, MPI_Comm comm, int keyval);
typedef int (QMPI_Comm_free_keyval_t) (QMPI_Context context, int tool_id, int *comm_keyval);
typedef int (QMPI_Keyval_free_t) (QMPI_Context context, int tool_id, int *keyval);
typedef int (QMPI_Comm_get_attr_t) (QMPI_Context context, int tool_id, MPI_Comm comm,
             int comm_keyval, void *attribute_val, int *flag);
typedef int (QMPI_Attr_get_t) (QMPI_Context context, int tool_id, MPI_Comm comm, int keyval,
             void *attribute_val, int *flag);
typedef int (QMPI_Comm_set_attr_t) (QMPI_Context context, int tool_id, MPI_Comm comm,
             int comm_keyval, void *attribute_val);
typedef int (QMPI_Attr_put_t) (QMPI_Context context, int tool_id, MPI_Comm comm, int keyval,
             void *attribute_val);
typedef int (QMPI_Type_create_keyval_t) (QMPI_Context context, int tool_id,
             MPI_Type_copy_attr_function *type_copy_attr_fn,
             MPI_Type_delete_attr_function *type_delete_attr_fn, int *type_keyval,
             void *extra_state);
typedef int (QMPI_Type_delete_attr_t) (QMPI_Context context, int tool_id, MPI_Datatype datatype,
             int type_keyval);
typedef int (QMPI_Type_free_keyval_t) (QMPI_Context context, int tool_id, int *type_keyval);
typedef int (QMPI_Type_get_attr_t) (QMPI_Context context, int tool_id, MPI_Datatype datatype,
             int type_keyval, void *attribute_val, int *flag);
typedef int (QMPI_Type_set_attr_t) (QMPI_Context context, int tool_id, MPI_Datatype datatype,
             int type_keyval, void *attribute_val);
typedef int (QMPI_Win_create_keyval_t) (QMPI_Context context, int tool_id,
             MPI_Win_copy_attr_function *win_copy_attr_fn,
             MPI_Win_delete_attr_function *win_delete_attr_fn, int *win_keyval, void *extra_state);
typedef int (QMPI_Win_delete_attr_t) (QMPI_Context context, int tool_id, MPI_Win win,
             int win_keyval);
typedef int (QMPI_Win_free_keyval_t) (QMPI_Context context, int tool_id, int *win_keyval);
typedef int (QMPI_Win_get_attr_t) (QMPI_Context context, int tool_id, MPI_Win win, int win_keyval,
             void *attribute_val, int *flag);
typedef int (QMPI_Win_set_attr_t) (QMPI_Context context, int tool_id, MPI_Win win, int win_keyval,
             void *attribute_val);
typedef int (QMPIX_Op_create_x_t) (QMPI_Context context, int tool_id,
             MPIX_User_function_x *user_fn_x, MPIX_Destructor_function *destructor_fn, int commute,
             void *extra_state, MPI_Op *op);
typedef int (QMPIX_Comm_create_errhandler_x_t) (QMPI_Context context, int tool_id,
             MPIX_Comm_errhandler_function_x *comm_errhandler_fn_x,
             MPIX_Destructor_function *destructor_fn, void *extra_state,
             MPI_Errhandler *errhandler);
typedef int (QMPIX_Win_create_errhandler_x_t) (QMPI_Context context, int tool_id,
             MPIX_Win_errhandler_function_x *comm_errhandler_fn_x,
             MPIX_Destructor_function *destructor_fn, void *extra_state,
             MPI_Errhandler *errhandler);
typedef int (QMPIX_File_create_errhandler_x_t) (QMPI_Context context, int tool_id,
             MPIX_File_errhandler_function_x *comm_errhandler_fn_x,
             MPIX_Destructor_function *destructor_fn, void *extra_state,
             MPI_Errhandler *errhandler);
typedef int (QMPIX_Session_create_errhandler_x_t) (QMPI_Context context, int tool_id,
             MPIX_Session_errhandler_function_x *comm_errhandler_fn_x,
             MPIX_Destructor_function *destructor_fn, void *extra_state,
             MPI_Errhandler *errhandler);
typedef int (QMPI_Allgather_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount,
             MPI_Datatype recvtype, MPI_Comm comm);
typedef int (QMPI_Allgather_c_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             MPI_Count sendcount, MPI_Datatype sendtype, void *recvbuf, MPI_Count recvcount,
             MPI_Datatype recvtype, MPI_Comm comm);
typedef int (QMPI_Allgather_init_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount,
             MPI_Datatype recvtype, MPI_Comm comm, MPI_Info info, MPI_Request *request);
typedef int (QMPI_Allgather_init_c_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             MPI_Count sendcount, MPI_Datatype sendtype, void *recvbuf, MPI_Count recvcount,
             MPI_Datatype recvtype, MPI_Comm comm, MPI_Info info, MPI_Request *request);
typedef int (QMPI_Allgatherv_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             int sendcount, MPI_Datatype sendtype, void *recvbuf, const int recvcounts[],
             const int displs[], MPI_Datatype recvtype, MPI_Comm comm);
typedef int (QMPI_Allgatherv_c_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             MPI_Count sendcount, MPI_Datatype sendtype, void *recvbuf,
             const MPI_Count recvcounts[], const MPI_Aint displs[], MPI_Datatype recvtype,
             MPI_Comm comm);
typedef int (QMPI_Allgatherv_init_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             int sendcount, MPI_Datatype sendtype, void *recvbuf, const int recvcounts[],
             const int displs[], MPI_Datatype recvtype, MPI_Comm comm, MPI_Info info,
             MPI_Request *request);
typedef int (QMPI_Allgatherv_init_c_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             MPI_Count sendcount, MPI_Datatype sendtype, void *recvbuf,
             const MPI_Count recvcounts[], const MPI_Aint displs[], MPI_Datatype recvtype,
             MPI_Comm comm, MPI_Info info, MPI_Request *request);
typedef int (QMPI_Allreduce_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm);
typedef int (QMPI_Allreduce_c_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             void *recvbuf, MPI_Count count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm);
typedef int (QMPI_Allreduce_init_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm,
             MPI_Info info, MPI_Request *request);
typedef int (QMPI_Allreduce_init_c_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             void *recvbuf, MPI_Count count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm,
             MPI_Info info, MPI_Request *request);
typedef int (QMPI_Alltoall_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount,
             MPI_Datatype recvtype, MPI_Comm comm);
typedef int (QMPI_Alltoall_c_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             MPI_Count sendcount, MPI_Datatype sendtype, void *recvbuf, MPI_Count recvcount,
             MPI_Datatype recvtype, MPI_Comm comm);
typedef int (QMPI_Alltoall_init_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount,
             MPI_Datatype recvtype, MPI_Comm comm, MPI_Info info, MPI_Request *request);
typedef int (QMPI_Alltoall_init_c_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             MPI_Count sendcount, MPI_Datatype sendtype, void *recvbuf, MPI_Count recvcount,
             MPI_Datatype recvtype, MPI_Comm comm, MPI_Info info, MPI_Request *request);
typedef int (QMPI_Alltoallv_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             const int sendcounts[], const int sdispls[], MPI_Datatype sendtype, void *recvbuf,
             const int recvcounts[], const int rdispls[], MPI_Datatype recvtype, MPI_Comm comm);
typedef int (QMPI_Alltoallv_c_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             const MPI_Count sendcounts[], const MPI_Aint sdispls[], MPI_Datatype sendtype,
             void *recvbuf, const MPI_Count recvcounts[], const MPI_Aint rdispls[],
             MPI_Datatype recvtype, MPI_Comm comm);
typedef int (QMPI_Alltoallv_init_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             const int sendcounts[], const int sdispls[], MPI_Datatype sendtype, void *recvbuf,
             const int recvcounts[], const int rdispls[], MPI_Datatype recvtype, MPI_Comm comm,
             MPI_Info info, MPI_Request *request);
typedef int (QMPI_Alltoallv_init_c_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             const MPI_Count sendcounts[], const MPI_Aint sdispls[], MPI_Datatype sendtype,
             void *recvbuf, const MPI_Count recvcounts[], const MPI_Aint rdispls[],
             MPI_Datatype recvtype, MPI_Comm comm, MPI_Info info, MPI_Request *request);
typedef int (QMPI_Alltoallw_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             const int sendcounts[], const int sdispls[], const MPI_Datatype sendtypes[],
             void *recvbuf, const int recvcounts[], const int rdispls[],
             const MPI_Datatype recvtypes[], MPI_Comm comm);
typedef int (QMPI_Alltoallw_c_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             const MPI_Count sendcounts[], const MPI_Aint sdispls[], const MPI_Datatype sendtypes[],
             void *recvbuf, const MPI_Count recvcounts[], const MPI_Aint rdispls[],
             const MPI_Datatype recvtypes[], MPI_Comm comm);
typedef int (QMPI_Alltoallw_init_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             const int sendcounts[], const int sdispls[], const MPI_Datatype sendtypes[],
             void *recvbuf, const int recvcounts[], const int rdispls[],
             const MPI_Datatype recvtypes[], MPI_Comm comm, MPI_Info info, MPI_Request *request);
typedef int (QMPI_Alltoallw_init_c_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             const MPI_Count sendcounts[], const MPI_Aint sdispls[], const MPI_Datatype sendtypes[],
             void *recvbuf, const MPI_Count recvcounts[], const MPI_Aint rdispls[],
             const MPI_Datatype recvtypes[], MPI_Comm comm, MPI_Info info, MPI_Request *request);
typedef int (QMPI_Barrier_t) (QMPI_Context context, int tool_id, MPI_Comm comm);
typedef int (QMPI_Barrier_init_t) (QMPI_Context context, int tool_id, MPI_Comm comm, MPI_Info info,
             MPI_Request *request);
typedef int (QMPI_Bcast_t) (QMPI_Context context, int tool_id, void *buffer, int count,
             MPI_Datatype datatype, int root, MPI_Comm comm);
typedef int (QMPI_Bcast_c_t) (QMPI_Context context, int tool_id, void *buffer, MPI_Count count,
             MPI_Datatype datatype, int root, MPI_Comm comm);
typedef int (QMPI_Bcast_init_t) (QMPI_Context context, int tool_id, void *buffer, int count,
             MPI_Datatype datatype, int root, MPI_Comm comm, MPI_Info info, MPI_Request *request);
typedef int (QMPI_Bcast_init_c_t) (QMPI_Context context, int tool_id, void *buffer, MPI_Count count,
             MPI_Datatype datatype, int root, MPI_Comm comm, MPI_Info info, MPI_Request *request);
typedef int (QMPI_Exscan_t) (QMPI_Context context, int tool_id, const void *sendbuf, void *recvbuf,
             int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm);
typedef int (QMPI_Exscan_c_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             void *recvbuf, MPI_Count count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm);
typedef int (QMPI_Exscan_init_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm,
             MPI_Info info, MPI_Request *request);
typedef int (QMPI_Exscan_init_c_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             void *recvbuf, MPI_Count count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm,
             MPI_Info info, MPI_Request *request);
typedef int (QMPI_Gather_t) (QMPI_Context context, int tool_id, const void *sendbuf, int sendcount,
             MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype, int root,
             MPI_Comm comm);
typedef int (QMPI_Gather_c_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             MPI_Count sendcount, MPI_Datatype sendtype, void *recvbuf, MPI_Count recvcount,
             MPI_Datatype recvtype, int root, MPI_Comm comm);
typedef int (QMPI_Gather_init_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount,
             MPI_Datatype recvtype, int root, MPI_Comm comm, MPI_Info info, MPI_Request *request);
typedef int (QMPI_Gather_init_c_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             MPI_Count sendcount, MPI_Datatype sendtype, void *recvbuf, MPI_Count recvcount,
             MPI_Datatype recvtype, int root, MPI_Comm comm, MPI_Info info, MPI_Request *request);
typedef int (QMPI_Gatherv_t) (QMPI_Context context, int tool_id, const void *sendbuf, int sendcount,
             MPI_Datatype sendtype, void *recvbuf, const int recvcounts[], const int displs[],
             MPI_Datatype recvtype, int root, MPI_Comm comm);
typedef int (QMPI_Gatherv_c_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             MPI_Count sendcount, MPI_Datatype sendtype, void *recvbuf,
             const MPI_Count recvcounts[], const MPI_Aint displs[], MPI_Datatype recvtype, int root,
             MPI_Comm comm);
typedef int (QMPI_Gatherv_init_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             int sendcount, MPI_Datatype sendtype, void *recvbuf, const int recvcounts[],
             const int displs[], MPI_Datatype recvtype, int root, MPI_Comm comm, MPI_Info info,
             MPI_Request *request);
typedef int (QMPI_Gatherv_init_c_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             MPI_Count sendcount, MPI_Datatype sendtype, void *recvbuf,
             const MPI_Count recvcounts[], const MPI_Aint displs[], MPI_Datatype recvtype, int root,
             MPI_Comm comm, MPI_Info info, MPI_Request *request);
typedef int (QMPI_Iallgather_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount,
             MPI_Datatype recvtype, MPI_Comm comm, MPI_Request *request);
typedef int (QMPI_Iallgather_c_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             MPI_Count sendcount, MPI_Datatype sendtype, void *recvbuf, MPI_Count recvcount,
             MPI_Datatype recvtype, MPI_Comm comm, MPI_Request *request);
typedef int (QMPI_Iallgatherv_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             int sendcount, MPI_Datatype sendtype, void *recvbuf, const int recvcounts[],
             const int displs[], MPI_Datatype recvtype, MPI_Comm comm, MPI_Request *request);
typedef int (QMPI_Iallgatherv_c_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             MPI_Count sendcount, MPI_Datatype sendtype, void *recvbuf,
             const MPI_Count recvcounts[], const MPI_Aint displs[], MPI_Datatype recvtype,
             MPI_Comm comm, MPI_Request *request);
typedef int (QMPI_Iallreduce_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm,
             MPI_Request *request);
typedef int (QMPI_Iallreduce_c_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             void *recvbuf, MPI_Count count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm,
             MPI_Request *request);
typedef int (QMPI_Ialltoall_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount,
             MPI_Datatype recvtype, MPI_Comm comm, MPI_Request *request);
typedef int (QMPI_Ialltoall_c_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             MPI_Count sendcount, MPI_Datatype sendtype, void *recvbuf, MPI_Count recvcount,
             MPI_Datatype recvtype, MPI_Comm comm, MPI_Request *request);
typedef int (QMPI_Ialltoallv_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             const int sendcounts[], const int sdispls[], MPI_Datatype sendtype, void *recvbuf,
             const int recvcounts[], const int rdispls[], MPI_Datatype recvtype, MPI_Comm comm,
             MPI_Request *request);
typedef int (QMPI_Ialltoallv_c_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             const MPI_Count sendcounts[], const MPI_Aint sdispls[], MPI_Datatype sendtype,
             void *recvbuf, const MPI_Count recvcounts[], const MPI_Aint rdispls[],
             MPI_Datatype recvtype, MPI_Comm comm, MPI_Request *request);
typedef int (QMPI_Ialltoallw_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             const int sendcounts[], const int sdispls[], const MPI_Datatype sendtypes[],
             void *recvbuf, const int recvcounts[], const int rdispls[],
             const MPI_Datatype recvtypes[], MPI_Comm comm, MPI_Request *request);
typedef int (QMPI_Ialltoallw_c_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             const MPI_Count sendcounts[], const MPI_Aint sdispls[], const MPI_Datatype sendtypes[],
             void *recvbuf, const MPI_Count recvcounts[], const MPI_Aint rdispls[],
             const MPI_Datatype recvtypes[], MPI_Comm comm, MPI_Request *request);
typedef int (QMPI_Ibarrier_t) (QMPI_Context context, int tool_id, MPI_Comm comm,
             MPI_Request *request);
typedef int (QMPI_Ibcast_t) (QMPI_Context context, int tool_id, void *buffer, int count,
             MPI_Datatype datatype, int root, MPI_Comm comm, MPI_Request *request);
typedef int (QMPI_Ibcast_c_t) (QMPI_Context context, int tool_id, void *buffer, MPI_Count count,
             MPI_Datatype datatype, int root, MPI_Comm comm, MPI_Request *request);
typedef int (QMPI_Iexscan_t) (QMPI_Context context, int tool_id, const void *sendbuf, void *recvbuf,
             int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm, MPI_Request *request);
typedef int (QMPI_Iexscan_c_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             void *recvbuf, MPI_Count count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm,
             MPI_Request *request);
typedef int (QMPI_Igather_t) (QMPI_Context context, int tool_id, const void *sendbuf, int sendcount,
             MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype, int root,
             MPI_Comm comm, MPI_Request *request);
typedef int (QMPI_Igather_c_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             MPI_Count sendcount, MPI_Datatype sendtype, void *recvbuf, MPI_Count recvcount,
             MPI_Datatype recvtype, int root, MPI_Comm comm, MPI_Request *request);
typedef int (QMPI_Igatherv_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             int sendcount, MPI_Datatype sendtype, void *recvbuf, const int recvcounts[],
             const int displs[], MPI_Datatype recvtype, int root, MPI_Comm comm,
             MPI_Request *request);
typedef int (QMPI_Igatherv_c_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             MPI_Count sendcount, MPI_Datatype sendtype, void *recvbuf,
             const MPI_Count recvcounts[], const MPI_Aint displs[], MPI_Datatype recvtype, int root,
             MPI_Comm comm, MPI_Request *request);
typedef int (QMPI_Ineighbor_allgather_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount,
             MPI_Datatype recvtype, MPI_Comm comm, MPI_Request *request);
typedef int (QMPI_Ineighbor_allgather_c_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             MPI_Count sendcount, MPI_Datatype sendtype, void *recvbuf, MPI_Count recvcount,
             MPI_Datatype recvtype, MPI_Comm comm, MPI_Request *request);
typedef int (QMPI_Ineighbor_allgatherv_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             int sendcount, MPI_Datatype sendtype, void *recvbuf, const int recvcounts[],
             const int displs[], MPI_Datatype recvtype, MPI_Comm comm, MPI_Request *request);
typedef int (QMPI_Ineighbor_allgatherv_c_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             MPI_Count sendcount, MPI_Datatype sendtype, void *recvbuf,
             const MPI_Count recvcounts[], const MPI_Aint displs[], MPI_Datatype recvtype,
             MPI_Comm comm, MPI_Request *request);
typedef int (QMPI_Ineighbor_alltoall_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount,
             MPI_Datatype recvtype, MPI_Comm comm, MPI_Request *request);
typedef int (QMPI_Ineighbor_alltoall_c_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             MPI_Count sendcount, MPI_Datatype sendtype, void *recvbuf, MPI_Count recvcount,
             MPI_Datatype recvtype, MPI_Comm comm, MPI_Request *request);
typedef int (QMPI_Ineighbor_alltoallv_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             const int sendcounts[], const int sdispls[], MPI_Datatype sendtype, void *recvbuf,
             const int recvcounts[], const int rdispls[], MPI_Datatype recvtype, MPI_Comm comm,
             MPI_Request *request);
typedef int (QMPI_Ineighbor_alltoallv_c_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             const MPI_Count sendcounts[], const MPI_Aint sdispls[], MPI_Datatype sendtype,
             void *recvbuf, const MPI_Count recvcounts[], const MPI_Aint rdispls[],
             MPI_Datatype recvtype, MPI_Comm comm, MPI_Request *request);
typedef int (QMPI_Ineighbor_alltoallw_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             const int sendcounts[], const MPI_Aint sdispls[], const MPI_Datatype sendtypes[],
             void *recvbuf, const int recvcounts[], const MPI_Aint rdispls[],
             const MPI_Datatype recvtypes[], MPI_Comm comm, MPI_Request *request);
typedef int (QMPI_Ineighbor_alltoallw_c_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             const MPI_Count sendcounts[], const MPI_Aint sdispls[], const MPI_Datatype sendtypes[],
             void *recvbuf, const MPI_Count recvcounts[], const MPI_Aint rdispls[],
             const MPI_Datatype recvtypes[], MPI_Comm comm, MPI_Request *request);
typedef int (QMPI_Ireduce_t) (QMPI_Context context, int tool_id, const void *sendbuf, void *recvbuf,
             int count, MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm,
             MPI_Request *request);
typedef int (QMPI_Ireduce_c_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             void *recvbuf, MPI_Count count, MPI_Datatype datatype, MPI_Op op, int root,
             MPI_Comm comm, MPI_Request *request);
typedef int (QMPI_Ireduce_scatter_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             void *recvbuf, const int recvcounts[], MPI_Datatype datatype, MPI_Op op, MPI_Comm comm,
             MPI_Request *request);
typedef int (QMPI_Ireduce_scatter_c_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             void *recvbuf, const MPI_Count recvcounts[], MPI_Datatype datatype, MPI_Op op,
             MPI_Comm comm, MPI_Request *request);
typedef int (QMPI_Ireduce_scatter_block_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             void *recvbuf, int recvcount, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm,
             MPI_Request *request);
typedef int (QMPI_Ireduce_scatter_block_c_t) (QMPI_Context context, int tool_id,
             const void *sendbuf, void *recvbuf, MPI_Count recvcount, MPI_Datatype datatype,
             MPI_Op op, MPI_Comm comm, MPI_Request *request);
typedef int (QMPI_Iscan_t) (QMPI_Context context, int tool_id, const void *sendbuf, void *recvbuf,
             int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm, MPI_Request *request);
typedef int (QMPI_Iscan_c_t) (QMPI_Context context, int tool_id, const void *sendbuf, void *recvbuf,
             MPI_Count count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm,
             MPI_Request *request);
typedef int (QMPI_Iscatter_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount,
             MPI_Datatype recvtype, int root, MPI_Comm comm, MPI_Request *request);
typedef int (QMPI_Iscatter_c_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             MPI_Count sendcount, MPI_Datatype sendtype, void *recvbuf, MPI_Count recvcount,
             MPI_Datatype recvtype, int root, MPI_Comm comm, MPI_Request *request);
typedef int (QMPI_Iscatterv_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             const int sendcounts[], const int displs[], MPI_Datatype sendtype, void *recvbuf,
             int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm, MPI_Request *request);
typedef int (QMPI_Iscatterv_c_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             const MPI_Count sendcounts[], const MPI_Aint displs[], MPI_Datatype sendtype,
             void *recvbuf, MPI_Count recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm,
             MPI_Request *request);
typedef int (QMPI_Neighbor_allgather_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount,
             MPI_Datatype recvtype, MPI_Comm comm);
typedef int (QMPI_Neighbor_allgather_c_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             MPI_Count sendcount, MPI_Datatype sendtype, void *recvbuf, MPI_Count recvcount,
             MPI_Datatype recvtype, MPI_Comm comm);
typedef int (QMPI_Neighbor_allgather_init_t) (QMPI_Context context, int tool_id,
             const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
             int recvcount, MPI_Datatype recvtype, MPI_Comm comm, MPI_Info info,
             MPI_Request *request);
typedef int (QMPI_Neighbor_allgather_init_c_t) (QMPI_Context context, int tool_id,
             const void *sendbuf, MPI_Count sendcount, MPI_Datatype sendtype, void *recvbuf,
             MPI_Count recvcount, MPI_Datatype recvtype, MPI_Comm comm, MPI_Info info,
             MPI_Request *request);
typedef int (QMPI_Neighbor_allgatherv_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             int sendcount, MPI_Datatype sendtype, void *recvbuf, const int recvcounts[],
             const int displs[], MPI_Datatype recvtype, MPI_Comm comm);
typedef int (QMPI_Neighbor_allgatherv_c_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             MPI_Count sendcount, MPI_Datatype sendtype, void *recvbuf,
             const MPI_Count recvcounts[], const MPI_Aint displs[], MPI_Datatype recvtype,
             MPI_Comm comm);
typedef int (QMPI_Neighbor_allgatherv_init_t) (QMPI_Context context, int tool_id,
             const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
             const int recvcounts[], const int displs[], MPI_Datatype recvtype, MPI_Comm comm,
             MPI_Info info, MPI_Request *request);
typedef int (QMPI_Neighbor_allgatherv_init_c_t) (QMPI_Context context, int tool_id,
             const void *sendbuf, MPI_Count sendcount, MPI_Datatype sendtype, void *recvbuf,
             const MPI_Count recvcounts[], const MPI_Aint displs[], MPI_Datatype recvtype,
             MPI_Comm comm, MPI_Info info, MPI_Request *request);
typedef int (QMPI_Neighbor_alltoall_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount,
             MPI_Datatype recvtype, MPI_Comm comm);
typedef int (QMPI_Neighbor_alltoall_c_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             MPI_Count sendcount, MPI_Datatype sendtype, void *recvbuf, MPI_Count recvcount,
             MPI_Datatype recvtype, MPI_Comm comm);
typedef int (QMPI_Neighbor_alltoall_init_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount,
             MPI_Datatype recvtype, MPI_Comm comm, MPI_Info info, MPI_Request *request);
typedef int (QMPI_Neighbor_alltoall_init_c_t) (QMPI_Context context, int tool_id,
             const void *sendbuf, MPI_Count sendcount, MPI_Datatype sendtype, void *recvbuf,
             MPI_Count recvcount, MPI_Datatype recvtype, MPI_Comm comm, MPI_Info info,
             MPI_Request *request);
typedef int (QMPI_Neighbor_alltoallv_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             const int sendcounts[], const int sdispls[], MPI_Datatype sendtype, void *recvbuf,
             const int recvcounts[], const int rdispls[], MPI_Datatype recvtype, MPI_Comm comm);
typedef int (QMPI_Neighbor_alltoallv_c_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             const MPI_Count sendcounts[], const MPI_Aint sdispls[], MPI_Datatype sendtype,
             void *recvbuf, const MPI_Count recvcounts[], const MPI_Aint rdispls[],
             MPI_Datatype recvtype, MPI_Comm comm);
typedef int (QMPI_Neighbor_alltoallv_init_t) (QMPI_Context context, int tool_id,
             const void *sendbuf, const int sendcounts[], const int sdispls[],
             MPI_Datatype sendtype, void *recvbuf, const int recvcounts[], const int rdispls[],
             MPI_Datatype recvtype, MPI_Comm comm, MPI_Info info, MPI_Request *request);
typedef int (QMPI_Neighbor_alltoallv_init_c_t) (QMPI_Context context, int tool_id,
             const void *sendbuf, const MPI_Count sendcounts[], const MPI_Aint sdispls[],
             MPI_Datatype sendtype, void *recvbuf, const MPI_Count recvcounts[],
             const MPI_Aint rdispls[], MPI_Datatype recvtype, MPI_Comm comm, MPI_Info info,
             MPI_Request *request);
typedef int (QMPI_Neighbor_alltoallw_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             const int sendcounts[], const MPI_Aint sdispls[], const MPI_Datatype sendtypes[],
             void *recvbuf, const int recvcounts[], const MPI_Aint rdispls[],
             const MPI_Datatype recvtypes[], MPI_Comm comm);
typedef int (QMPI_Neighbor_alltoallw_c_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             const MPI_Count sendcounts[], const MPI_Aint sdispls[], const MPI_Datatype sendtypes[],
             void *recvbuf, const MPI_Count recvcounts[], const MPI_Aint rdispls[],
             const MPI_Datatype recvtypes[], MPI_Comm comm);
typedef int (QMPI_Neighbor_alltoallw_init_t) (QMPI_Context context, int tool_id,
             const void *sendbuf, const int sendcounts[], const MPI_Aint sdispls[],
             const MPI_Datatype sendtypes[], void *recvbuf, const int recvcounts[],
             const MPI_Aint rdispls[], const MPI_Datatype recvtypes[], MPI_Comm comm, MPI_Info info,
             MPI_Request *request);
typedef int (QMPI_Neighbor_alltoallw_init_c_t) (QMPI_Context context, int tool_id,
             const void *sendbuf, const MPI_Count sendcounts[], const MPI_Aint sdispls[],
             const MPI_Datatype sendtypes[], void *recvbuf, const MPI_Count recvcounts[],
             const MPI_Aint rdispls[], const MPI_Datatype recvtypes[], MPI_Comm comm, MPI_Info info,
             MPI_Request *request);
typedef int (QMPI_Reduce_t) (QMPI_Context context, int tool_id, const void *sendbuf, void *recvbuf,
             int count, MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm);
typedef int (QMPI_Reduce_c_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             void *recvbuf, MPI_Count count, MPI_Datatype datatype, MPI_Op op, int root,
             MPI_Comm comm);
typedef int (QMPI_Reduce_init_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm,
             MPI_Info info, MPI_Request *request);
typedef int (QMPI_Reduce_init_c_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             void *recvbuf, MPI_Count count, MPI_Datatype datatype, MPI_Op op, int root,
             MPI_Comm comm, MPI_Info info, MPI_Request *request);
typedef int (QMPI_Reduce_local_t) (QMPI_Context context, int tool_id, const void *inbuf,
             void *inoutbuf, int count, MPI_Datatype datatype, MPI_Op op);
typedef int (QMPI_Reduce_local_c_t) (QMPI_Context context, int tool_id, const void *inbuf,
             void *inoutbuf, MPI_Count count, MPI_Datatype datatype, MPI_Op op);
typedef int (QMPI_Reduce_scatter_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             void *recvbuf, const int recvcounts[], MPI_Datatype datatype, MPI_Op op,
             MPI_Comm comm);
typedef int (QMPI_Reduce_scatter_c_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             void *recvbuf, const MPI_Count recvcounts[], MPI_Datatype datatype, MPI_Op op,
             MPI_Comm comm);
typedef int (QMPI_Reduce_scatter_block_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             void *recvbuf, int recvcount, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm);
typedef int (QMPI_Reduce_scatter_block_c_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             void *recvbuf, MPI_Count recvcount, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm);
typedef int (QMPI_Reduce_scatter_block_init_t) (QMPI_Context context, int tool_id,
             const void *sendbuf, void *recvbuf, int recvcount, MPI_Datatype datatype, MPI_Op op,
             MPI_Comm comm, MPI_Info info, MPI_Request *request);
typedef int (QMPI_Reduce_scatter_block_init_c_t) (QMPI_Context context, int tool_id,
             const void *sendbuf, void *recvbuf, MPI_Count recvcount, MPI_Datatype datatype,
             MPI_Op op, MPI_Comm comm, MPI_Info info, MPI_Request *request);
typedef int (QMPI_Reduce_scatter_init_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             void *recvbuf, const int recvcounts[], MPI_Datatype datatype, MPI_Op op, MPI_Comm comm,
             MPI_Info info, MPI_Request *request);
typedef int (QMPI_Reduce_scatter_init_c_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             void *recvbuf, const MPI_Count recvcounts[], MPI_Datatype datatype, MPI_Op op,
             MPI_Comm comm, MPI_Info info, MPI_Request *request);
typedef int (QMPI_Scan_t) (QMPI_Context context, int tool_id, const void *sendbuf, void *recvbuf,
             int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm);
typedef int (QMPI_Scan_c_t) (QMPI_Context context, int tool_id, const void *sendbuf, void *recvbuf,
             MPI_Count count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm);
typedef int (QMPI_Scan_init_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm,
             MPI_Info info, MPI_Request *request);
typedef int (QMPI_Scan_init_c_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             void *recvbuf, MPI_Count count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm,
             MPI_Info info, MPI_Request *request);
typedef int (QMPI_Scatter_t) (QMPI_Context context, int tool_id, const void *sendbuf, int sendcount,
             MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype, int root,
             MPI_Comm comm);
typedef int (QMPI_Scatter_c_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             MPI_Count sendcount, MPI_Datatype sendtype, void *recvbuf, MPI_Count recvcount,
             MPI_Datatype recvtype, int root, MPI_Comm comm);
typedef int (QMPI_Scatter_init_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount,
             MPI_Datatype recvtype, int root, MPI_Comm comm, MPI_Info info, MPI_Request *request);
typedef int (QMPI_Scatter_init_c_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             MPI_Count sendcount, MPI_Datatype sendtype, void *recvbuf, MPI_Count recvcount,
             MPI_Datatype recvtype, int root, MPI_Comm comm, MPI_Info info, MPI_Request *request);
typedef int (QMPI_Scatterv_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             const int sendcounts[], const int displs[], MPI_Datatype sendtype, void *recvbuf,
             int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm);
typedef int (QMPI_Scatterv_c_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             const MPI_Count sendcounts[], const MPI_Aint displs[], MPI_Datatype sendtype,
             void *recvbuf, MPI_Count recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm);
typedef int (QMPI_Scatterv_init_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             const int sendcounts[], const int displs[], MPI_Datatype sendtype, void *recvbuf,
             int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm, MPI_Info info,
             MPI_Request *request);
typedef int (QMPI_Scatterv_init_c_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             const MPI_Count sendcounts[], const MPI_Aint displs[], MPI_Datatype sendtype,
             void *recvbuf, MPI_Count recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm,
             MPI_Info info, MPI_Request *request);
typedef int (QMPI_Comm_compare_t) (QMPI_Context context, int tool_id, MPI_Comm comm1,
             MPI_Comm comm2, int *result);
typedef int (QMPI_Comm_create_t) (QMPI_Context context, int tool_id, MPI_Comm comm, MPI_Group group,
             MPI_Comm *newcomm);
typedef int (QMPI_Comm_create_group_t) (QMPI_Context context, int tool_id, MPI_Comm comm,
             MPI_Group group, int tag, MPI_Comm *newcomm);
typedef int (QMPI_Comm_dup_t) (QMPI_Context context, int tool_id, MPI_Comm comm,
             MPI_Comm *newcomm);
typedef int (QMPI_Comm_dup_with_info_t) (QMPI_Context context, int tool_id, MPI_Comm comm,
             MPI_Info info, MPI_Comm *newcomm);
typedef int (QMPI_Comm_free_t) (QMPI_Context context, int tool_id, MPI_Comm *comm);
typedef int (QMPI_Comm_get_info_t) (QMPI_Context context, int tool_id, MPI_Comm comm,
             MPI_Info *info_used);
typedef int (QMPI_Comm_get_name_t) (QMPI_Context context, int tool_id, MPI_Comm comm,
             char *comm_name, int *resultlen);
typedef int (QMPI_Comm_group_t) (QMPI_Context context, int tool_id, MPI_Comm comm,
             MPI_Group *group);
typedef int (QMPI_Comm_idup_t) (QMPI_Context context, int tool_id, MPI_Comm comm, MPI_Comm *newcomm,
             MPI_Request *request);
typedef int (QMPI_Comm_idup_with_info_t) (QMPI_Context context, int tool_id, MPI_Comm comm,
             MPI_Info info, MPI_Comm *newcomm, MPI_Request *request);
typedef int (QMPI_Comm_rank_t) (QMPI_Context context, int tool_id, MPI_Comm comm, int *rank);
typedef int (QMPI_Comm_remote_group_t) (QMPI_Context context, int tool_id, MPI_Comm comm,
             MPI_Group *group);
typedef int (QMPI_Comm_remote_size_t) (QMPI_Context context, int tool_id, MPI_Comm comm,
             int *size);
typedef int (QMPI_Comm_set_info_t) (QMPI_Context context, int tool_id, MPI_Comm comm,
             MPI_Info info);
typedef int (QMPI_Comm_set_name_t) (QMPI_Context context, int tool_id, MPI_Comm comm,
             const char *comm_name);
typedef int (QMPI_Comm_size_t) (QMPI_Context context, int tool_id, MPI_Comm comm, int *size);
typedef int (QMPI_Comm_split_t) (QMPI_Context context, int tool_id, MPI_Comm comm, int color,
             int key, MPI_Comm *newcomm);
typedef int (QMPI_Comm_split_type_t) (QMPI_Context context, int tool_id, MPI_Comm comm,
             int split_type, int key, MPI_Info info, MPI_Comm *newcomm);
typedef int (QMPI_Comm_test_inter_t) (QMPI_Context context, int tool_id, MPI_Comm comm, int *flag);
typedef int (QMPI_Intercomm_create_t) (QMPI_Context context, int tool_id, MPI_Comm local_comm,
             int local_leader, MPI_Comm peer_comm, int remote_leader, int tag,
             MPI_Comm *newintercomm);
typedef int (QMPI_Intercomm_create_from_groups_t) (QMPI_Context context, int tool_id,
             MPI_Group local_group, int local_leader, MPI_Group remote_group, int remote_leader,
             const char *stringtag, MPI_Info info, MPI_Errhandler errhandler,
             MPI_Comm *newintercomm);
typedef int (QMPI_Intercomm_merge_t) (QMPI_Context context, int tool_id, MPI_Comm intercomm,
             int high, MPI_Comm *newintracomm);
typedef int (QMPIX_Comm_test_threadcomm_t) (QMPI_Context context, int tool_id, MPI_Comm comm,
             int *flag);
typedef int (QMPIX_Comm_revoke_t) (QMPI_Context context, int tool_id, MPI_Comm comm);
typedef int (QMPIX_Comm_shrink_t) (QMPI_Context context, int tool_id, MPI_Comm comm,
             MPI_Comm *newcomm);
typedef int (QMPIX_Comm_failure_ack_t) (QMPI_Context context, int tool_id, MPI_Comm comm);
typedef int (QMPIX_Comm_failure_get_acked_t) (QMPI_Context context, int tool_id, MPI_Comm comm,
             MPI_Group *failedgrp);
typedef int (QMPIX_Comm_agree_t) (QMPI_Context context, int tool_id, MPI_Comm comm, int *flag);
typedef int (QMPIX_Comm_get_failed_t) (QMPI_Context context, int tool_id, MPI_Comm comm,
             MPI_Group *failedgrp);
typedef int (QMPI_Get_address_t) (QMPI_Context context, int tool_id, const void *location,
             MPI_Aint *address);
typedef int (QMPI_Address_t) (QMPI_Context context, int tool_id, void *location,
             MPI_Aint *address);
typedef int (QMPI_Get_count_t) (QMPI_Context context, int tool_id, const MPI_Status *status,
             MPI_Datatype datatype, int *count);
typedef int (QMPI_Get_count_c_t) (QMPI_Context context, int tool_id, const MPI_Status *status,
             MPI_Datatype datatype, MPI_Count *count);
typedef int (QMPI_Get_elements_t) (QMPI_Context context, int tool_id, const MPI_Status *status,
             MPI_Datatype datatype, int *count);
typedef int (QMPI_Get_elements_c_t) (QMPI_Context context, int tool_id, const MPI_Status *status,
             MPI_Datatype datatype, MPI_Count *count);
typedef int (QMPI_Get_elements_x_t) (QMPI_Context context, int tool_id, const MPI_Status *status,
             MPI_Datatype datatype, MPI_Count *count);
typedef int (QMPI_Pack_t) (QMPI_Context context, int tool_id, const void *inbuf, int incount,
             MPI_Datatype datatype, void *outbuf, int outsize, int *position, MPI_Comm comm);
typedef int (QMPI_Pack_c_t) (QMPI_Context context, int tool_id, const void *inbuf,
             MPI_Count incount, MPI_Datatype datatype, void *outbuf, MPI_Count outsize,
             MPI_Count *position, MPI_Comm comm);
typedef int (QMPI_Pack_external_t) (QMPI_Context context, int tool_id, const char *datarep,
             const void *inbuf, int incount, MPI_Datatype datatype, void *outbuf, MPI_Aint outsize,
             MPI_Aint *position);
typedef int (QMPI_Pack_external_c_t) (QMPI_Context context, int tool_id, const char *datarep,
             const void *inbuf, MPI_Count incount, MPI_Datatype datatype, void *outbuf,
             MPI_Count outsize, MPI_Count *position);
typedef int (QMPI_Pack_external_size_t) (QMPI_Context context, int tool_id, const char *datarep,
             int incount, MPI_Datatype datatype, MPI_Aint *size);
typedef int (QMPI_Pack_external_size_c_t) (QMPI_Context context, int tool_id, const char *datarep,
             MPI_Count incount, MPI_Datatype datatype, MPI_Count *size);
typedef int (QMPI_Pack_size_t) (QMPI_Context context, int tool_id, int incount,
             MPI_Datatype datatype, MPI_Comm comm, int *size);
typedef int (QMPI_Pack_size_c_t) (QMPI_Context context, int tool_id, MPI_Count incount,
             MPI_Datatype datatype, MPI_Comm comm, MPI_Count *size);
typedef int (QMPI_Status_set_elements_t) (QMPI_Context context, int tool_id, MPI_Status *status,
             MPI_Datatype datatype, int count);
typedef int (QMPI_Status_set_elements_c_t) (QMPI_Context context, int tool_id, MPI_Status *status,
             MPI_Datatype datatype, MPI_Count count);
typedef int (QMPI_Status_set_elements_x_t) (QMPI_Context context, int tool_id, MPI_Status *status,
             MPI_Datatype datatype, MPI_Count count);
typedef int (QMPI_Type_commit_t) (QMPI_Context context, int tool_id, MPI_Datatype *datatype);
typedef int (QMPI_Type_contiguous_t) (QMPI_Context context, int tool_id, int count,
             MPI_Datatype oldtype, MPI_Datatype *newtype);
typedef int (QMPI_Type_contiguous_c_t) (QMPI_Context context, int tool_id, MPI_Count count,
             MPI_Datatype oldtype, MPI_Datatype *newtype);
typedef int (QMPI_Type_create_darray_t) (QMPI_Context context, int tool_id, int size, int rank,
             int ndims, const int array_of_gsizes[], const int array_of_distribs[],
             const int array_of_dargs[], const int array_of_psizes[], int order,
             MPI_Datatype oldtype, MPI_Datatype *newtype);
typedef int (QMPI_Type_create_darray_c_t) (QMPI_Context context, int tool_id, int size, int rank,
             int ndims, const MPI_Count array_of_gsizes[], const int array_of_distribs[],
             const int array_of_dargs[], const int array_of_psizes[], int order,
             MPI_Datatype oldtype, MPI_Datatype *newtype);
typedef int (QMPI_Type_create_f90_complex_t) (QMPI_Context context, int tool_id, int p, int r,
             MPI_Datatype *newtype);
typedef int (QMPI_Type_create_f90_integer_t) (QMPI_Context context, int tool_id, int r,
             MPI_Datatype *newtype);
typedef int (QMPI_Type_create_f90_real_t) (QMPI_Context context, int tool_id, int p, int r,
             MPI_Datatype *newtype);
typedef int (QMPI_Type_create_hindexed_t) (QMPI_Context context, int tool_id, int count,
             const int array_of_blocklengths[], const MPI_Aint array_of_displacements[],
             MPI_Datatype oldtype, MPI_Datatype *newtype);
typedef int (QMPI_Type_create_hindexed_c_t) (QMPI_Context context, int tool_id, MPI_Count count,
             const MPI_Count array_of_blocklengths[], const MPI_Count array_of_displacements[],
             MPI_Datatype oldtype, MPI_Datatype *newtype);
typedef int (QMPI_Type_hindexed_t) (QMPI_Context context, int tool_id, int count,
             int array_of_blocklengths[], MPI_Aint array_of_displacements[], MPI_Datatype oldtype,
             MPI_Datatype *newtype);
typedef int (QMPI_Type_create_hindexed_block_t) (QMPI_Context context, int tool_id, int count,
             int blocklength, const MPI_Aint array_of_displacements[], MPI_Datatype oldtype,
             MPI_Datatype *newtype);
typedef int (QMPI_Type_create_hindexed_block_c_t) (QMPI_Context context, int tool_id,
             MPI_Count count, MPI_Count blocklength, const MPI_Count array_of_displacements[],
             MPI_Datatype oldtype, MPI_Datatype *newtype);
typedef int (QMPI_Type_create_hvector_t) (QMPI_Context context, int tool_id, int count,
             int blocklength, MPI_Aint stride, MPI_Datatype oldtype, MPI_Datatype *newtype);
typedef int (QMPI_Type_create_hvector_c_t) (QMPI_Context context, int tool_id, MPI_Count count,
             MPI_Count blocklength, MPI_Count stride, MPI_Datatype oldtype, MPI_Datatype *newtype);
typedef int (QMPI_Type_hvector_t) (QMPI_Context context, int tool_id, int count, int blocklength,
             MPI_Aint stride, MPI_Datatype oldtype, MPI_Datatype *newtype);
typedef int (QMPI_Type_create_indexed_block_t) (QMPI_Context context, int tool_id, int count,
             int blocklength, const int array_of_displacements[], MPI_Datatype oldtype,
             MPI_Datatype *newtype);
typedef int (QMPI_Type_create_indexed_block_c_t) (QMPI_Context context, int tool_id,
             MPI_Count count, MPI_Count blocklength, const MPI_Count array_of_displacements[],
             MPI_Datatype oldtype, MPI_Datatype *newtype);
typedef int (QMPI_Type_create_resized_t) (QMPI_Context context, int tool_id, MPI_Datatype oldtype,
             MPI_Aint lb, MPI_Aint extent, MPI_Datatype *newtype);
typedef int (QMPI_Type_create_resized_c_t) (QMPI_Context context, int tool_id, MPI_Datatype oldtype,
             MPI_Count lb, MPI_Count extent, MPI_Datatype *newtype);
typedef int (QMPI_Type_create_struct_t) (QMPI_Context context, int tool_id, int count,
             const int array_of_blocklengths[], const MPI_Aint array_of_displacements[],
             const MPI_Datatype array_of_types[], MPI_Datatype *newtype);
typedef int (QMPI_Type_create_struct_c_t) (QMPI_Context context, int tool_id, MPI_Count count,
             const MPI_Count array_of_blocklengths[], const MPI_Count array_of_displacements[],
             const MPI_Datatype array_of_types[], MPI_Datatype *newtype);
typedef int (QMPI_Type_struct_t) (QMPI_Context context, int tool_id, int count,
             int array_of_blocklengths[], MPI_Aint array_of_displacements[],
             MPI_Datatype array_of_types[], MPI_Datatype *newtype);
typedef int (QMPI_Type_create_subarray_t) (QMPI_Context context, int tool_id, int ndims,
             const int array_of_sizes[], const int array_of_subsizes[], const int array_of_starts[],
             int order, MPI_Datatype oldtype, MPI_Datatype *newtype);
typedef int (QMPI_Type_create_subarray_c_t) (QMPI_Context context, int tool_id, int ndims,
             const MPI_Count array_of_sizes[], const MPI_Count array_of_subsizes[],
             const MPI_Count array_of_starts[], int order, MPI_Datatype oldtype,
             MPI_Datatype *newtype);
typedef int (QMPI_Type_dup_t) (QMPI_Context context, int tool_id, MPI_Datatype oldtype,
             MPI_Datatype *newtype);
typedef int (QMPI_Type_free_t) (QMPI_Context context, int tool_id, MPI_Datatype *datatype);
typedef int (QMPI_Type_get_contents_t) (QMPI_Context context, int tool_id, MPI_Datatype datatype,
             int max_integers, int max_addresses, int max_datatypes, int array_of_integers[],
             MPI_Aint array_of_addresses[], MPI_Datatype array_of_datatypes[]);
typedef int (QMPI_Type_get_contents_c_t) (QMPI_Context context, int tool_id, MPI_Datatype datatype,
             MPI_Count max_integers, MPI_Count max_addresses, MPI_Count max_large_counts,
             MPI_Count max_datatypes, int array_of_integers[], MPI_Aint array_of_addresses[],
             MPI_Count array_of_large_counts[], MPI_Datatype array_of_datatypes[]);
typedef int (QMPI_Type_get_envelope_t) (QMPI_Context context, int tool_id, MPI_Datatype datatype,
             int *num_integers, int *num_addresses, int *num_datatypes, int *combiner);
typedef int (QMPI_Type_get_envelope_c_t) (QMPI_Context context, int tool_id, MPI_Datatype datatype,
             MPI_Count *num_integers, MPI_Count *num_addresses, MPI_Count *num_large_counts,
             MPI_Count *num_datatypes, int *combiner);
typedef int (QMPI_Type_get_extent_t) (QMPI_Context context, int tool_id, MPI_Datatype datatype,
             MPI_Aint *lb, MPI_Aint *extent);
typedef int (QMPI_Type_get_extent_c_t) (QMPI_Context context, int tool_id, MPI_Datatype datatype,
             MPI_Count *lb, MPI_Count *extent);
typedef int (QMPI_Type_get_extent_x_t) (QMPI_Context context, int tool_id, MPI_Datatype datatype,
             MPI_Count *lb, MPI_Count *extent);
typedef int (QMPI_Type_get_name_t) (QMPI_Context context, int tool_id, MPI_Datatype datatype,
             char *type_name, int *resultlen);
typedef int (QMPI_Type_get_true_extent_t) (QMPI_Context context, int tool_id, MPI_Datatype datatype,
             MPI_Aint *true_lb, MPI_Aint *true_extent);
typedef int (QMPI_Type_get_true_extent_c_t) (QMPI_Context context, int tool_id,
             MPI_Datatype datatype, MPI_Count *true_lb, MPI_Count *true_extent);
typedef int (QMPI_Type_get_true_extent_x_t) (QMPI_Context context, int tool_id,
             MPI_Datatype datatype, MPI_Count *true_lb, MPI_Count *true_extent);
typedef int (QMPI_Type_get_value_index_t) (QMPI_Context context, int tool_id,
             MPI_Datatype value_type, MPI_Datatype index_type, MPI_Datatype *pair_type);
typedef int (QMPI_Type_indexed_t) (QMPI_Context context, int tool_id, int count,
             const int array_of_blocklengths[], const int array_of_displacements[],
             MPI_Datatype oldtype, MPI_Datatype *newtype);
typedef int (QMPI_Type_indexed_c_t) (QMPI_Context context, int tool_id, MPI_Count count,
             const MPI_Count array_of_blocklengths[], const MPI_Count array_of_displacements[],
             MPI_Datatype oldtype, MPI_Datatype *newtype);
typedef int (QMPI_Type_match_size_t) (QMPI_Context context, int tool_id, int typeclass, int size,
             MPI_Datatype *datatype);
typedef int (QMPI_Type_set_name_t) (QMPI_Context context, int tool_id, MPI_Datatype datatype,
             const char *type_name);
typedef int (QMPI_Type_size_t) (QMPI_Context context, int tool_id, MPI_Datatype datatype,
             int *size);
typedef int (QMPI_Type_size_c_t) (QMPI_Context context, int tool_id, MPI_Datatype datatype,
             MPI_Count *size);
typedef int (QMPI_Type_size_x_t) (QMPI_Context context, int tool_id, MPI_Datatype datatype,
             MPI_Count *size);
typedef int (QMPI_Type_vector_t) (QMPI_Context context, int tool_id, int count, int blocklength,
             int stride, MPI_Datatype oldtype, MPI_Datatype *newtype);
typedef int (QMPI_Type_vector_c_t) (QMPI_Context context, int tool_id, MPI_Count count,
             MPI_Count blocklength, MPI_Count stride, MPI_Datatype oldtype, MPI_Datatype *newtype);
typedef int (QMPI_Unpack_t) (QMPI_Context context, int tool_id, const void *inbuf, int insize,
             int *position, void *outbuf, int outcount, MPI_Datatype datatype, MPI_Comm comm);
typedef int (QMPI_Unpack_c_t) (QMPI_Context context, int tool_id, const void *inbuf,
             MPI_Count insize, MPI_Count *position, void *outbuf, MPI_Count outcount,
             MPI_Datatype datatype, MPI_Comm comm);
typedef int (QMPI_Unpack_external_t) (QMPI_Context context, int tool_id, const char datarep[],
             const void *inbuf, MPI_Aint insize, MPI_Aint *position, void *outbuf, int outcount,
             MPI_Datatype datatype);
typedef int (QMPI_Unpack_external_c_t) (QMPI_Context context, int tool_id, const char datarep[],
             const void *inbuf, MPI_Count insize, MPI_Count *position, void *outbuf,
             MPI_Count outcount, MPI_Datatype datatype);
typedef int (QMPI_Type_extent_t) (QMPI_Context context, int tool_id, MPI_Datatype datatype,
             MPI_Aint *extent);
typedef int (QMPI_Type_lb_t) (QMPI_Context context, int tool_id, MPI_Datatype datatype,
             MPI_Aint *displacement);
typedef int (QMPI_Type_ub_t) (QMPI_Context context, int tool_id, MPI_Datatype datatype,
             MPI_Aint *displacement);
typedef int (QMPIX_Type_iov_len_t) (QMPI_Context context, int tool_id, MPI_Datatype datatype,
             MPI_Count max_iov_bytes, MPI_Count *iov_len, MPI_Count *actual_iov_bytes);
typedef int (QMPIX_Type_iov_t) (QMPI_Context context, int tool_id, MPI_Datatype datatype,
             MPI_Count iov_offset, MPIX_Iov *iov, MPI_Count max_iov_len,
             MPI_Count *actual_iov_len);
typedef int (QMPI_Add_error_class_t) (QMPI_Context context, int tool_id, int *errorclass);
typedef int (QMPI_Add_error_code_t) (QMPI_Context context, int tool_id, int errorclass,
             int *errorcode);
typedef int (QMPI_Add_error_string_t) (QMPI_Context context, int tool_id, int errorcode,
             const char *string);
typedef int (QMPI_Comm_call_errhandler_t) (QMPI_Context context, int tool_id, MPI_Comm comm,
             int errorcode);
typedef int (QMPI_Comm_create_errhandler_t) (QMPI_Context context, int tool_id,
             MPI_Comm_errhandler_function *comm_errhandler_fn, MPI_Errhandler *errhandler);
typedef int (QMPI_Errhandler_create_t) (QMPI_Context context, int tool_id,
             MPI_Comm_errhandler_function *comm_errhandler_fn, MPI_Errhandler *errhandler);
typedef int (QMPI_Comm_get_errhandler_t) (QMPI_Context context, int tool_id, MPI_Comm comm,
             MPI_Errhandler *errhandler);
typedef int (QMPI_Errhandler_get_t) (QMPI_Context context, int tool_id, MPI_Comm comm,
             MPI_Errhandler *errhandler);
typedef int (QMPI_Comm_set_errhandler_t) (QMPI_Context context, int tool_id, MPI_Comm comm,
             MPI_Errhandler errhandler);
typedef int (QMPI_Errhandler_set_t) (QMPI_Context context, int tool_id, MPI_Comm comm,
             MPI_Errhandler errhandler);
typedef int (QMPI_Errhandler_free_t) (QMPI_Context context, int tool_id,
             MPI_Errhandler *errhandler);
typedef int (QMPI_Error_class_t) (QMPI_Context context, int tool_id, int errorcode,
             int *errorclass);
typedef int (QMPI_Error_string_t) (QMPI_Context context, int tool_id, int errorcode, char *string,
             int *resultlen);
typedef int (QMPI_File_call_errhandler_t) (QMPI_Context context, int tool_id, MPI_File fh,
             int errorcode);
typedef int (QMPI_File_create_errhandler_t) (QMPI_Context context, int tool_id,
             MPI_File_errhandler_function *file_errhandler_fn, MPI_Errhandler *errhandler);
typedef int (QMPI_File_get_errhandler_t) (QMPI_Context context, int tool_id, MPI_File file,
             MPI_Errhandler *errhandler);
typedef int (QMPI_File_set_errhandler_t) (QMPI_Context context, int tool_id, MPI_File file,
             MPI_Errhandler errhandler);
typedef int (QMPI_Remove_error_class_t) (QMPI_Context context, int tool_id, int errorclass);
typedef int (QMPI_Remove_error_code_t) (QMPI_Context context, int tool_id, int errorcode);
typedef int (QMPI_Remove_error_string_t) (QMPI_Context context, int tool_id, int errorcode);
typedef int (QMPI_Session_call_errhandler_t) (QMPI_Context context, int tool_id,
             MPI_Session session, int errorcode);
typedef int (QMPI_Session_create_errhandler_t) (QMPI_Context context, int tool_id,
             MPI_Session_errhandler_function *session_errhandler_fn, MPI_Errhandler *errhandler);
typedef int (QMPI_Session_get_errhandler_t) (QMPI_Context context, int tool_id, MPI_Session session,
             MPI_Errhandler *errhandler);
typedef int (QMPI_Session_set_errhandler_t) (QMPI_Context context, int tool_id, MPI_Session session,
             MPI_Errhandler errhandler);
typedef int (QMPI_Win_call_errhandler_t) (QMPI_Context context, int tool_id, MPI_Win win,
             int errorcode);
typedef int (QMPI_Win_create_errhandler_t) (QMPI_Context context, int tool_id,
             MPI_Win_errhandler_function *win_errhandler_fn, MPI_Errhandler *errhandler);
typedef int (QMPI_Win_get_errhandler_t) (QMPI_Context context, int tool_id, MPI_Win win,
             MPI_Errhandler *errhandler);
typedef int (QMPI_Win_set_errhandler_t) (QMPI_Context context, int tool_id, MPI_Win win,
             MPI_Errhandler errhandler);
typedef MPI_Fint (QMPI_Comm_c2f_t) (QMPI_Context context, int tool_id, MPI_Comm comm);
typedef MPI_Comm (QMPI_Comm_f2c_t) (QMPI_Context context, int tool_id, MPI_Fint comm);
typedef MPI_Fint (QMPI_Errhandler_c2f_t) (QMPI_Context context, int tool_id,
                  MPI_Errhandler errhandler);
typedef MPI_Errhandler (QMPI_Errhandler_f2c_t) (QMPI_Context context, int tool_id,
                        MPI_Fint errhandler);
typedef MPI_Fint (QMPI_Group_c2f_t) (QMPI_Context context, int tool_id, MPI_Group group);
typedef MPI_Group (QMPI_Group_f2c_t) (QMPI_Context context, int tool_id, MPI_Fint group);
typedef MPI_Fint (QMPI_Info_c2f_t) (QMPI_Context context, int tool_id, MPI_Info info);
typedef MPI_Info (QMPI_Info_f2c_t) (QMPI_Context context, int tool_id, MPI_Fint info);
typedef MPI_Fint (QMPI_Message_c2f_t) (QMPI_Context context, int tool_id, MPI_Message message);
typedef MPI_Message (QMPI_Message_f2c_t) (QMPI_Context context, int tool_id, MPI_Fint message);
typedef MPI_Fint (QMPI_Op_c2f_t) (QMPI_Context context, int tool_id, MPI_Op op);
typedef MPI_Op (QMPI_Op_f2c_t) (QMPI_Context context, int tool_id, MPI_Fint op);
typedef MPI_Fint (QMPI_Request_c2f_t) (QMPI_Context context, int tool_id, MPI_Request request);
typedef MPI_Request (QMPI_Request_f2c_t) (QMPI_Context context, int tool_id, MPI_Fint request);
typedef MPI_Fint (QMPI_Session_c2f_t) (QMPI_Context context, int tool_id, MPI_Session session);
typedef MPI_Session (QMPI_Session_f2c_t) (QMPI_Context context, int tool_id, MPI_Fint session);
typedef int (QMPI_Status_c2f_t) (QMPI_Context context, int tool_id, const MPI_Status *c_status,
             MPI_Fint *f_status);
typedef int (QMPI_Status_c2f08_t) (QMPI_Context context, int tool_id, const MPI_Status *c_status,
             MPI_F08_status *f08_status);
typedef int (QMPI_Status_f082c_t) (QMPI_Context context, int tool_id,
             const MPI_F08_status *f08_status, MPI_Status *c_status);
typedef int (QMPI_Status_f082f_t) (QMPI_Context context, int tool_id,
             const MPI_F08_status *f08_status, MPI_Fint *f_status);
typedef int (QMPI_Status_f2c_t) (QMPI_Context context, int tool_id, const MPI_Fint *f_status,
             MPI_Status *c_status);
typedef int (QMPI_Status_f2f08_t) (QMPI_Context context, int tool_id, const MPI_Fint *f_status,
             MPI_F08_status *f08_status);
typedef MPI_Fint (QMPI_Type_c2f_t) (QMPI_Context context, int tool_id, MPI_Datatype datatype);
typedef MPI_Datatype (QMPI_Type_f2c_t) (QMPI_Context context, int tool_id, MPI_Fint datatype);
typedef MPI_Fint (QMPI_Win_c2f_t) (QMPI_Context context, int tool_id, MPI_Win win);
typedef MPI_Win (QMPI_Win_f2c_t) (QMPI_Context context, int tool_id, MPI_Fint win);
typedef int (QMPI_Group_compare_t) (QMPI_Context context, int tool_id, MPI_Group group1,
             MPI_Group group2, int *result);
typedef int (QMPI_Group_difference_t) (QMPI_Context context, int tool_id, MPI_Group group1,
             MPI_Group group2, MPI_Group *newgroup);
typedef int (QMPI_Group_excl_t) (QMPI_Context context, int tool_id, MPI_Group group, int n,
             const int ranks[], MPI_Group *newgroup);
typedef int (QMPI_Group_free_t) (QMPI_Context context, int tool_id, MPI_Group *group);
typedef int (QMPI_Group_incl_t) (QMPI_Context context, int tool_id, MPI_Group group, int n,
             const int ranks[], MPI_Group *newgroup);
typedef int (QMPI_Group_intersection_t) (QMPI_Context context, int tool_id, MPI_Group group1,
             MPI_Group group2, MPI_Group *newgroup);
typedef int (QMPI_Group_range_excl_t) (QMPI_Context context, int tool_id, MPI_Group group, int n,
             int ranges[][3], MPI_Group *newgroup);
typedef int (QMPI_Group_range_incl_t) (QMPI_Context context, int tool_id, MPI_Group group, int n,
             int ranges[][3], MPI_Group *newgroup);
typedef int (QMPI_Group_rank_t) (QMPI_Context context, int tool_id, MPI_Group group, int *rank);
typedef int (QMPI_Group_size_t) (QMPI_Context context, int tool_id, MPI_Group group, int *size);
typedef int (QMPI_Group_translate_ranks_t) (QMPI_Context context, int tool_id, MPI_Group group1,
             int n, const int ranks1[], MPI_Group group2, int ranks2[]);
typedef int (QMPI_Group_union_t) (QMPI_Context context, int tool_id, MPI_Group group1,
             MPI_Group group2, MPI_Group *newgroup);
typedef int (QMPI_Info_create_t) (QMPI_Context context, int tool_id, MPI_Info *info);
typedef int (QMPI_Info_create_env_t) (QMPI_Context context, int tool_id, int argc, char *argv[],
             MPI_Info *info);
typedef int (QMPI_Info_delete_t) (QMPI_Context context, int tool_id, MPI_Info info,
             const char *key);
typedef int (QMPI_Info_dup_t) (QMPI_Context context, int tool_id, MPI_Info info,
             MPI_Info *newinfo);
typedef int (QMPI_Info_free_t) (QMPI_Context context, int tool_id, MPI_Info *info);
typedef int (QMPI_Info_get_t) (QMPI_Context context, int tool_id, MPI_Info info, const char *key,
             int valuelen, char *value, int *flag);
typedef int (QMPI_Info_get_nkeys_t) (QMPI_Context context, int tool_id, MPI_Info info, int *nkeys);
typedef int (QMPI_Info_get_nthkey_t) (QMPI_Context context, int tool_id, MPI_Info info, int n,
             char *key);
typedef int (QMPI_Info_get_string_t) (QMPI_Context context, int tool_id, MPI_Info info,
             const char *key, int *buflen, char *value, int *flag);
typedef int (QMPI_Info_get_valuelen_t) (QMPI_Context context, int tool_id, MPI_Info info,
             const char *key, int *valuelen, int *flag);
typedef int (QMPI_Info_set_t) (QMPI_Context context, int tool_id, MPI_Info info, const char *key,
             const char *value);
typedef int (QMPIX_Info_set_hex_t) (QMPI_Context context, int tool_id, MPI_Info info,
             const char *key, const void *value, int value_size);
typedef int (QMPI_Abort_t) (QMPI_Context context, int tool_id, MPI_Comm comm, int errorcode);
typedef int (QMPI_Comm_create_from_group_t) (QMPI_Context context, int tool_id, MPI_Group group,
             const char *stringtag, MPI_Info info, MPI_Errhandler errhandler, MPI_Comm *newcomm);
typedef int (QMPI_Finalize_t) (QMPI_Context context, int tool_id);
typedef int (QMPI_Finalized_t) (QMPI_Context context, int tool_id, int *flag);
typedef int (QMPI_Group_from_session_pset_t) (QMPI_Context context, int tool_id,
             MPI_Session session, const char *pset_name, MPI_Group *newgroup);
typedef int (QMPI_Init_t) (QMPI_Context context, int tool_id, int *argc, char ***argv);
typedef int (QMPI_Init_thread_t) (QMPI_Context context, int tool_id, int *argc, char ***argv,
             int required, int *provided);
typedef int (QMPI_Initialized_t) (QMPI_Context context, int tool_id, int *flag);
typedef int (QMPI_Is_thread_main_t) (QMPI_Context context, int tool_id, int *flag);
typedef int (QMPI_Query_thread_t) (QMPI_Context context, int tool_id, int *provided);
typedef int (QMPI_Session_finalize_t) (QMPI_Context context, int tool_id, MPI_Session *session);
typedef int (QMPI_Session_get_info_t) (QMPI_Context context, int tool_id, MPI_Session session,
             MPI_Info *info_used);
typedef int (QMPI_Session_get_nth_pset_t) (QMPI_Context context, int tool_id, MPI_Session session,
             MPI_Info info, int n, int *pset_len, char *pset_name);
typedef int (QMPI_Session_get_num_psets_t) (QMPI_Context context, int tool_id, MPI_Session session,
             MPI_Info info, int *npset_names);
typedef int (QMPI_Session_get_pset_info_t) (QMPI_Context context, int tool_id, MPI_Session session,
             const char *pset_name, MPI_Info *info);
typedef int (QMPI_Session_init_t) (QMPI_Context context, int tool_id, MPI_Info info,
             MPI_Errhandler errhandler, MPI_Session *session);
typedef MPI_Aint (QMPI_Aint_add_t) (QMPI_Context context, int tool_id, MPI_Aint base,
                  MPI_Aint disp);
typedef MPI_Aint (QMPI_Aint_diff_t) (QMPI_Context context, int tool_id, MPI_Aint addr1,
                  MPI_Aint addr2);
typedef int (QMPI_Get_library_version_t) (QMPI_Context context, int tool_id, char *version,
             int *resultlen);
typedef int (QMPI_Get_processor_name_t) (QMPI_Context context, int tool_id, char *name,
             int *resultlen);
typedef int (QMPI_Get_version_t) (QMPI_Context context, int tool_id, int *version,
             int *subversion);
typedef int (QMPI_Pcontrol_t) (QMPI_Context context, int tool_id, const int level, ...);
typedef int (QMPIX_GPU_query_support_t) (QMPI_Context context, int tool_id, int gpu_type,
             int *is_supported);
typedef int (QMPIX_Query_cuda_support_t) (QMPI_Context context, int tool_id);
typedef int (QMPIX_Query_ze_support_t) (QMPI_Context context, int tool_id);
typedef int (QMPIX_Query_hip_support_t) (QMPI_Context context, int tool_id);
typedef int (QMPI_T_category_changed_t) (QMPI_Context context, int tool_id, int *update_number);
typedef int (QMPI_T_category_get_categories_t) (QMPI_Context context, int tool_id, int cat_index,
             int len, int indices[]);
typedef int (QMPI_T_category_get_cvars_t) (QMPI_Context context, int tool_id, int cat_index,
             int len, int indices[]);
typedef int (QMPI_T_category_get_events_t) (QMPI_Context context, int tool_id, int cat_index,
             int len, int indices[]);
typedef int (QMPI_T_category_get_index_t) (QMPI_Context context, int tool_id, const char *name,
             int *cat_index);
typedef int (QMPI_T_category_get_info_t) (QMPI_Context context, int tool_id, int cat_index,
             char *name, int *name_len, char *desc, int *desc_len, int *num_cvars, int *num_pvars,
             int *num_categories);
typedef int (QMPI_T_category_get_num_t) (QMPI_Context context, int tool_id, int *num_cat);
typedef int (QMPI_T_category_get_num_events_t) (QMPI_Context context, int tool_id, int cat_index,
             int *num_events);
typedef int (QMPI_T_category_get_pvars_t) (QMPI_Context context, int tool_id, int cat_index,
             int len, int indices[]);
typedef int (QMPI_T_cvar_get_index_t) (QMPI_Context context, int tool_id, const char *name,
             int *cvar_index);
typedef int (QMPI_T_cvar_get_info_t) (QMPI_Context context, int tool_id, int cvar_index, char *name,
             int *name_len, int *verbosity, MPI_Datatype *datatype, MPI_T_enum *enumtype,
             char *desc, int *desc_len, int *bind, int *scope);
typedef int (QMPI_T_cvar_get_num_t) (QMPI_Context context, int tool_id, int *num_cvar);
typedef int (QMPI_T_cvar_handle_alloc_t) (QMPI_Context context, int tool_id, int cvar_index,
             void *obj_handle, MPI_T_cvar_handle *handle, int *count);
typedef int (QMPI_T_cvar_handle_free_t) (QMPI_Context context, int tool_id,
             MPI_T_cvar_handle *handle);
typedef int (QMPI_T_cvar_read_t) (QMPI_Context context, int tool_id, MPI_T_cvar_handle handle,
             void *buf);
typedef int (QMPI_T_cvar_write_t) (QMPI_Context context, int tool_id, MPI_T_cvar_handle handle,
             const void *buf);
typedef int (QMPI_T_enum_get_info_t) (QMPI_Context context, int tool_id, MPI_T_enum enumtype,
             int *num, char *name, int *name_len);
typedef int (QMPI_T_enum_get_item_t) (QMPI_Context context, int tool_id, MPI_T_enum enumtype,
             int indx, int *value, char *name, int *name_len);
typedef int (QMPI_T_event_callback_get_info_t) (QMPI_Context context, int tool_id,
             MPI_T_event_registration event_registration, MPI_T_cb_safety cb_safety,
             MPI_Info *info_used);
typedef int (QMPI_T_event_callback_set_info_t) (QMPI_Context context, int tool_id,
             MPI_T_event_registration event_registration, MPI_T_cb_safety cb_safety,
             MPI_Info info);
typedef int (QMPI_T_event_copy_t) (QMPI_Context context, int tool_id,
             MPI_T_event_instance event_instance, void *buffer);
typedef int (QMPI_T_event_get_index_t) (QMPI_Context context, int tool_id, const char *name,
             int *event_index);
typedef int (QMPI_T_event_get_info_t) (QMPI_Context context, int tool_id, int event_index,
             char *name, int *name_len, int *verbosity, MPI_Datatype array_of_datatypes[],
             MPI_Aint array_of_displacements[], int *num_elements, MPI_T_enum *enumtype,
             MPI_Info *info, char *desc, int *desc_len, int *bind);
typedef int (QMPI_T_event_get_num_t) (QMPI_Context context, int tool_id, int *num_events);
typedef int (QMPI_T_event_get_source_t) (QMPI_Context context, int tool_id,
             MPI_T_event_instance event_instance, int *source_index);
typedef int (QMPI_T_event_get_timestamp_t) (QMPI_Context context, int tool_id,
             MPI_T_event_instance event_instance, MPI_Count *event_timestamp);
typedef int (QMPI_T_event_handle_alloc_t) (QMPI_Context context, int tool_id, int event_index,
             void *obj_handle, MPI_Info info, MPI_T_event_registration *event_registration);
typedef int (QMPI_T_event_handle_free_t) (QMPI_Context context, int tool_id,
             MPI_T_event_registration event_registration, void *user_data,
             MPI_T_event_free_cb_function free_cb_function);
typedef int (QMPI_T_event_handle_get_info_t) (QMPI_Context context, int tool_id,
             MPI_T_event_registration event_registration, MPI_Info *info_used);
typedef int (QMPI_T_event_handle_set_info_t) (QMPI_Context context, int tool_id,
             MPI_T_event_registration event_registration, MPI_Info info);
typedef int (QMPI_T_event_read_t) (QMPI_Context context, int tool_id,
             MPI_T_event_instance event_instance, int element_index, void *buffer);
typedef int (QMPI_T_event_register_callback_t) (QMPI_Context context, int tool_id,
             MPI_T_event_registration event_registration, MPI_T_cb_safety cb_safety, MPI_Info info,
             void *user_data, MPI_T_event_cb_function event_cb_function);
typedef int (QMPI_T_event_set_dropped_handler_t) (QMPI_Context context, int tool_id,
             MPI_T_event_registration event_registration,
             MPI_T_event_dropped_cb_function dropped_cb_function);
typedef int (QMPI_T_finalize_t) (QMPI_Context context, int tool_id);
typedef int (QMPI_T_init_thread_t) (QMPI_Context context, int tool_id, int required,
             int *provided);
typedef int (QMPI_T_pvar_get_index_t) (QMPI_Context context, int tool_id, const char *name,
             int var_class, int *pvar_index);
typedef int (QMPI_T_pvar_get_info_t) (QMPI_Context context, int tool_id, int pvar_index, char *name,
             int *name_len, int *verbosity, int *var_class, MPI_Datatype *datatype,
             MPI_T_enum *enumtype, char *desc, int *desc_len, int *bind, int *readonly,
             int *continuous, int *atomic);
typedef int (QMPI_T_pvar_get_num_t) (QMPI_Context context, int tool_id, int *num_pvar);
typedef int (QMPI_T_pvar_handle_alloc_t) (QMPI_Context context, int tool_id,
             MPI_T_pvar_session session, int pvar_index, void *obj_handle,
             MPI_T_pvar_handle *handle, int *count);
typedef int (QMPI_T_pvar_handle_free_t) (QMPI_Context context, int tool_id,
             MPI_T_pvar_session session, MPI_T_pvar_handle *handle);
typedef int (QMPI_T_pvar_read_t) (QMPI_Context context, int tool_id, MPI_T_pvar_session session,
             MPI_T_pvar_handle handle, void *buf);
typedef int (QMPI_T_pvar_readreset_t) (QMPI_Context context, int tool_id,
             MPI_T_pvar_session session, MPI_T_pvar_handle handle, void *buf);
typedef int (QMPI_T_pvar_reset_t) (QMPI_Context context, int tool_id, MPI_T_pvar_session session,
             MPI_T_pvar_handle handle);
typedef int (QMPI_T_pvar_session_create_t) (QMPI_Context context, int tool_id,
             MPI_T_pvar_session *session);
typedef int (QMPI_T_pvar_session_free_t) (QMPI_Context context, int tool_id,
             MPI_T_pvar_session *session);
typedef int (QMPI_T_pvar_start_t) (QMPI_Context context, int tool_id, MPI_T_pvar_session session,
             MPI_T_pvar_handle handle);
typedef int (QMPI_T_pvar_stop_t) (QMPI_Context context, int tool_id, MPI_T_pvar_session session,
             MPI_T_pvar_handle handle);
typedef int (QMPI_T_pvar_write_t) (QMPI_Context context, int tool_id, MPI_T_pvar_session session,
             MPI_T_pvar_handle handle, const void *buf);
typedef int (QMPI_T_source_get_info_t) (QMPI_Context context, int tool_id, int source_index,
             char *name, int *name_len, char *desc, int *desc_len, MPI_T_source_order *ordering,
             MPI_Count *ticks_per_second, MPI_Count *max_ticks, MPI_Info *info);
typedef int (QMPI_T_source_get_num_t) (QMPI_Context context, int tool_id, int *num_sources);
typedef int (QMPI_T_source_get_timestamp_t) (QMPI_Context context, int tool_id, int source_index,
             MPI_Count *timestamp);
typedef int (QMPI_Op_commutative_t) (QMPI_Context context, int tool_id, MPI_Op op, int *commute);
typedef int (QMPI_Op_create_t) (QMPI_Context context, int tool_id, MPI_User_function *user_fn,
             int commute, MPI_Op *op);
typedef int (QMPI_Op_create_c_t) (QMPI_Context context, int tool_id, MPI_User_function_c *user_fn,
             int commute, MPI_Op *op);
typedef int (QMPI_Op_free_t) (QMPI_Context context, int tool_id, MPI_Op *op);
typedef int (QMPI_Parrived_t) (QMPI_Context context, int tool_id, MPI_Request request,
             int partition, int *flag);
typedef int (QMPI_Pready_t) (QMPI_Context context, int tool_id, int partition,
             MPI_Request request);
typedef int (QMPI_Pready_list_t) (QMPI_Context context, int tool_id, int length,
             const int array_of_partitions[], MPI_Request request);
typedef int (QMPI_Pready_range_t) (QMPI_Context context, int tool_id, int partition_low,
             int partition_high, MPI_Request request);
typedef int (QMPI_Precv_init_t) (QMPI_Context context, int tool_id, void *buf, int partitions,
             MPI_Count count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm,
             MPI_Info info, MPI_Request *request);
typedef int (QMPI_Psend_init_t) (QMPI_Context context, int tool_id, const void *buf, int partitions,
             MPI_Count count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm,
             MPI_Info info, MPI_Request *request);
typedef int (QMPI_Bsend_t) (QMPI_Context context, int tool_id, const void *buf, int count,
             MPI_Datatype datatype, int dest, int tag, MPI_Comm comm);
typedef int (QMPI_Bsend_c_t) (QMPI_Context context, int tool_id, const void *buf, MPI_Count count,
             MPI_Datatype datatype, int dest, int tag, MPI_Comm comm);
typedef int (QMPI_Bsend_init_t) (QMPI_Context context, int tool_id, const void *buf, int count,
             MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Request *request);
typedef int (QMPI_Bsend_init_c_t) (QMPI_Context context, int tool_id, const void *buf,
             MPI_Count count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm,
             MPI_Request *request);
typedef int (QMPI_Buffer_attach_t) (QMPI_Context context, int tool_id, void *buffer, int size);
typedef int (QMPI_Buffer_attach_c_t) (QMPI_Context context, int tool_id, void *buffer,
             MPI_Count size);
typedef int (QMPI_Buffer_detach_t) (QMPI_Context context, int tool_id, void *buffer_addr,
             int *size);
typedef int (QMPI_Buffer_detach_c_t) (QMPI_Context context, int tool_id, void *buffer_addr,
             MPI_Count *size);
typedef int (QMPI_Buffer_flush_t) (QMPI_Context context, int tool_id);
typedef int (QMPI_Buffer_iflush_t) (QMPI_Context context, int tool_id, MPI_Request *request);
typedef int (QMPI_Comm_attach_buffer_t) (QMPI_Context context, int tool_id, MPI_Comm comm,
             void *buffer, int size);
typedef int (QMPI_Comm_attach_buffer_c_t) (QMPI_Context context, int tool_id, MPI_Comm comm,
             void *buffer, MPI_Count size);
typedef int (QMPI_Comm_detach_buffer_t) (QMPI_Context context, int tool_id, MPI_Comm comm,
             void *buffer_addr, int *size);
typedef int (QMPI_Comm_detach_buffer_c_t) (QMPI_Context context, int tool_id, MPI_Comm comm,
             void *buffer_addr, MPI_Count *size);
typedef int (QMPI_Comm_flush_buffer_t) (QMPI_Context context, int tool_id, MPI_Comm comm);
typedef int (QMPI_Comm_iflush_buffer_t) (QMPI_Context context, int tool_id, MPI_Comm comm,
             MPI_Request *request);
typedef int (QMPI_Ibsend_t) (QMPI_Context context, int tool_id, const void *buf, int count,
             MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Request *request);
typedef int (QMPI_Ibsend_c_t) (QMPI_Context context, int tool_id, const void *buf, MPI_Count count,
             MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Request *request);
typedef int (QMPI_Improbe_t) (QMPI_Context context, int tool_id, int source, int tag, MPI_Comm comm,
             int *flag, MPI_Message *message, MPI_Status *status);
typedef int (QMPI_Imrecv_t) (QMPI_Context context, int tool_id, void *buf, int count,
             MPI_Datatype datatype, MPI_Message *message, MPI_Request *request);
typedef int (QMPI_Imrecv_c_t) (QMPI_Context context, int tool_id, void *buf, MPI_Count count,
             MPI_Datatype datatype, MPI_Message *message, MPI_Request *request);
typedef int (QMPI_Iprobe_t) (QMPI_Context context, int tool_id, int source, int tag, MPI_Comm comm,
             int *flag, MPI_Status *status);
typedef int (QMPI_Irecv_t) (QMPI_Context context, int tool_id, void *buf, int count,
             MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Request *request);
typedef int (QMPI_Irecv_c_t) (QMPI_Context context, int tool_id, void *buf, MPI_Count count,
             MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Request *request);
typedef int (QMPI_Irsend_t) (QMPI_Context context, int tool_id, const void *buf, int count,
             MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Request *request);
typedef int (QMPI_Irsend_c_t) (QMPI_Context context, int tool_id, const void *buf, MPI_Count count,
             MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Request *request);
typedef int (QMPI_Isend_t) (QMPI_Context context, int tool_id, const void *buf, int count,
             MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Request *request);
typedef int (QMPI_Isend_c_t) (QMPI_Context context, int tool_id, const void *buf, MPI_Count count,
             MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Request *request);
typedef int (QMPI_Isendrecv_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             int sendcount, MPI_Datatype sendtype, int dest, int sendtag, void *recvbuf,
             int recvcount, MPI_Datatype recvtype, int source, int recvtag, MPI_Comm comm,
             MPI_Request *request);
typedef int (QMPI_Isendrecv_c_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             MPI_Count sendcount, MPI_Datatype sendtype, int dest, int sendtag, void *recvbuf,
             MPI_Count recvcount, MPI_Datatype recvtype, int source, int recvtag, MPI_Comm comm,
             MPI_Request *request);
typedef int (QMPI_Isendrecv_replace_t) (QMPI_Context context, int tool_id, void *buf, int count,
             MPI_Datatype datatype, int dest, int sendtag, int source, int recvtag, MPI_Comm comm,
             MPI_Request *request);
typedef int (QMPI_Isendrecv_replace_c_t) (QMPI_Context context, int tool_id, void *buf,
             MPI_Count count, MPI_Datatype datatype, int dest, int sendtag, int source, int recvtag,
             MPI_Comm comm, MPI_Request *request);
typedef int (QMPI_Issend_t) (QMPI_Context context, int tool_id, const void *buf, int count,
             MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Request *request);
typedef int (QMPI_Issend_c_t) (QMPI_Context context, int tool_id, const void *buf, MPI_Count count,
             MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Request *request);
typedef int (QMPI_Mprobe_t) (QMPI_Context context, int tool_id, int source, int tag, MPI_Comm comm,
             MPI_Message *message, MPI_Status *status);
typedef int (QMPI_Mrecv_t) (QMPI_Context context, int tool_id, void *buf, int count,
             MPI_Datatype datatype, MPI_Message *message, MPI_Status *status);
typedef int (QMPI_Mrecv_c_t) (QMPI_Context context, int tool_id, void *buf, MPI_Count count,
             MPI_Datatype datatype, MPI_Message *message, MPI_Status *status);
typedef int (QMPI_Probe_t) (QMPI_Context context, int tool_id, int source, int tag, MPI_Comm comm,
             MPI_Status *status);
typedef int (QMPI_Recv_t) (QMPI_Context context, int tool_id, void *buf, int count,
             MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Status *status);
typedef int (QMPI_Recv_c_t) (QMPI_Context context, int tool_id, void *buf, MPI_Count count,
             MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Status *status);
typedef int (QMPI_Recv_init_t) (QMPI_Context context, int tool_id, void *buf, int count,
             MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Request *request);
typedef int (QMPI_Recv_init_c_t) (QMPI_Context context, int tool_id, void *buf, MPI_Count count,
             MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Request *request);
typedef int (QMPI_Rsend_t) (QMPI_Context context, int tool_id, const void *buf, int count,
             MPI_Datatype datatype, int dest, int tag, MPI_Comm comm);
typedef int (QMPI_Rsend_c_t) (QMPI_Context context, int tool_id, const void *buf, MPI_Count count,
             MPI_Datatype datatype, int dest, int tag, MPI_Comm comm);
typedef int (QMPI_Rsend_init_t) (QMPI_Context context, int tool_id, const void *buf, int count,
             MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Request *request);
typedef int (QMPI_Rsend_init_c_t) (QMPI_Context context, int tool_id, const void *buf,
             MPI_Count count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm,
             MPI_Request *request);
typedef int (QMPI_Send_t) (QMPI_Context context, int tool_id, const void *buf, int count,
             MPI_Datatype datatype, int dest, int tag, MPI_Comm comm);
typedef int (QMPI_Send_c_t) (QMPI_Context context, int tool_id, const void *buf, MPI_Count count,
             MPI_Datatype datatype, int dest, int tag, MPI_Comm comm);
typedef int (QMPI_Send_init_t) (QMPI_Context context, int tool_id, const void *buf, int count,
             MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Request *request);
typedef int (QMPI_Send_init_c_t) (QMPI_Context context, int tool_id, const void *buf,
             MPI_Count count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm,
             MPI_Request *request);
typedef int (QMPI_Sendrecv_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             int sendcount, MPI_Datatype sendtype, int dest, int sendtag, void *recvbuf,
             int recvcount, MPI_Datatype recvtype, int source, int recvtag, MPI_Comm comm,
             MPI_Status *status);
typedef int (QMPI_Sendrecv_c_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             MPI_Count sendcount, MPI_Datatype sendtype, int dest, int sendtag, void *recvbuf,
             MPI_Count recvcount, MPI_Datatype recvtype, int source, int recvtag, MPI_Comm comm,
             MPI_Status *status);
typedef int (QMPI_Sendrecv_replace_t) (QMPI_Context context, int tool_id, void *buf, int count,
             MPI_Datatype datatype, int dest, int sendtag, int source, int recvtag, MPI_Comm comm,
             MPI_Status *status);
typedef int (QMPI_Sendrecv_replace_c_t) (QMPI_Context context, int tool_id, void *buf,
             MPI_Count count, MPI_Datatype datatype, int dest, int sendtag, int source, int recvtag,
             MPI_Comm comm, MPI_Status *status);
typedef int (QMPI_Session_attach_buffer_t) (QMPI_Context context, int tool_id, MPI_Session session,
             void *buffer, int size);
typedef int (QMPI_Session_attach_buffer_c_t) (QMPI_Context context, int tool_id,
             MPI_Session session, void *buffer, MPI_Count size);
typedef int (QMPI_Session_detach_buffer_t) (QMPI_Context context, int tool_id, MPI_Session session,
             void *buffer_addr, int *size);
typedef int (QMPI_Session_detach_buffer_c_t) (QMPI_Context context, int tool_id,
             MPI_Session session, void *buffer_addr, MPI_Count *size);
typedef int (QMPI_Session_flush_buffer_t) (QMPI_Context context, int tool_id, MPI_Session session);
typedef int (QMPI_Session_iflush_buffer_t) (QMPI_Context context, int tool_id, MPI_Session session,
             MPI_Request *request);
typedef int (QMPI_Ssend_t) (QMPI_Context context, int tool_id, const void *buf, int count,
             MPI_Datatype datatype, int dest, int tag, MPI_Comm comm);
typedef int (QMPI_Ssend_c_t) (QMPI_Context context, int tool_id, const void *buf, MPI_Count count,
             MPI_Datatype datatype, int dest, int tag, MPI_Comm comm);
typedef int (QMPI_Ssend_init_t) (QMPI_Context context, int tool_id, const void *buf, int count,
             MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Request *request);
typedef int (QMPI_Ssend_init_c_t) (QMPI_Context context, int tool_id, const void *buf,
             MPI_Count count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm,
             MPI_Request *request);
typedef int (QMPI_Cancel_t) (QMPI_Context context, int tool_id, MPI_Request *request);
typedef int (QMPI_Grequest_complete_t) (QMPI_Context context, int tool_id, MPI_Request request);
typedef int (QMPI_Grequest_start_t) (QMPI_Context context, int tool_id,
             MPI_Grequest_query_function *query_fn, MPI_Grequest_free_function *free_fn,
             MPI_Grequest_cancel_function *cancel_fn, void *extra_state, MPI_Request *request);
typedef int (QMPI_Request_free_t) (QMPI_Context context, int tool_id, MPI_Request *request);
typedef int (QMPI_Request_get_status_t) (QMPI_Context context, int tool_id, MPI_Request request,
             int *flag, MPI_Status *status);
typedef int (QMPI_Request_get_status_all_t) (QMPI_Context context, int tool_id, int count,
             const MPI_Request array_of_requests[], int *flag, MPI_Status *array_of_statuses);
typedef int (QMPI_Request_get_status_any_t) (QMPI_Context context, int tool_id, int count,
             const MPI_Request array_of_requests[], int *indx, int *flag, MPI_Status *status);
typedef int (QMPI_Request_get_status_some_t) (QMPI_Context context, int tool_id, int incount,
             const MPI_Request array_of_requests[], int *outcount, int array_of_indices[],
             MPI_Status *array_of_statuses);
typedef int (QMPI_Start_t) (QMPI_Context context, int tool_id, MPI_Request *request);
typedef int (QMPI_Startall_t) (QMPI_Context context, int tool_id, int count,
             MPI_Request array_of_requests[]);
typedef int (QMPI_Status_get_error_t) (QMPI_Context context, int tool_id, const MPI_Status *status,
             int *error);
typedef int (QMPI_Status_get_source_t) (QMPI_Context context, int tool_id, const MPI_Status *status,
             int *source);
typedef int (QMPI_Status_get_tag_t) (QMPI_Context context, int tool_id, const MPI_Status *status,
             int *tag);
typedef int (QMPI_Status_set_error_t) (QMPI_Context context, int tool_id, MPI_Status *status,
             int error);
typedef int (QMPI_Status_set_source_t) (QMPI_Context context, int tool_id, MPI_Status *status,
             int source);
typedef int (QMPI_Status_set_tag_t) (QMPI_Context context, int tool_id, MPI_Status *status,
             int tag);
typedef int (QMPI_Status_set_cancelled_t) (QMPI_Context context, int tool_id, MPI_Status *status,
             int flag);
typedef int (QMPI_Test_t) (QMPI_Context context, int tool_id, MPI_Request *request, int *flag,
             MPI_Status *status);
typedef int (QMPI_Test_cancelled_t) (QMPI_Context context, int tool_id, const MPI_Status *status,
             int *flag);
typedef int (QMPI_Testall_t) (QMPI_Context context, int tool_id, int count,
             MPI_Request array_of_requests[], int *flag, MPI_Status *array_of_statuses);
typedef int (QMPI_Testany_t) (QMPI_Context context, int tool_id, int count,
             MPI_Request array_of_requests[], int *indx, int *flag, MPI_Status *status);
typedef int (QMPI_Testsome_t) (QMPI_Context context, int tool_id, int incount,
             MPI_Request array_of_requests[], int *outcount, int array_of_indices[],
             MPI_Status *array_of_statuses);
typedef int (QMPI_Wait_t) (QMPI_Context context, int tool_id, MPI_Request *request,
             MPI_Status *status);
typedef int (QMPI_Waitall_t) (QMPI_Context context, int tool_id, int count,
             MPI_Request array_of_requests[], MPI_Status *array_of_statuses);
typedef int (QMPI_Waitany_t) (QMPI_Context context, int tool_id, int count,
             MPI_Request array_of_requests[], int *indx, MPI_Status *status);
typedef int (QMPI_Waitsome_t) (QMPI_Context context, int tool_id, int incount,
             MPI_Request array_of_requests[], int *outcount, int array_of_indices[],
             MPI_Status *array_of_statuses);
typedef int (QMPIX_Grequest_start_t) (QMPI_Context context, int tool_id,
             MPI_Grequest_query_function *query_fn, MPI_Grequest_free_function *free_fn,
             MPI_Grequest_cancel_function *cancel_fn, MPIX_Grequest_poll_function *poll_fn,
             MPIX_Grequest_wait_function *wait_fn, void *extra_state, MPI_Request *request);
typedef int (QMPIX_Grequest_class_create_t) (QMPI_Context context, int tool_id,
             MPI_Grequest_query_function *query_fn, MPI_Grequest_free_function *free_fn,
             MPI_Grequest_cancel_function *cancel_fn, MPIX_Grequest_poll_function *poll_fn,
             MPIX_Grequest_wait_function *wait_fn, MPIX_Grequest_class *greq_class);
typedef int (QMPIX_Grequest_class_allocate_t) (QMPI_Context context, int tool_id,
             MPIX_Grequest_class greq_class, void *extra_state, MPI_Request *request);
typedef int (QMPIX_Request_is_complete_t) (QMPI_Context context, int tool_id, MPI_Request request);
typedef int (QMPI_Accumulate_t) (QMPI_Context context, int tool_id, const void *origin_addr,
             int origin_count, MPI_Datatype origin_datatype, int target_rank, MPI_Aint target_disp,
             int target_count, MPI_Datatype target_datatype, MPI_Op op, MPI_Win win);
typedef int (QMPI_Accumulate_c_t) (QMPI_Context context, int tool_id, const void *origin_addr,
             MPI_Count origin_count, MPI_Datatype origin_datatype, int target_rank,
             MPI_Aint target_disp, MPI_Count target_count, MPI_Datatype target_datatype, MPI_Op op,
             MPI_Win win);
typedef int (QMPI_Alloc_mem_t) (QMPI_Context context, int tool_id, MPI_Aint size, MPI_Info info,
             void *baseptr);
typedef int (QMPI_Compare_and_swap_t) (QMPI_Context context, int tool_id, const void *origin_addr,
             const void *compare_addr, void *result_addr, MPI_Datatype datatype, int target_rank,
             MPI_Aint target_disp, MPI_Win win);
typedef int (QMPI_Fetch_and_op_t) (QMPI_Context context, int tool_id, const void *origin_addr,
             void *result_addr, MPI_Datatype datatype, int target_rank, MPI_Aint target_disp,
             MPI_Op op, MPI_Win win);
typedef int (QMPI_Free_mem_t) (QMPI_Context context, int tool_id, void *base);
typedef int (QMPI_Get_t) (QMPI_Context context, int tool_id, void *origin_addr, int origin_count,
             MPI_Datatype origin_datatype, int target_rank, MPI_Aint target_disp, int target_count,
             MPI_Datatype target_datatype, MPI_Win win);
typedef int (QMPI_Get_c_t) (QMPI_Context context, int tool_id, void *origin_addr,
             MPI_Count origin_count, MPI_Datatype origin_datatype, int target_rank,
             MPI_Aint target_disp, MPI_Count target_count, MPI_Datatype target_datatype,
             MPI_Win win);
typedef int (QMPI_Get_accumulate_t) (QMPI_Context context, int tool_id, const void *origin_addr,
             int origin_count, MPI_Datatype origin_datatype, void *result_addr, int result_count,
             MPI_Datatype result_datatype, int target_rank, MPI_Aint target_disp, int target_count,
             MPI_Datatype target_datatype, MPI_Op op, MPI_Win win);
typedef int (QMPI_Get_accumulate_c_t) (QMPI_Context context, int tool_id, const void *origin_addr,
             MPI_Count origin_count, MPI_Datatype origin_datatype, void *result_addr,
             MPI_Count result_count, MPI_Datatype result_datatype, int target_rank,
             MPI_Aint target_disp, MPI_Count target_count, MPI_Datatype target_datatype, MPI_Op op,
             MPI_Win win);
typedef int (QMPI_Put_t) (QMPI_Context context, int tool_id, const void *origin_addr,
             int origin_count, MPI_Datatype origin_datatype, int target_rank, MPI_Aint target_disp,
             int target_count, MPI_Datatype target_datatype, MPI_Win win);
typedef int (QMPI_Put_c_t) (QMPI_Context context, int tool_id, const void *origin_addr,
             MPI_Count origin_count, MPI_Datatype origin_datatype, int target_rank,
             MPI_Aint target_disp, MPI_Count target_count, MPI_Datatype target_datatype,
             MPI_Win win);
typedef int (QMPI_Raccumulate_t) (QMPI_Context context, int tool_id, const void *origin_addr,
             int origin_count, MPI_Datatype origin_datatype, int target_rank, MPI_Aint target_disp,
             int target_count, MPI_Datatype target_datatype, MPI_Op op, MPI_Win win,
             MPI_Request *request);
typedef int (QMPI_Raccumulate_c_t) (QMPI_Context context, int tool_id, const void *origin_addr,
             MPI_Count origin_count, MPI_Datatype origin_datatype, int target_rank,
             MPI_Aint target_disp, MPI_Count target_count, MPI_Datatype target_datatype, MPI_Op op,
             MPI_Win win, MPI_Request *request);
typedef int (QMPI_Rget_t) (QMPI_Context context, int tool_id, void *origin_addr, int origin_count,
             MPI_Datatype origin_datatype, int target_rank, MPI_Aint target_disp, int target_count,
             MPI_Datatype target_datatype, MPI_Win win, MPI_Request *request);
typedef int (QMPI_Rget_c_t) (QMPI_Context context, int tool_id, void *origin_addr,
             MPI_Count origin_count, MPI_Datatype origin_datatype, int target_rank,
             MPI_Aint target_disp, MPI_Count target_count, MPI_Datatype target_datatype,
             MPI_Win win, MPI_Request *request);
typedef int (QMPI_Rget_accumulate_t) (QMPI_Context context, int tool_id, const void *origin_addr,
             int origin_count, MPI_Datatype origin_datatype, void *result_addr, int result_count,
             MPI_Datatype result_datatype, int target_rank, MPI_Aint target_disp, int target_count,
             MPI_Datatype target_datatype, MPI_Op op, MPI_Win win, MPI_Request *request);
typedef int (QMPI_Rget_accumulate_c_t) (QMPI_Context context, int tool_id, const void *origin_addr,
             MPI_Count origin_count, MPI_Datatype origin_datatype, void *result_addr,
             MPI_Count result_count, MPI_Datatype result_datatype, int target_rank,
             MPI_Aint target_disp, MPI_Count target_count, MPI_Datatype target_datatype, MPI_Op op,
             MPI_Win win, MPI_Request *request);
typedef int (QMPI_Rput_t) (QMPI_Context context, int tool_id, const void *origin_addr,
             int origin_count, MPI_Datatype origin_datatype, int target_rank, MPI_Aint target_disp,
             int target_count, MPI_Datatype target_datatype, MPI_Win win, MPI_Request *request);
typedef int (QMPI_Rput_c_t) (QMPI_Context context, int tool_id, const void *origin_addr,
             MPI_Count origin_count, MPI_Datatype origin_datatype, int target_rank,
             MPI_Aint target_disp, MPI_Count target_count, MPI_Datatype target_datatype,
             MPI_Win win, MPI_Request *request);
typedef int (QMPI_Win_allocate_t) (QMPI_Context context, int tool_id, MPI_Aint size, int disp_unit,
             MPI_Info info, MPI_Comm comm, void *baseptr, MPI_Win *win);
typedef int (QMPI_Win_allocate_c_t) (QMPI_Context context, int tool_id, MPI_Aint size,
             MPI_Aint disp_unit, MPI_Info info, MPI_Comm comm, void *baseptr, MPI_Win *win);
typedef int (QMPI_Win_allocate_shared_t) (QMPI_Context context, int tool_id, MPI_Aint size,
             int disp_unit, MPI_Info info, MPI_Comm comm, void *baseptr, MPI_Win *win);
typedef int (QMPI_Win_allocate_shared_c_t) (QMPI_Context context, int tool_id, MPI_Aint size,
             MPI_Aint disp_unit, MPI_Info info, MPI_Comm comm, void *baseptr, MPI_Win *win);
typedef int (QMPI_Win_attach_t) (QMPI_Context context, int tool_id, MPI_Win win, void *base,
             MPI_Aint size);
typedef int (QMPI_Win_complete_t) (QMPI_Context context, int tool_id, MPI_Win win);
typedef int (QMPI_Win_create_t) (QMPI_Context context, int tool_id, void *base, MPI_Aint size,
             int disp_unit, MPI_Info info, MPI_Comm comm, MPI_Win *win);
typedef int (QMPI_Win_create_c_t) (QMPI_Context context, int tool_id, void *base, MPI_Aint size,
             MPI_Aint disp_unit, MPI_Info info, MPI_Comm comm, MPI_Win *win);
typedef int (QMPI_Win_create_dynamic_t) (QMPI_Context context, int tool_id, MPI_Info info,
             MPI_Comm comm, MPI_Win *win);
typedef int (QMPI_Win_detach_t) (QMPI_Context context, int tool_id, MPI_Win win, const void *base);
typedef int (QMPI_Win_fence_t) (QMPI_Context context, int tool_id, int assert, MPI_Win win);
typedef int (QMPI_Win_flush_t) (QMPI_Context context, int tool_id, int rank, MPI_Win win);
typedef int (QMPI_Win_flush_all_t) (QMPI_Context context, int tool_id, MPI_Win win);
typedef int (QMPI_Win_flush_local_t) (QMPI_Context context, int tool_id, int rank, MPI_Win win);
typedef int (QMPI_Win_flush_local_all_t) (QMPI_Context context, int tool_id, MPI_Win win);
typedef int (QMPI_Win_free_t) (QMPI_Context context, int tool_id, MPI_Win *win);
typedef int (QMPI_Win_get_group_t) (QMPI_Context context, int tool_id, MPI_Win win,
             MPI_Group *group);
typedef int (QMPI_Win_get_info_t) (QMPI_Context context, int tool_id, MPI_Win win,
             MPI_Info *info_used);
typedef int (QMPI_Win_get_name_t) (QMPI_Context context, int tool_id, MPI_Win win, char *win_name,
             int *resultlen);
typedef int (QMPI_Win_lock_t) (QMPI_Context context, int tool_id, int lock_type, int rank,
             int assert, MPI_Win win);
typedef int (QMPI_Win_lock_all_t) (QMPI_Context context, int tool_id, int assert, MPI_Win win);
typedef int (QMPI_Win_post_t) (QMPI_Context context, int tool_id, MPI_Group group, int assert,
             MPI_Win win);
typedef int (QMPI_Win_set_info_t) (QMPI_Context context, int tool_id, MPI_Win win, MPI_Info info);
typedef int (QMPI_Win_set_name_t) (QMPI_Context context, int tool_id, MPI_Win win,
             const char *win_name);
typedef int (QMPI_Win_shared_query_t) (QMPI_Context context, int tool_id, MPI_Win win, int rank,
             MPI_Aint *size, int *disp_unit, void *baseptr);
typedef int (QMPI_Win_shared_query_c_t) (QMPI_Context context, int tool_id, MPI_Win win, int rank,
             MPI_Aint *size, MPI_Aint *disp_unit, void *baseptr);
typedef int (QMPI_Win_start_t) (QMPI_Context context, int tool_id, MPI_Group group, int assert,
             MPI_Win win);
typedef int (QMPI_Win_sync_t) (QMPI_Context context, int tool_id, MPI_Win win);
typedef int (QMPI_Win_test_t) (QMPI_Context context, int tool_id, MPI_Win win, int *flag);
typedef int (QMPI_Win_unlock_t) (QMPI_Context context, int tool_id, int rank, MPI_Win win);
typedef int (QMPI_Win_unlock_all_t) (QMPI_Context context, int tool_id, MPI_Win win);
typedef int (QMPI_Win_wait_t) (QMPI_Context context, int tool_id, MPI_Win win);
typedef int (QMPI_Close_port_t) (QMPI_Context context, int tool_id, const char *port_name);
typedef int (QMPI_Comm_accept_t) (QMPI_Context context, int tool_id, const char *port_name,
             MPI_Info info, int root, MPI_Comm comm, MPI_Comm *newcomm);
typedef int (QMPI_Comm_connect_t) (QMPI_Context context, int tool_id, const char *port_name,
             MPI_Info info, int root, MPI_Comm comm, MPI_Comm *newcomm);
typedef int (QMPI_Comm_disconnect_t) (QMPI_Context context, int tool_id, MPI_Comm *comm);
typedef int (QMPI_Comm_get_parent_t) (QMPI_Context context, int tool_id, MPI_Comm *parent);
typedef int (QMPI_Comm_join_t) (QMPI_Context context, int tool_id, int fd, MPI_Comm *intercomm);
typedef int (QMPI_Comm_spawn_t) (QMPI_Context context, int tool_id, const char *command,
             char *argv[], int maxprocs, MPI_Info info, int root, MPI_Comm comm,
             MPI_Comm *intercomm, int array_of_errcodes[]);
typedef int (QMPI_Comm_spawn_multiple_t) (QMPI_Context context, int tool_id, int count,
             char *array_of_commands[], char **array_of_argv[], const int array_of_maxprocs[],
             const MPI_Info array_of_info[], int root, MPI_Comm comm, MPI_Comm *intercomm,
             int array_of_errcodes[]);
typedef int (QMPI_Lookup_name_t) (QMPI_Context context, int tool_id, const char *service_name,
             MPI_Info info, char *port_name);
typedef int (QMPI_Open_port_t) (QMPI_Context context, int tool_id, MPI_Info info, char *port_name);
typedef int (QMPI_Publish_name_t) (QMPI_Context context, int tool_id, const char *service_name,
             MPI_Info info, const char *port_name);
typedef int (QMPI_Unpublish_name_t) (QMPI_Context context, int tool_id, const char *service_name,
             MPI_Info info, const char *port_name);
typedef int (QMPIX_Stream_create_t) (QMPI_Context context, int tool_id, MPI_Info info,
             MPIX_Stream *stream);
typedef int (QMPIX_Stream_free_t) (QMPI_Context context, int tool_id, MPIX_Stream *stream);
typedef int (QMPIX_Stream_comm_create_t) (QMPI_Context context, int tool_id, MPI_Comm comm,
             MPIX_Stream stream, MPI_Comm *newcomm);
typedef int (QMPIX_Stream_comm_create_multiplex_t) (QMPI_Context context, int tool_id,
             MPI_Comm comm, int count, MPIX_Stream array_of_streams[], MPI_Comm *newcomm);
typedef int (QMPIX_Comm_get_stream_t) (QMPI_Context context, int tool_id, MPI_Comm comm, int idx,
             MPIX_Stream *stream);
typedef int (QMPIX_Stream_progress_t) (QMPI_Context context, int tool_id, MPIX_Stream stream);
typedef int (QMPIX_Start_progress_thread_t) (QMPI_Context context, int tool_id,
             MPIX_Stream stream);
typedef int (QMPIX_Stop_progress_thread_t) (QMPI_Context context, int tool_id, MPIX_Stream stream);
typedef int (QMPIX_Stream_send_t) (QMPI_Context context, int tool_id, const void *buf, int count,
             MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, int source_stream_index,
             int dest_stream_index);
typedef int (QMPIX_Stream_send_c_t) (QMPI_Context context, int tool_id, const void *buf,
             MPI_Count count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm,
             int source_stream_index, int dest_stream_index);
typedef int (QMPIX_Stream_isend_t) (QMPI_Context context, int tool_id, const void *buf, int count,
             MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, int source_stream_index,
             int dest_stream_index, MPI_Request *request);
typedef int (QMPIX_Stream_isend_c_t) (QMPI_Context context, int tool_id, const void *buf,
             MPI_Count count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm,
             int source_stream_index, int dest_stream_index, MPI_Request *request);
typedef int (QMPIX_Stream_recv_t) (QMPI_Context context, int tool_id, void *buf, int count,
             MPI_Datatype datatype, int source, int tag, MPI_Comm comm, int source_stream_index,
             int dest_stream_index, MPI_Status *status);
typedef int (QMPIX_Stream_recv_c_t) (QMPI_Context context, int tool_id, void *buf, MPI_Count count,
             MPI_Datatype datatype, int source, int tag, MPI_Comm comm, int source_stream_index,
             int dest_stream_index, MPI_Status *status);
typedef int (QMPIX_Stream_irecv_t) (QMPI_Context context, int tool_id, void *buf, int count,
             MPI_Datatype datatype, int source, int tag, MPI_Comm comm, int source_stream_index,
             int dest_stream_index, MPI_Request *request);
typedef int (QMPIX_Stream_irecv_c_t) (QMPI_Context context, int tool_id, void *buf, MPI_Count count,
             MPI_Datatype datatype, int source, int tag, MPI_Comm comm, int source_stream_index,
             int dest_stream_index, MPI_Request *request);
typedef int (QMPIX_Send_enqueue_t) (QMPI_Context context, int tool_id, const void *buf, int count,
             MPI_Datatype datatype, int dest, int tag, MPI_Comm comm);
typedef int (QMPIX_Send_enqueue_c_t) (QMPI_Context context, int tool_id, const void *buf,
             MPI_Count count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm);
typedef int (QMPIX_Recv_enqueue_t) (QMPI_Context context, int tool_id, void *buf, int count,
             MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Status *status);
typedef int (QMPIX_Recv_enqueue_c_t) (QMPI_Context context, int tool_id, void *buf, MPI_Count count,
             MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Status *status);
typedef int (QMPIX_Isend_enqueue_t) (QMPI_Context context, int tool_id, const void *buf, int count,
             MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Request *request);
typedef int (QMPIX_Isend_enqueue_c_t) (QMPI_Context context, int tool_id, const void *buf,
             MPI_Count count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm,
             MPI_Request *request);
typedef int (QMPIX_Irecv_enqueue_t) (QMPI_Context context, int tool_id, void *buf, int count,
             MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Request *request);
typedef int (QMPIX_Irecv_enqueue_c_t) (QMPI_Context context, int tool_id, void *buf,
             MPI_Count count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm,
             MPI_Request *request);
typedef int (QMPIX_Wait_enqueue_t) (QMPI_Context context, int tool_id, MPI_Request *request,
             MPI_Status *status);
typedef int (QMPIX_Waitall_enqueue_t) (QMPI_Context context, int tool_id, int count,
             MPI_Request array_of_requests[], MPI_Status *array_of_statuses);
typedef int (QMPIX_Allreduce_enqueue_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm);
typedef int (QMPIX_Allreduce_enqueue_c_t) (QMPI_Context context, int tool_id, const void *sendbuf,
             void *recvbuf, MPI_Count count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm);
typedef int (QMPIX_Threadcomm_init_t) (QMPI_Context context, int tool_id, MPI_Comm comm,
             int num_threads, MPI_Comm *newthreadcomm);
typedef int (QMPIX_Threadcomm_free_t) (QMPI_Context context, int tool_id, MPI_Comm *threadcomm);
typedef int (QMPIX_Threadcomm_start_t) (QMPI_Context context, int tool_id, MPI_Comm threadcomm);
typedef int (QMPIX_Threadcomm_finish_t) (QMPI_Context context, int tool_id, MPI_Comm threadcomm);
typedef double (QMPI_Wtick_t) (QMPI_Context context, int tool_id);
typedef double (QMPI_Wtime_t) (QMPI_Context context, int tool_id);
typedef int (QMPI_Cart_coords_t) (QMPI_Context context, int tool_id, MPI_Comm comm, int rank,
             int maxdims, int coords[]);
typedef int (QMPI_Cart_create_t) (QMPI_Context context, int tool_id, MPI_Comm comm_old, int ndims,
             const int dims[], const int periods[], int reorder, MPI_Comm *comm_cart);
typedef int (QMPI_Cart_get_t) (QMPI_Context context, int tool_id, MPI_Comm comm, int maxdims,
             int dims[], int periods[], int coords[]);
typedef int (QMPI_Cart_map_t) (QMPI_Context context, int tool_id, MPI_Comm comm, int ndims,
             const int dims[], const int periods[], int *newrank);
typedef int (QMPI_Cart_rank_t) (QMPI_Context context, int tool_id, MPI_Comm comm,
             const int coords[], int *rank);
typedef int (QMPI_Cart_shift_t) (QMPI_Context context, int tool_id, MPI_Comm comm, int direction,
             int disp, int *rank_source, int *rank_dest);
typedef int (QMPI_Cart_sub_t) (QMPI_Context context, int tool_id, MPI_Comm comm,
             const int remain_dims[], MPI_Comm *newcomm);
typedef int (QMPI_Cartdim_get_t) (QMPI_Context context, int tool_id, MPI_Comm comm, int *ndims);
typedef int (QMPI_Dims_create_t) (QMPI_Context context, int tool_id, int nnodes, int ndims,
             int dims[]);
typedef int (QMPI_Dist_graph_create_t) (QMPI_Context context, int tool_id, MPI_Comm comm_old, int n,
             const int sources[], const int degrees[], const int destinations[],
             const int weights[], MPI_Info info, int reorder, MPI_Comm *comm_dist_graph);
typedef int (QMPI_Dist_graph_create_adjacent_t) (QMPI_Context context, int tool_id,
             MPI_Comm comm_old, int indegree, const int sources[], const int sourceweights[],
             int outdegree, const int destinations[], const int destweights[], MPI_Info info,
             int reorder, MPI_Comm *comm_dist_graph);
typedef int (QMPI_Dist_graph_neighbors_t) (QMPI_Context context, int tool_id, MPI_Comm comm,
             int maxindegree, int sources[], int sourceweights[], int maxoutdegree,
             int destinations[], int destweights[]);
typedef int (QMPI_Dist_graph_neighbors_count_t) (QMPI_Context context, int tool_id, MPI_Comm comm,
             int *indegree, int *outdegree, int *weighted);
typedef int (QMPI_Get_hw_resource_info_t) (QMPI_Context context, int tool_id, MPI_Info *hw_info);
typedef int (QMPI_Graph_create_t) (QMPI_Context context, int tool_id, MPI_Comm comm_old, int nnodes,
             const int indx[], const int edges[], int reorder, MPI_Comm *comm_graph);
typedef int (QMPI_Graph_get_t) (QMPI_Context context, int tool_id, MPI_Comm comm, int maxindex,
             int maxedges, int indx[], int edges[]);
typedef int (QMPI_Graph_map_t) (QMPI_Context context, int tool_id, MPI_Comm comm, int nnodes,
             const int indx[], const int edges[], int *newrank);
typedef int (QMPI_Graph_neighbors_t) (QMPI_Context context, int tool_id, MPI_Comm comm, int rank,
             int maxneighbors, int neighbors[]);
typedef int (QMPI_Graph_neighbors_count_t) (QMPI_Context context, int tool_id, MPI_Comm comm,
             int rank, int *nneighbors);
typedef int (QMPI_Graphdims_get_t) (QMPI_Context context, int tool_id, MPI_Comm comm, int *nnodes,
             int *nedges);
typedef int (QMPI_Topo_test_t) (QMPI_Context context, int tool_id, MPI_Comm comm, int *status);
typedef MPI_Fint (QMPI_File_c2f_t) (QMPI_Context context, int tool_id, MPI_File file);
typedef int (QMPI_File_close_t) (QMPI_Context context, int tool_id, MPI_File *fh);
typedef int (QMPI_File_delete_t) (QMPI_Context context, int tool_id, const char *filename,
             MPI_Info info);
typedef MPI_File (QMPI_File_f2c_t) (QMPI_Context context, int tool_id, MPI_Fint file);
typedef int (QMPI_File_get_amode_t) (QMPI_Context context, int tool_id, MPI_File fh, int *amode);
typedef int (QMPI_File_get_atomicity_t) (QMPI_Context context, int tool_id, MPI_File fh,
             int *flag);
typedef int (QMPI_File_get_byte_offset_t) (QMPI_Context context, int tool_id, MPI_File fh,
             MPI_Offset offset, MPI_Offset *disp);
typedef int (QMPI_File_get_group_t) (QMPI_Context context, int tool_id, MPI_File fh,
             MPI_Group *group);
typedef int (QMPI_File_get_info_t) (QMPI_Context context, int tool_id, MPI_File fh,
             MPI_Info *info_used);
typedef int (QMPI_File_get_position_t) (QMPI_Context context, int tool_id, MPI_File fh,
             MPI_Offset *offset);
typedef int (QMPI_File_get_position_shared_t) (QMPI_Context context, int tool_id, MPI_File fh,
             MPI_Offset *offset);
typedef int (QMPI_File_get_size_t) (QMPI_Context context, int tool_id, MPI_File fh,
             MPI_Offset *size);
typedef int (QMPI_File_get_type_extent_t) (QMPI_Context context, int tool_id, MPI_File fh,
             MPI_Datatype datatype, MPI_Aint *extent);
typedef int (QMPI_File_get_type_extent_c_t) (QMPI_Context context, int tool_id, MPI_File fh,
             MPI_Datatype datatype, MPI_Count *extent);
typedef int (QMPI_File_get_view_t) (QMPI_Context context, int tool_id, MPI_File fh,
             MPI_Offset *disp, MPI_Datatype *etype, MPI_Datatype *filetype, char *datarep);
typedef int (QMPI_File_iread_t) (QMPI_Context context, int tool_id, MPI_File fh, void *buf,
             int count, MPI_Datatype datatype, MPI_Request *request);
typedef int (QMPI_File_iread_c_t) (QMPI_Context context, int tool_id, MPI_File fh, void *buf,
             MPI_Count count, MPI_Datatype datatype, MPI_Request *request);
typedef int (QMPI_File_iread_all_t) (QMPI_Context context, int tool_id, MPI_File fh, void *buf,
             int count, MPI_Datatype datatype, MPI_Request *request);
typedef int (QMPI_File_iread_all_c_t) (QMPI_Context context, int tool_id, MPI_File fh, void *buf,
             MPI_Count count, MPI_Datatype datatype, MPI_Request *request);
typedef int (QMPI_File_iread_at_t) (QMPI_Context context, int tool_id, MPI_File fh,
             MPI_Offset offset, void *buf, int count, MPI_Datatype datatype, MPI_Request *request);
typedef int (QMPI_File_iread_at_c_t) (QMPI_Context context, int tool_id, MPI_File fh,
             MPI_Offset offset, void *buf, MPI_Count count, MPI_Datatype datatype,
             MPI_Request *request);
typedef int (QMPI_File_iread_at_all_t) (QMPI_Context context, int tool_id, MPI_File fh,
             MPI_Offset offset, void *buf, int count, MPI_Datatype datatype, MPI_Request *request);
typedef int (QMPI_File_iread_at_all_c_t) (QMPI_Context context, int tool_id, MPI_File fh,
             MPI_Offset offset, void *buf, MPI_Count count, MPI_Datatype datatype,
             MPI_Request *request);
typedef int (QMPI_File_iread_shared_t) (QMPI_Context context, int tool_id, MPI_File fh, void *buf,
             int count, MPI_Datatype datatype, MPI_Request *request);
typedef int (QMPI_File_iread_shared_c_t) (QMPI_Context context, int tool_id, MPI_File fh, void *buf,
             MPI_Count count, MPI_Datatype datatype, MPI_Request *request);
typedef int (QMPI_File_iwrite_t) (QMPI_Context context, int tool_id, MPI_File fh, const void *buf,
             int count, MPI_Datatype datatype, MPI_Request *request);
typedef int (QMPI_File_iwrite_c_t) (QMPI_Context context, int tool_id, MPI_File fh, const void *buf,
             MPI_Count count, MPI_Datatype datatype, MPI_Request *request);
typedef int (QMPI_File_iwrite_all_t) (QMPI_Context context, int tool_id, MPI_File fh,
             const void *buf, int count, MPI_Datatype datatype, MPI_Request *request);
typedef int (QMPI_File_iwrite_all_c_t) (QMPI_Context context, int tool_id, MPI_File fh,
             const void *buf, MPI_Count count, MPI_Datatype datatype, MPI_Request *request);
typedef int (QMPI_File_iwrite_at_t) (QMPI_Context context, int tool_id, MPI_File fh,
             MPI_Offset offset, const void *buf, int count, MPI_Datatype datatype,
             MPI_Request *request);
typedef int (QMPI_File_iwrite_at_c_t) (QMPI_Context context, int tool_id, MPI_File fh,
             MPI_Offset offset, const void *buf, MPI_Count count, MPI_Datatype datatype,
             MPI_Request *request);
typedef int (QMPI_File_iwrite_at_all_t) (QMPI_Context context, int tool_id, MPI_File fh,
             MPI_Offset offset, const void *buf, int count, MPI_Datatype datatype,
             MPI_Request *request);
typedef int (QMPI_File_iwrite_at_all_c_t) (QMPI_Context context, int tool_id, MPI_File fh,
             MPI_Offset offset, const void *buf, MPI_Count count, MPI_Datatype datatype,
             MPI_Request *request);
typedef int (QMPI_File_iwrite_shared_t) (QMPI_Context context, int tool_id, MPI_File fh,
             const void *buf, int count, MPI_Datatype datatype, MPI_Request *request);
typedef int (QMPI_File_iwrite_shared_c_t) (QMPI_Context context, int tool_id, MPI_File fh,
             const void *buf, MPI_Count count, MPI_Datatype datatype, MPI_Request *request);
typedef int (QMPI_File_open_t) (QMPI_Context context, int tool_id, MPI_Comm comm,
             const char *filename, int amode, MPI_Info info, MPI_File *fh);
typedef int (QMPI_File_preallocate_t) (QMPI_Context context, int tool_id, MPI_File fh,
             MPI_Offset size);
typedef int (QMPI_File_read_t) (QMPI_Context context, int tool_id, MPI_File fh, void *buf,
             int count, MPI_Datatype datatype, MPI_Status *status);
typedef int (QMPI_File_read_c_t) (QMPI_Context context, int tool_id, MPI_File fh, void *buf,
             MPI_Count count, MPI_Datatype datatype, MPI_Status *status);
typedef int (QMPI_File_read_all_t) (QMPI_Context context, int tool_id, MPI_File fh, void *buf,
             int count, MPI_Datatype datatype, MPI_Status *status);
typedef int (QMPI_File_read_all_c_t) (QMPI_Context context, int tool_id, MPI_File fh, void *buf,
             MPI_Count count, MPI_Datatype datatype, MPI_Status *status);
typedef int (QMPI_File_read_all_begin_t) (QMPI_Context context, int tool_id, MPI_File fh, void *buf,
             int count, MPI_Datatype datatype);
typedef int (QMPI_File_read_all_begin_c_t) (QMPI_Context context, int tool_id, MPI_File fh,
             void *buf, MPI_Count count, MPI_Datatype datatype);
typedef int (QMPI_File_read_all_end_t) (QMPI_Context context, int tool_id, MPI_File fh, void *buf,
             MPI_Status *status);
typedef int (QMPI_File_read_at_t) (QMPI_Context context, int tool_id, MPI_File fh,
             MPI_Offset offset, void *buf, int count, MPI_Datatype datatype, MPI_Status *status);
typedef int (QMPI_File_read_at_c_t) (QMPI_Context context, int tool_id, MPI_File fh,
             MPI_Offset offset, void *buf, MPI_Count count, MPI_Datatype datatype,
             MPI_Status *status);
typedef int (QMPI_File_read_at_all_t) (QMPI_Context context, int tool_id, MPI_File fh,
             MPI_Offset offset, void *buf, int count, MPI_Datatype datatype, MPI_Status *status);
typedef int (QMPI_File_read_at_all_c_t) (QMPI_Context context, int tool_id, MPI_File fh,
             MPI_Offset offset, void *buf, MPI_Count count, MPI_Datatype datatype,
             MPI_Status *status);
typedef int (QMPI_File_read_at_all_begin_t) (QMPI_Context context, int tool_id, MPI_File fh,
             MPI_Offset offset, void *buf, int count, MPI_Datatype datatype);
typedef int (QMPI_File_read_at_all_begin_c_t) (QMPI_Context context, int tool_id, MPI_File fh,
             MPI_Offset offset, void *buf, MPI_Count count, MPI_Datatype datatype);
typedef int (QMPI_File_read_at_all_end_t) (QMPI_Context context, int tool_id, MPI_File fh,
             void *buf, MPI_Status *status);
typedef int (QMPI_File_read_ordered_t) (QMPI_Context context, int tool_id, MPI_File fh, void *buf,
             int count, MPI_Datatype datatype, MPI_Status *status);
typedef int (QMPI_File_read_ordered_c_t) (QMPI_Context context, int tool_id, MPI_File fh, void *buf,
             MPI_Count count, MPI_Datatype datatype, MPI_Status *status);
typedef int (QMPI_File_read_ordered_begin_t) (QMPI_Context context, int tool_id, MPI_File fh,
             void *buf, int count, MPI_Datatype datatype);
typedef int (QMPI_File_read_ordered_begin_c_t) (QMPI_Context context, int tool_id, MPI_File fh,
             void *buf, MPI_Count count, MPI_Datatype datatype);
typedef int (QMPI_File_read_ordered_end_t) (QMPI_Context context, int tool_id, MPI_File fh,
             void *buf, MPI_Status *status);
typedef int (QMPI_File_read_shared_t) (QMPI_Context context, int tool_id, MPI_File fh, void *buf,
             int count, MPI_Datatype datatype, MPI_Status *status);
typedef int (QMPI_File_read_shared_c_t) (QMPI_Context context, int tool_id, MPI_File fh, void *buf,
             MPI_Count count, MPI_Datatype datatype, MPI_Status *status);
typedef int (QMPI_File_seek_t) (QMPI_Context context, int tool_id, MPI_File fh, MPI_Offset offset,
             int whence);
typedef int (QMPI_File_seek_shared_t) (QMPI_Context context, int tool_id, MPI_File fh,
             MPI_Offset offset, int whence);
typedef int (QMPI_File_set_atomicity_t) (QMPI_Context context, int tool_id, MPI_File fh, int flag);
typedef int (QMPI_File_set_info_t) (QMPI_Context context, int tool_id, MPI_File fh, MPI_Info info);
typedef int (QMPI_File_set_size_t) (QMPI_Context context, int tool_id, MPI_File fh,
             MPI_Offset size);
typedef int (QMPI_File_set_view_t) (QMPI_Context context, int tool_id, MPI_File fh, MPI_Offset disp,
             MPI_Datatype etype, MPI_Datatype filetype, const char *datarep, MPI_Info info);
typedef int (QMPI_File_sync_t) (QMPI_Context context, int tool_id, MPI_File fh);
typedef int (QMPI_File_write_t) (QMPI_Context context, int tool_id, MPI_File fh, const void *buf,
             int count, MPI_Datatype datatype, MPI_Status *status);
typedef int (QMPI_File_write_c_t) (QMPI_Context context, int tool_id, MPI_File fh, const void *buf,
             MPI_Count count, MPI_Datatype datatype, MPI_Status *status);
typedef int (QMPI_File_write_all_t) (QMPI_Context context, int tool_id, MPI_File fh,
             const void *buf, int count, MPI_Datatype datatype, MPI_Status *status);
typedef int (QMPI_File_write_all_c_t) (QMPI_Context context, int tool_id, MPI_File fh,
             const void *buf, MPI_Count count, MPI_Datatype datatype, MPI_Status *status);
typedef int (QMPI_File_write_all_begin_t) (QMPI_Context context, int tool_id, MPI_File fh,
             const void *buf, int count, MPI_Datatype datatype);
typedef int (QMPI_File_write_all_begin_c_t) (QMPI_Context context, int tool_id, MPI_File fh,
             const void *buf, MPI_Count count, MPI_Datatype datatype);
typedef int (QMPI_File_write_all_end_t) (QMPI_Context context, int tool_id, MPI_File fh,
             const void *buf, MPI_Status *status);
typedef int (QMPI_File_write_at_t) (QMPI_Context context, int tool_id, MPI_File fh,
             MPI_Offset offset, const void *buf, int count, MPI_Datatype datatype,
             MPI_Status *status);
typedef int (QMPI_File_write_at_c_t) (QMPI_Context context, int tool_id, MPI_File fh,
             MPI_Offset offset, const void *buf, MPI_Count count, MPI_Datatype datatype,
             MPI_Status *status);
typedef int (QMPI_File_write_at_all_t) (QMPI_Context context, int tool_id, MPI_File fh,
             MPI_Offset offset, const void *buf, int count, MPI_Datatype datatype,
             MPI_Status *status);
typedef int (QMPI_File_write_at_all_c_t) (QMPI_Context context, int tool_id, MPI_File fh,
             MPI_Offset offset, const void *buf, MPI_Count count, MPI_Datatype datatype,
             MPI_Status *status);
typedef int (QMPI_File_write_at_all_begin_t) (QMPI_Context context, int tool_id, MPI_File fh,
             MPI_Offset offset, const void *buf, int count, MPI_Datatype datatype);
typedef int (QMPI_File_write_at_all_begin_c_t) (QMPI_Context context, int tool_id, MPI_File fh,
             MPI_Offset offset, const void *buf, MPI_Count count, MPI_Datatype datatype);
typedef int (QMPI_File_write_at_all_end_t) (QMPI_Context context, int tool_id, MPI_File fh,
             const void *buf, MPI_Status *status);
typedef int (QMPI_File_write_ordered_t) (QMPI_Context context, int tool_id, MPI_File fh,
             const void *buf, int count, MPI_Datatype datatype, MPI_Status *status);
typedef int (QMPI_File_write_ordered_c_t) (QMPI_Context context, int tool_id, MPI_File fh,
             const void *buf, MPI_Count count, MPI_Datatype datatype, MPI_Status *status);
typedef int (QMPI_File_write_ordered_begin_t) (QMPI_Context context, int tool_id, MPI_File fh,
             const void *buf, int count, MPI_Datatype datatype);
typedef int (QMPI_File_write_ordered_begin_c_t) (QMPI_Context context, int tool_id, MPI_File fh,
             const void *buf, MPI_Count count, MPI_Datatype datatype);
typedef int (QMPI_File_write_ordered_end_t) (QMPI_Context context, int tool_id, MPI_File fh,
             const void *buf, MPI_Status *status);
typedef int (QMPI_File_write_shared_t) (QMPI_Context context, int tool_id, MPI_File fh,
             const void *buf, int count, MPI_Datatype datatype, MPI_Status *status);
typedef int (QMPI_File_write_shared_c_t) (QMPI_Context context, int tool_id, MPI_File fh,
             const void *buf, MPI_Count count, MPI_Datatype datatype, MPI_Status *status);
typedef int (QMPI_Register_datarep_t) (QMPI_Context context, int tool_id, const char *datarep,
             MPI_Datarep_conversion_function *read_conversion_fn,
             MPI_Datarep_conversion_function *write_conversion_fn,
             MPI_Datarep_extent_function *dtype_file_extent_fn, void *extra_state);
typedef int (QMPI_Register_datarep_c_t) (QMPI_Context context, int tool_id, const char *datarep,
             MPI_Datarep_conversion_function_c *read_conversion_fn,
             MPI_Datarep_conversion_function_c *write_conversion_fn,
             MPI_Datarep_extent_function *dtype_file_extent_fn, void *extra_state);
typedef int (QMPI_File_toint_t) (QMPI_Context context, int tool_id, MPI_File file);
typedef MPI_File (QMPI_File_fromint_t) (QMPI_Context context, int tool_id, int file);
typedef __darwin_clock_t clock_t;
typedef __darwin_time_t time_t;
struct timespec
{
 __darwin_time_t tv_sec;
 long tv_nsec;
};
struct tm {
 int tm_sec;
 int tm_min;
 int tm_hour;
 int tm_mday;
 int tm_mon;
 int tm_year;
 int tm_wday;
 int tm_yday;
 int tm_isdst;
 long tm_gmtoff;
 char * tm_zone;
};
extern char * tzname[];
extern int getdate_err;
extern long timezone __asm("_" "timezone" );
extern int daylight;
char * asctime(const struct tm *);
clock_t clock(void) __asm("_" "clock" );
char * ctime(const time_t *);
double difftime(time_t, time_t);
struct tm *getdate(const char *);
struct tm *gmtime(const time_t *);
struct tm *localtime(const time_t *);
time_t mktime(struct tm *) __asm("_" "mktime" );
size_t strftime(char * restrict, size_t __maxsize, const char * restrict, const struct tm * restrict) __asm("_" "strftime" );
char * strptime(const char * restrict, const char * restrict, struct tm * restrict) __asm("_" "strptime" );
time_t time(time_t *);
void tzset(void);
char * asctime_r(const struct tm * restrict, char * restrict );
char * ctime_r(const time_t *, char *);
struct tm *gmtime_r(const time_t * restrict, struct tm * restrict);
struct tm *localtime_r(const time_t * restrict, struct tm * restrict);
time_t posix2time(time_t);
void tzsetwall(void);
time_t time2posix(time_t);
time_t timelocal(struct tm * const);
time_t timegm(struct tm * const);
int nanosleep(const struct timespec *__rqtp, struct timespec *__rmtp) __asm("_" "nanosleep" );
typedef enum {
_CLOCK_REALTIME __attribute__((availability(macosx,introduced=10.12))) __attribute__((availability(ios,introduced=10.0))) __attribute__((availability(tvos,introduced=10.0))) __attribute__((availability(watchos,introduced=3.0))) = 0,
_CLOCK_MONOTONIC __attribute__((availability(macosx,introduced=10.12))) __attribute__((availability(ios,introduced=10.0))) __attribute__((availability(tvos,introduced=10.0))) __attribute__((availability(watchos,introduced=3.0))) = 6,
_CLOCK_MONOTONIC_RAW __attribute__((availability(macosx,introduced=10.12))) __attribute__((availability(ios,introduced=10.0))) __attribute__((availability(tvos,introduced=10.0))) __attribute__((availability(watchos,introduced=3.0))) = 4,
_CLOCK_MONOTONIC_RAW_APPROX __attribute__((availability(macosx,introduced=10.12))) __attribute__((availability(ios,introduced=10.0))) __attribute__((availability(tvos,introduced=10.0))) __attribute__((availability(watchos,introduced=3.0))) = 5,
_CLOCK_UPTIME_RAW __attribute__((availability(macosx,introduced=10.12))) __attribute__((availability(ios,introduced=10.0))) __attribute__((availability(tvos,introduced=10.0))) __attribute__((availability(watchos,introduced=3.0))) = 8,
_CLOCK_UPTIME_RAW_APPROX __attribute__((availability(macosx,introduced=10.12))) __attribute__((availability(ios,introduced=10.0))) __attribute__((availability(tvos,introduced=10.0))) __attribute__((availability(watchos,introduced=3.0))) = 9,
_CLOCK_PROCESS_CPUTIME_ID __attribute__((availability(macosx,introduced=10.12))) __attribute__((availability(ios,introduced=10.0))) __attribute__((availability(tvos,introduced=10.0))) __attribute__((availability(watchos,introduced=3.0))) = 12,
_CLOCK_THREAD_CPUTIME_ID __attribute__((availability(macosx,introduced=10.12))) __attribute__((availability(ios,introduced=10.0))) __attribute__((availability(tvos,introduced=10.0))) __attribute__((availability(watchos,introduced=3.0))) = 16
} clockid_t;
__attribute__((availability(macosx,introduced=10.12))) __attribute__((availability(ios,introduced=10.0))) __attribute__((availability(tvos,introduced=10.0))) __attribute__((availability(watchos,introduced=3.0)))
int clock_getres(clockid_t __clock_id, struct timespec *__res);
__attribute__((availability(macosx,introduced=10.12))) __attribute__((availability(ios,introduced=10.0))) __attribute__((availability(tvos,introduced=10.0))) __attribute__((availability(watchos,introduced=3.0)))
int clock_gettime(clockid_t __clock_id, struct timespec *__tp);
__attribute__((availability(macosx,introduced=10.12))) __attribute__((availability(ios,introduced=10.0))) __attribute__((availability(tvos,introduced=10.0))) __attribute__((availability(watchos,introduced=3.0)))
__uint64_t clock_gettime_nsec_np(clockid_t __clock_id);
__attribute__((availability(macosx,introduced=10.12))) __attribute__((availability(ios,unavailable)))
__attribute__((availability(tvos,unavailable))) __attribute__((availability(watchos,unavailable)))
int clock_settime(clockid_t __clock_id, const struct timespec *__tp);
__attribute__((availability(macos,introduced=10.15))) __attribute__((availability(ios,introduced=13.0))) __attribute__((availability(tvos,introduced=13.0))) __attribute__((availability(watchos,introduced=6.0)))
int timespec_get(struct timespec *ts, int base);
struct _43d5iqv
{ long _gnzei3z;
    double _zle6kwh;
    float _wd1xmxb;
    float _rg2gzr7;
    float _x5qjicj;
    float _5bwikot;
    float _plozjcn;
    float _g87oi0e;
    float _7yt2m6t;
    float _fqewaou;
    float _8drc11e;
    float _z0yncxl;
    unsigned int _jlteiyz;
};
struct _gzui1b0
{ float _0dtirc9;
    float _txcb48y;
    float _hez74la;
    float _hj3muur;
    float _iq2fa39;
    float _p8ijl2z;
    unsigned int _otn9uze;
    unsigned int _y25m8kb;
    unsigned int _1nlc8jp;
    unsigned int _d7x48wz;
    unsigned int _mygty3a;
    unsigned short int _w2maiet;
    unsigned short int _k00ilu6;
    unsigned short int _25qbrmq;
};
struct _rcbpljd
{ float _d31k9qm;
    float _m7pyxdw;
    float _31yrujo;
    float _8kzruin;
    float _glqj5sx;
    float _9iy4x41;
    float _ckeh5sw;
    float _ktakpux;
    float _qphr62g;
    float _dswm31y;
    unsigned int _cs4kek4;
    unsigned int _7yvdqbs;
    unsigned short int _cw2lbm2;
    unsigned short int _q8hke22;
    unsigned short int _p7bd4w1;
    unsigned short int _t9nm3v2;
    unsigned short int _xo100vh;
    unsigned short int _fhv1xg1;
};
extern void SubtractVect3(float *restrict _9yo18pa, const float *restrict _0ljch2v, const float *restrict _6zs62yw);
extern float VectorLength(const float *restrict _2ljk4q5);
extern void NrmlzeVect(float *restrict _2ljk4q5);
extern void GetStkAndDipVect(const float *restrict _ljefy3z, float *restrict _wrb7f6f, float *restrict _jssuuqo);
extern void RotStressG2L(float *restrict _yj9boq0, const float *restrict _v3z6bz8, const float *restrict _e1g983e);
extern void ScaleStress6(float *restrict _kibu4fm, float _088o98b);
extern float ran0(long *_1poxrte);
extern long int GetByteSize(int _eoz0kz4, int _tajw253, int _0fqldb2);
extern void ReAssignElemVals(struct _43d5iqv *restrict _j3trp0j, const struct _gzui1b0 *restrict _s3kkv9r, unsigned int i, unsigned int *restrict _mgvu4is, float *restrict _vwvm4w3, float *restrict _wyx9trf, unsigned int _gvz9aaw);
extern void StrainHS_Nikkhoo(float _e8f8409[6], float _wyzx4w6[6], float X, float Y, float Z, float _0lzjjeg[3], float _9skdgig[3], float _5r5uqpg[3], float _q0oxfq8, float _tb94874, float _3jgzlma, const float _9wszv83, const float _ufwbucj);
extern void StrainFS_Nikkhoo(float _e8f8409[6], float _wyzx4w6[6], float X, float Y, float Z, float _0lzjjeg[3], float _9skdgig[3], float _5r5uqpg[3], float _q0oxfq8, float _tb94874, float _3jgzlma, const float _9wszv83, const float _ufwbucj);
void _18hafey(float _kb4ehou, unsigned int _o31yhxs, const float *restrict _v78hu7f, unsigned int _mda9oek, const unsigned int *restrict _qw4c6lb, const float *restrict _nbhnal3, unsigned int ***restrict _l9ue421, unsigned int ***restrict _t4q9bc2, unsigned int *restrict _vfua9bq);
void _h75lxo6(unsigned int _mda9oek, const unsigned int *restrict _qw4c6lb, const float *restrict _nbhnal3, const float *restrict _r1dwwo8, unsigned int **restrict _2h2p39y, unsigned int _vfua9bq, float ***restrict _vlnknd6);
void _j4t0sek(const char *restrict _slwf90k, unsigned int _2vpgstq, const float *restrict _q01tee1, unsigned int **restrict _qc6dr08, unsigned int _qsq0p82, float **restrict _zh13cd5, unsigned int _04k0lth, const float *restrict _o6l1f90, unsigned int **restrict _1zlbori, unsigned int _6parsg7, float **restrict _s3jonou);
void MakeNew_Khmatrix(int _b6gs8g8, int _tk63n2q, const char *restrict _slwf90k, struct _43d5iqv *_j3trp0j, struct _gzui1b0 *restrict _s3kkv9r, const struct _rcbpljd *restrict _84u89yj, const int *restrict _6zx1tv1, const int *restrict _emruyp0, const unsigned int *restrict _5izolnn, const unsigned int *restrict _laym0wb, const float *restrict _lqbmae3, const float *restrict _h55rbc7, const unsigned int *restrict _egok01t, const unsigned int *restrict _qydbpgv, const float *restrict _3n3bzt9, const float *restrict _s9bp4m1, const float *restrict _n8la100, const float *restrict _t8qzpky, unsigned int *restrict _mgvu4is, float *restrict _vwvm4w3, float *restrict _wyx9trf, float *restrict _dntcagr, unsigned short int *restrict _hthsjr0, unsigned short int *restrict _pmlvak9, unsigned short int *restrict _vk0ucc1, unsigned short int *restrict _vshkoec, unsigned short int ***restrict _k2wbxan, unsigned short int ***restrict _n14eb87, unsigned short int ***restrict _va1nkbt, unsigned short int ***restrict _eurxhw0, float ***restrict _0isb5iq, float ***restrict _zzf8s2z, float ***restrict _3jczucq, float ***restrict _sh3tow4, float ***restrict _jbc2sws, float ***restrict _e9ftaxp, float ***restrict _6fokzix, float ***restrict _z39pecn, float ***restrict _t997e4u, float ***restrict _d1vn8lf, float ***restrict _4xkss0k)
{
    struct timespec _zn4regb, _ofn1r4k;
    double _weehwjq;
    unsigned int _pwh0mx1 = 0u;
    unsigned int _kicwmpe = 0u;
    unsigned int **_f7fr91i = ((void *)0), **_wshmlfz = ((void *)0), **_jzr25md = ((void *)0), **_jlpefja = ((void *)0);
    float **_pqvl00e = ((void *)0), **_xf7jxjn = ((void *)0);
    clock_gettime(_CLOCK_REALTIME, &_zn4regb);
    _18hafey( _84u89yj->_d31k9qm, _84u89yj->_cw2lbm2, _lqbmae3, _s3kkv9r->_otn9uze, _egok01t, _3n3bzt9, &_f7fr91i, &_jzr25md, &_pwh0mx1);
    if (_s3kkv9r->_y25m8kb > 0u)
    { _18hafey(_84u89yj->_m7pyxdw, _84u89yj->_q8hke22, _h55rbc7, _s3kkv9r->_y25m8kb, _qydbpgv, _s9bp4m1, &_wshmlfz, &_jlpefja, &_kicwmpe);
    }
    MPI_Barrier( ((MPI_Comm)0x44000000) );
    clock_gettime(_CLOCK_REALTIME, &_ofn1r4k);
    _weehwjq = (_ofn1r4k.tv_sec - _zn4regb.tv_sec) + (_ofn1r4k.tv_nsec - _zn4regb.tv_nsec)/1.0E+9;
    if (_tk63n2q == 0) { fprintf(__stdoutp,"BranchCountF: %7u\nTime taken (in seconds): %lf\n", _pwh0mx1, _weehwjq); }
    _h75lxo6( _s3kkv9r->_otn9uze, _egok01t, _3n3bzt9, _n8la100, _f7fr91i, _pwh0mx1, &_pqvl00e);
    if (_s3kkv9r->_y25m8kb > 0u)
    { _h75lxo6(_s3kkv9r->_y25m8kb, _qydbpgv, _s9bp4m1, _t8qzpky, _wshmlfz, _kicwmpe, &_xf7jxjn);
    }
    MPI_Barrier( ((MPI_Comm)0x44000000) );
    clock_gettime(_CLOCK_REALTIME, &_zn4regb);
    { unsigned int i, j;
        unsigned int _xdrkz40, _fu7bfxq, _be3yqy5, _iua08ek, _vmefozt, _o7vbgr8, _jg2y8rg, _ybxkpm2, _1b07g2r;
        float _66cnolj, _q52wspo, _fv177ud, _8e7flqh, _1vsmado;
        float _uug95rl[3], _dan085r[3], _e3nwxuj[3], _3ac6pe6[3], _3lbqgu3[3], _obyghgk[3], _2ljk4q5[3];
        float _mt8vt43[6], _pkxvhag[6], _dtg5nqj[6], _pt5r5rc[6], _kc4lorx[6], _qijc1yd[6], _dhgzz0j[6], _v71448n[9], _2292mtr[9], _gwaj1g1[9];
        long int _050ui22;
        unsigned int *_rgb93wk = ((void *)0), *_kb700uc = ((void *)0), *_ye0qzsa = ((void *)0), *_0s1f537 = ((void *)0);
        float *_nf74t60 = ((void *)0), *_ge2j51b = ((void *)0), *_mrunjdu = ((void *)0), *_ycz04cd = ((void *)0), *_ozctwjl = ((void *)0);
        float *_0wfoy4d = ((void *)0), *_zareiya = ((void *)0), *_0dl5e11 = ((void *)0), *_ny181ac = ((void *)0), *_jdtvv2y = ((void *)0), *_b68rerv = ((void *)0);
        _rgb93wk = malloc(_s3kkv9r->_otn9uze *sizeof(unsigned int));
        _kb700uc = malloc(_pwh0mx1 *sizeof(unsigned int));
        _050ui22 = GetByteSize(64, 2*_s3kkv9r->_otn9uze, sizeof(float));
        _nf74t60 = aligned_alloc(64, _050ui22);
        _ge2j51b = aligned_alloc(64, _050ui22);
        _mrunjdu = aligned_alloc(64, _050ui22);
        if (_s3kkv9r->_y25m8kb > 0u)
        { _ye0qzsa = malloc(_s3kkv9r->_y25m8kb *sizeof(unsigned int));
            _0s1f537 = malloc(_kicwmpe *sizeof(unsigned int));
            _050ui22 = GetByteSize(64, 3*_s3kkv9r->_y25m8kb, sizeof(float));
            _ycz04cd = aligned_alloc(64, _050ui22);
            _ozctwjl = aligned_alloc(64, _050ui22);
        }
        _050ui22 = GetByteSize(64, _6zx1tv1[_tk63n2q], sizeof(unsigned short int *));
        *_k2wbxan = aligned_alloc(64, _050ui22);
        _050ui22 = GetByteSize(64, _s3kkv9r->_otn9uze, sizeof(unsigned short int));
        for (i = 0u; i < _6zx1tv1[_tk63n2q]; i++)
        { (*_k2wbxan)[i] = aligned_alloc(64, _050ui22); __builtin___memset_chk ((*_k2wbxan)[i], 0, _050ui22, __builtin_object_size ((*_k2wbxan)[i], 0));
        }
        if (_s3kkv9r->_y25m8kb > 0u)
        { _050ui22 = GetByteSize(64, _6zx1tv1[_tk63n2q], sizeof(unsigned short int *));
            *_n14eb87 = aligned_alloc(64, _050ui22);
            _050ui22 = GetByteSize(64, _s3kkv9r->_y25m8kb, sizeof(unsigned short int));
            for (i = 0u; i < _6zx1tv1[_tk63n2q]; i++)
            { (*_n14eb87)[i] = aligned_alloc(64, _050ui22); __builtin___memset_chk ((*_n14eb87)[i], 0, _050ui22, __builtin_object_size ((*_n14eb87)[i], 0));
            }
        }
        _050ui22 = GetByteSize(64, _6zx1tv1[_tk63n2q], sizeof(float *));
        *_0isb5iq = aligned_alloc(64, _050ui22);
        *_zzf8s2z = aligned_alloc(64, _050ui22);
        *_3jczucq = aligned_alloc(64, _050ui22);
        if (_s3kkv9r->_y25m8kb > 0u)
        { *_sh3tow4 = aligned_alloc(64, _050ui22);
            *_jbc2sws = aligned_alloc(64, _050ui22);
        }
        for (i = 0; i < _6zx1tv1[_tk63n2q]; i++)
        { _vmefozt = _5izolnn[i];
            _050ui22 = GetByteSize(64, 2*_s3kkv9r->_otn9uze, sizeof(float));
            __builtin___memset_chk (_nf74t60, 0, _050ui22, __builtin_object_size (_nf74t60, 0));
            __builtin___memset_chk (_ge2j51b, 0, _050ui22, __builtin_object_size (_ge2j51b, 0));
            __builtin___memset_chk (_mrunjdu, 0, _050ui22, __builtin_object_size (_mrunjdu, 0));
            __builtin___memcpy_chk (_uug95rl, &_3n3bzt9[_vmefozt*16 +0], 3*sizeof(float), __builtin_object_size (_uug95rl, 0));
            __builtin___memcpy_chk (_v71448n, &_3n3bzt9[_vmefozt*16 +6], 9*sizeof(float), __builtin_object_size (_v71448n, 0));
            __builtin___memset_chk (_rgb93wk, 0, _s3kkv9r->_otn9uze*sizeof(unsigned int), __builtin_object_size (_rgb93wk, 0));
            __builtin___memset_chk (_kb700uc, 0, _pwh0mx1*sizeof(unsigned int), __builtin_object_size (_kb700uc, 0));
            _xdrkz40 = 0u;
            while (_xdrkz40 < _s3kkv9r->_otn9uze)
            {
                for (j = 0u; j < _s3kkv9r->_otn9uze; j++) { if (_rgb93wk[j] == 0u) { _o7vbgr8 = j; break; } }
                __builtin___memcpy_chk (_dan085r, &_3n3bzt9[_o7vbgr8*16 +0], 3*sizeof(float), __builtin_object_size (_dan085r, 0));
                SubtractVect3(_2ljk4q5, _uug95rl, _dan085r); _1vsmado = VectorLength(_2ljk4q5);
                _8e7flqh = (_1vsmado >= 1.75f*_84u89yj->_d31k9qm)*1.0f + (_1vsmado < 1.75f*_84u89yj->_d31k9qm)*0.0f;
                _8e7flqh = (_egok01t[_vmefozt*8 +4] == _egok01t[_o7vbgr8*8 +4])*1.0f + (_egok01t[_vmefozt*8 +4] != _egok01t[_o7vbgr8*8 +4])*_8e7flqh;
                _jg2y8rg = _f7fr91i[_o7vbgr8][0];
                _fu7bfxq = 0;
                _ybxkpm2 = 1u;
                while (_ybxkpm2 == 1u)
                {
                    _1b07g2r = 0u;
                    while (_fu7bfxq < _jg2y8rg-1)
                    { _be3yqy5 = _f7fr91i[_o7vbgr8][16 -_jg2y8rg +_fu7bfxq];
                        __builtin___memcpy_chk (_dan085r, &_pqvl00e[_be3yqy5][0], 3*sizeof(float), __builtin_object_size (_dan085r, 0));
                        SubtractVect3(_2ljk4q5, _uug95rl, _dan085r); _1vsmado = VectorLength(_2ljk4q5);
                        _1b07g2r = (_kb700uc[_be3yqy5] == 0u)*1u + (_kb700uc[_be3yqy5] != 0u)*0u;
                        _1b07g2r = (_jzr25md[_be3yqy5][0] > 1u)*_1b07g2r + (_jzr25md[_be3yqy5][0] <= 1u)*0u;
                        _1b07g2r = (_f7fr91i[_vmefozt][16 -_jg2y8rg +_fu7bfxq] != _be3yqy5)*_1b07g2r + (_f7fr91i[_vmefozt][16 -_jg2y8rg +_fu7bfxq] == _be3yqy5)*0u;
                        _1b07g2r = (_1vsmado > (_pqvl00e[_be3yqy5][5] +_84u89yj->_ktakpux*_84u89yj->_d31k9qm))*_1b07g2r + (_1vsmado <= (_pqvl00e[_be3yqy5][5] +_84u89yj->_ktakpux*_84u89yj->_d31k9qm))*0u;
                        if (_1b07g2r == 1u) { break; }
                        _fu7bfxq += 1u;
                    }
                    if (_1b07g2r == 0u)
                    { __builtin___memcpy_chk (_e3nwxuj, &_n8la100[_egok01t[_o7vbgr8*8 +0]*4 +0], 3*sizeof(float), __builtin_object_size (_e3nwxuj, 0));
                        __builtin___memcpy_chk (_3ac6pe6, &_n8la100[_egok01t[_o7vbgr8*8 +1]*4 +0], 3*sizeof(float), __builtin_object_size (_3ac6pe6, 0));
                        __builtin___memcpy_chk (_3lbqgu3, &_n8la100[_egok01t[_o7vbgr8*8 +2]*4 +0], 3*sizeof(float), __builtin_object_size (_3lbqgu3, 0));
                        if (_84u89yj->_fhv1xg1 == 1u)
                        { StrainHS_Nikkhoo(_mt8vt43, _dhgzz0j, _uug95rl[0], _uug95rl[1], _uug95rl[2], _e3nwxuj, _3ac6pe6, _3lbqgu3, _84u89yj->_31yrujo, 0.0f, 0.0f, _j3trp0j->_x5qjicj, _j3trp0j->_plozjcn);
                            StrainHS_Nikkhoo(_pkxvhag, _dhgzz0j, _uug95rl[0], _uug95rl[1], _uug95rl[2], _e3nwxuj, _3ac6pe6, _3lbqgu3, 0.0f, _84u89yj->_31yrujo, 0.0f, _j3trp0j->_x5qjicj, _j3trp0j->_plozjcn);
                            StrainHS_Nikkhoo(_dtg5nqj, _dhgzz0j, _uug95rl[0], _uug95rl[1], _uug95rl[2], _e3nwxuj, _3ac6pe6, _3lbqgu3, 0.0f, 0.0f, _84u89yj->_31yrujo, _j3trp0j->_x5qjicj, _j3trp0j->_plozjcn);
                        }
                        else
                        { StrainFS_Nikkhoo(_mt8vt43, _dhgzz0j, _uug95rl[0], _uug95rl[1], _uug95rl[2], _e3nwxuj, _3ac6pe6, _3lbqgu3, _84u89yj->_31yrujo, 0.0f, 0.0f, _j3trp0j->_x5qjicj, _j3trp0j->_plozjcn);
                            StrainFS_Nikkhoo(_pkxvhag, _dhgzz0j, _uug95rl[0], _uug95rl[1], _uug95rl[2], _e3nwxuj, _3ac6pe6, _3lbqgu3, 0.0f, _84u89yj->_31yrujo, 0.0f, _j3trp0j->_x5qjicj, _j3trp0j->_plozjcn);
                            StrainFS_Nikkhoo(_dtg5nqj, _dhgzz0j, _uug95rl[0], _uug95rl[1], _uug95rl[2], _e3nwxuj, _3ac6pe6, _3lbqgu3, 0.0f, 0.0f, _84u89yj->_31yrujo, _j3trp0j->_x5qjicj, _j3trp0j->_plozjcn);
                        }
                        ScaleStress6(_mt8vt43, (_8e7flqh/(_84u89yj->_31yrujo*_3n3bzt9[_o7vbgr8*16 +15]))); RotStressG2L(&_2292mtr[0], _v71448n, _mt8vt43);
                        ScaleStress6(_pkxvhag, (_8e7flqh/(_84u89yj->_31yrujo*_3n3bzt9[_o7vbgr8*16 +15]))); RotStressG2L(&_2292mtr[3], _v71448n, _pkxvhag);
                        ScaleStress6(_dtg5nqj, (_8e7flqh/(_84u89yj->_31yrujo*_3n3bzt9[_o7vbgr8*16 +15]))); RotStressG2L(&_2292mtr[6], _v71448n, _dtg5nqj);
                        if (_vmefozt == _o7vbgr8)
                        { _vwvm4w3[i*16 +5] = _2292mtr[0]*_3n3bzt9[_o7vbgr8*16 +15];
                            _vwvm4w3[i*16 +6] = _2292mtr[4]*_3n3bzt9[_o7vbgr8*16 +15];
                            _vwvm4w3[i*16 +7] = _2292mtr[8]*_3n3bzt9[_o7vbgr8*16 +15];
                            _wyx9trf[i*16 +1] = ( ((fabsf(_vwvm4w3[i*16 +5])) > (fabsf(_vwvm4w3[i*16 +6]))) *(fabsf(_vwvm4w3[i*16 +5])) + ((fabsf(_vwvm4w3[i*16 +5])) <= (fabsf(_vwvm4w3[i*16 +6]))) *(fabsf(_vwvm4w3[i*16 +6])) );
                            _wyx9trf[i*16 +1] = ( ((_wyx9trf[i*16 +1]) > (fabsf(_vwvm4w3[i*16 +7]))) *(_wyx9trf[i*16 +1]) + ((_wyx9trf[i*16 +1]) <= (fabsf(_vwvm4w3[i*16 +7]))) *(fabsf(_vwvm4w3[i*16 +7])) );
                            ReAssignElemVals(_j3trp0j, _s3kkv9r, i, _mgvu4is, _vwvm4w3, _wyx9trf, 1);
                        }
                        _nf74t60[_hthsjr0[i]*2 +0] = _2292mtr[0]; _nf74t60[_hthsjr0[i]*2 +1] = _2292mtr[3];
                        _ge2j51b[_hthsjr0[i]*2 +0] = _2292mtr[1]; _ge2j51b[_hthsjr0[i]*2 +1] = _2292mtr[4];
                        _mrunjdu[_hthsjr0[i]*2 +0] = _2292mtr[2]; _mrunjdu[_hthsjr0[i]*2 +1] = _2292mtr[5];
                        (*_k2wbxan)[i][_o7vbgr8] = _hthsjr0[i];
                        _hthsjr0[i] += 1u;
                        _rgb93wk[_o7vbgr8] = 1u;
                        _xdrkz40 = 0u;
                        for (j = 0u; j < _s3kkv9r->_otn9uze; j++) { _xdrkz40 += _rgb93wk[j]; }
                        for (j = 1u; j < _jg2y8rg; j++) { _kb700uc[ _f7fr91i[_o7vbgr8][16 -_jg2y8rg +j] ] = 1u; }
                        _ybxkpm2 = 0u;
                    }
                    else
                    { _be3yqy5 = _f7fr91i[_o7vbgr8][16 -_jg2y8rg +_fu7bfxq];
                        __builtin___memcpy_chk (_e3nwxuj, &_pqvl00e[_be3yqy5][9], 3*sizeof(float), __builtin_object_size (_e3nwxuj, 0));
                        __builtin___memcpy_chk (_3ac6pe6, &_pqvl00e[_be3yqy5][12], 3*sizeof(float), __builtin_object_size (_3ac6pe6, 0));
                        __builtin___memcpy_chk (_3lbqgu3, &_pqvl00e[_be3yqy5][15], 3*sizeof(float), __builtin_object_size (_3lbqgu3, 0));
                        __builtin___memcpy_chk (_obyghgk, &_pqvl00e[_be3yqy5][18], 3*sizeof(float), __builtin_object_size (_obyghgk, 0));
                        if (_84u89yj->_fhv1xg1 == 1u)
                        { StrainHS_Nikkhoo(_pt5r5rc, _dhgzz0j, _uug95rl[0], _uug95rl[1], _uug95rl[2], _e3nwxuj, _3ac6pe6, _obyghgk, _84u89yj->_31yrujo, 0.0f, 0.0f, _j3trp0j->_x5qjicj, _j3trp0j->_plozjcn);
                            StrainHS_Nikkhoo(_kc4lorx, _dhgzz0j, _uug95rl[0], _uug95rl[1], _uug95rl[2], _e3nwxuj, _3ac6pe6, _obyghgk, 0.0f, _84u89yj->_31yrujo, 0.0f, _j3trp0j->_x5qjicj, _j3trp0j->_plozjcn);
                            StrainHS_Nikkhoo(_mt8vt43, _dhgzz0j, _uug95rl[0], _uug95rl[1], _uug95rl[2], _3ac6pe6, _3lbqgu3, _obyghgk, _84u89yj->_31yrujo, 0.0f, 0.0f, _j3trp0j->_x5qjicj, _j3trp0j->_plozjcn);
                            StrainHS_Nikkhoo(_pkxvhag, _dhgzz0j, _uug95rl[0], _uug95rl[1], _uug95rl[2], _3ac6pe6, _3lbqgu3, _obyghgk, 0.0f, _84u89yj->_31yrujo, 0.0f, _j3trp0j->_x5qjicj, _j3trp0j->_plozjcn);
                            _mt8vt43[0] += _pt5r5rc[0]; _mt8vt43[1] += _pt5r5rc[1]; _mt8vt43[2] += _pt5r5rc[2];
                            _mt8vt43[3] += _pt5r5rc[3]; _mt8vt43[4] += _pt5r5rc[4]; _mt8vt43[5] += _pt5r5rc[5];
                            _pkxvhag[0] += _kc4lorx[0]; _pkxvhag[1] += _kc4lorx[1]; _pkxvhag[2] += _kc4lorx[2];
                            _pkxvhag[3] += _kc4lorx[3]; _pkxvhag[4] += _kc4lorx[4]; _pkxvhag[5] += _kc4lorx[5];
                        }
                        else
                        { StrainFS_Nikkhoo(_pt5r5rc, _dhgzz0j, _uug95rl[0], _uug95rl[1], _uug95rl[2], _e3nwxuj, _3ac6pe6, _obyghgk, _84u89yj->_31yrujo, 0.0f, 0.0f, _j3trp0j->_x5qjicj, _j3trp0j->_plozjcn);
                            StrainFS_Nikkhoo(_kc4lorx, _dhgzz0j, _uug95rl[0], _uug95rl[1], _uug95rl[2], _e3nwxuj, _3ac6pe6, _obyghgk, 0.0f, _84u89yj->_31yrujo, 0.0f, _j3trp0j->_x5qjicj, _j3trp0j->_plozjcn);
                            StrainFS_Nikkhoo(_mt8vt43, _dhgzz0j, _uug95rl[0], _uug95rl[1], _uug95rl[2], _3ac6pe6, _3lbqgu3, _obyghgk, _84u89yj->_31yrujo, 0.0f, 0.0f, _j3trp0j->_x5qjicj, _j3trp0j->_plozjcn);
                            StrainFS_Nikkhoo(_pkxvhag, _dhgzz0j, _uug95rl[0], _uug95rl[1], _uug95rl[2], _3ac6pe6, _3lbqgu3, _obyghgk, 0.0f, _84u89yj->_31yrujo, 0.0f, _j3trp0j->_x5qjicj, _j3trp0j->_plozjcn);
                            _mt8vt43[0] += _pt5r5rc[0]; _mt8vt43[1] += _pt5r5rc[1]; _mt8vt43[2] += _pt5r5rc[2];
                            _mt8vt43[3] += _pt5r5rc[3]; _mt8vt43[4] += _pt5r5rc[4]; _mt8vt43[5] += _pt5r5rc[5];
                            _pkxvhag[0] += _kc4lorx[0]; _pkxvhag[1] += _kc4lorx[1]; _pkxvhag[2] += _kc4lorx[2];
                            _pkxvhag[3] += _kc4lorx[3]; _pkxvhag[4] += _kc4lorx[4]; _pkxvhag[5] += _kc4lorx[5];
                        }
                        ScaleStress6(_mt8vt43, (_8e7flqh/(_84u89yj->_31yrujo*_pqvl00e[_be3yqy5][3]))); RotStressG2L(&_2292mtr[0], _v71448n, _mt8vt43);
                        ScaleStress6(_pkxvhag, (_8e7flqh/(_84u89yj->_31yrujo*_pqvl00e[_be3yqy5][3]))); RotStressG2L(&_2292mtr[3], _v71448n, _pkxvhag);
                        _1b07g2r = 1u;
                        for (j = 0; j < _jzr25md[_be3yqy5][0]; j++)
                        { _iua08ek = _jzr25md[_be3yqy5][(j+1)];
                            __builtin___memcpy_chk (_e3nwxuj, &_pqvl00e[_iua08ek][9], 3*sizeof(float), __builtin_object_size (_e3nwxuj, 0));
                            __builtin___memcpy_chk (_3ac6pe6, &_pqvl00e[_iua08ek][12], 3*sizeof(float), __builtin_object_size (_3ac6pe6, 0));
                            __builtin___memcpy_chk (_3lbqgu3, &_pqvl00e[_iua08ek][15], 3*sizeof(float), __builtin_object_size (_3lbqgu3, 0));
                            __builtin___memcpy_chk (_obyghgk, &_pqvl00e[_iua08ek][18], 3*sizeof(float), __builtin_object_size (_obyghgk, 0));
                            if (_84u89yj->_fhv1xg1 == 1u)
                            { StrainHS_Nikkhoo(_pt5r5rc, _dhgzz0j, _uug95rl[0], _uug95rl[1], _uug95rl[2], _e3nwxuj, _3ac6pe6, _obyghgk, _84u89yj->_31yrujo, 0.0f, 0.0f, _j3trp0j->_x5qjicj, _j3trp0j->_plozjcn);
                                StrainHS_Nikkhoo(_kc4lorx, _dhgzz0j, _uug95rl[0], _uug95rl[1], _uug95rl[2], _e3nwxuj, _3ac6pe6, _obyghgk, 0.0f, _84u89yj->_31yrujo, 0.0f, _j3trp0j->_x5qjicj, _j3trp0j->_plozjcn);
                                StrainHS_Nikkhoo(_mt8vt43, _dhgzz0j, _uug95rl[0], _uug95rl[1], _uug95rl[2], _3ac6pe6, _3lbqgu3, _obyghgk, _84u89yj->_31yrujo, 0.0f, 0.0f, _j3trp0j->_x5qjicj, _j3trp0j->_plozjcn);
                                StrainHS_Nikkhoo(_pkxvhag, _dhgzz0j, _uug95rl[0], _uug95rl[1], _uug95rl[2], _3ac6pe6, _3lbqgu3, _obyghgk, 0.0f, _84u89yj->_31yrujo, 0.0f, _j3trp0j->_x5qjicj, _j3trp0j->_plozjcn);
                                _mt8vt43[0] += _pt5r5rc[0]; _mt8vt43[1] += _pt5r5rc[1]; _mt8vt43[2] += _pt5r5rc[2];
                                _mt8vt43[3] += _pt5r5rc[3]; _mt8vt43[4] += _pt5r5rc[4]; _mt8vt43[5] += _pt5r5rc[5];
                                _pkxvhag[0] += _kc4lorx[0]; _pkxvhag[1] += _kc4lorx[1]; _pkxvhag[2] += _kc4lorx[2];
                                _pkxvhag[3] += _kc4lorx[3]; _pkxvhag[4] += _kc4lorx[4]; _pkxvhag[5] += _kc4lorx[5];
                            }
                            else
                            { StrainFS_Nikkhoo(_pt5r5rc, _dhgzz0j, _uug95rl[0], _uug95rl[1], _uug95rl[2], _e3nwxuj, _3ac6pe6, _obyghgk, _84u89yj->_31yrujo, 0.0f, 0.0f, _j3trp0j->_x5qjicj, _j3trp0j->_plozjcn);
                                StrainFS_Nikkhoo(_kc4lorx, _dhgzz0j, _uug95rl[0], _uug95rl[1], _uug95rl[2], _e3nwxuj, _3ac6pe6, _obyghgk, 0.0f, _84u89yj->_31yrujo, 0.0f, _j3trp0j->_x5qjicj, _j3trp0j->_plozjcn);
                                StrainFS_Nikkhoo(_mt8vt43, _dhgzz0j, _uug95rl[0], _uug95rl[1], _uug95rl[2], _3ac6pe6, _3lbqgu3, _obyghgk, _84u89yj->_31yrujo, 0.0f, 0.0f, _j3trp0j->_x5qjicj, _j3trp0j->_plozjcn);
                                StrainFS_Nikkhoo(_pkxvhag, _dhgzz0j, _uug95rl[0], _uug95rl[1], _uug95rl[2], _3ac6pe6, _3lbqgu3, _obyghgk, 0.0f, _84u89yj->_31yrujo, 0.0f, _j3trp0j->_x5qjicj, _j3trp0j->_plozjcn);
                                _mt8vt43[0] += _pt5r5rc[0]; _mt8vt43[1] += _pt5r5rc[1]; _mt8vt43[2] += _pt5r5rc[2];
                                _mt8vt43[3] += _pt5r5rc[3]; _mt8vt43[4] += _pt5r5rc[4]; _mt8vt43[5] += _pt5r5rc[5];
                                _pkxvhag[0] += _kc4lorx[0]; _pkxvhag[1] += _kc4lorx[1]; _pkxvhag[2] += _kc4lorx[2];
                                _pkxvhag[3] += _kc4lorx[3]; _pkxvhag[4] += _kc4lorx[4]; _pkxvhag[5] += _kc4lorx[5];
                            }
                            ScaleStress6(_mt8vt43, (_8e7flqh/(_84u89yj->_31yrujo*_pqvl00e[_iua08ek][3]))); RotStressG2L(&_gwaj1g1[0], _v71448n, _mt8vt43);
                            ScaleStress6(_pkxvhag, (_8e7flqh/(_84u89yj->_31yrujo*_pqvl00e[_iua08ek][3]))); RotStressG2L(&_gwaj1g1[3], _v71448n, _pkxvhag);
                            __builtin___memcpy_chk (_dan085r, &_pqvl00e[_iua08ek][0], 3*sizeof(float), __builtin_object_size (_dan085r, 0));
                            SubtractVect3(_2ljk4q5, _uug95rl, _dan085r); _1vsmado = VectorLength(_2ljk4q5);
                            _fv177ud = (_1vsmado/_84u89yj->_d31k9qm - _84u89yj->_ktakpux)/(_84u89yj->_qphr62g - _84u89yj->_ktakpux);
                            _fv177ud = (_fv177ud >= 0.0f)*_fv177ud + (_fv177ud < 0.0f)*0.0f;
                            _fv177ud = (_fv177ud <= 1.0f)*_fv177ud + (_fv177ud > 1.0f)*1.0f;
                            _fv177ud = _84u89yj->_dswm31y*_fv177ud;
                            _66cnolj = sqrtf( (_gwaj1g1[0]-_2292mtr[0])*(_gwaj1g1[0]-_2292mtr[0]) + (_gwaj1g1[1]-_2292mtr[1])*(_gwaj1g1[1]-_2292mtr[1]) + (_gwaj1g1[2]-_2292mtr[2])*(_gwaj1g1[2]-_2292mtr[2]) + (_gwaj1g1[3]-_2292mtr[3])*(_gwaj1g1[3]-_2292mtr[3]) + (_gwaj1g1[4]-_2292mtr[4])*(_gwaj1g1[4]-_2292mtr[4]) + (_gwaj1g1[5]-_2292mtr[5])*(_gwaj1g1[5]-_2292mtr[5]) );
                            _q52wspo = sqrtf( _2292mtr[0] *_2292mtr[0] + _2292mtr[1] *_2292mtr[1] + _2292mtr[2] *_2292mtr[2] + _2292mtr[3] *_2292mtr[3] + _2292mtr[4] *_2292mtr[4] + _2292mtr[5] *_2292mtr[5] );
                            _1b07g2r = ((_66cnolj/_q52wspo) <= _fv177ud)*_1b07g2r + ((_66cnolj/_q52wspo) > _fv177ud)*0u;
                        }
                        if (_1b07g2r == 1u)
                        { _nf74t60[_hthsjr0[i]*2 +0] = _2292mtr[0]; _nf74t60[_hthsjr0[i]*2 +1] = _2292mtr[3];
                            _ge2j51b[_hthsjr0[i]*2 +0] = _2292mtr[1]; _ge2j51b[_hthsjr0[i]*2 +1] = _2292mtr[4];
                            _mrunjdu[_hthsjr0[i]*2 +0] = _2292mtr[2]; _mrunjdu[_hthsjr0[i]*2 +1] = _2292mtr[5];
                            _xdrkz40 = 0u;
                            for (j = 0u; j < _s3kkv9r->_otn9uze; j++)
                            { _rgb93wk[j] = (_f7fr91i[j][16 -_jg2y8rg +_fu7bfxq] == _be3yqy5)*1u + (_f7fr91i[j][16 -_jg2y8rg +_fu7bfxq] != _be3yqy5)*_rgb93wk[j];
                                (*_k2wbxan)[i][j] = (_f7fr91i[j][16 -_jg2y8rg +_fu7bfxq] == _be3yqy5)*_hthsjr0[i] + (_f7fr91i[j][16 -_jg2y8rg +_fu7bfxq] != _be3yqy5)*(*_k2wbxan)[i][j];
                                _xdrkz40 += _rgb93wk[j];
                            }
                            _hthsjr0[i] += 1u;
                            for (j = 0u; j < _jg2y8rg; j++) { _kb700uc[ _f7fr91i[_o7vbgr8][16 -_jg2y8rg +j] ] = 1u; }
                            _ybxkpm2 = 0u;
                        }
                    _fu7bfxq += 1u;
            } } }
            _050ui22 = GetByteSize(64, 2*_hthsjr0[i], sizeof(float));
            (*_0isb5iq)[i] = aligned_alloc(64, _050ui22); __builtin___memcpy_chk ((*_0isb5iq)[i], _nf74t60, 2*_hthsjr0[i]* sizeof(float), __builtin_object_size ((*_0isb5iq)[i], 0));
            (*_zzf8s2z)[i] = aligned_alloc(64, _050ui22); __builtin___memcpy_chk ((*_zzf8s2z)[i], _ge2j51b, 2*_hthsjr0[i]* sizeof(float), __builtin_object_size ((*_zzf8s2z)[i], 0));
            (*_3jczucq)[i] = aligned_alloc(64, _050ui22); __builtin___memcpy_chk ((*_3jczucq)[i], _mrunjdu, 2*_hthsjr0[i]* sizeof(float), __builtin_object_size ((*_3jczucq)[i], 0));
            if ((_s3kkv9r->_k00ilu6 == 1)&&(_tk63n2q == 0)) { fprintf(__stdoutp,"\nElement %d of %d    with     %10u resolved sources", i, _6zx1tv1[_tk63n2q], _hthsjr0[i]); }
            if (_s3kkv9r->_y25m8kb > 0u)
            { _050ui22 = GetByteSize(64, 3*_s3kkv9r->_y25m8kb, sizeof(float));
                __builtin___memset_chk (_ycz04cd, 0, _050ui22, __builtin_object_size (_ycz04cd, 0));
                __builtin___memset_chk (_ozctwjl, 0, _050ui22, __builtin_object_size (_ozctwjl, 0));
                __builtin___memset_chk (_ye0qzsa, 0, _s3kkv9r->_y25m8kb*sizeof(unsigned int), __builtin_object_size (_ye0qzsa, 0));
                __builtin___memset_chk (_0s1f537, 0, _kicwmpe*sizeof(unsigned int), __builtin_object_size (_0s1f537, 0));
                _xdrkz40 = 0u;
                while (_xdrkz40 < _s3kkv9r->_y25m8kb)
                {
                    for (j = 0u; j < _s3kkv9r->_y25m8kb; j++) { if (_ye0qzsa[j] == 0u) { _o7vbgr8 = j; break; } }
                    __builtin___memcpy_chk (_dan085r, &_s9bp4m1[_o7vbgr8*16 +0], 3u*sizeof(float), __builtin_object_size (_dan085r, 0));
                    SubtractVect3(_2ljk4q5, _uug95rl, _dan085r); _1vsmado = VectorLength(_2ljk4q5);
                    _8e7flqh = (_1vsmado >= 1.75f*_84u89yj->_m7pyxdw)*1.0f + (_1vsmado < 1.75f*_84u89yj->_m7pyxdw)*0.0f;
                    _jg2y8rg = _wshmlfz[_o7vbgr8][0];
                    _fu7bfxq = 0;
                    _ybxkpm2 = 1u;
                    while (_ybxkpm2 == 1u)
                    {
                        _1b07g2r = 0u;
                        while (_fu7bfxq < _jg2y8rg-1)
                        { _be3yqy5 = _wshmlfz[_o7vbgr8][16 -_jg2y8rg +_fu7bfxq];
                            __builtin___memcpy_chk (_dan085r, &_xf7jxjn[_be3yqy5][0], 3*sizeof(float), __builtin_object_size (_dan085r, 0));
                            SubtractVect3(_2ljk4q5, _uug95rl, _dan085r); _1vsmado = VectorLength(_2ljk4q5);
                            _1b07g2r = (_0s1f537[_be3yqy5] == 0u)*1u + (_0s1f537[_be3yqy5] != 0u)*0u;
                            _1b07g2r = (_jlpefja[_be3yqy5][0] > 1u)*_1b07g2r + (_jlpefja[_be3yqy5][0] <= 1u)*0u;
                            _1b07g2r = (_1vsmado > (_xf7jxjn[_be3yqy5][5] +_84u89yj->_ktakpux*_84u89yj->_m7pyxdw))*_1b07g2r + (_1vsmado <= (_xf7jxjn[_be3yqy5][5] +_84u89yj->_ktakpux*_84u89yj->_m7pyxdw))*0u;
                            if (_1b07g2r == 1u) { break; }
                            _fu7bfxq += 1u;
                        }
                        if (_1b07g2r == 0u)
                        { __builtin___memcpy_chk (_e3nwxuj, &_t8qzpky[_qydbpgv[_o7vbgr8*8 +0]*4 +0], 3*sizeof(float), __builtin_object_size (_e3nwxuj, 0));
                            __builtin___memcpy_chk (_3ac6pe6, &_t8qzpky[_qydbpgv[_o7vbgr8*8 +1]*4 +0], 3*sizeof(float), __builtin_object_size (_3ac6pe6, 0));
                            __builtin___memcpy_chk (_3lbqgu3, &_t8qzpky[_qydbpgv[_o7vbgr8*8 +2]*4 +0], 3*sizeof(float), __builtin_object_size (_3lbqgu3, 0));
                            if (_84u89yj->_fhv1xg1 == 1u)
                            { StrainHS_Nikkhoo(_mt8vt43, _dhgzz0j, _uug95rl[0], _uug95rl[1], _uug95rl[2], _e3nwxuj, _3ac6pe6, _3lbqgu3, _84u89yj->_8kzruin, 0.0f, 0.0f, _j3trp0j->_x5qjicj, _j3trp0j->_plozjcn);
                                StrainHS_Nikkhoo(_pkxvhag, _dhgzz0j, _uug95rl[0], _uug95rl[1], _uug95rl[2], _e3nwxuj, _3ac6pe6, _3lbqgu3, 0.0f, _84u89yj->_8kzruin, 0.0f, _j3trp0j->_x5qjicj, _j3trp0j->_plozjcn);
                                StrainHS_Nikkhoo(_dtg5nqj, _dhgzz0j, _uug95rl[0], _uug95rl[1], _uug95rl[2], _e3nwxuj, _3ac6pe6, _3lbqgu3, 0.0f, 0.0f, _84u89yj->_8kzruin, _j3trp0j->_x5qjicj, _j3trp0j->_plozjcn);
                            }
                            else
                            { StrainFS_Nikkhoo(_mt8vt43, _dhgzz0j, _uug95rl[0], _uug95rl[1], _uug95rl[2], _e3nwxuj, _3ac6pe6, _3lbqgu3, _84u89yj->_8kzruin, 0.0f, 0.0f, _j3trp0j->_x5qjicj, _j3trp0j->_plozjcn);
                                StrainFS_Nikkhoo(_pkxvhag, _dhgzz0j, _uug95rl[0], _uug95rl[1], _uug95rl[2], _e3nwxuj, _3ac6pe6, _3lbqgu3, 0.0f, _84u89yj->_8kzruin, 0.0f, _j3trp0j->_x5qjicj, _j3trp0j->_plozjcn);
                                StrainFS_Nikkhoo(_dtg5nqj, _dhgzz0j, _uug95rl[0], _uug95rl[1], _uug95rl[2], _e3nwxuj, _3ac6pe6, _3lbqgu3, 0.0f, 0.0f, _84u89yj->_8kzruin, _j3trp0j->_x5qjicj, _j3trp0j->_plozjcn);
                            }
                            ScaleStress6(_mt8vt43, (_8e7flqh/(_84u89yj->_8kzruin*_s9bp4m1[_o7vbgr8*16 +15]))); RotStressG2L(&_2292mtr[0], _v71448n, _mt8vt43);
                            ScaleStress6(_pkxvhag, (_8e7flqh/(_84u89yj->_8kzruin*_s9bp4m1[_o7vbgr8*16 +15]))); RotStressG2L(&_2292mtr[3], _v71448n, _pkxvhag);
                            ScaleStress6(_dtg5nqj, (_8e7flqh/(_84u89yj->_8kzruin*_s9bp4m1[_o7vbgr8*16 +15]))); RotStressG2L(&_2292mtr[6], _v71448n, _dtg5nqj);
                            _ycz04cd[_pmlvak9[i]*3 +0] = _2292mtr[0]; _ycz04cd[_pmlvak9[i]*3 +1] = _2292mtr[3]; _ycz04cd[_pmlvak9[i]*3 +2] = _2292mtr[6];
                            _ozctwjl[_pmlvak9[i]*3 +0] = _2292mtr[1]; _ozctwjl[_pmlvak9[i]*3 +1] = _2292mtr[4]; _ozctwjl[_pmlvak9[i]*3 +2] = _2292mtr[7];
                            (*_n14eb87)[i][_o7vbgr8] = _pmlvak9[i];
                            _pmlvak9[i] += 1u;
                            _ye0qzsa[_o7vbgr8] = 1u;
                            _xdrkz40 = 0u;
                            for (j = 0u; j < _s3kkv9r->_y25m8kb; j++) { _xdrkz40 += _ye0qzsa[j]; }
                            for (j = 0u; j < _jg2y8rg; j++) { _0s1f537[ _wshmlfz[_o7vbgr8][16 -_jg2y8rg +j] ] = 1u; }
                            _ybxkpm2 = 0u;
                        }
                        else
                        { _be3yqy5 = _wshmlfz[_o7vbgr8][16 -_jg2y8rg +_fu7bfxq];
                            __builtin___memcpy_chk (_e3nwxuj, &_xf7jxjn[_be3yqy5][9], 3*sizeof(float), __builtin_object_size (_e3nwxuj, 0));
                            __builtin___memcpy_chk (_3ac6pe6, &_xf7jxjn[_be3yqy5][12], 3*sizeof(float), __builtin_object_size (_3ac6pe6, 0));
                            __builtin___memcpy_chk (_3lbqgu3, &_xf7jxjn[_be3yqy5][15], 3*sizeof(float), __builtin_object_size (_3lbqgu3, 0));
                            __builtin___memcpy_chk (_obyghgk, &_xf7jxjn[_be3yqy5][18], 3*sizeof(float), __builtin_object_size (_obyghgk, 0));
                            if (_84u89yj->_fhv1xg1 == 1u)
                            { StrainHS_Nikkhoo(_pt5r5rc, _dhgzz0j, _uug95rl[0], _uug95rl[1], _uug95rl[2], _e3nwxuj, _3ac6pe6, _obyghgk, _84u89yj->_8kzruin, 0.0f, 0.0f, _j3trp0j->_x5qjicj, _j3trp0j->_plozjcn);
                                StrainHS_Nikkhoo(_kc4lorx, _dhgzz0j, _uug95rl[0], _uug95rl[1], _uug95rl[2], _e3nwxuj, _3ac6pe6, _obyghgk, 0.0f, _84u89yj->_8kzruin, 0.0f, _j3trp0j->_x5qjicj, _j3trp0j->_plozjcn);
                                StrainHS_Nikkhoo(_qijc1yd, _dhgzz0j, _uug95rl[0], _uug95rl[1], _uug95rl[2], _e3nwxuj, _3ac6pe6, _obyghgk, 0.0f, 0.0f, _84u89yj->_8kzruin, _j3trp0j->_x5qjicj, _j3trp0j->_plozjcn);
                                StrainHS_Nikkhoo(_mt8vt43, _dhgzz0j, _uug95rl[0], _uug95rl[1], _uug95rl[2], _3ac6pe6, _3lbqgu3, _obyghgk, _84u89yj->_8kzruin, 0.0f, 0.0f, _j3trp0j->_x5qjicj, _j3trp0j->_plozjcn);
                                StrainHS_Nikkhoo(_pkxvhag, _dhgzz0j, _uug95rl[0], _uug95rl[1], _uug95rl[2], _3ac6pe6, _3lbqgu3, _obyghgk, 0.0f, _84u89yj->_8kzruin, 0.0f, _j3trp0j->_x5qjicj, _j3trp0j->_plozjcn);
                                StrainHS_Nikkhoo(_dtg5nqj, _dhgzz0j, _uug95rl[0], _uug95rl[1], _uug95rl[2], _3ac6pe6, _3lbqgu3, _obyghgk, 0.0f, 0.0f, _84u89yj->_8kzruin, _j3trp0j->_x5qjicj, _j3trp0j->_plozjcn);
                                _mt8vt43[0] += _pt5r5rc[0]; _mt8vt43[1] += _pt5r5rc[1]; _mt8vt43[2] += _pt5r5rc[2];
                                _mt8vt43[3] += _pt5r5rc[3]; _mt8vt43[4] += _pt5r5rc[4]; _mt8vt43[5] += _pt5r5rc[5];
                                _pkxvhag[0] += _kc4lorx[0]; _pkxvhag[1] += _kc4lorx[1]; _pkxvhag[2] += _kc4lorx[2];
                                _pkxvhag[3] += _kc4lorx[3]; _pkxvhag[4] += _kc4lorx[4]; _pkxvhag[5] += _kc4lorx[5];
                                _dtg5nqj[0] += _qijc1yd[0]; _dtg5nqj[1] += _qijc1yd[1]; _dtg5nqj[2] += _qijc1yd[2];
                                _dtg5nqj[3] += _qijc1yd[3]; _dtg5nqj[4] += _qijc1yd[4]; _dtg5nqj[5] += _qijc1yd[5];
                            }
                            else
                            { StrainFS_Nikkhoo(_pt5r5rc, _dhgzz0j, _uug95rl[0], _uug95rl[1], _uug95rl[2], _e3nwxuj, _3ac6pe6, _obyghgk, _84u89yj->_8kzruin, 0.0f, 0.0f, _j3trp0j->_x5qjicj, _j3trp0j->_plozjcn);
                                StrainFS_Nikkhoo(_kc4lorx, _dhgzz0j, _uug95rl[0], _uug95rl[1], _uug95rl[2], _e3nwxuj, _3ac6pe6, _obyghgk, 0.0f, _84u89yj->_8kzruin, 0.0f, _j3trp0j->_x5qjicj, _j3trp0j->_plozjcn);
                                StrainFS_Nikkhoo(_qijc1yd, _dhgzz0j, _uug95rl[0], _uug95rl[1], _uug95rl[2], _e3nwxuj, _3ac6pe6, _obyghgk, 0.0f, 0.0f, _84u89yj->_8kzruin, _j3trp0j->_x5qjicj, _j3trp0j->_plozjcn);
                                StrainFS_Nikkhoo(_mt8vt43, _dhgzz0j, _uug95rl[0], _uug95rl[1], _uug95rl[2], _3ac6pe6, _3lbqgu3, _obyghgk, _84u89yj->_8kzruin, 0.0f, 0.0f, _j3trp0j->_x5qjicj, _j3trp0j->_plozjcn);
                                StrainFS_Nikkhoo(_pkxvhag, _dhgzz0j, _uug95rl[0], _uug95rl[1], _uug95rl[2], _3ac6pe6, _3lbqgu3, _obyghgk, 0.0f, _84u89yj->_8kzruin, 0.0f, _j3trp0j->_x5qjicj, _j3trp0j->_plozjcn);
                                StrainFS_Nikkhoo(_dtg5nqj, _dhgzz0j, _uug95rl[0], _uug95rl[1], _uug95rl[2], _3ac6pe6, _3lbqgu3, _obyghgk, 0.0f, 0.0f, _84u89yj->_8kzruin, _j3trp0j->_x5qjicj, _j3trp0j->_plozjcn);
                                _mt8vt43[0] += _pt5r5rc[0]; _mt8vt43[1] += _pt5r5rc[1]; _mt8vt43[2] += _pt5r5rc[2];
                                _mt8vt43[3] += _pt5r5rc[3]; _mt8vt43[4] += _pt5r5rc[4]; _mt8vt43[5] += _pt5r5rc[5];
                                _pkxvhag[0] += _kc4lorx[0]; _pkxvhag[1] += _kc4lorx[1]; _pkxvhag[2] += _kc4lorx[2];
                                _pkxvhag[3] += _kc4lorx[3]; _pkxvhag[4] += _kc4lorx[4]; _pkxvhag[5] += _kc4lorx[5];
                                _dtg5nqj[0] += _qijc1yd[0]; _dtg5nqj[1] += _qijc1yd[1]; _dtg5nqj[2] += _qijc1yd[2];
                                _dtg5nqj[3] += _qijc1yd[3]; _dtg5nqj[4] += _qijc1yd[4]; _dtg5nqj[5] += _qijc1yd[5];
                            }
                            ScaleStress6(_mt8vt43, (_8e7flqh/(_84u89yj->_8kzruin*_xf7jxjn[_be3yqy5][3]))); RotStressG2L(&_2292mtr[0], _v71448n, _mt8vt43);
                            ScaleStress6(_pkxvhag, (_8e7flqh/(_84u89yj->_8kzruin*_xf7jxjn[_be3yqy5][3]))); RotStressG2L(&_2292mtr[3], _v71448n, _pkxvhag);
                            ScaleStress6(_dtg5nqj, (_8e7flqh/(_84u89yj->_8kzruin*_xf7jxjn[_be3yqy5][3]))); RotStressG2L(&_2292mtr[6], _v71448n, _dtg5nqj);
                            _1b07g2r = 1u;
                            for (j = 0; j < _jlpefja[_be3yqy5][0]; j++)
                            { _iua08ek = _jlpefja[_be3yqy5][(j+1)];
                                __builtin___memcpy_chk (_e3nwxuj, &_xf7jxjn[_iua08ek][9], 3*sizeof(float), __builtin_object_size (_e3nwxuj, 0));
                                __builtin___memcpy_chk (_3ac6pe6, &_xf7jxjn[_iua08ek][12], 3*sizeof(float), __builtin_object_size (_3ac6pe6, 0));
                                __builtin___memcpy_chk (_3lbqgu3, &_xf7jxjn[_iua08ek][15], 3*sizeof(float), __builtin_object_size (_3lbqgu3, 0));
                                __builtin___memcpy_chk (_obyghgk, &_xf7jxjn[_iua08ek][18], 3*sizeof(float), __builtin_object_size (_obyghgk, 0));
                                if (_84u89yj->_fhv1xg1 == 1u)
                                { StrainHS_Nikkhoo(_pt5r5rc, _dhgzz0j, _uug95rl[0], _uug95rl[1], _uug95rl[2], _e3nwxuj, _3ac6pe6, _obyghgk, _84u89yj->_8kzruin, 0.0f, 0.0f, _j3trp0j->_x5qjicj, _j3trp0j->_plozjcn);
                                    StrainHS_Nikkhoo(_kc4lorx, _dhgzz0j, _uug95rl[0], _uug95rl[1], _uug95rl[2], _e3nwxuj, _3ac6pe6, _obyghgk, 0.0f, _84u89yj->_8kzruin, 0.0f, _j3trp0j->_x5qjicj, _j3trp0j->_plozjcn);
                                    StrainHS_Nikkhoo(_qijc1yd, _dhgzz0j, _uug95rl[0], _uug95rl[1], _uug95rl[2], _e3nwxuj, _3ac6pe6, _obyghgk, 0.0f, 0.0f, _84u89yj->_8kzruin, _j3trp0j->_x5qjicj, _j3trp0j->_plozjcn);
                                    StrainHS_Nikkhoo(_mt8vt43, _dhgzz0j, _uug95rl[0], _uug95rl[1], _uug95rl[2], _3ac6pe6, _3lbqgu3, _obyghgk, _84u89yj->_8kzruin, 0.0f, 0.0f, _j3trp0j->_x5qjicj, _j3trp0j->_plozjcn);
                                    StrainHS_Nikkhoo(_pkxvhag, _dhgzz0j, _uug95rl[0], _uug95rl[1], _uug95rl[2], _3ac6pe6, _3lbqgu3, _obyghgk, 0.0f, _84u89yj->_8kzruin, 0.0f, _j3trp0j->_x5qjicj, _j3trp0j->_plozjcn);
                                    StrainHS_Nikkhoo(_dtg5nqj, _dhgzz0j, _uug95rl[0], _uug95rl[1], _uug95rl[2], _3ac6pe6, _3lbqgu3, _obyghgk, 0.0f, 0.0f, _84u89yj->_8kzruin, _j3trp0j->_x5qjicj, _j3trp0j->_plozjcn);
                                    _mt8vt43[0] += _pt5r5rc[0]; _mt8vt43[1] += _pt5r5rc[1]; _mt8vt43[2] += _pt5r5rc[2];
                                    _mt8vt43[3] += _pt5r5rc[3]; _mt8vt43[4] += _pt5r5rc[4]; _mt8vt43[5] += _pt5r5rc[5];
                                    _pkxvhag[0] += _kc4lorx[0]; _pkxvhag[1] += _kc4lorx[1]; _pkxvhag[2] += _kc4lorx[2];
                                    _pkxvhag[3] += _kc4lorx[3]; _pkxvhag[4] += _kc4lorx[4]; _pkxvhag[5] += _kc4lorx[5];
                                    _dtg5nqj[0] += _qijc1yd[0]; _dtg5nqj[1] += _qijc1yd[1]; _dtg5nqj[2] += _qijc1yd[2];
                                    _dtg5nqj[3] += _qijc1yd[3]; _dtg5nqj[4] += _qijc1yd[4]; _dtg5nqj[5] += _qijc1yd[5];
                                }
                                else
                                { StrainFS_Nikkhoo(_pt5r5rc, _dhgzz0j, _uug95rl[0], _uug95rl[1], _uug95rl[2], _e3nwxuj, _3ac6pe6, _obyghgk, _84u89yj->_8kzruin, 0.0f, 0.0f, _j3trp0j->_x5qjicj, _j3trp0j->_plozjcn);
                                    StrainFS_Nikkhoo(_kc4lorx, _dhgzz0j, _uug95rl[0], _uug95rl[1], _uug95rl[2], _e3nwxuj, _3ac6pe6, _obyghgk, 0.0f, _84u89yj->_8kzruin, 0.0f, _j3trp0j->_x5qjicj, _j3trp0j->_plozjcn);
                                    StrainFS_Nikkhoo(_qijc1yd, _dhgzz0j, _uug95rl[0], _uug95rl[1], _uug95rl[2], _e3nwxuj, _3ac6pe6, _obyghgk, 0.0f, 0.0f, _84u89yj->_8kzruin, _j3trp0j->_x5qjicj, _j3trp0j->_plozjcn);
                                    StrainFS_Nikkhoo(_mt8vt43, _dhgzz0j, _uug95rl[0], _uug95rl[1], _uug95rl[2], _3ac6pe6, _3lbqgu3, _obyghgk, _84u89yj->_8kzruin, 0.0f, 0.0f, _j3trp0j->_x5qjicj, _j3trp0j->_plozjcn);
                                    StrainFS_Nikkhoo(_pkxvhag, _dhgzz0j, _uug95rl[0], _uug95rl[1], _uug95rl[2], _3ac6pe6, _3lbqgu3, _obyghgk, 0.0f, _84u89yj->_8kzruin, 0.0f, _j3trp0j->_x5qjicj, _j3trp0j->_plozjcn);
                                    StrainFS_Nikkhoo(_dtg5nqj, _dhgzz0j, _uug95rl[0], _uug95rl[1], _uug95rl[2], _3ac6pe6, _3lbqgu3, _obyghgk, 0.0f, 0.0f, _84u89yj->_8kzruin, _j3trp0j->_x5qjicj, _j3trp0j->_plozjcn);
                                    _mt8vt43[0] += _pt5r5rc[0]; _mt8vt43[1] += _pt5r5rc[1]; _mt8vt43[2] += _pt5r5rc[2];
                                    _mt8vt43[3] += _pt5r5rc[3]; _mt8vt43[4] += _pt5r5rc[4]; _mt8vt43[5] += _pt5r5rc[5];
                                    _pkxvhag[0] += _kc4lorx[0]; _pkxvhag[1] += _kc4lorx[1]; _pkxvhag[2] += _kc4lorx[2];
                                    _pkxvhag[3] += _kc4lorx[3]; _pkxvhag[4] += _kc4lorx[4]; _pkxvhag[5] += _kc4lorx[5];
                                    _dtg5nqj[0] += _qijc1yd[0]; _dtg5nqj[1] += _qijc1yd[1]; _dtg5nqj[2] += _qijc1yd[2];
                                    _dtg5nqj[3] += _qijc1yd[3]; _dtg5nqj[4] += _qijc1yd[4]; _dtg5nqj[5] += _qijc1yd[5];
                                }
                                ScaleStress6(_mt8vt43, (_8e7flqh/(_84u89yj->_8kzruin*_xf7jxjn[_iua08ek][3]))); RotStressG2L(&_gwaj1g1[0], _v71448n, _mt8vt43);
                                ScaleStress6(_pkxvhag, (_8e7flqh/(_84u89yj->_8kzruin*_xf7jxjn[_iua08ek][3]))); RotStressG2L(&_gwaj1g1[3], _v71448n, _pkxvhag);
                                ScaleStress6(_dtg5nqj, (_8e7flqh/(_84u89yj->_8kzruin*_xf7jxjn[_iua08ek][3]))); RotStressG2L(&_gwaj1g1[6], _v71448n, _dtg5nqj);
                                __builtin___memcpy_chk (_dan085r, &_xf7jxjn[_iua08ek][0], 3*sizeof(float), __builtin_object_size (_dan085r, 0));
                                SubtractVect3(_2ljk4q5, _uug95rl, _dan085r); _1vsmado = VectorLength(_2ljk4q5);
                                _fv177ud = (_1vsmado/_84u89yj->_m7pyxdw - _84u89yj->_ktakpux)/(_84u89yj->_qphr62g - _84u89yj->_ktakpux);
                                _fv177ud = (_fv177ud >= 0.0f)*_fv177ud + (_fv177ud < 0.0f)*0.0f;
                                _fv177ud = (_fv177ud <= 1.0f)*_fv177ud + (_fv177ud > 1.0f)*1.0f;
                                _fv177ud = _84u89yj->_dswm31y*_fv177ud;
                                _66cnolj = sqrtf( (_gwaj1g1[0]-_2292mtr[0])*(_gwaj1g1[0]-_2292mtr[0]) + (_gwaj1g1[1]-_2292mtr[1])*(_gwaj1g1[1]-_2292mtr[1]) + (_gwaj1g1[3]-_2292mtr[3])*(_gwaj1g1[3]-_2292mtr[3]) + (_gwaj1g1[4]-_2292mtr[4])*(_gwaj1g1[4]-_2292mtr[4]) + (_gwaj1g1[6]-_2292mtr[6])*(_gwaj1g1[6]-_2292mtr[6]) + (_gwaj1g1[7]-_2292mtr[7])*(_gwaj1g1[7]-_2292mtr[7]) );
                                _q52wspo = sqrtf( _2292mtr[0] *_2292mtr[0] + _2292mtr[1] *_2292mtr[1] + _2292mtr[3] *_2292mtr[3] + _2292mtr[4] *_2292mtr[4] + _2292mtr[6] *_2292mtr[6] + _2292mtr[7] *_2292mtr[7] );
                                _1b07g2r = ((_66cnolj/_q52wspo) <= _fv177ud)*_1b07g2r + ((_66cnolj/_q52wspo) > _fv177ud)*0u;
                            }
                            if (_1b07g2r == 1u)
                            { _ycz04cd[_pmlvak9[i]*3 +0] = _2292mtr[0]; _ycz04cd[_pmlvak9[i]*3 +1] = _2292mtr[3]; _ycz04cd[_pmlvak9[i]*3 +2] = _2292mtr[6];
                                _ozctwjl[_pmlvak9[i]*3 +0] = _2292mtr[1]; _ozctwjl[_pmlvak9[i]*3 +1] = _2292mtr[4]; _ozctwjl[_pmlvak9[i]*3 +2] = _2292mtr[7];
                                _xdrkz40 = 0u;
                                for (j = 0u; j < _s3kkv9r->_y25m8kb; j++)
                                { _ye0qzsa[j] = (_wshmlfz[j][16 -_jg2y8rg +_fu7bfxq] == _be3yqy5)*1u + (_wshmlfz[j][16 -_jg2y8rg +_fu7bfxq] != _be3yqy5)*_ye0qzsa[j];
                                    (*_n14eb87)[i][j] = (_wshmlfz[j][16 -_jg2y8rg +_fu7bfxq] == _be3yqy5)*_pmlvak9[i] + (_wshmlfz[j][16 -_jg2y8rg +_fu7bfxq] != _be3yqy5)*(*_n14eb87)[i][j];
                                    _xdrkz40 += _ye0qzsa[j];
                                }
                                _pmlvak9[i] += 1u;
                                for (j = 0u; j < _jg2y8rg; j++) { _0s1f537[ _wshmlfz[_o7vbgr8][16 -_jg2y8rg +j] ] = 1u; }
                                _ybxkpm2 = 0u;
                            }
                            _fu7bfxq += 1u;
                } } }
                _050ui22 = GetByteSize(64, 3*_pmlvak9[i], sizeof(float));
                (*_sh3tow4)[i] = aligned_alloc(64, _050ui22); __builtin___memcpy_chk ((*_sh3tow4)[i], _ycz04cd, 3*_pmlvak9[i]* sizeof(float), __builtin_object_size ((*_sh3tow4)[i], 0));
                (*_jbc2sws)[i] = aligned_alloc(64, _050ui22); __builtin___memcpy_chk ((*_jbc2sws)[i], _ozctwjl, 3*_pmlvak9[i]* sizeof(float), __builtin_object_size ((*_jbc2sws)[i], 0));
            }
            else
            {
            }
        }
        free(_nf74t60); free(_ge2j51b); free(_mrunjdu); free(_ycz04cd); free(_ozctwjl);
        if (_s3kkv9r->_y25m8kb > 0u)
        { _050ui22 = GetByteSize(64, 3*_s3kkv9r->_y25m8kb, sizeof(float));
            _0wfoy4d = aligned_alloc(64, _050ui22);
            _zareiya = aligned_alloc(64, _050ui22);
            _0dl5e11 = aligned_alloc(64, _050ui22);
            _050ui22 = GetByteSize(64, 2*_s3kkv9r->_otn9uze, sizeof(float));
            _ny181ac = aligned_alloc(64, _050ui22);
            _jdtvv2y = aligned_alloc(64, _050ui22);
            _b68rerv = aligned_alloc(64, _050ui22);
            _050ui22 = GetByteSize(64, _emruyp0[_tk63n2q], sizeof(unsigned short int *));
            *_va1nkbt = aligned_alloc(64, _050ui22);
            _050ui22 = GetByteSize(64, _s3kkv9r->_y25m8kb, sizeof(unsigned short int));
            for (i = 0u; i < _emruyp0[_tk63n2q]; i++)
            { (*_va1nkbt)[i] = aligned_alloc(64, _050ui22); __builtin___memset_chk ((*_va1nkbt)[i], 0, _050ui22, __builtin_object_size ((*_va1nkbt)[i], 0));
            }
            _050ui22 = GetByteSize(64, _emruyp0[_tk63n2q], sizeof(unsigned short int *));
            *_eurxhw0 = aligned_alloc(64, _050ui22);
            _050ui22 = GetByteSize(64, _s3kkv9r->_otn9uze, sizeof(unsigned short int));
            for (i = 0u; i < _emruyp0[_tk63n2q]; i++)
            { (*_eurxhw0)[i] = aligned_alloc(64, _050ui22); __builtin___memset_chk ((*_eurxhw0)[i], 0, _050ui22, __builtin_object_size ((*_eurxhw0)[i], 0));
            }
            _050ui22 = GetByteSize(64, _emruyp0[_tk63n2q], sizeof(float *));
            *_e9ftaxp = aligned_alloc(64, _050ui22);
            *_6fokzix = aligned_alloc(64, _050ui22);
            *_z39pecn = aligned_alloc(64, _050ui22);
            *_t997e4u = aligned_alloc(64, _050ui22);
            *_d1vn8lf = aligned_alloc(64, _050ui22);
            *_4xkss0k = aligned_alloc(64, _050ui22);
            for (i = 0u; i < _emruyp0[_tk63n2q]; i++)
            { _vmefozt = _laym0wb[i];
                _050ui22 = GetByteSize(64, 3*_s3kkv9r->_y25m8kb, sizeof(float));
                __builtin___memset_chk (_0wfoy4d, 0, _050ui22, __builtin_object_size (_0wfoy4d, 0));
                __builtin___memset_chk (_zareiya, 0, _050ui22, __builtin_object_size (_zareiya, 0));
                __builtin___memset_chk (_0dl5e11, 0, _050ui22, __builtin_object_size (_0dl5e11, 0));
                __builtin___memcpy_chk (_uug95rl, &_s9bp4m1[_vmefozt*16 +0], 3*sizeof(float), __builtin_object_size (_uug95rl, 0));
                __builtin___memcpy_chk (_v71448n, &_s9bp4m1[_vmefozt*16 +6], 9*sizeof(float), __builtin_object_size (_v71448n, 0));
                __builtin___memset_chk (_ye0qzsa, 0, _s3kkv9r->_y25m8kb*sizeof(unsigned int), __builtin_object_size (_ye0qzsa, 0));
                __builtin___memset_chk (_0s1f537, 0, _kicwmpe*sizeof(unsigned int), __builtin_object_size (_0s1f537, 0));
                _xdrkz40 = 0u;
                while (_xdrkz40 < _s3kkv9r->_y25m8kb)
                {
                    for (j = 0u; j < _s3kkv9r->_y25m8kb; j++) { if (_ye0qzsa[j] == 0u) { _o7vbgr8 = j; break; } }
                    __builtin___memcpy_chk (_dan085r, &_s9bp4m1[_o7vbgr8*16 +0], 3*sizeof(float), __builtin_object_size (_dan085r, 0));
                    SubtractVect3(_2ljk4q5, _uug95rl, _dan085r); _1vsmado = VectorLength(_2ljk4q5);
                    _8e7flqh = (_1vsmado >= 1.75f*_84u89yj->_m7pyxdw)*1.0f + (_1vsmado < 1.75f*_84u89yj->_m7pyxdw)*0.0f;
                    _8e7flqh = (_qydbpgv[_vmefozt*8 +3] == _qydbpgv[_o7vbgr8*8 +3])*1.0f + (_qydbpgv[_vmefozt*8 +3] != _qydbpgv[_o7vbgr8*8 +3])*_8e7flqh;
                    _jg2y8rg = _wshmlfz[_o7vbgr8][0];
                    _fu7bfxq = 0;
                    _ybxkpm2 = 1u;
                    while (_ybxkpm2 == 1u)
                    {
                        _1b07g2r = 0u;
                        while (_fu7bfxq < _jg2y8rg-1)
                        { _be3yqy5 = _wshmlfz[_o7vbgr8][16 -_jg2y8rg +_fu7bfxq];
                            __builtin___memcpy_chk (_dan085r, &_xf7jxjn[_be3yqy5][0], 3*sizeof(float), __builtin_object_size (_dan085r, 0));
                            SubtractVect3(_2ljk4q5, _uug95rl, _dan085r); _1vsmado = VectorLength(_2ljk4q5);
                            _1b07g2r = (_0s1f537[_be3yqy5] == 0u)*1u + (_0s1f537[_be3yqy5] != 0u)*0u;
                            _1b07g2r = (_jlpefja[_be3yqy5][0] > 1u)*_1b07g2r + (_jlpefja[_be3yqy5][0] <= 1u)*0u;
                            _1b07g2r = (_wshmlfz[_vmefozt][16 -_jg2y8rg +_fu7bfxq] != _be3yqy5)*_1b07g2r + (_wshmlfz[_vmefozt][16 -_jg2y8rg +_fu7bfxq] == _be3yqy5)*0u;
                            _1b07g2r = (_1vsmado > (_xf7jxjn[_be3yqy5][5] +_84u89yj->_ktakpux*_84u89yj->_m7pyxdw))*_1b07g2r + (_1vsmado <= (_xf7jxjn[_be3yqy5][5] +_84u89yj->_ktakpux*_84u89yj->_m7pyxdw))*0u;
                            if (_1b07g2r == 1u) { break; }
                            _fu7bfxq += 1u;
                        }
                        if (_1b07g2r == 0u)
                        { __builtin___memcpy_chk (_e3nwxuj, &_t8qzpky[_qydbpgv[_o7vbgr8*8 +0]*4 +0], 3*sizeof(float), __builtin_object_size (_e3nwxuj, 0));
                            __builtin___memcpy_chk (_3ac6pe6, &_t8qzpky[_qydbpgv[_o7vbgr8*8 +1]*4 +0], 3*sizeof(float), __builtin_object_size (_3ac6pe6, 0));
                            __builtin___memcpy_chk (_3lbqgu3, &_t8qzpky[_qydbpgv[_o7vbgr8*8 +2]*4 +0], 3*sizeof(float), __builtin_object_size (_3lbqgu3, 0));
                            if (_84u89yj->_fhv1xg1 == 1u)
                            { StrainHS_Nikkhoo(_mt8vt43, _dhgzz0j, _uug95rl[0], _uug95rl[1], _uug95rl[2], _e3nwxuj, _3ac6pe6, _3lbqgu3, _84u89yj->_8kzruin, 0.0f, 0.0f, _j3trp0j->_x5qjicj, _j3trp0j->_plozjcn);
                                StrainHS_Nikkhoo(_pkxvhag, _dhgzz0j, _uug95rl[0], _uug95rl[1], _uug95rl[2], _e3nwxuj, _3ac6pe6, _3lbqgu3, 0.0f, _84u89yj->_8kzruin, 0.0f, _j3trp0j->_x5qjicj, _j3trp0j->_plozjcn);
                                StrainHS_Nikkhoo(_dtg5nqj, _dhgzz0j, _uug95rl[0], _uug95rl[1], _uug95rl[2], _e3nwxuj, _3ac6pe6, _3lbqgu3, 0.0f, 0.0f, _84u89yj->_8kzruin, _j3trp0j->_x5qjicj, _j3trp0j->_plozjcn);
                            }
                            else
                            { StrainFS_Nikkhoo(_mt8vt43, _dhgzz0j, _uug95rl[0], _uug95rl[1], _uug95rl[2], _e3nwxuj, _3ac6pe6, _3lbqgu3, _84u89yj->_8kzruin, 0.0f, 0.0f, _j3trp0j->_x5qjicj, _j3trp0j->_plozjcn);
                                StrainFS_Nikkhoo(_pkxvhag, _dhgzz0j, _uug95rl[0], _uug95rl[1], _uug95rl[2], _e3nwxuj, _3ac6pe6, _3lbqgu3, 0.0f, _84u89yj->_8kzruin, 0.0f, _j3trp0j->_x5qjicj, _j3trp0j->_plozjcn);
                                StrainFS_Nikkhoo(_dtg5nqj, _dhgzz0j, _uug95rl[0], _uug95rl[1], _uug95rl[2], _e3nwxuj, _3ac6pe6, _3lbqgu3, 0.0f, 0.0f, _84u89yj->_8kzruin, _j3trp0j->_x5qjicj, _j3trp0j->_plozjcn);
                            }
                            ScaleStress6(_mt8vt43, (_8e7flqh/(_84u89yj->_8kzruin*_s9bp4m1[_o7vbgr8*16 +15]))); RotStressG2L(&_2292mtr[0], _v71448n, _mt8vt43);
                            ScaleStress6(_pkxvhag, (_8e7flqh/(_84u89yj->_8kzruin*_s9bp4m1[_o7vbgr8*16 +15]))); RotStressG2L(&_2292mtr[3], _v71448n, _pkxvhag);
                            ScaleStress6(_dtg5nqj, (_8e7flqh/(_84u89yj->_8kzruin*_s9bp4m1[_o7vbgr8*16 +15]))); RotStressG2L(&_2292mtr[6], _v71448n, _dtg5nqj);
                            if (_vmefozt == _o7vbgr8)
                            {
                                _dntcagr[i*16 +5] = _2292mtr[0]*_s9bp4m1[_o7vbgr8*16 +15];
                                _dntcagr[i*16 +6] = _2292mtr[4]*_s9bp4m1[_o7vbgr8*16 +15];
                                _dntcagr[i*16 +7] = _2292mtr[8]*_s9bp4m1[_o7vbgr8*16 +15];
                            }
                            _0wfoy4d[_vk0ucc1[i]*3 +0] = _2292mtr[0]; _0wfoy4d[_vk0ucc1[i]*3 +1] = _2292mtr[3]; _0wfoy4d[_vk0ucc1[i]*3 +2] = _2292mtr[6];
                            _zareiya[_vk0ucc1[i]*3 +0] = _2292mtr[1]; _zareiya[_vk0ucc1[i]*3 +1] = _2292mtr[4]; _zareiya[_vk0ucc1[i]*3 +2] = _2292mtr[7];
                            _0dl5e11[_vk0ucc1[i]*3 +0] = _2292mtr[2]; _0dl5e11[_vk0ucc1[i]*3 +1] = _2292mtr[5]; _0dl5e11[_vk0ucc1[i]*3 +2] = _2292mtr[8];
                            (*_va1nkbt)[i][_o7vbgr8] = _vk0ucc1[i];
                            _vk0ucc1[i] += 1u;
                            _ye0qzsa[_o7vbgr8] = 1u;
                            _xdrkz40 = 0u;
                            for (j = 0u; j < _s3kkv9r->_y25m8kb; j++) { _xdrkz40 += _ye0qzsa[j]; }
                            for (j = 0u; j < _jg2y8rg; j++) { _0s1f537[ _wshmlfz[_o7vbgr8][16 -_jg2y8rg +j] ] = 1u; }
                            _ybxkpm2 = 0u;
                        }
                        else
                        { _be3yqy5 = _wshmlfz[_o7vbgr8][16 -_jg2y8rg +_fu7bfxq];
                            __builtin___memcpy_chk (_e3nwxuj, &_xf7jxjn[_be3yqy5][9], 3*sizeof(float), __builtin_object_size (_e3nwxuj, 0));
                            __builtin___memcpy_chk (_3ac6pe6, &_xf7jxjn[_be3yqy5][12], 3*sizeof(float), __builtin_object_size (_3ac6pe6, 0));
                            __builtin___memcpy_chk (_3lbqgu3, &_xf7jxjn[_be3yqy5][15], 3*sizeof(float), __builtin_object_size (_3lbqgu3, 0));
                            __builtin___memcpy_chk (_obyghgk, &_xf7jxjn[_be3yqy5][18], 3*sizeof(float), __builtin_object_size (_obyghgk, 0));
                            if (_84u89yj->_fhv1xg1 == 1u)
                            { StrainHS_Nikkhoo(_pt5r5rc, _dhgzz0j, _uug95rl[0], _uug95rl[1], _uug95rl[2], _e3nwxuj, _3ac6pe6, _obyghgk, _84u89yj->_8kzruin, 0.0f, 0.0f, _j3trp0j->_x5qjicj, _j3trp0j->_plozjcn);
                                StrainHS_Nikkhoo(_kc4lorx, _dhgzz0j, _uug95rl[0], _uug95rl[1], _uug95rl[2], _e3nwxuj, _3ac6pe6, _obyghgk, 0.0f, _84u89yj->_8kzruin, 0.0f, _j3trp0j->_x5qjicj, _j3trp0j->_plozjcn);
                                StrainHS_Nikkhoo(_qijc1yd, _dhgzz0j, _uug95rl[0], _uug95rl[1], _uug95rl[2], _e3nwxuj, _3ac6pe6, _obyghgk, 0.0f, 0.0f, _84u89yj->_8kzruin, _j3trp0j->_x5qjicj, _j3trp0j->_plozjcn);
                                StrainHS_Nikkhoo(_mt8vt43, _dhgzz0j, _uug95rl[0], _uug95rl[1], _uug95rl[2], _3ac6pe6, _3lbqgu3, _obyghgk, _84u89yj->_8kzruin, 0.0f, 0.0f, _j3trp0j->_x5qjicj, _j3trp0j->_plozjcn);
                                StrainHS_Nikkhoo(_pkxvhag, _dhgzz0j, _uug95rl[0], _uug95rl[1], _uug95rl[2], _3ac6pe6, _3lbqgu3, _obyghgk, 0.0f, _84u89yj->_8kzruin, 0.0f, _j3trp0j->_x5qjicj, _j3trp0j->_plozjcn);
                                StrainHS_Nikkhoo(_dtg5nqj, _dhgzz0j, _uug95rl[0], _uug95rl[1], _uug95rl[2], _3ac6pe6, _3lbqgu3, _obyghgk, 0.0f, 0.0f, _84u89yj->_8kzruin, _j3trp0j->_x5qjicj, _j3trp0j->_plozjcn);
                                _mt8vt43[0] += _pt5r5rc[0]; _mt8vt43[1] += _pt5r5rc[1]; _mt8vt43[2] += _pt5r5rc[2];
                                _mt8vt43[3] += _pt5r5rc[3]; _mt8vt43[4] += _pt5r5rc[4]; _mt8vt43[5] += _pt5r5rc[5];
                                _pkxvhag[0] += _kc4lorx[0]; _pkxvhag[1] += _kc4lorx[1]; _pkxvhag[2] += _kc4lorx[2];
                                _pkxvhag[3] += _kc4lorx[3]; _pkxvhag[4] += _kc4lorx[4]; _pkxvhag[5] += _kc4lorx[5];
                                _dtg5nqj[0] += _qijc1yd[0]; _dtg5nqj[1] += _qijc1yd[1]; _dtg5nqj[2] += _qijc1yd[2];
                                _dtg5nqj[3] += _qijc1yd[3]; _dtg5nqj[4] += _qijc1yd[4]; _dtg5nqj[5] += _qijc1yd[5];
                            }
                            else
                            { StrainFS_Nikkhoo(_pt5r5rc, _dhgzz0j, _uug95rl[0], _uug95rl[1], _uug95rl[2], _e3nwxuj, _3ac6pe6, _obyghgk, _84u89yj->_8kzruin, 0.0f, 0.0f, _j3trp0j->_x5qjicj, _j3trp0j->_plozjcn);
                                StrainFS_Nikkhoo(_kc4lorx, _dhgzz0j, _uug95rl[0], _uug95rl[1], _uug95rl[2], _e3nwxuj, _3ac6pe6, _obyghgk, 0.0f, _84u89yj->_8kzruin, 0.0f, _j3trp0j->_x5qjicj, _j3trp0j->_plozjcn);
                                StrainFS_Nikkhoo(_qijc1yd, _dhgzz0j, _uug95rl[0], _uug95rl[1], _uug95rl[2], _e3nwxuj, _3ac6pe6, _obyghgk, 0.0f, 0.0f, _84u89yj->_8kzruin, _j3trp0j->_x5qjicj, _j3trp0j->_plozjcn);
                                StrainFS_Nikkhoo(_mt8vt43, _dhgzz0j, _uug95rl[0], _uug95rl[1], _uug95rl[2], _3ac6pe6, _3lbqgu3, _obyghgk, _84u89yj->_8kzruin, 0.0f, 0.0f, _j3trp0j->_x5qjicj, _j3trp0j->_plozjcn);
                                StrainFS_Nikkhoo(_pkxvhag, _dhgzz0j, _uug95rl[0], _uug95rl[1], _uug95rl[2], _3ac6pe6, _3lbqgu3, _obyghgk, 0.0f, _84u89yj->_8kzruin, 0.0f, _j3trp0j->_x5qjicj, _j3trp0j->_plozjcn);
                                StrainFS_Nikkhoo(_dtg5nqj, _dhgzz0j, _uug95rl[0], _uug95rl[1], _uug95rl[2], _3ac6pe6, _3lbqgu3, _obyghgk, 0.0f, 0.0f, _84u89yj->_8kzruin, _j3trp0j->_x5qjicj, _j3trp0j->_plozjcn);
                                _mt8vt43[0] += _pt5r5rc[0]; _mt8vt43[1] += _pt5r5rc[1]; _mt8vt43[2] += _pt5r5rc[2];
                                _mt8vt43[3] += _pt5r5rc[3]; _mt8vt43[4] += _pt5r5rc[4]; _mt8vt43[5] += _pt5r5rc[5];
                                _pkxvhag[0] += _kc4lorx[0]; _pkxvhag[1] += _kc4lorx[1]; _pkxvhag[2] += _kc4lorx[2];
                                _pkxvhag[3] += _kc4lorx[3]; _pkxvhag[4] += _kc4lorx[4]; _pkxvhag[5] += _kc4lorx[5];
                                _dtg5nqj[0] += _qijc1yd[0]; _dtg5nqj[1] += _qijc1yd[1]; _dtg5nqj[2] += _qijc1yd[2];
                                _dtg5nqj[3] += _qijc1yd[3]; _dtg5nqj[4] += _qijc1yd[4]; _dtg5nqj[5] += _qijc1yd[5];
                            }
                            ScaleStress6(_mt8vt43, (_8e7flqh/(_84u89yj->_8kzruin*_xf7jxjn[_be3yqy5][3]))); RotStressG2L(&_2292mtr[0], _v71448n, _mt8vt43);
                            ScaleStress6(_pkxvhag, (_8e7flqh/(_84u89yj->_8kzruin*_xf7jxjn[_be3yqy5][3]))); RotStressG2L(&_2292mtr[3], _v71448n, _pkxvhag);
                            ScaleStress6(_dtg5nqj, (_8e7flqh/(_84u89yj->_8kzruin*_xf7jxjn[_be3yqy5][3]))); RotStressG2L(&_2292mtr[6], _v71448n, _dtg5nqj);
                            _1b07g2r = 1u;
                            for (j = 0; j < _jlpefja[_be3yqy5][0]; j++)
                            { _iua08ek = _jlpefja[_be3yqy5][(j+1)];
                                __builtin___memcpy_chk (_e3nwxuj, &_xf7jxjn[_iua08ek][9], 3*sizeof(float), __builtin_object_size (_e3nwxuj, 0));
                                __builtin___memcpy_chk (_3ac6pe6, &_xf7jxjn[_iua08ek][12], 3*sizeof(float), __builtin_object_size (_3ac6pe6, 0));
                                __builtin___memcpy_chk (_3lbqgu3, &_xf7jxjn[_iua08ek][15], 3*sizeof(float), __builtin_object_size (_3lbqgu3, 0));
                                __builtin___memcpy_chk (_obyghgk, &_xf7jxjn[_iua08ek][18], 3*sizeof(float), __builtin_object_size (_obyghgk, 0));
                                if (_84u89yj->_fhv1xg1 == 1u)
                                { StrainHS_Nikkhoo(_pt5r5rc, _dhgzz0j, _uug95rl[0], _uug95rl[1], _uug95rl[2], _e3nwxuj, _3ac6pe6, _obyghgk, _84u89yj->_8kzruin, 0.0f, 0.0f, _j3trp0j->_x5qjicj, _j3trp0j->_plozjcn);
                                    StrainHS_Nikkhoo(_kc4lorx, _dhgzz0j, _uug95rl[0], _uug95rl[1], _uug95rl[2], _e3nwxuj, _3ac6pe6, _obyghgk, 0.0f, _84u89yj->_8kzruin, 0.0f, _j3trp0j->_x5qjicj, _j3trp0j->_plozjcn);
                                    StrainHS_Nikkhoo(_qijc1yd, _dhgzz0j, _uug95rl[0], _uug95rl[1], _uug95rl[2], _e3nwxuj, _3ac6pe6, _obyghgk, 0.0f, 0.0f, _84u89yj->_8kzruin, _j3trp0j->_x5qjicj, _j3trp0j->_plozjcn);
                                    StrainHS_Nikkhoo(_mt8vt43, _dhgzz0j, _uug95rl[0], _uug95rl[1], _uug95rl[2], _3ac6pe6, _3lbqgu3, _obyghgk, _84u89yj->_8kzruin, 0.0f, 0.0f, _j3trp0j->_x5qjicj, _j3trp0j->_plozjcn);
                                    StrainHS_Nikkhoo(_pkxvhag, _dhgzz0j, _uug95rl[0], _uug95rl[1], _uug95rl[2], _3ac6pe6, _3lbqgu3, _obyghgk, 0.0f, _84u89yj->_8kzruin, 0.0f, _j3trp0j->_x5qjicj, _j3trp0j->_plozjcn);
                                    StrainHS_Nikkhoo(_dtg5nqj, _dhgzz0j, _uug95rl[0], _uug95rl[1], _uug95rl[2], _3ac6pe6, _3lbqgu3, _obyghgk, 0.0f, 0.0f, _84u89yj->_8kzruin, _j3trp0j->_x5qjicj, _j3trp0j->_plozjcn);
                                    _mt8vt43[0] += _pt5r5rc[0]; _mt8vt43[1] += _pt5r5rc[1]; _mt8vt43[2] += _pt5r5rc[2];
                                    _mt8vt43[3] += _pt5r5rc[3]; _mt8vt43[4] += _pt5r5rc[4]; _mt8vt43[5] += _pt5r5rc[5];
                                    _pkxvhag[0] += _kc4lorx[0]; _pkxvhag[1] += _kc4lorx[1]; _pkxvhag[2] += _kc4lorx[2];
                                    _pkxvhag[3] += _kc4lorx[3]; _pkxvhag[4] += _kc4lorx[4]; _pkxvhag[5] += _kc4lorx[5];
                                    _dtg5nqj[0] += _qijc1yd[0]; _dtg5nqj[1] += _qijc1yd[1]; _dtg5nqj[2] += _qijc1yd[2];
                                    _dtg5nqj[3] += _qijc1yd[3]; _dtg5nqj[4] += _qijc1yd[4]; _dtg5nqj[5] += _qijc1yd[5];
                                }
                                else
                                { StrainFS_Nikkhoo(_pt5r5rc, _dhgzz0j, _uug95rl[0], _uug95rl[1], _uug95rl[2], _e3nwxuj, _3ac6pe6, _obyghgk, _84u89yj->_8kzruin, 0.0f, 0.0f, _j3trp0j->_x5qjicj, _j3trp0j->_plozjcn);
                                    StrainFS_Nikkhoo(_kc4lorx, _dhgzz0j, _uug95rl[0], _uug95rl[1], _uug95rl[2], _e3nwxuj, _3ac6pe6, _obyghgk, 0.0f, _84u89yj->_8kzruin, 0.0f, _j3trp0j->_x5qjicj, _j3trp0j->_plozjcn);
                                    StrainFS_Nikkhoo(_qijc1yd, _dhgzz0j, _uug95rl[0], _uug95rl[1], _uug95rl[2], _e3nwxuj, _3ac6pe6, _obyghgk, 0.0f, 0.0f, _84u89yj->_8kzruin, _j3trp0j->_x5qjicj, _j3trp0j->_plozjcn);
                                    StrainFS_Nikkhoo(_mt8vt43, _dhgzz0j, _uug95rl[0], _uug95rl[1], _uug95rl[2], _3ac6pe6, _3lbqgu3, _obyghgk, _84u89yj->_8kzruin, 0.0f, 0.0f, _j3trp0j->_x5qjicj, _j3trp0j->_plozjcn);
                                    StrainFS_Nikkhoo(_pkxvhag, _dhgzz0j, _uug95rl[0], _uug95rl[1], _uug95rl[2], _3ac6pe6, _3lbqgu3, _obyghgk, 0.0f, _84u89yj->_8kzruin, 0.0f, _j3trp0j->_x5qjicj, _j3trp0j->_plozjcn);
                                    StrainFS_Nikkhoo(_dtg5nqj, _dhgzz0j, _uug95rl[0], _uug95rl[1], _uug95rl[2], _3ac6pe6, _3lbqgu3, _obyghgk, 0.0f, 0.0f, _84u89yj->_8kzruin, _j3trp0j->_x5qjicj, _j3trp0j->_plozjcn);
                                    _mt8vt43[0] += _pt5r5rc[0]; _mt8vt43[1] += _pt5r5rc[1]; _mt8vt43[2] += _pt5r5rc[2];
                                    _mt8vt43[3] += _pt5r5rc[3]; _mt8vt43[4] += _pt5r5rc[4]; _mt8vt43[5] += _pt5r5rc[5];
                                    _pkxvhag[0] += _kc4lorx[0]; _pkxvhag[1] += _kc4lorx[1]; _pkxvhag[2] += _kc4lorx[2];
                                    _pkxvhag[3] += _kc4lorx[3]; _pkxvhag[4] += _kc4lorx[4]; _pkxvhag[5] += _kc4lorx[5];
                                    _dtg5nqj[0] += _qijc1yd[0]; _dtg5nqj[1] += _qijc1yd[1]; _dtg5nqj[2] += _qijc1yd[2];
                                    _dtg5nqj[3] += _qijc1yd[3]; _dtg5nqj[4] += _qijc1yd[4]; _dtg5nqj[5] += _qijc1yd[5];
                                }
                                ScaleStress6(_mt8vt43, (_8e7flqh/(_84u89yj->_8kzruin*_xf7jxjn[_iua08ek][3]))); RotStressG2L(&_gwaj1g1[0], _v71448n, _mt8vt43);
                                ScaleStress6(_pkxvhag, (_8e7flqh/(_84u89yj->_8kzruin*_xf7jxjn[_iua08ek][3]))); RotStressG2L(&_gwaj1g1[3], _v71448n, _pkxvhag);
                                ScaleStress6(_dtg5nqj, (_8e7flqh/(_84u89yj->_8kzruin*_xf7jxjn[_iua08ek][3]))); RotStressG2L(&_gwaj1g1[6], _v71448n, _dtg5nqj);
                                __builtin___memcpy_chk (_dan085r, &_xf7jxjn[_iua08ek][0], 3*sizeof(float), __builtin_object_size (_dan085r, 0));
                                SubtractVect3(_2ljk4q5, _uug95rl, _dan085r); _1vsmado = VectorLength(_2ljk4q5);
                                _fv177ud = (_1vsmado/_84u89yj->_m7pyxdw - _84u89yj->_ktakpux)/(_84u89yj->_qphr62g - _84u89yj->_ktakpux);
                                _fv177ud = (_fv177ud >= 0.0f)*_fv177ud + (_fv177ud < 0.0f)*0.0f;
                                _fv177ud = (_fv177ud <= 1.0f)*_fv177ud + (_fv177ud > 1.0f)*1.0f;
                                _fv177ud = _84u89yj->_dswm31y*_fv177ud;
                                _66cnolj = sqrtf( (_gwaj1g1[0]-_2292mtr[0])*(_gwaj1g1[0]-_2292mtr[0]) + (_gwaj1g1[1]-_2292mtr[1])*(_gwaj1g1[1]-_2292mtr[1]) + (_gwaj1g1[2]-_2292mtr[2])*(_gwaj1g1[2]-_2292mtr[2]) + (_gwaj1g1[3]-_2292mtr[3])*(_gwaj1g1[3]-_2292mtr[3]) + (_gwaj1g1[4]-_2292mtr[4])*(_gwaj1g1[4]-_2292mtr[4]) + (_gwaj1g1[5]-_2292mtr[5])*(_gwaj1g1[5]-_2292mtr[5]) + (_gwaj1g1[6]-_2292mtr[6])*(_gwaj1g1[6]-_2292mtr[6]) + (_gwaj1g1[7]-_2292mtr[7])*(_gwaj1g1[7]-_2292mtr[7]) + (_gwaj1g1[8]-_2292mtr[8])*(_gwaj1g1[8]-_2292mtr[8]));
                                _q52wspo = sqrtf( _2292mtr[0] *_2292mtr[0] + _2292mtr[1] *_2292mtr[1] + _2292mtr[2] *_2292mtr[2] + _2292mtr[3] *_2292mtr[3] + _2292mtr[4] *_2292mtr[4] + _2292mtr[5] *_2292mtr[5] + _2292mtr[6] *_2292mtr[6] + _2292mtr[7] *_2292mtr[7] + _2292mtr[8] *_2292mtr[8]);
                                _1b07g2r = ((_66cnolj/_q52wspo) <= _fv177ud)*_1b07g2r + ((_66cnolj/_q52wspo) > _fv177ud)*0u;
                            }
                            if (_1b07g2r == 1u)
                            { _0wfoy4d[_vk0ucc1[i]*3 +0] = _2292mtr[0]; _0wfoy4d[_vk0ucc1[i]*3 +1] = _2292mtr[3]; _0wfoy4d[_vk0ucc1[i]*3 +2] = _2292mtr[6];
                                _zareiya[_vk0ucc1[i]*3 +0] = _2292mtr[1]; _zareiya[_vk0ucc1[i]*3 +1] = _2292mtr[4]; _zareiya[_vk0ucc1[i]*3 +2] = _2292mtr[7];
                                _0dl5e11[_vk0ucc1[i]*3 +0] = _2292mtr[2]; _0dl5e11[_vk0ucc1[i]*3 +1] = _2292mtr[5]; _0dl5e11[_vk0ucc1[i]*3 +2] = _2292mtr[8];
                                _xdrkz40 = 0u;
                                for (j = 0u; j < _s3kkv9r->_y25m8kb; j++)
                                { _ye0qzsa[j] = (_wshmlfz[j][16 -_jg2y8rg +_fu7bfxq] == _be3yqy5)*1u + (_wshmlfz[j][16 -_jg2y8rg +_fu7bfxq] != _be3yqy5)*_ye0qzsa[j];
                                    (*_va1nkbt)[i][j] = (_wshmlfz[j][16 -_jg2y8rg +_fu7bfxq] == _be3yqy5)*_vk0ucc1[i] + (_wshmlfz[j][16 -_jg2y8rg +_fu7bfxq] != _be3yqy5)*(*_va1nkbt)[i][j];
                                    _xdrkz40 += _ye0qzsa[j];
                                }
                                _vk0ucc1[i] += 1u;
                                for (j = 0u; j < _jg2y8rg; j++) { _0s1f537[ _wshmlfz[_o7vbgr8][16 -_jg2y8rg +j] ] = 1u; }
                                _ybxkpm2 = 0u;
                            }
                            _fu7bfxq += 1u;
                } } }
                _050ui22 = GetByteSize(64, 3*_vk0ucc1[i], sizeof(float));
                (*_e9ftaxp)[i] = aligned_alloc(64, _050ui22); __builtin___memcpy_chk ((*_e9ftaxp)[i], _0wfoy4d, 3*_vk0ucc1[i]* sizeof(float), __builtin_object_size ((*_e9ftaxp)[i], 0));
                (*_6fokzix)[i] = aligned_alloc(64, _050ui22); __builtin___memcpy_chk ((*_6fokzix)[i], _zareiya, 3*_vk0ucc1[i]* sizeof(float), __builtin_object_size ((*_6fokzix)[i], 0));
                (*_z39pecn)[i] = aligned_alloc(64, _050ui22); __builtin___memcpy_chk ((*_z39pecn)[i], _0dl5e11, 3*_vk0ucc1[i]* sizeof(float), __builtin_object_size ((*_z39pecn)[i], 0));
                _050ui22 = GetByteSize(64, 2*_s3kkv9r->_otn9uze, sizeof(float));
                __builtin___memset_chk (_ny181ac, 0, _050ui22, __builtin_object_size (_ny181ac, 0));
                __builtin___memset_chk (_jdtvv2y, 0, _050ui22, __builtin_object_size (_jdtvv2y, 0));
                __builtin___memset_chk (_b68rerv, 0, _050ui22, __builtin_object_size (_b68rerv, 0));
                __builtin___memset_chk (_rgb93wk, 0, _s3kkv9r->_otn9uze*sizeof(unsigned int), __builtin_object_size (_rgb93wk, 0));
                __builtin___memset_chk (_kb700uc, 0, _pwh0mx1*sizeof(unsigned int), __builtin_object_size (_kb700uc, 0));
                _xdrkz40 = 0u;
                while (_xdrkz40 < _s3kkv9r->_otn9uze)
                {
                    for (j = 0u; j < _s3kkv9r->_otn9uze; j++) { if (_rgb93wk[j] == 0u) { _o7vbgr8 = j; break; } }
                    __builtin___memcpy_chk (_dan085r, &_3n3bzt9[_o7vbgr8*16 +0], 3*sizeof(float), __builtin_object_size (_dan085r, 0));
                    SubtractVect3(_2ljk4q5, _uug95rl, _dan085r); _1vsmado = VectorLength(_2ljk4q5);
                    _8e7flqh = (_1vsmado >= 1.75f*_84u89yj->_d31k9qm)*1.0f + (_1vsmado < 1.75f*_84u89yj->_d31k9qm)*0.0f;
                    _jg2y8rg = _f7fr91i[_o7vbgr8][0];
                    _fu7bfxq = 0;
                    _ybxkpm2 = 1u;
                    while (_ybxkpm2 == 1u)
                    {
                        _1b07g2r = 0u;
                        while (_fu7bfxq < _jg2y8rg-1)
                        { _be3yqy5 = _f7fr91i[_o7vbgr8][16 -_jg2y8rg +_fu7bfxq];
                            __builtin___memcpy_chk (_dan085r, &_pqvl00e[_be3yqy5][0], 3*sizeof(float), __builtin_object_size (_dan085r, 0));
                            SubtractVect3(_2ljk4q5, _uug95rl, _dan085r); _1vsmado = VectorLength(_2ljk4q5);
                            _1b07g2r = (_kb700uc[_be3yqy5] == 0u)*1u + (_kb700uc[_be3yqy5] != 0u)*0u;
                            _1b07g2r = (_jzr25md[_be3yqy5][0] > 1u)*_1b07g2r + (_jzr25md[_be3yqy5][0] <= 1u)*0u;
                            _1b07g2r = (_1vsmado > (_pqvl00e[_be3yqy5][5] +_84u89yj->_ktakpux*_84u89yj->_d31k9qm))*_1b07g2r + (_1vsmado <= (_pqvl00e[_be3yqy5][5] +_84u89yj->_ktakpux*_84u89yj->_d31k9qm))*0u;
                            if (_1b07g2r == 1u) { break; }
                            _fu7bfxq += 1u;
                        }
                        if (_1b07g2r == 0u)
                        { __builtin___memcpy_chk (_e3nwxuj, &_n8la100[_egok01t[_o7vbgr8*8 +0]*4 +0], 3*sizeof(float), __builtin_object_size (_e3nwxuj, 0));
                            __builtin___memcpy_chk (_3ac6pe6, &_n8la100[_egok01t[_o7vbgr8*8 +1]*4 +0], 3*sizeof(float), __builtin_object_size (_3ac6pe6, 0));
                            __builtin___memcpy_chk (_3lbqgu3, &_n8la100[_egok01t[_o7vbgr8*8 +2]*4 +0], 3*sizeof(float), __builtin_object_size (_3lbqgu3, 0));
                            if (_84u89yj->_fhv1xg1 == 1u)
                            { StrainHS_Nikkhoo(_mt8vt43, _dhgzz0j, _uug95rl[0], _uug95rl[1], _uug95rl[2], _e3nwxuj, _3ac6pe6, _3lbqgu3, _84u89yj->_31yrujo, 0.0f, 0.0f, _j3trp0j->_x5qjicj, _j3trp0j->_plozjcn);
                                StrainHS_Nikkhoo(_pkxvhag, _dhgzz0j, _uug95rl[0], _uug95rl[1], _uug95rl[2], _e3nwxuj, _3ac6pe6, _3lbqgu3, 0.0f, _84u89yj->_31yrujo, 0.0f, _j3trp0j->_x5qjicj, _j3trp0j->_plozjcn);
                           }
                            else
                            { StrainFS_Nikkhoo(_mt8vt43, _dhgzz0j, _uug95rl[0], _uug95rl[1], _uug95rl[2], _e3nwxuj, _3ac6pe6, _3lbqgu3, _84u89yj->_31yrujo, 0.0f, 0.0f, _j3trp0j->_x5qjicj, _j3trp0j->_plozjcn);
                                StrainFS_Nikkhoo(_pkxvhag, _dhgzz0j, _uug95rl[0], _uug95rl[1], _uug95rl[2], _e3nwxuj, _3ac6pe6, _3lbqgu3, 0.0f, _84u89yj->_31yrujo, 0.0f, _j3trp0j->_x5qjicj, _j3trp0j->_plozjcn);
                            }
                            ScaleStress6(_mt8vt43, (_8e7flqh/(_84u89yj->_31yrujo*_3n3bzt9[_o7vbgr8*16 +15]))); RotStressG2L(&_2292mtr[0], _v71448n, _mt8vt43);
                            ScaleStress6(_pkxvhag, (_8e7flqh/(_84u89yj->_31yrujo*_3n3bzt9[_o7vbgr8*16 +15]))); RotStressG2L(&_2292mtr[3], _v71448n, _pkxvhag);
                            _ny181ac[_vshkoec[i]*2 +0] = _2292mtr[0]; _ny181ac[_vshkoec[i]*2 +1] = _2292mtr[3];
                            _jdtvv2y[_vshkoec[i]*2 +0] = _2292mtr[1]; _jdtvv2y[_vshkoec[i]*2 +1] = _2292mtr[4];
                            _b68rerv[_vshkoec[i]*2 +0] = _2292mtr[2]; _b68rerv[_vshkoec[i]*2 +1] = _2292mtr[5];
                            (*_eurxhw0)[i][_o7vbgr8] = _vshkoec[i];
                            _vshkoec[i] += 1u;
                            _rgb93wk[_o7vbgr8] = 1u;
                            _xdrkz40 = 0u;
                            for (j = 0u; j < _s3kkv9r->_otn9uze; j++) { _xdrkz40 += _rgb93wk[j]; }
                            for (j = 0u; j < _jg2y8rg; j++) { _kb700uc[ _f7fr91i[_o7vbgr8][16 -_jg2y8rg +j] ] = 1u; }
                            _ybxkpm2 = 0u;
                        }
                        else
                        { _be3yqy5 = _f7fr91i[_o7vbgr8][16 -_jg2y8rg +_fu7bfxq];
                            __builtin___memcpy_chk (_e3nwxuj, &_pqvl00e[_be3yqy5][9], 3*sizeof(float), __builtin_object_size (_e3nwxuj, 0));
                            __builtin___memcpy_chk (_3ac6pe6, &_pqvl00e[_be3yqy5][12], 3*sizeof(float), __builtin_object_size (_3ac6pe6, 0));
                            __builtin___memcpy_chk (_3lbqgu3, &_pqvl00e[_be3yqy5][15], 3*sizeof(float), __builtin_object_size (_3lbqgu3, 0));
                            __builtin___memcpy_chk (_obyghgk, &_pqvl00e[_be3yqy5][18], 3*sizeof(float), __builtin_object_size (_obyghgk, 0));
                            if (_84u89yj->_fhv1xg1 == 1u)
                            { StrainHS_Nikkhoo(_pt5r5rc, _dhgzz0j, _uug95rl[0], _uug95rl[1], _uug95rl[2], _e3nwxuj, _3ac6pe6, _obyghgk, _84u89yj->_31yrujo, 0.0f, 0.0f, _j3trp0j->_x5qjicj, _j3trp0j->_plozjcn);
                                StrainHS_Nikkhoo(_kc4lorx, _dhgzz0j, _uug95rl[0], _uug95rl[1], _uug95rl[2], _e3nwxuj, _3ac6pe6, _obyghgk, 0.0f, _84u89yj->_31yrujo, 0.0f, _j3trp0j->_x5qjicj, _j3trp0j->_plozjcn);
                                StrainHS_Nikkhoo(_mt8vt43, _dhgzz0j, _uug95rl[0], _uug95rl[1], _uug95rl[2], _3ac6pe6, _3lbqgu3, _obyghgk, _84u89yj->_31yrujo, 0.0f, 0.0f, _j3trp0j->_x5qjicj, _j3trp0j->_plozjcn);
                                StrainHS_Nikkhoo(_pkxvhag, _dhgzz0j, _uug95rl[0], _uug95rl[1], _uug95rl[2], _3ac6pe6, _3lbqgu3, _obyghgk, 0.0f, _84u89yj->_31yrujo, 0.0f, _j3trp0j->_x5qjicj, _j3trp0j->_plozjcn);
                                _mt8vt43[0] += _pt5r5rc[0]; _mt8vt43[1] += _pt5r5rc[1]; _mt8vt43[2] += _pt5r5rc[2];
                                _mt8vt43[3] += _pt5r5rc[3]; _mt8vt43[4] += _pt5r5rc[4]; _mt8vt43[5] += _pt5r5rc[5];
                                _pkxvhag[0] += _kc4lorx[0]; _pkxvhag[1] += _kc4lorx[1]; _pkxvhag[2] += _kc4lorx[2];
                                _pkxvhag[3] += _kc4lorx[3]; _pkxvhag[4] += _kc4lorx[4]; _pkxvhag[5] += _kc4lorx[5];
                            }
                            else
                            { StrainFS_Nikkhoo(_pt5r5rc, _dhgzz0j, _uug95rl[0], _uug95rl[1], _uug95rl[2], _e3nwxuj, _3ac6pe6, _obyghgk, _84u89yj->_31yrujo, 0.0f, 0.0f, _j3trp0j->_x5qjicj, _j3trp0j->_plozjcn);
                                StrainFS_Nikkhoo(_kc4lorx, _dhgzz0j, _uug95rl[0], _uug95rl[1], _uug95rl[2], _e3nwxuj, _3ac6pe6, _obyghgk, 0.0f, _84u89yj->_31yrujo, 0.0f, _j3trp0j->_x5qjicj, _j3trp0j->_plozjcn);
                                StrainFS_Nikkhoo(_mt8vt43, _dhgzz0j, _uug95rl[0], _uug95rl[1], _uug95rl[2], _3ac6pe6, _3lbqgu3, _obyghgk, _84u89yj->_31yrujo, 0.0f, 0.0f, _j3trp0j->_x5qjicj, _j3trp0j->_plozjcn);
                                StrainFS_Nikkhoo(_pkxvhag, _dhgzz0j, _uug95rl[0], _uug95rl[1], _uug95rl[2], _3ac6pe6, _3lbqgu3, _obyghgk, 0.0f, _84u89yj->_31yrujo, 0.0f, _j3trp0j->_x5qjicj, _j3trp0j->_plozjcn);
                                _mt8vt43[0] += _pt5r5rc[0]; _mt8vt43[1] += _pt5r5rc[1]; _mt8vt43[2] += _pt5r5rc[2];
                                _mt8vt43[3] += _pt5r5rc[3]; _mt8vt43[4] += _pt5r5rc[4]; _mt8vt43[5] += _pt5r5rc[5];
                                _pkxvhag[0] += _kc4lorx[0]; _pkxvhag[1] += _kc4lorx[1]; _pkxvhag[2] += _kc4lorx[2];
                                _pkxvhag[3] += _kc4lorx[3]; _pkxvhag[4] += _kc4lorx[4]; _pkxvhag[5] += _kc4lorx[5];
                            }
                            ScaleStress6(_mt8vt43, (_8e7flqh/(_84u89yj->_31yrujo*_pqvl00e[_be3yqy5][3]))); RotStressG2L(&_2292mtr[0], _v71448n, _mt8vt43);
                            ScaleStress6(_pkxvhag, (_8e7flqh/(_84u89yj->_31yrujo*_pqvl00e[_be3yqy5][3]))); RotStressG2L(&_2292mtr[3], _v71448n, _pkxvhag);
                            _1b07g2r = 1u;
                            for (j = 0; j < _jzr25md[_be3yqy5][0]; j++)
                            { _iua08ek = _jzr25md[_be3yqy5][(j+1)];
                                __builtin___memcpy_chk (_e3nwxuj, &_pqvl00e[_iua08ek][9], 3*sizeof(float), __builtin_object_size (_e3nwxuj, 0));
                                __builtin___memcpy_chk (_3ac6pe6, &_pqvl00e[_iua08ek][12], 3*sizeof(float), __builtin_object_size (_3ac6pe6, 0));
                                __builtin___memcpy_chk (_3lbqgu3, &_pqvl00e[_iua08ek][15], 3*sizeof(float), __builtin_object_size (_3lbqgu3, 0));
                                __builtin___memcpy_chk (_obyghgk, &_pqvl00e[_iua08ek][18], 3*sizeof(float), __builtin_object_size (_obyghgk, 0));
                                if (_84u89yj->_fhv1xg1 == 1u)
                                { StrainHS_Nikkhoo(_pt5r5rc, _dhgzz0j, _uug95rl[0], _uug95rl[1], _uug95rl[2], _e3nwxuj, _3ac6pe6, _obyghgk, _84u89yj->_31yrujo, 0.0f, 0.0f, _j3trp0j->_x5qjicj, _j3trp0j->_plozjcn);
                                    StrainHS_Nikkhoo(_kc4lorx, _dhgzz0j, _uug95rl[0], _uug95rl[1], _uug95rl[2], _e3nwxuj, _3ac6pe6, _obyghgk, 0.0f, _84u89yj->_31yrujo, 0.0f, _j3trp0j->_x5qjicj, _j3trp0j->_plozjcn);
                                    StrainHS_Nikkhoo(_mt8vt43, _dhgzz0j, _uug95rl[0], _uug95rl[1], _uug95rl[2], _3ac6pe6, _3lbqgu3, _obyghgk, _84u89yj->_31yrujo, 0.0f, 0.0f, _j3trp0j->_x5qjicj, _j3trp0j->_plozjcn);
                                    StrainHS_Nikkhoo(_pkxvhag, _dhgzz0j, _uug95rl[0], _uug95rl[1], _uug95rl[2], _3ac6pe6, _3lbqgu3, _obyghgk, 0.0f, _84u89yj->_31yrujo, 0.0f, _j3trp0j->_x5qjicj, _j3trp0j->_plozjcn);
                                    _mt8vt43[0] += _pt5r5rc[0]; _mt8vt43[1] += _pt5r5rc[1]; _mt8vt43[2] += _pt5r5rc[2];
                                    _mt8vt43[3] += _pt5r5rc[3]; _mt8vt43[4] += _pt5r5rc[4]; _mt8vt43[5] += _pt5r5rc[5];
                                    _pkxvhag[0] += _kc4lorx[0]; _pkxvhag[1] += _kc4lorx[1]; _pkxvhag[2] += _kc4lorx[2];
                                    _pkxvhag[3] += _kc4lorx[3]; _pkxvhag[4] += _kc4lorx[4]; _pkxvhag[5] += _kc4lorx[5];
                                }
                                else
                                { StrainFS_Nikkhoo(_pt5r5rc, _dhgzz0j, _uug95rl[0], _uug95rl[1], _uug95rl[2], _e3nwxuj, _3ac6pe6, _obyghgk, _84u89yj->_31yrujo, 0.0f, 0.0f, _j3trp0j->_x5qjicj, _j3trp0j->_plozjcn);
                                    StrainFS_Nikkhoo(_kc4lorx, _dhgzz0j, _uug95rl[0], _uug95rl[1], _uug95rl[2], _e3nwxuj, _3ac6pe6, _obyghgk, 0.0f, _84u89yj->_31yrujo, 0.0f, _j3trp0j->_x5qjicj, _j3trp0j->_plozjcn);
                                    StrainFS_Nikkhoo(_mt8vt43, _dhgzz0j, _uug95rl[0], _uug95rl[1], _uug95rl[2], _3ac6pe6, _3lbqgu3, _obyghgk, _84u89yj->_31yrujo, 0.0f, 0.0f, _j3trp0j->_x5qjicj, _j3trp0j->_plozjcn);
                                    StrainFS_Nikkhoo(_pkxvhag, _dhgzz0j, _uug95rl[0], _uug95rl[1], _uug95rl[2], _3ac6pe6, _3lbqgu3, _obyghgk, 0.0f, _84u89yj->_31yrujo, 0.0f, _j3trp0j->_x5qjicj, _j3trp0j->_plozjcn);
                                    _mt8vt43[0] += _pt5r5rc[0]; _mt8vt43[1] += _pt5r5rc[1]; _mt8vt43[2] += _pt5r5rc[2];
                                    _mt8vt43[3] += _pt5r5rc[3]; _mt8vt43[4] += _pt5r5rc[4]; _mt8vt43[5] += _pt5r5rc[5];
                                    _pkxvhag[0] += _kc4lorx[0]; _pkxvhag[1] += _kc4lorx[1]; _pkxvhag[2] += _kc4lorx[2];
                                    _pkxvhag[3] += _kc4lorx[3]; _pkxvhag[4] += _kc4lorx[4]; _pkxvhag[5] += _kc4lorx[5];
                                }
                                ScaleStress6(_mt8vt43, (_8e7flqh/(_84u89yj->_31yrujo*_pqvl00e[_iua08ek][3]))); RotStressG2L(&_gwaj1g1[0], _v71448n, _mt8vt43);
                                ScaleStress6(_pkxvhag, (_8e7flqh/(_84u89yj->_31yrujo*_pqvl00e[_iua08ek][3]))); RotStressG2L(&_gwaj1g1[3], _v71448n, _pkxvhag);
                                __builtin___memcpy_chk (_dan085r, &_pqvl00e[_iua08ek][0], 3*sizeof(float), __builtin_object_size (_dan085r, 0));
                                SubtractVect3(_2ljk4q5, _uug95rl, _dan085r); _1vsmado = VectorLength(_2ljk4q5);
                                _fv177ud = (_1vsmado/_84u89yj->_d31k9qm - _84u89yj->_ktakpux)/(_84u89yj->_qphr62g - _84u89yj->_ktakpux);
                                _fv177ud = (_fv177ud >= 0.0f)*_fv177ud + (_fv177ud < 0.0f)*0.0f;
                                _fv177ud = (_fv177ud <= 1.0f)*_fv177ud + (_fv177ud > 1.0f)*1.0f;
                                _fv177ud = _84u89yj->_dswm31y*_fv177ud;
                                _66cnolj = sqrtf( (_gwaj1g1[0]-_2292mtr[0])*(_gwaj1g1[0]-_2292mtr[0]) + (_gwaj1g1[1]-_2292mtr[1])*(_gwaj1g1[1]-_2292mtr[1]) + (_gwaj1g1[2]-_2292mtr[2])*(_gwaj1g1[2]-_2292mtr[2]) + (_gwaj1g1[3]-_2292mtr[3])*(_gwaj1g1[3]-_2292mtr[3]) + (_gwaj1g1[4]-_2292mtr[4])*(_gwaj1g1[4]-_2292mtr[4]) + (_gwaj1g1[5]-_2292mtr[5])*(_gwaj1g1[5]-_2292mtr[5]) );
                                _q52wspo = sqrtf( _2292mtr[0] *_2292mtr[0] + _2292mtr[1] *_2292mtr[1] + _2292mtr[2] *_2292mtr[2] + _2292mtr[3] *_2292mtr[3] + _2292mtr[4] *_2292mtr[4] + _2292mtr[5] *_2292mtr[5] );
                                _1b07g2r = ((_66cnolj/_q52wspo) <= _fv177ud)*_1b07g2r + ((_66cnolj/_q52wspo) > _fv177ud)*0u;
                            }
                            if (_1b07g2r == 1u)
                            { _ny181ac[_vshkoec[i]*2 +0] = _2292mtr[0]; _ny181ac[_vshkoec[i]*2 +1] = _2292mtr[3];
                                _jdtvv2y[_vshkoec[i]*2 +0] = _2292mtr[1]; _jdtvv2y[_vshkoec[i]*2 +1] = _2292mtr[4];
                                _b68rerv[_vshkoec[i]*2 +0] = _2292mtr[2]; _b68rerv[_vshkoec[i]*2 +1] = _2292mtr[5];
                                _xdrkz40 = 0u;
                                for (j = 0u; j < _s3kkv9r->_otn9uze; j++)
                                { _rgb93wk[j] = (_f7fr91i[j][16 -_jg2y8rg +_fu7bfxq] == _be3yqy5)*1u + (_f7fr91i[j][16 -_jg2y8rg +_fu7bfxq] != _be3yqy5)*_rgb93wk[j];
                                    (*_eurxhw0)[i][j] = (_f7fr91i[j][16 -_jg2y8rg +_fu7bfxq] == _be3yqy5)*_vshkoec[i] + (_f7fr91i[j][16 -_jg2y8rg +_fu7bfxq] != _be3yqy5)*(*_eurxhw0)[i][j];
                                    _xdrkz40 += _rgb93wk[j];
                                }
                                _vshkoec[i] += 1u;
                                for (j = 0u; j < _jg2y8rg; j++) { _kb700uc[ _f7fr91i[_o7vbgr8][16 -_jg2y8rg +j] ] = 1u; }
                                _ybxkpm2 = 0u;
                            }
                            _fu7bfxq += 1u;
                } } }
                _050ui22 = GetByteSize(64, 2*_vshkoec[i], sizeof(float));
                (*_t997e4u)[i] = aligned_alloc(64, _050ui22); __builtin___memcpy_chk ((*_t997e4u)[i], _ny181ac, 2*_vshkoec[i]* sizeof(float), __builtin_object_size ((*_t997e4u)[i], 0));
                (*_d1vn8lf)[i] = aligned_alloc(64, _050ui22); __builtin___memcpy_chk ((*_d1vn8lf)[i], _jdtvv2y, 2*_vshkoec[i]* sizeof(float), __builtin_object_size ((*_d1vn8lf)[i], 0));
                (*_4xkss0k)[i] = aligned_alloc(64, _050ui22); __builtin___memcpy_chk ((*_4xkss0k)[i], _b68rerv, 2*_vshkoec[i]* sizeof(float), __builtin_object_size ((*_4xkss0k)[i], 0));
            }
        }
        free(_rgb93wk); free(_kb700uc); free(_ye0qzsa); free(_0s1f537);
        free(_0wfoy4d); free(_zareiya); free(_0dl5e11); free(_ny181ac); free(_jdtvv2y); free(_b68rerv);
    }
    MPI_Barrier( ((MPI_Comm)0x44000000) );
    clock_gettime(_CLOCK_REALTIME, &_ofn1r4k);
    _weehwjq = (_ofn1r4k.tv_sec - _zn4regb.tv_sec) + (_ofn1r4k.tv_nsec - _zn4regb.tv_nsec)/1.0E+9;
    if (_tk63n2q == 0) { fprintf(__stdoutp,"\nTime taken for Kh_matrix(in seconds): %lf\n", _weehwjq); }
    for (int i = 0; i < _pwh0mx1; i++) { free(_pqvl00e[i]); }
    for (int i = 0; i < _s3kkv9r->_otn9uze; i++) { free(_f7fr91i[i]); }
    for (int i = 0; i < _s3kkv9r->_otn9uze*4; i++) { free(_jzr25md[i]); }
    for (int i = 0; i < _kicwmpe; i++) { free(_xf7jxjn[i]); }
    for (int i = 0; i < _s3kkv9r->_y25m8kb; i++) { free(_wshmlfz[i]); }
    for (int i = 0; i < _s3kkv9r->_y25m8kb*4; i++) { free(_jlpefja[i]); }
    free(_pqvl00e); free(_xf7jxjn); free(_f7fr91i); free(_wshmlfz); free(_jzr25md); free(_jlpefja);
    unsigned int i, _jp61lfc[4], _r7vgrf1[4], _6zjq7w7[4], _e4ex6u8[4];
    __builtin___memset_chk (_jp61lfc, 0, 4*sizeof(unsigned int), __builtin_object_size (_jp61lfc, 0));
    __builtin___memset_chk (_r7vgrf1, 0, 4*sizeof(unsigned int), __builtin_object_size (_r7vgrf1, 0));
    __builtin___memset_chk (_6zjq7w7, 0, 4*sizeof(unsigned int), __builtin_object_size (_6zjq7w7, 0));
    __builtin___memset_chk (_e4ex6u8, 0, 4*sizeof(unsigned int), __builtin_object_size (_e4ex6u8, 0));
    for (i = 0u; i < _6zx1tv1[_tk63n2q]; i++)
    { _jp61lfc[0] += _hthsjr0[i]; _e4ex6u8[0] = (_e4ex6u8[0] >= _hthsjr0[i])*_e4ex6u8[0] + (_e4ex6u8[0] < _hthsjr0[i])*_hthsjr0[i];
    }
    if (_s3kkv9r->_y25m8kb > 0u)
    { for (i = 0u; i < _6zx1tv1[_tk63n2q]; i++)
        { _jp61lfc[1] += _pmlvak9[i]; _e4ex6u8[1] = (_e4ex6u8[1] >= _pmlvak9[i])*_e4ex6u8[1] + (_e4ex6u8[1] < _pmlvak9[i])*_pmlvak9[i];
        }
    }
    for (i = 0u; i < _emruyp0[_tk63n2q]; i++)
    { _jp61lfc[2] += _vk0ucc1[i]; _e4ex6u8[2] = (_e4ex6u8[2] >= _vk0ucc1[i])*_e4ex6u8[2] + (_e4ex6u8[2] < _vk0ucc1[i])*_vk0ucc1[i];
        _jp61lfc[3] += _vshkoec[i]; _e4ex6u8[3] = (_e4ex6u8[3] >= _vshkoec[i])*_e4ex6u8[3] + (_e4ex6u8[3] < _vshkoec[i])*_vshkoec[i];
    }
    MPI_Allreduce((void *) -1, _e4ex6u8, 4, ((MPI_Datatype)0x4c000406), (MPI_Op)(0x58000001), ((MPI_Comm)0x44000000));
    _s3kkv9r->_d7x48wz = ( ((_e4ex6u8[0]) > (_e4ex6u8[3])) *(_e4ex6u8[0]) + ((_e4ex6u8[0]) <= (_e4ex6u8[3])) *(_e4ex6u8[3]) );
    _s3kkv9r->_mygty3a = ( ((_e4ex6u8[1]) > (_e4ex6u8[2])) *(_e4ex6u8[1]) + ((_e4ex6u8[1]) <= (_e4ex6u8[2])) *(_e4ex6u8[2]) );
    _jp61lfc[0] = _jp61lfc[0]/_6zx1tv1[_tk63n2q];
    _jp61lfc[1] = _jp61lfc[1]/_6zx1tv1[_tk63n2q];
    _jp61lfc[2] = (_emruyp0[_tk63n2q] > 0u) ? (_jp61lfc[2]/_emruyp0[_tk63n2q]) : 0u;
    _jp61lfc[3] = (_emruyp0[_tk63n2q] > 0u) ? (_jp61lfc[3]/_emruyp0[_tk63n2q]) : 0u;
    __builtin___memcpy_chk (_r7vgrf1, _jp61lfc, 4*sizeof(unsigned int), __builtin_object_size (_r7vgrf1, 0));
    __builtin___memcpy_chk (_6zjq7w7, _jp61lfc, 4*sizeof(unsigned int), __builtin_object_size (_6zjq7w7, 0));
    MPI_Allreduce((void *) -1, _r7vgrf1, 4, ((MPI_Datatype)0x4c000406), (MPI_Op)(0x58000002), ((MPI_Comm)0x44000000));
    MPI_Allreduce((void *) -1, _6zjq7w7, 4, ((MPI_Datatype)0x4c000406), (MPI_Op)(0x58000001), ((MPI_Comm)0x44000000));
    MPI_Allreduce((void *) -1, _jp61lfc, 4, ((MPI_Datatype)0x4c000406), (MPI_Op)(0x58000003), ((MPI_Comm)0x44000000));
    _jp61lfc[0] /= _b6gs8g8; _jp61lfc[1] /= _b6gs8g8; _jp61lfc[2] /= _b6gs8g8; _jp61lfc[3] /= _b6gs8g8;
    if (_tk63n2q == 0)
    { fprintf(__stdoutp,"Kh_lengths (average and range) across RANKs\n");
        fprintf(__stdoutp,"Minimum:         %u\n", _r7vgrf1[0]);
        fprintf(__stdoutp,"Average:         %u\n", _jp61lfc[0]);
        fprintf(__stdoutp,"Maximum:         %u\n", _6zjq7w7[0]);
    }
    MPI_Barrier( ((MPI_Comm)0x44000000) );
    return;
}
void Save_Khmatrix(int _tk63n2q, const char *restrict _xqt4p50, const struct _gzui1b0 *restrict _s3kkv9r, const struct _rcbpljd *restrict _84u89yj, const int *restrict _npuevjl, const int *restrict _yh1fkj4, const int *restrict _6zx1tv1, const int *restrict _emruyp0, const unsigned int *restrict _5izolnn, const unsigned int *restrict _laym0wb, const float *restrict _vwvm4w3, const float *restrict _dntcagr, unsigned short int *restrict _hthsjr0, unsigned short int *restrict _pmlvak9, unsigned short int *restrict _vk0ucc1, unsigned short int *restrict _vshkoec, unsigned short int **restrict _k2wbxan, unsigned short int **restrict _n14eb87, unsigned short int **restrict _va1nkbt, unsigned short int **restrict _eurxhw0, float **restrict _0isb5iq, float **restrict _zzf8s2z, float **restrict _3jczucq, float **restrict _sh3tow4, float **restrict _jbc2sws, float **restrict _e9ftaxp, float **restrict _6fokzix, float **restrict _z39pecn, float **restrict _t997e4u, float **restrict _d1vn8lf, float **restrict _4xkss0k)
{
    MPI_Status _myx444n;
    MPI_Offset _r1yxr2d;
    MPI_File _f445cri;
    unsigned int i;
    unsigned long long _9lzs34j;
    unsigned short int *_02g3yu6 = ((void *)0), *_ui3do70 = ((void *)0), *_h4a3yjt = ((void *)0), *_u7fbmen = ((void *)0);
    unsigned long long *_bh6nqqb = ((void *)0), *_0rmhrcn = ((void *)0), *_sosx43w = ((void *)0), *_fsqg8x8 = ((void *)0);
    float *_8sryo4r = ((void *)0), *_6ub7l8v = ((void *)0);
    _02g3yu6 = malloc(_s3kkv9r->_otn9uze *sizeof(unsigned short int)); __builtin___memset_chk (_02g3yu6, 0, _s3kkv9r->_otn9uze *sizeof(unsigned short int), __builtin_object_size (_02g3yu6, 0));
    _bh6nqqb = malloc(_s3kkv9r->_otn9uze *sizeof(unsigned long long)); __builtin___memset_chk (_bh6nqqb, 0, _s3kkv9r->_otn9uze *sizeof(unsigned long long), __builtin_object_size (_bh6nqqb, 0));
    _8sryo4r = malloc((3*_s3kkv9r->_otn9uze) *sizeof(float)); __builtin___memset_chk (_8sryo4r, 0,(3*_s3kkv9r->_otn9uze) *sizeof(float), __builtin_object_size (_8sryo4r, 0));
    if (_s3kkv9r->_y25m8kb > 0u)
    { _ui3do70 = malloc(_s3kkv9r->_otn9uze *sizeof(unsigned short int)); __builtin___memset_chk (_ui3do70, 0, _s3kkv9r->_otn9uze *sizeof(unsigned short int), __builtin_object_size (_ui3do70, 0));
        _0rmhrcn = malloc(_s3kkv9r->_otn9uze *sizeof(unsigned long long)); __builtin___memset_chk (_0rmhrcn, 0, _s3kkv9r->_otn9uze *sizeof(unsigned long long), __builtin_object_size (_0rmhrcn, 0));
        _h4a3yjt = malloc(_s3kkv9r->_y25m8kb *sizeof(unsigned short int)); __builtin___memset_chk (_h4a3yjt, 0, _s3kkv9r->_y25m8kb *sizeof(unsigned short int), __builtin_object_size (_h4a3yjt, 0));
        _sosx43w = malloc(_s3kkv9r->_y25m8kb *sizeof(unsigned long long)); __builtin___memset_chk (_sosx43w, 0, _s3kkv9r->_y25m8kb *sizeof(unsigned long long), __builtin_object_size (_sosx43w, 0));
        _u7fbmen = malloc(_s3kkv9r->_y25m8kb *sizeof(unsigned short int)); __builtin___memset_chk (_u7fbmen, 0, _s3kkv9r->_y25m8kb *sizeof(unsigned short int), __builtin_object_size (_u7fbmen, 0));
        _fsqg8x8 = malloc(_s3kkv9r->_y25m8kb *sizeof(unsigned long long)); __builtin___memset_chk (_fsqg8x8, 0, _s3kkv9r->_y25m8kb *sizeof(unsigned long long), __builtin_object_size (_fsqg8x8, 0));
        _6ub7l8v = malloc((3*_s3kkv9r->_y25m8kb) *sizeof(float)); __builtin___memset_chk (_6ub7l8v, 0,(3*_s3kkv9r->_y25m8kb) *sizeof(float), __builtin_object_size (_6ub7l8v, 0));
    }
    for (i = 0u; i < _6zx1tv1[_tk63n2q]; i++) { _02g3yu6[_5izolnn[i]] = _hthsjr0[i]; }
    MPI_Allreduce((void *) -1, _02g3yu6, _s3kkv9r->_otn9uze, ((MPI_Datatype)0x4c000204), (MPI_Op)(0x58000003), ((MPI_Comm)0x44000000));
    for (i = 1u; i < _s3kkv9r->_otn9uze; i++)
    { _bh6nqqb[i] = _bh6nqqb[i-1] + (_s3kkv9r->_otn9uze*sizeof(unsigned short int) + _02g3yu6[i-1]*6*sizeof(float));
    }
    for (i = 0u; i < _6zx1tv1[_tk63n2q]; i++) { _8sryo4r[_5izolnn[i] +0*_s3kkv9r->_otn9uze] = _vwvm4w3[i*16 +5]; _8sryo4r[_5izolnn[i] +1*_s3kkv9r->_otn9uze] = _vwvm4w3[i*16 +6]; _8sryo4r[_5izolnn[i] +2*_s3kkv9r->_otn9uze] = _vwvm4w3[i*16 +7]; }
    MPI_Allreduce((void *) -1, _8sryo4r, (3*_s3kkv9r->_otn9uze), ((MPI_Datatype)0x4c00040a), (MPI_Op)(0x58000003), ((MPI_Comm)0x44000000));
    if (_s3kkv9r->_y25m8kb > 0u)
    { for (i = 0u; i < _6zx1tv1[_tk63n2q]; i++) { _ui3do70[_5izolnn[i]] = _pmlvak9[i]; }
        for (i = 0u; i < _emruyp0[_tk63n2q]; i++) { _h4a3yjt[_laym0wb[i]] = _vk0ucc1[i]; _u7fbmen[_laym0wb[i]] = _vshkoec[i]; }
        MPI_Allreduce((void *) -1, _ui3do70, _s3kkv9r->_otn9uze, ((MPI_Datatype)0x4c000204), (MPI_Op)(0x58000003), ((MPI_Comm)0x44000000));
        MPI_Allreduce((void *) -1, _h4a3yjt, _s3kkv9r->_y25m8kb, ((MPI_Datatype)0x4c000204), (MPI_Op)(0x58000003), ((MPI_Comm)0x44000000));
        MPI_Allreduce((void *) -1, _u7fbmen, _s3kkv9r->_y25m8kb, ((MPI_Datatype)0x4c000204), (MPI_Op)(0x58000003), ((MPI_Comm)0x44000000));
        for (i = 1u; i < _s3kkv9r->_otn9uze; i++)
        { _0rmhrcn[i] = _0rmhrcn[i-1] + (_s3kkv9r->_y25m8kb*sizeof(unsigned short int) + _ui3do70[i-1]*6*sizeof(float));
        }
        for (i = 1u; i < _s3kkv9r->_y25m8kb; i++)
        { _sosx43w[i] = _sosx43w[i-1] + (_s3kkv9r->_y25m8kb*sizeof(unsigned short int) + _h4a3yjt[i-1]*9*sizeof(float));
            _fsqg8x8[i] = _fsqg8x8[i-1] + (_s3kkv9r->_otn9uze*sizeof(unsigned short int) + _u7fbmen[i-1]*6*sizeof(float));
        }
        for (i = 0u; i < _emruyp0[_tk63n2q]; i++) { _6ub7l8v[_laym0wb[i] +0*_s3kkv9r->_y25m8kb] = _dntcagr[i*16 +5]; _6ub7l8v[_laym0wb[i] +1*_s3kkv9r->_y25m8kb] = _dntcagr[i*16 +6]; _6ub7l8v[_laym0wb[i] +2*_s3kkv9r->_y25m8kb] = _dntcagr[i*16 +7]; }
        MPI_Allreduce((void *) -1, _6ub7l8v, (3*_s3kkv9r->_y25m8kb), ((MPI_Datatype)0x4c00040a), (MPI_Op)(0x58000003), ((MPI_Comm)0x44000000));
    }
    MPI_File_delete(_xqt4p50, ((MPI_Info)0x1c000000));
    MPI_File_open(((MPI_Comm)0x44000000), _xqt4p50, 1|4, ((MPI_Info)0x1c000000), &_f445cri);
    if (_tk63n2q == 0)
    { MPI_File_write(_f445cri, &_s3kkv9r->_otn9uze, 1, ((MPI_Datatype)0x4c000406), &_myx444n); MPI_File_write(_f445cri, &_84u89yj->_cs4kek4, 1, ((MPI_Datatype)0x4c000406), &_myx444n);
        MPI_File_write(_f445cri, &_s3kkv9r->_y25m8kb, 1, ((MPI_Datatype)0x4c000406), &_myx444n); MPI_File_write(_f445cri, &_84u89yj->_7yvdqbs, 1, ((MPI_Datatype)0x4c000406), &_myx444n);
    }
    _r1yxr2d = (4 * sizeof(unsigned int));
    MPI_File_write_at(_f445cri, (_r1yxr2d + _npuevjl[_tk63n2q]*sizeof(unsigned short int)), &_02g3yu6[_npuevjl[_tk63n2q]], _6zx1tv1[_tk63n2q], ((MPI_Datatype)0x4c000204), &_myx444n); _r1yxr2d += (_s3kkv9r->_otn9uze * sizeof(unsigned short int));
    if (_s3kkv9r->_y25m8kb > 0u)
    { MPI_File_write_at(_f445cri, (_r1yxr2d + _npuevjl[_tk63n2q]*sizeof(unsigned short int)), &_ui3do70[_npuevjl[_tk63n2q]], _6zx1tv1[_tk63n2q], ((MPI_Datatype)0x4c000204), &_myx444n); _r1yxr2d += (_s3kkv9r->_otn9uze * sizeof(unsigned short int));
        MPI_File_write_at(_f445cri, (_r1yxr2d + _yh1fkj4[_tk63n2q]*sizeof(unsigned short int)), &_h4a3yjt[_yh1fkj4[_tk63n2q]], _emruyp0[_tk63n2q], ((MPI_Datatype)0x4c000204), &_myx444n); _r1yxr2d += (_s3kkv9r->_y25m8kb * sizeof(unsigned short int));
        MPI_File_write_at(_f445cri, (_r1yxr2d + _yh1fkj4[_tk63n2q]*sizeof(unsigned short int)), &_u7fbmen[_yh1fkj4[_tk63n2q]], _emruyp0[_tk63n2q], ((MPI_Datatype)0x4c000204), &_myx444n); _r1yxr2d += (_s3kkv9r->_y25m8kb * sizeof(unsigned short int));
    }
    MPI_File_write_at(_f445cri, (_r1yxr2d + _npuevjl[_tk63n2q]*sizeof(float)), &_8sryo4r[0*_s3kkv9r->_otn9uze +_npuevjl[_tk63n2q]], _6zx1tv1[_tk63n2q], ((MPI_Datatype)0x4c00040a), &_myx444n); _r1yxr2d += (_s3kkv9r->_otn9uze * sizeof(float));
    MPI_File_write_at(_f445cri, (_r1yxr2d + _npuevjl[_tk63n2q]*sizeof(float)), &_8sryo4r[1*_s3kkv9r->_otn9uze +_npuevjl[_tk63n2q]], _6zx1tv1[_tk63n2q], ((MPI_Datatype)0x4c00040a), &_myx444n); _r1yxr2d += (_s3kkv9r->_otn9uze * sizeof(float));
    MPI_File_write_at(_f445cri, (_r1yxr2d + _npuevjl[_tk63n2q]*sizeof(float)), &_8sryo4r[2*_s3kkv9r->_otn9uze +_npuevjl[_tk63n2q]], _6zx1tv1[_tk63n2q], ((MPI_Datatype)0x4c00040a), &_myx444n); _r1yxr2d += (_s3kkv9r->_otn9uze * sizeof(float));
    if (_s3kkv9r->_y25m8kb > 0u)
    { MPI_File_write_at(_f445cri, (_r1yxr2d + _yh1fkj4[_tk63n2q]*sizeof(float)), &_6ub7l8v[0*_s3kkv9r->_y25m8kb +_yh1fkj4[_tk63n2q]], _emruyp0[_tk63n2q], ((MPI_Datatype)0x4c00040a), &_myx444n); _r1yxr2d += (_s3kkv9r->_y25m8kb * sizeof(float));
        MPI_File_write_at(_f445cri, (_r1yxr2d + _yh1fkj4[_tk63n2q]*sizeof(float)), &_6ub7l8v[1*_s3kkv9r->_y25m8kb +_yh1fkj4[_tk63n2q]], _emruyp0[_tk63n2q], ((MPI_Datatype)0x4c00040a), &_myx444n); _r1yxr2d += (_s3kkv9r->_y25m8kb * sizeof(float));
        MPI_File_write_at(_f445cri, (_r1yxr2d + _yh1fkj4[_tk63n2q]*sizeof(float)), &_6ub7l8v[2*_s3kkv9r->_y25m8kb +_yh1fkj4[_tk63n2q]], _emruyp0[_tk63n2q], ((MPI_Datatype)0x4c00040a), &_myx444n); _r1yxr2d += (_s3kkv9r->_y25m8kb * sizeof(float));
    }
    for (i = 0u; i < _6zx1tv1[_tk63n2q]; i++)
    { _9lzs34j = 0;
        MPI_File_write_at(_f445cri, (_r1yxr2d +_bh6nqqb[_5izolnn[i]] +_9lzs34j), _k2wbxan[i], _s3kkv9r->_otn9uze, ((MPI_Datatype)0x4c000204), &_myx444n); _9lzs34j += _s3kkv9r->_otn9uze *sizeof(unsigned short int);
        MPI_File_write_at(_f445cri, (_r1yxr2d +_bh6nqqb[_5izolnn[i]] +_9lzs34j), _0isb5iq[i], 2*_hthsjr0[i], ((MPI_Datatype)0x4c00040a), &_myx444n); _9lzs34j += (_hthsjr0[i]*2*sizeof(float));
        MPI_File_write_at(_f445cri, (_r1yxr2d +_bh6nqqb[_5izolnn[i]] +_9lzs34j), _zzf8s2z[i], 2*_hthsjr0[i], ((MPI_Datatype)0x4c00040a), &_myx444n); _9lzs34j += (_hthsjr0[i]*2*sizeof(float));
        MPI_File_write_at(_f445cri, (_r1yxr2d +_bh6nqqb[_5izolnn[i]] +_9lzs34j), _3jczucq[i], 2*_hthsjr0[i], ((MPI_Datatype)0x4c00040a), &_myx444n); _9lzs34j += (_hthsjr0[i]*2*sizeof(float));
    }
    _r1yxr2d += (_bh6nqqb[_s3kkv9r->_otn9uze-1] + _s3kkv9r->_otn9uze*sizeof(unsigned short int) + 6*_02g3yu6[_s3kkv9r->_otn9uze-1]*sizeof(float) );
    if (_s3kkv9r->_y25m8kb > 0u)
    { for (i = 0u; i < _6zx1tv1[_tk63n2q]; i++)
        { _9lzs34j = 0;
            MPI_File_write_at(_f445cri, (_r1yxr2d +_0rmhrcn[_5izolnn[i]] +_9lzs34j), _n14eb87[i], _s3kkv9r->_y25m8kb, ((MPI_Datatype)0x4c000204), &_myx444n); _9lzs34j += _s3kkv9r->_y25m8kb *sizeof(unsigned short int);
            MPI_File_write_at(_f445cri, (_r1yxr2d +_0rmhrcn[_5izolnn[i]] +_9lzs34j), _sh3tow4[i], 3*_pmlvak9[i], ((MPI_Datatype)0x4c00040a), &_myx444n); _9lzs34j += (_pmlvak9[i]*3*sizeof(float));
            MPI_File_write_at(_f445cri, (_r1yxr2d +_0rmhrcn[_5izolnn[i]] +_9lzs34j), _jbc2sws[i], 3*_pmlvak9[i], ((MPI_Datatype)0x4c00040a), &_myx444n); _9lzs34j += (_pmlvak9[i]*3*sizeof(float));
        }
        _r1yxr2d += (_0rmhrcn[_s3kkv9r->_otn9uze-1] + _s3kkv9r->_y25m8kb*sizeof(unsigned short int) + 6*_ui3do70[_s3kkv9r->_otn9uze-1]*sizeof(float) );
        if (_s3kkv9r->_y25m8kb > 0)
        { for (i = 0u; i < _emruyp0[_tk63n2q]; i++)
            { _9lzs34j = 0;
                MPI_File_write_at(_f445cri, (_r1yxr2d +_sosx43w[_laym0wb[i]] +_9lzs34j), _va1nkbt[i], _s3kkv9r->_y25m8kb, ((MPI_Datatype)0x4c000204), &_myx444n); _9lzs34j += _s3kkv9r->_y25m8kb *sizeof(unsigned short int);
                MPI_File_write_at(_f445cri, (_r1yxr2d +_sosx43w[_laym0wb[i]] +_9lzs34j), _e9ftaxp[i], 3*_vk0ucc1[i], ((MPI_Datatype)0x4c00040a), &_myx444n); _9lzs34j += (_vk0ucc1[i]*3*sizeof(float));
                MPI_File_write_at(_f445cri, (_r1yxr2d +_sosx43w[_laym0wb[i]] +_9lzs34j), _6fokzix[i], 3*_vk0ucc1[i], ((MPI_Datatype)0x4c00040a), &_myx444n); _9lzs34j += (_vk0ucc1[i]*3*sizeof(float));
                MPI_File_write_at(_f445cri, (_r1yxr2d +_sosx43w[_laym0wb[i]] +_9lzs34j), _z39pecn[i], 3*_vk0ucc1[i], ((MPI_Datatype)0x4c00040a), &_myx444n); _9lzs34j += (_vk0ucc1[i]*3*sizeof(float));
            }
            _r1yxr2d += (_sosx43w[_s3kkv9r->_y25m8kb-1] + _s3kkv9r->_y25m8kb*sizeof(unsigned short int) + 9*_h4a3yjt[_s3kkv9r->_y25m8kb-1]*sizeof(float) );
            for (i = 0u; i < _emruyp0[_tk63n2q]; i++)
            { _9lzs34j = 0;
                MPI_File_write_at(_f445cri, (_r1yxr2d +_fsqg8x8[_laym0wb[i]] +_9lzs34j), _eurxhw0[i], _s3kkv9r->_otn9uze, ((MPI_Datatype)0x4c000204), &_myx444n); _9lzs34j += _s3kkv9r->_otn9uze *sizeof(unsigned short int);
                MPI_File_write_at(_f445cri, (_r1yxr2d +_fsqg8x8[_laym0wb[i]] +_9lzs34j), _t997e4u[i], 2*_vshkoec[i], ((MPI_Datatype)0x4c00040a), &_myx444n); _9lzs34j += (_vshkoec[i]*2*sizeof(float));
                MPI_File_write_at(_f445cri, (_r1yxr2d +_fsqg8x8[_laym0wb[i]] +_9lzs34j), _d1vn8lf[i], 2*_vshkoec[i], ((MPI_Datatype)0x4c00040a), &_myx444n); _9lzs34j += (_vshkoec[i]*2*sizeof(float));
                MPI_File_write_at(_f445cri, (_r1yxr2d +_fsqg8x8[_laym0wb[i]] +_9lzs34j), _4xkss0k[i], 2*_vshkoec[i], ((MPI_Datatype)0x4c00040a), &_myx444n); _9lzs34j += (_vshkoec[i]*2*sizeof(float));
            }
        }
    }
    MPI_Barrier( ((MPI_Comm)0x44000000) );
    MPI_File_close(&_f445cri);
    free(_02g3yu6); free(_ui3do70); free(_h4a3yjt); free(_u7fbmen); free(_bh6nqqb); free(_0rmhrcn); free(_sosx43w); free(_fsqg8x8); free(_8sryo4r); free(_6ub7l8v);
    return;
}
void Load_Khmatrix(int _b6gs8g8, int _tk63n2q, const char *restrict _xqt4p50, struct _43d5iqv *restrict _j3trp0j, struct _gzui1b0 *restrict _s3kkv9r, const struct _rcbpljd *restrict _84u89yj, const int *restrict _npuevjl, const int *restrict _yh1fkj4, const int *restrict _6zx1tv1, const int *restrict _emruyp0, const unsigned int *restrict _5izolnn, unsigned int *restrict _mgvu4is, float *restrict _vwvm4w3, float *restrict _wyx9trf, const unsigned int *restrict _laym0wb, float *restrict _dntcagr, unsigned short int *restrict _hthsjr0, unsigned short int *restrict _pmlvak9, unsigned short int *restrict _vk0ucc1, unsigned short int *restrict _vshkoec, unsigned short int ***restrict _k2wbxan, unsigned short int ***restrict _n14eb87, unsigned short int ***restrict _va1nkbt, unsigned short int ***restrict _eurxhw0, float ***restrict _0isb5iq, float ***restrict _zzf8s2z, float ***restrict _3jczucq, float ***restrict _sh3tow4, float ***restrict _jbc2sws, float ***restrict _e9ftaxp, float ***restrict _6fokzix, float ***restrict _z39pecn, float ***restrict _t997e4u, float ***restrict _d1vn8lf, float ***restrict _4xkss0k)
{
    unsigned int i;
    unsigned int _ryf0992, _9w4uov5, _aegvfq3, _xlul6lg;
    long int _050ui22;
    MPI_Status _myx444n;
    MPI_Offset _r1yxr2d;
    MPI_File _f445cri;
    MPI_File_open(((MPI_Comm)0x44000000), _xqt4p50, 2, ((MPI_Info)0x1c000000), &_f445cri);
    MPI_File_read(_f445cri, &_ryf0992, 1, ((MPI_Datatype)0x4c000406), &_myx444n); MPI_File_read(_f445cri, &_aegvfq3, 1, ((MPI_Datatype)0x4c000406), &_myx444n);
    MPI_File_read(_f445cri, &_9w4uov5, 1, ((MPI_Datatype)0x4c000406), &_myx444n); MPI_File_read(_f445cri, &_xlul6lg, 1, ((MPI_Datatype)0x4c000406), &_myx444n);
    if ((_ryf0992 != _s3kkv9r->_otn9uze) || (_9w4uov5 != _s3kkv9r->_y25m8kb) || (_aegvfq3 != _84u89yj->_cs4kek4) || (_xlul6lg != _84u89yj->_7yvdqbs))
    { perror("Kh_matrix file was openend, but fault/boundary element numbers not matching => exit\n"); exit(10); }
    _050ui22 = GetByteSize(64, _6zx1tv1[_tk63n2q], sizeof(unsigned short int *));
    *_k2wbxan = aligned_alloc(64, _050ui22);
    _050ui22 = GetByteSize(64, _6zx1tv1[_tk63n2q], sizeof(float *));
    *_0isb5iq = aligned_alloc(64, _050ui22);
    *_zzf8s2z = aligned_alloc(64, _050ui22);
    *_3jczucq = aligned_alloc(64, _050ui22);
    if (_s3kkv9r->_y25m8kb > 0u)
    { _050ui22 = GetByteSize(64, _6zx1tv1[_tk63n2q], sizeof(unsigned short int *));
        *_n14eb87 = aligned_alloc(64, _050ui22);
        _050ui22 = GetByteSize(64, _emruyp0[_tk63n2q], sizeof(unsigned short int *));
        *_va1nkbt = aligned_alloc(64, _050ui22);
        *_eurxhw0 = aligned_alloc(64, _050ui22);
        _050ui22 = GetByteSize(64, _6zx1tv1[_tk63n2q], sizeof(float *));
        *_sh3tow4 = aligned_alloc(64, _050ui22);
        *_jbc2sws = aligned_alloc(64, _050ui22);
        _050ui22 = GetByteSize(64, _emruyp0[_tk63n2q], sizeof(float *));
        *_e9ftaxp = aligned_alloc(64, _050ui22);
        *_6fokzix = aligned_alloc(64, _050ui22);
        *_z39pecn = aligned_alloc(64, _050ui22);
        *_t997e4u = aligned_alloc(64, _050ui22);
        *_d1vn8lf = aligned_alloc(64, _050ui22);
        *_4xkss0k = aligned_alloc(64, _050ui22);
    }
    _r1yxr2d = (4 * sizeof(unsigned int)) ;
    unsigned long long _9lzs34j;
    unsigned short int *_02g3yu6 = ((void *)0), *_ui3do70 = ((void *)0), *_h4a3yjt = ((void *)0), *_u7fbmen = ((void *)0);
    unsigned long long *_bh6nqqb = ((void *)0), *_0rmhrcn = ((void *)0), *_sosx43w = ((void *)0), *_fsqg8x8 = ((void *)0);
    float *_8sryo4r = ((void *)0), *_6ub7l8v = ((void *)0);
    _02g3yu6 = malloc(_s3kkv9r->_otn9uze *sizeof(unsigned short int)); __builtin___memset_chk (_02g3yu6, 0, _s3kkv9r->_otn9uze *sizeof(unsigned short int), __builtin_object_size (_02g3yu6, 0));
    _bh6nqqb = malloc(_s3kkv9r->_otn9uze *sizeof(unsigned long long)); __builtin___memset_chk (_bh6nqqb, 0, _s3kkv9r->_otn9uze *sizeof(unsigned long long), __builtin_object_size (_bh6nqqb, 0));
    _8sryo4r = malloc((3*_s3kkv9r->_otn9uze) *sizeof(float)); __builtin___memset_chk (_8sryo4r, 0,(3*_s3kkv9r->_otn9uze) *sizeof(float), __builtin_object_size (_8sryo4r, 0));
    if (_s3kkv9r->_y25m8kb > 0u)
    { _ui3do70 = malloc(_s3kkv9r->_otn9uze *sizeof(unsigned short int)); __builtin___memset_chk (_ui3do70, 0, _s3kkv9r->_otn9uze *sizeof(unsigned short int), __builtin_object_size (_ui3do70, 0));
        _0rmhrcn = malloc(_s3kkv9r->_otn9uze *sizeof(unsigned long long)); __builtin___memset_chk (_0rmhrcn, 0, _s3kkv9r->_otn9uze *sizeof(unsigned long long), __builtin_object_size (_0rmhrcn, 0));
        _h4a3yjt = malloc(_s3kkv9r->_y25m8kb *sizeof(unsigned short int)); __builtin___memset_chk (_h4a3yjt, 0, _s3kkv9r->_y25m8kb *sizeof(unsigned short int), __builtin_object_size (_h4a3yjt, 0));
        _sosx43w = malloc(_s3kkv9r->_y25m8kb *sizeof(unsigned long long)); __builtin___memset_chk (_sosx43w, 0, _s3kkv9r->_y25m8kb *sizeof(unsigned long long), __builtin_object_size (_sosx43w, 0));
        _u7fbmen = malloc(_s3kkv9r->_y25m8kb *sizeof(unsigned short int)); __builtin___memset_chk (_u7fbmen, 0, _s3kkv9r->_y25m8kb *sizeof(unsigned short int), __builtin_object_size (_u7fbmen, 0));
        _fsqg8x8 = malloc(_s3kkv9r->_y25m8kb *sizeof(unsigned long long)); __builtin___memset_chk (_fsqg8x8, 0, _s3kkv9r->_y25m8kb *sizeof(unsigned long long), __builtin_object_size (_fsqg8x8, 0));
        _6ub7l8v = malloc((3*_s3kkv9r->_y25m8kb) *sizeof(float)); __builtin___memset_chk (_6ub7l8v, 0,(3*_s3kkv9r->_y25m8kb) *sizeof(float), __builtin_object_size (_6ub7l8v, 0));
    }
    MPI_File_read_at(_f445cri, (_r1yxr2d + _npuevjl[_tk63n2q]*sizeof(unsigned short int)), &_02g3yu6[_npuevjl[_tk63n2q]], _6zx1tv1[_tk63n2q], ((MPI_Datatype)0x4c000204), &_myx444n); _r1yxr2d += (_s3kkv9r->_otn9uze * sizeof(unsigned short int));
    if (_s3kkv9r->_y25m8kb > 0u)
    { MPI_File_read_at(_f445cri, (_r1yxr2d + _npuevjl[_tk63n2q]*sizeof(unsigned short int)), &_ui3do70[_npuevjl[_tk63n2q]], _6zx1tv1[_tk63n2q], ((MPI_Datatype)0x4c000204), &_myx444n); _r1yxr2d += (_s3kkv9r->_otn9uze * sizeof(unsigned short int));
        MPI_File_read_at(_f445cri, (_r1yxr2d + _yh1fkj4[_tk63n2q]*sizeof(unsigned short int)), &_h4a3yjt[_yh1fkj4[_tk63n2q]], _emruyp0[_tk63n2q], ((MPI_Datatype)0x4c000204), &_myx444n); _r1yxr2d += (_s3kkv9r->_y25m8kb * sizeof(unsigned short int));
        MPI_File_read_at(_f445cri, (_r1yxr2d + _yh1fkj4[_tk63n2q]*sizeof(unsigned short int)), &_u7fbmen[_yh1fkj4[_tk63n2q]], _emruyp0[_tk63n2q], ((MPI_Datatype)0x4c000204), &_myx444n); _r1yxr2d += (_s3kkv9r->_y25m8kb * sizeof(unsigned short int));
    }
    MPI_Allreduce((void *) -1, _02g3yu6, _s3kkv9r->_otn9uze, ((MPI_Datatype)0x4c000204), (MPI_Op)(0x58000003), ((MPI_Comm)0x44000000));
    if (_s3kkv9r->_y25m8kb > 0u)
    { MPI_Allreduce((void *) -1, _ui3do70, _s3kkv9r->_otn9uze, ((MPI_Datatype)0x4c000204), (MPI_Op)(0x58000003), ((MPI_Comm)0x44000000));
        MPI_Allreduce((void *) -1, _h4a3yjt, _s3kkv9r->_y25m8kb, ((MPI_Datatype)0x4c000204), (MPI_Op)(0x58000003), ((MPI_Comm)0x44000000));
        MPI_Allreduce((void *) -1, _u7fbmen, _s3kkv9r->_y25m8kb, ((MPI_Datatype)0x4c000204), (MPI_Op)(0x58000003), ((MPI_Comm)0x44000000));
    }
    for (i = 1u; i < _s3kkv9r->_otn9uze; i++)
    { _bh6nqqb[i] = _bh6nqqb[i-1] + (_s3kkv9r->_otn9uze*sizeof(unsigned short int) + _02g3yu6[i-1]*6*sizeof(float));
    }
    if (_s3kkv9r->_y25m8kb > 0u)
    { for (i = 1u; i < _s3kkv9r->_otn9uze; i++)
        { _0rmhrcn[i] = _0rmhrcn[i-1] + (_s3kkv9r->_y25m8kb*sizeof(unsigned short int) + _ui3do70[i-1]*6*sizeof(float));
        }
        for (i = 1u; i < _s3kkv9r->_y25m8kb; i++)
        { _sosx43w[i] = _sosx43w[i-1] + (_s3kkv9r->_y25m8kb*sizeof(unsigned short int) + _h4a3yjt[i-1]*9*sizeof(float));
            _fsqg8x8[i] = _fsqg8x8[i-1] + (_s3kkv9r->_otn9uze*sizeof(unsigned short int) + _u7fbmen[i-1]*6*sizeof(float));
        }
    }
    MPI_File_read_at(_f445cri, (_r1yxr2d + _npuevjl[_tk63n2q]*sizeof(float)), &_8sryo4r[0*_s3kkv9r->_otn9uze +_npuevjl[_tk63n2q]], _6zx1tv1[_tk63n2q], ((MPI_Datatype)0x4c00040a), &_myx444n); _r1yxr2d += (_s3kkv9r->_otn9uze * sizeof(float));
    MPI_File_read_at(_f445cri, (_r1yxr2d + _npuevjl[_tk63n2q]*sizeof(float)), &_8sryo4r[1*_s3kkv9r->_otn9uze +_npuevjl[_tk63n2q]], _6zx1tv1[_tk63n2q], ((MPI_Datatype)0x4c00040a), &_myx444n); _r1yxr2d += (_s3kkv9r->_otn9uze * sizeof(float));
    MPI_File_read_at(_f445cri, (_r1yxr2d + _npuevjl[_tk63n2q]*sizeof(float)), &_8sryo4r[2*_s3kkv9r->_otn9uze +_npuevjl[_tk63n2q]], _6zx1tv1[_tk63n2q], ((MPI_Datatype)0x4c00040a), &_myx444n); _r1yxr2d += (_s3kkv9r->_otn9uze * sizeof(float));
    MPI_Allreduce((void *) -1, _8sryo4r, (3*_s3kkv9r->_otn9uze), ((MPI_Datatype)0x4c00040a), (MPI_Op)(0x58000003), ((MPI_Comm)0x44000000));
    if (_s3kkv9r->_y25m8kb > 0u)
    { MPI_File_read_at(_f445cri, (_r1yxr2d + _yh1fkj4[_tk63n2q]*sizeof(float)), &_6ub7l8v[0*_s3kkv9r->_y25m8kb +_yh1fkj4[_tk63n2q]], _emruyp0[_tk63n2q], ((MPI_Datatype)0x4c00040a), &_myx444n); _r1yxr2d += (_s3kkv9r->_y25m8kb * sizeof(float));
        MPI_File_read_at(_f445cri, (_r1yxr2d + _yh1fkj4[_tk63n2q]*sizeof(float)), &_6ub7l8v[1*_s3kkv9r->_y25m8kb +_yh1fkj4[_tk63n2q]], _emruyp0[_tk63n2q], ((MPI_Datatype)0x4c00040a), &_myx444n); _r1yxr2d += (_s3kkv9r->_y25m8kb * sizeof(float));
        MPI_File_read_at(_f445cri, (_r1yxr2d + _yh1fkj4[_tk63n2q]*sizeof(float)), &_6ub7l8v[2*_s3kkv9r->_y25m8kb +_yh1fkj4[_tk63n2q]], _emruyp0[_tk63n2q], ((MPI_Datatype)0x4c00040a), &_myx444n); _r1yxr2d += (_s3kkv9r->_y25m8kb * sizeof(float));
        MPI_Allreduce((void *) -1, _6ub7l8v, (3*_s3kkv9r->_y25m8kb), ((MPI_Datatype)0x4c00040a), (MPI_Op)(0x58000003), ((MPI_Comm)0x44000000));
    }
    for (i = 0u; i < _6zx1tv1[_tk63n2q]; i++)
    { _hthsjr0[i] = _02g3yu6[_5izolnn[i]];
        _vwvm4w3[i*16 +5] = _8sryo4r[_5izolnn[i] +0*_s3kkv9r->_otn9uze]; _vwvm4w3[i*16 +6] = _8sryo4r[_5izolnn[i] +1*_s3kkv9r->_otn9uze]; _vwvm4w3[i*16 +7] = _8sryo4r[_5izolnn[i] +2*_s3kkv9r->_otn9uze];
    }
    if (_s3kkv9r->_y25m8kb > 0u)
    { for (i = 0u; i < _6zx1tv1[_tk63n2q]; i++)
        { _pmlvak9[i] = _ui3do70[_5izolnn[i]]; }
        for (i = 0u; i < _emruyp0[_tk63n2q]; i++)
        { _vk0ucc1[i] = _h4a3yjt[_laym0wb[i]]; _vshkoec[i] = _u7fbmen[_laym0wb[i]];
            _dntcagr[i*16 +5] = _6ub7l8v[_laym0wb[i] +0*_s3kkv9r->_y25m8kb]; _dntcagr[i*16 +6] = _6ub7l8v[_laym0wb[i] +1*_s3kkv9r->_y25m8kb]; _dntcagr[i*16 +7] = _6ub7l8v[_laym0wb[i] +2*_s3kkv9r->_y25m8kb];
        }
    }
    for (i = 0u; i < _6zx1tv1[_tk63n2q]; i++)
    { _9lzs34j = 0;
        _050ui22 = GetByteSize(64, _s3kkv9r->_otn9uze, sizeof(unsigned short int));
        (*_k2wbxan)[i] = aligned_alloc(64, _050ui22);
        _050ui22 = GetByteSize(64, 2*_hthsjr0[i], sizeof(float));
        (*_0isb5iq)[i] = aligned_alloc(64, _050ui22);
        (*_zzf8s2z)[i] = aligned_alloc(64, _050ui22);
        (*_3jczucq)[i] = aligned_alloc(64, _050ui22);
        MPI_File_read_at(_f445cri, (_r1yxr2d +_bh6nqqb[_5izolnn[i]] +_9lzs34j), (*_k2wbxan)[i], _s3kkv9r->_otn9uze, ((MPI_Datatype)0x4c000204), &_myx444n); _9lzs34j += _s3kkv9r->_otn9uze *sizeof(unsigned short int);
        MPI_File_read_at(_f445cri, (_r1yxr2d +_bh6nqqb[_5izolnn[i]] +_9lzs34j), (*_0isb5iq)[i], 2*_hthsjr0[i], ((MPI_Datatype)0x4c00040a), &_myx444n); _9lzs34j += (_hthsjr0[i]*2*sizeof(float));
        MPI_File_read_at(_f445cri, (_r1yxr2d +_bh6nqqb[_5izolnn[i]] +_9lzs34j), (*_zzf8s2z)[i], 2*_hthsjr0[i], ((MPI_Datatype)0x4c00040a), &_myx444n); _9lzs34j += (_hthsjr0[i]*2*sizeof(float));
        MPI_File_read_at(_f445cri, (_r1yxr2d +_bh6nqqb[_5izolnn[i]] +_9lzs34j), (*_3jczucq)[i], 2*_hthsjr0[i], ((MPI_Datatype)0x4c00040a), &_myx444n); _9lzs34j += (_hthsjr0[i]*2*sizeof(float));
    }
    _r1yxr2d += (_bh6nqqb[_s3kkv9r->_otn9uze-1] + _s3kkv9r->_otn9uze*sizeof(unsigned short int) + 6*_02g3yu6[_s3kkv9r->_otn9uze-1]*sizeof(float) );
    if (_s3kkv9r->_y25m8kb > 0u)
    { for (i = 0u; i < _6zx1tv1[_tk63n2q]; i++)
        { _9lzs34j = 0;
            _050ui22 = GetByteSize(64, _s3kkv9r->_y25m8kb, sizeof(unsigned short int));
            (*_n14eb87)[i] = aligned_alloc(64, _050ui22);
            _050ui22 = GetByteSize(64, 3*_pmlvak9[i], sizeof(float));
            (*_sh3tow4)[i] = aligned_alloc(64, _050ui22);
            (*_jbc2sws)[i] = aligned_alloc(64, _050ui22);
            MPI_File_read_at(_f445cri, (_r1yxr2d +_0rmhrcn[_5izolnn[i]] +_9lzs34j), (*_n14eb87)[i], _s3kkv9r->_y25m8kb, ((MPI_Datatype)0x4c000204), &_myx444n); _9lzs34j += _s3kkv9r->_y25m8kb *sizeof(unsigned short int);
            MPI_File_read_at(_f445cri, (_r1yxr2d +_0rmhrcn[_5izolnn[i]] +_9lzs34j), (*_sh3tow4)[i], 3*_pmlvak9[i], ((MPI_Datatype)0x4c00040a), &_myx444n); _9lzs34j += (_pmlvak9[i]*3*sizeof(float));
            MPI_File_read_at(_f445cri, (_r1yxr2d +_0rmhrcn[_5izolnn[i]] +_9lzs34j), (*_jbc2sws)[i], 3*_pmlvak9[i], ((MPI_Datatype)0x4c00040a), &_myx444n); _9lzs34j += (_pmlvak9[i]*3*sizeof(float));
        }
        _r1yxr2d += (_0rmhrcn[_s3kkv9r->_otn9uze-1] + _s3kkv9r->_y25m8kb*sizeof(unsigned short int) + 6*_ui3do70[_s3kkv9r->_otn9uze-1]*sizeof(float) );
        for (i = 0u; i < _emruyp0[_tk63n2q]; i++)
        { _9lzs34j = 0;
            _050ui22 = GetByteSize(64, _s3kkv9r->_y25m8kb, sizeof(unsigned short int));
            (*_va1nkbt)[i] = aligned_alloc(64, _050ui22);
            _050ui22 = GetByteSize(64, 3*_vk0ucc1[i], sizeof(float));
            (*_e9ftaxp)[i] = aligned_alloc(64, _050ui22);
            (*_6fokzix)[i] = aligned_alloc(64, _050ui22);
            (*_z39pecn)[i] = aligned_alloc(64, _050ui22);
            MPI_File_read_at(_f445cri, (_r1yxr2d +_sosx43w[_laym0wb[i]] +_9lzs34j), (*_va1nkbt)[i], _s3kkv9r->_y25m8kb, ((MPI_Datatype)0x4c000204), &_myx444n); _9lzs34j += _s3kkv9r->_y25m8kb *sizeof(unsigned short int);
            MPI_File_read_at(_f445cri, (_r1yxr2d +_sosx43w[_laym0wb[i]] +_9lzs34j), (*_e9ftaxp)[i], 3*_vk0ucc1[i], ((MPI_Datatype)0x4c00040a), &_myx444n); _9lzs34j += (_vk0ucc1[i]*3*sizeof(float));
            MPI_File_read_at(_f445cri, (_r1yxr2d +_sosx43w[_laym0wb[i]] +_9lzs34j), (*_6fokzix)[i], 3*_vk0ucc1[i], ((MPI_Datatype)0x4c00040a), &_myx444n); _9lzs34j += (_vk0ucc1[i]*3*sizeof(float));
            MPI_File_read_at(_f445cri, (_r1yxr2d +_sosx43w[_laym0wb[i]] +_9lzs34j), (*_z39pecn)[i], 3*_vk0ucc1[i], ((MPI_Datatype)0x4c00040a), &_myx444n); _9lzs34j += (_vk0ucc1[i]*3*sizeof(float));
        }
        _r1yxr2d += (_sosx43w[_s3kkv9r->_y25m8kb-1] + _s3kkv9r->_y25m8kb*sizeof(unsigned short int) + 9*_h4a3yjt[_s3kkv9r->_y25m8kb-1]*sizeof(float) );
        for (i = 0u; i < _emruyp0[_tk63n2q]; i++)
        { _9lzs34j = 0;
            _050ui22 = GetByteSize(64, _s3kkv9r->_otn9uze, sizeof(unsigned short int));
            (*_eurxhw0)[i] = aligned_alloc(64, _050ui22);
            _050ui22 = GetByteSize(64, 2*_vshkoec[i], sizeof(float));
            (*_t997e4u)[i] = aligned_alloc(64, _050ui22);
            (*_d1vn8lf)[i] = aligned_alloc(64, _050ui22);
            (*_4xkss0k)[i] = aligned_alloc(64, _050ui22);
            MPI_File_read_at(_f445cri, (_r1yxr2d +_fsqg8x8[_laym0wb[i]] +_9lzs34j), (*_eurxhw0)[i], _s3kkv9r->_otn9uze, ((MPI_Datatype)0x4c000204), &_myx444n); _9lzs34j += _s3kkv9r->_otn9uze *sizeof(unsigned short int);
            MPI_File_read_at(_f445cri, (_r1yxr2d +_fsqg8x8[_laym0wb[i]] +_9lzs34j), (*_t997e4u)[i], 2*_vshkoec[i], ((MPI_Datatype)0x4c00040a), &_myx444n); _9lzs34j += (_vshkoec[i]*2*sizeof(float));
            MPI_File_read_at(_f445cri, (_r1yxr2d +_fsqg8x8[_laym0wb[i]] +_9lzs34j), (*_d1vn8lf)[i], 2*_vshkoec[i], ((MPI_Datatype)0x4c00040a), &_myx444n); _9lzs34j += (_vshkoec[i]*2*sizeof(float));
            MPI_File_read_at(_f445cri, (_r1yxr2d +_fsqg8x8[_laym0wb[i]] +_9lzs34j), (*_4xkss0k)[i], 2*_vshkoec[i], ((MPI_Datatype)0x4c00040a), &_myx444n); _9lzs34j += (_vshkoec[i]*2*sizeof(float));
        }
    }
    MPI_Barrier( ((MPI_Comm)0x44000000) );
    MPI_File_close(&_f445cri);
    free(_02g3yu6); free(_ui3do70); free(_h4a3yjt); free(_u7fbmen); free(_bh6nqqb); free(_0rmhrcn); free(_sosx43w); free(_fsqg8x8); free(_8sryo4r); free(_6ub7l8v);
    for (i = 0u; i < _6zx1tv1[_tk63n2q]; i++)
    {
        _wyx9trf[i*16 +1] = ( ((fabsf(_vwvm4w3[i*16 +5])) > (fabsf(_vwvm4w3[i*16 +6]))) *(fabsf(_vwvm4w3[i*16 +5])) + ((fabsf(_vwvm4w3[i*16 +5])) <= (fabsf(_vwvm4w3[i*16 +6]))) *(fabsf(_vwvm4w3[i*16 +6])) );
        _wyx9trf[i*16 +1] = ( ((_wyx9trf[i*16 +1]) > (fabsf(_vwvm4w3[i*16 +7]))) *(_wyx9trf[i*16 +1]) + ((_wyx9trf[i*16 +1]) <= (fabsf(_vwvm4w3[i*16 +7]))) *(fabsf(_vwvm4w3[i*16 +7])) );
        ReAssignElemVals(_j3trp0j, _s3kkv9r, i, _mgvu4is, _vwvm4w3, _wyx9trf, 1);
    }
    unsigned int _jp61lfc[4], _r7vgrf1[4], _6zjq7w7[4], _e4ex6u8[4];
    __builtin___memset_chk (_jp61lfc, 0, 4*sizeof(unsigned int), __builtin_object_size (_jp61lfc, 0));
    __builtin___memset_chk (_r7vgrf1, 0, 4*sizeof(unsigned int), __builtin_object_size (_r7vgrf1, 0));
    __builtin___memset_chk (_6zjq7w7, 0, 4*sizeof(unsigned int), __builtin_object_size (_6zjq7w7, 0));
    __builtin___memset_chk (_e4ex6u8, 0, 4*sizeof(unsigned int), __builtin_object_size (_e4ex6u8, 0));
    for (i = 0u; i < _6zx1tv1[_tk63n2q]; i++)
    { _jp61lfc[0] += _hthsjr0[i]; _e4ex6u8[0] = (_e4ex6u8[0] >= _hthsjr0[i])*_e4ex6u8[0] + (_e4ex6u8[0] < _hthsjr0[i])*_hthsjr0[i];
    }
    if (_s3kkv9r->_y25m8kb > 0u)
    { for (i = 0u; i < _6zx1tv1[_tk63n2q]; i++)
        { _jp61lfc[1] += _pmlvak9[i]; _e4ex6u8[1] = (_e4ex6u8[1] >= _pmlvak9[i])*_e4ex6u8[1] + (_e4ex6u8[1] < _pmlvak9[i])*_pmlvak9[i];
        }
        for (i = 0u; i < _emruyp0[_tk63n2q]; i++)
        { _jp61lfc[2] += _vk0ucc1[i]; _e4ex6u8[2] = (_e4ex6u8[2] >= _vk0ucc1[i])*_e4ex6u8[2] + (_e4ex6u8[2] < _vk0ucc1[i])*_vk0ucc1[i];
            _jp61lfc[3] += _vshkoec[i]; _e4ex6u8[3] = (_e4ex6u8[3] >= _vshkoec[i])*_e4ex6u8[3] + (_e4ex6u8[3] < _vshkoec[i])*_vshkoec[i];
        }
    }
    MPI_Allreduce((void *) -1, _e4ex6u8, 4, ((MPI_Datatype)0x4c000406), (MPI_Op)(0x58000001), ((MPI_Comm)0x44000000));
    _s3kkv9r->_d7x48wz = ( ((_e4ex6u8[0]) > (_e4ex6u8[3])) *(_e4ex6u8[0]) + ((_e4ex6u8[0]) <= (_e4ex6u8[3])) *(_e4ex6u8[3]) );
    _s3kkv9r->_mygty3a = ( ((_e4ex6u8[1]) > (_e4ex6u8[2])) *(_e4ex6u8[1]) + ((_e4ex6u8[1]) <= (_e4ex6u8[2])) *(_e4ex6u8[2]) );
    _jp61lfc[0] = _jp61lfc[0]/_6zx1tv1[_tk63n2q];
    _jp61lfc[1] = _jp61lfc[1]/_6zx1tv1[_tk63n2q];
    _jp61lfc[2] = (_emruyp0[_tk63n2q] > 0u) ? (_jp61lfc[2]/_emruyp0[_tk63n2q]) : 0u;
    _jp61lfc[3] = (_emruyp0[_tk63n2q] > 0u) ? (_jp61lfc[3]/_emruyp0[_tk63n2q]) : 0u;
    __builtin___memcpy_chk (_r7vgrf1, _jp61lfc, 4*sizeof(unsigned int), __builtin_object_size (_r7vgrf1, 0));
    __builtin___memcpy_chk (_6zjq7w7, _jp61lfc, 4*sizeof(unsigned int), __builtin_object_size (_6zjq7w7, 0));
    MPI_Allreduce((void *) -1, _r7vgrf1, 4, ((MPI_Datatype)0x4c000406), (MPI_Op)(0x58000002), ((MPI_Comm)0x44000000));
    MPI_Allreduce((void *) -1, _6zjq7w7, 4, ((MPI_Datatype)0x4c000406), (MPI_Op)(0x58000001), ((MPI_Comm)0x44000000));
    MPI_Allreduce((void *) -1, _jp61lfc, 4, ((MPI_Datatype)0x4c000406), (MPI_Op)(0x58000003), ((MPI_Comm)0x44000000));
    _jp61lfc[0] /= _b6gs8g8; _jp61lfc[1] /= _b6gs8g8; _jp61lfc[2] /= _b6gs8g8; _jp61lfc[3] /= _b6gs8g8;
    if (_tk63n2q == 0)
    { fprintf(__stdoutp,"Kh_lengths (average and range) across RANKs\n");
        fprintf(__stdoutp,"Minimum:         %u\n", _r7vgrf1[0]);
        fprintf(__stdoutp,"Average:         %u\n", _jp61lfc[0]);
        fprintf(__stdoutp,"Maximum:         %u\n", _6zjq7w7[0]);
    }
    return;
}
void _18hafey(float _kb4ehou, unsigned int _o31yhxs, const float *restrict _v78hu7f, unsigned int _mda9oek, const unsigned int *restrict _qw4c6lb, const float *restrict _nbhnal3, unsigned int ***restrict _l9ue421, unsigned int ***restrict _t4q9bc2, unsigned int *restrict _vfua9bq)
{
    float _gqtoamd = 3.0f;
    float _nwywrra = 0.666f;
    float _tatj9a3 = 1.0f/_nwywrra;
    long int _050ui22;
    _050ui22 = GetByteSize(64, _mda9oek, sizeof(unsigned int *));
    *_l9ue421 = aligned_alloc(64, _050ui22);
    _050ui22 = GetByteSize(64, 16, sizeof(unsigned int));
    for (int i = 0u; i < _mda9oek; i++)
    { (*_l9ue421)[i] = aligned_alloc(64, _050ui22); __builtin___memset_chk ((*_l9ue421)[i], 0, _050ui22, __builtin_object_size ((*_l9ue421)[i], 0));
    }
    _050ui22 = GetByteSize(64, _mda9oek*4, sizeof(unsigned int *));
    *_t4q9bc2 = aligned_alloc(64, _050ui22);
    _050ui22 = GetByteSize(64, 8, sizeof(unsigned int));
    for (int i = 0u; i < (_mda9oek*4); i++)
    { (*_t4q9bc2)[i] = aligned_alloc(64, _050ui22); __builtin___memset_chk ((*_t4q9bc2)[i], 0, _050ui22, __builtin_object_size ((*_t4q9bc2)[i], 0));
    }
    unsigned int i, j, k, _i4w6ozx, _7qqyzqi, _50b8y5d, _fftouhy, _v4zwl3z, _jg2y8rg;
    float _2zcfzw0, _y36ddfg;
    unsigned int _hx6zdcu[2][4];
    float _ljefy3z[3], _wrb7f6f[3], _jssuuqo[3], _h61iaoi[6];
    unsigned int *_lu6k7p9 = malloc(_o31yhxs *sizeof(unsigned int));
    unsigned int **_xd356wa = malloc(_o31yhxs *sizeof(unsigned int *));
    __builtin___memset_chk (_lu6k7p9, 0, _o31yhxs*sizeof(unsigned int), __builtin_object_size (_lu6k7p9, 0));
    for (i = 0u; i < _mda9oek; i++) { _lu6k7p9[_qw4c6lb[i*8 +3]] += 1u; }
    for (i = 0u; i < _o31yhxs; i++) { _xd356wa[i] = malloc(_lu6k7p9[i]*sizeof(unsigned int)); }
    __builtin___memset_chk (_lu6k7p9, 0, _o31yhxs*sizeof(unsigned int), __builtin_object_size (_lu6k7p9, 0));
    for (i = 0u; i < _mda9oek; i++)
    { _xd356wa[_qw4c6lb[i*8 +3]][_lu6k7p9[_qw4c6lb[i*8 +3]]] = i;
        _lu6k7p9[_qw4c6lb[i*8 +3]] += 1u;
    }
    for (i = 0u; i < _o31yhxs; i++)
    { float _0i9g5dz[3], _ukwqvr1[4] = {3.40282347e+38F, -3.40282347e+38F, 3.40282347e+38F, -3.40282347e+38F};
        unsigned int *_yt91qe3 = malloc(_lu6k7p9[i]* sizeof(unsigned int)); __builtin___memset_chk (_yt91qe3, 0, _lu6k7p9[i]* sizeof(unsigned int), __builtin_object_size (_yt91qe3, 0));
        unsigned int *_pj8t22c = malloc(_lu6k7p9[i]* sizeof(unsigned int)); __builtin___memset_chk (_pj8t22c, 0, _lu6k7p9[i]* sizeof(unsigned int), __builtin_object_size (_pj8t22c, 0));
        unsigned int *_jcmrsez = malloc(_lu6k7p9[i]* sizeof(unsigned int)); __builtin___memset_chk (_jcmrsez, 0, _lu6k7p9[i]* sizeof(unsigned int), __builtin_object_size (_jcmrsez, 0));
        unsigned int **_2h2p39y = malloc(_lu6k7p9[i] *sizeof(unsigned int *));
        float **_ldalmje = malloc(_lu6k7p9[i] *sizeof(float *));
        for (j = 0u; j < _lu6k7p9[i]; j++)
        { _050ui22 = GetByteSize(64, 16, sizeof(unsigned int));
            _2h2p39y[j] = aligned_alloc(64, _050ui22); __builtin___memset_chk (_2h2p39y[j], 0, _050ui22, __builtin_object_size (_2h2p39y[j], 0));
            _050ui22 = GetByteSize(64, 16, sizeof(float));
            _ldalmje[j] = aligned_alloc(64, _050ui22); __builtin___memset_chk (_ldalmje[j], 0, _050ui22, __builtin_object_size (_ldalmje[j], 0));
        }
        __builtin___memcpy_chk (_ljefy3z, &_v78hu7f[i*8 +3], 3*sizeof(float), __builtin_object_size (_ljefy3z, 0));
        GetStkAndDipVect(_ljefy3z, _wrb7f6f, _jssuuqo);
        _h61iaoi[0] = 3.40282347e+38F; _h61iaoi[1] = 0.0f; _h61iaoi[2] = -3.40282347e+38F;
        _h61iaoi[3] = 3.40282347e+38F; _h61iaoi[4] = 0.0f; _h61iaoi[5] = -3.40282347e+38F;
        for (j = 0; j < _lu6k7p9[i]; j++)
        { _0i9g5dz[0] = _nbhnal3[ _xd356wa[i][j]*16 +0] - _v78hu7f[i*8 +0];
            _0i9g5dz[1] = _nbhnal3[ _xd356wa[i][j]*16 +1] - _v78hu7f[i*8 +1];
            _0i9g5dz[2] = _nbhnal3[ _xd356wa[i][j]*16 +2] - _v78hu7f[i*8 +2];
            _ldalmje[j][8] = _ljefy3z[0]*_0i9g5dz[0] + _ljefy3z[1]*_0i9g5dz[1] + _ljefy3z[2]*_0i9g5dz[2];
            _ldalmje[j][9] = _wrb7f6f[0]*_0i9g5dz[0] + _wrb7f6f[1]*_0i9g5dz[1] + _wrb7f6f[2]*_0i9g5dz[2];
            _ldalmje[j][10]= _jssuuqo[0]*_0i9g5dz[0] + _jssuuqo[1]*_0i9g5dz[1] + _jssuuqo[2]*_0i9g5dz[2];
            _ukwqvr1[0] = ( ((_ukwqvr1[0]) < (_ldalmje[j][9])) *(_ukwqvr1[0]) + ((_ukwqvr1[0]) >= (_ldalmje[j][9])) *(_ldalmje[j][9]) ); _ukwqvr1[1] = ( ((_ukwqvr1[1]) > (_ldalmje[j][9])) *(_ukwqvr1[1]) + ((_ukwqvr1[1]) <= (_ldalmje[j][9])) *(_ldalmje[j][9]) );
            _ukwqvr1[2] = ( ((_ukwqvr1[2]) < (_ldalmje[j][10])) *(_ukwqvr1[2]) + ((_ukwqvr1[2]) >= (_ldalmje[j][10])) *(_ldalmje[j][10]) ); _ukwqvr1[3] = ( ((_ukwqvr1[3]) > (_ldalmje[j][10])) *(_ukwqvr1[3]) + ((_ukwqvr1[3]) <= (_ldalmje[j][10])) *(_ldalmje[j][10]) );
        }
        _h61iaoi[0] = ( ((_h61iaoi[0]) < (_ukwqvr1[0])) *(_h61iaoi[0]) + ((_h61iaoi[0]) >= (_ukwqvr1[0])) *(_ukwqvr1[0]) ); _h61iaoi[2] = ( ((_h61iaoi[2]) > (_ukwqvr1[1])) *(_h61iaoi[2]) + ((_h61iaoi[2]) <= (_ukwqvr1[1])) *(_ukwqvr1[1]) );
        _h61iaoi[3] = ( ((_h61iaoi[3]) < (_ukwqvr1[2])) *(_h61iaoi[3]) + ((_h61iaoi[3]) >= (_ukwqvr1[2])) *(_ukwqvr1[2]) ); _h61iaoi[5] = ( ((_h61iaoi[5]) > (_ukwqvr1[3])) *(_h61iaoi[5]) + ((_h61iaoi[5]) <= (_ukwqvr1[3])) *(_ukwqvr1[3]) );
        _2zcfzw0 = ( (((_h61iaoi[2] -_h61iaoi[0])) > ((_h61iaoi[5] -_h61iaoi[3]))) *((_h61iaoi[2] -_h61iaoi[0])) + (((_h61iaoi[2] -_h61iaoi[0])) <= ((_h61iaoi[5] -_h61iaoi[3]))) *((_h61iaoi[5] -_h61iaoi[3])) );
        _ldalmje[0][0] = _h61iaoi[0]; _ldalmje[0][1] = _h61iaoi[2];
        _ldalmje[0][2] = _h61iaoi[3]; _ldalmje[0][3] = _h61iaoi[5];
        _7qqyzqi = 0u;
        _v4zwl3z = 0u;
        _i4w6ozx = 1u;
        while (0.5f*_2zcfzw0 >= _gqtoamd*_kb4ehou)
        { _7qqyzqi += 1u;
            _2zcfzw0 *= 0.5f;
            _50b8y5d = 0u;
            for (j = 0u; j < _i4w6ozx; j++)
            { _h61iaoi[0] = _ldalmje[j][0]; _h61iaoi[2] = _ldalmje[j][1]; _h61iaoi[1] = 0.5f*(_h61iaoi[0] +_h61iaoi[2]);
                _h61iaoi[3] = _ldalmje[j][2]; _h61iaoi[5] = _ldalmje[j][3]; _h61iaoi[4] = 0.5f*(_h61iaoi[3] +_h61iaoi[5]);
                _fftouhy = 0u;
                __builtin___memset_chk (_jcmrsez, 0, _lu6k7p9[i]* sizeof(unsigned int), __builtin_object_size (_jcmrsez, 0));
                _hx6zdcu[0][0] = 0u; _hx6zdcu[0][1] = 0u; _hx6zdcu[0][2] = 0u; _hx6zdcu[0][3] = 0u;
                _hx6zdcu[1][0] = 0u; _hx6zdcu[1][1] = 0u; _hx6zdcu[1][2] = 0u; _hx6zdcu[1][3] = 0u;
                for (k = 0u; k < _lu6k7p9[i]; k++)
                { _jg2y8rg = (_2h2p39y[k][_7qqyzqi-1] == _yt91qe3[j])*1u + (_2h2p39y[k][_7qqyzqi-1] != _yt91qe3[j])*0u;
                    _jcmrsez[_fftouhy] = (_jg2y8rg == 1u)*k + (_jg2y8rg != 1u)*_jcmrsez[_fftouhy];
                    _fftouhy = (_jg2y8rg == 1u)*(_fftouhy+1u) + (_jg2y8rg != 1u)*_fftouhy;
                }
                _y36ddfg = (_h61iaoi[2]-_h61iaoi[0]+_kb4ehou)/(_h61iaoi[5]-_h61iaoi[3]+_kb4ehou);
                if (_y36ddfg >= _tatj9a3)
                { for (k = 0u; k < _fftouhy; k++)
                    { if (_ldalmje[_jcmrsez[k]][9] < _h61iaoi[1])
                        { if (_hx6zdcu[0][0] == 0u)
                            { _v4zwl3z += 1u;
                                _hx6zdcu[0][0] = 1u; _hx6zdcu[1][0] = _v4zwl3z;
                                _ldalmje[_50b8y5d][4] = _h61iaoi[0]; _ldalmje[_50b8y5d][5] = _h61iaoi[1];
                                _ldalmje[_50b8y5d][6] = _h61iaoi[3]; _ldalmje[_50b8y5d][7] = _h61iaoi[5];
                                _pj8t22c[_50b8y5d] = _v4zwl3z; _50b8y5d += 1u;
                                (*_t4q9bc2)[(_yt91qe3[j]+ *_vfua9bq)][0] += 1u;
                                (*_t4q9bc2)[(_yt91qe3[j]+ *_vfua9bq)][(*_t4q9bc2)[(_yt91qe3[j]+ *_vfua9bq)][0]] = (_v4zwl3z + *_vfua9bq);
                            }
                            _2h2p39y[_jcmrsez[k]][_7qqyzqi] = _hx6zdcu[1][0];
                        }
                        else if (_ldalmje[_jcmrsez[k]][9] >= _h61iaoi[1])
                        { if (_hx6zdcu[0][1] == 0u)
                            { _v4zwl3z += 1u;
                                _hx6zdcu[0][1] = 1u; _hx6zdcu[1][1] = _v4zwl3z;
                                _ldalmje[_50b8y5d][4] = _h61iaoi[1]; _ldalmje[_50b8y5d][5] = _h61iaoi[2];
                                _ldalmje[_50b8y5d][6] = _h61iaoi[3]; _ldalmje[_50b8y5d][7] = _h61iaoi[5];
                                _pj8t22c[_50b8y5d] = _v4zwl3z; _50b8y5d += 1u;
                                (*_t4q9bc2)[(_yt91qe3[j]+ *_vfua9bq)][0] += 1u;
                                (*_t4q9bc2)[(_yt91qe3[j]+ *_vfua9bq)][(*_t4q9bc2)[(_yt91qe3[j]+ *_vfua9bq)][0]] = (_v4zwl3z + *_vfua9bq);
                            }
                            _2h2p39y[_jcmrsez[k]][_7qqyzqi] = _hx6zdcu[1][1];
                } } }
                else if (_y36ddfg <= _nwywrra)
                { for (k = 0u; k < _fftouhy; k++)
                    { if (_ldalmje[_jcmrsez[k]][10] < _h61iaoi[4])
                        { if (_hx6zdcu[0][0] == 0u)
                            { _v4zwl3z += 1u;
                                _hx6zdcu[0][0] = 1u; _hx6zdcu[1][0] = _v4zwl3z;
                                _ldalmje[_50b8y5d][4] = _h61iaoi[0]; _ldalmje[_50b8y5d][5] = _h61iaoi[2];
                                _ldalmje[_50b8y5d][6] = _h61iaoi[3]; _ldalmje[_50b8y5d][7] = _h61iaoi[4];
                                _pj8t22c[_50b8y5d] = _v4zwl3z; _50b8y5d += 1u;
                                (*_t4q9bc2)[(_yt91qe3[j]+ *_vfua9bq)][0] += 1u;
                                (*_t4q9bc2)[(_yt91qe3[j]+ *_vfua9bq)][(*_t4q9bc2)[(_yt91qe3[j]+ *_vfua9bq)][0]] = (_v4zwl3z + *_vfua9bq);
                                                                   }
                            _2h2p39y[_jcmrsez[k]][_7qqyzqi] = _hx6zdcu[1][0];
                        }
                        else if (_ldalmje[_jcmrsez[k]][10] >= _h61iaoi[4])
                        { if (_hx6zdcu[0][1] == 0u)
                            { _v4zwl3z += 1u;
                                _hx6zdcu[0][1] = 1u; _hx6zdcu[1][1] = _v4zwl3z;
                                _ldalmje[_50b8y5d][4] = _h61iaoi[0]; _ldalmje[_50b8y5d][5] = _h61iaoi[2];
                                _ldalmje[_50b8y5d][6] = _h61iaoi[4]; _ldalmje[_50b8y5d][7] = _h61iaoi[5];
                                _pj8t22c[_50b8y5d] = _v4zwl3z; _50b8y5d += 1u;
                                (*_t4q9bc2)[(_yt91qe3[j]+ *_vfua9bq)][0] += 1u;
                                (*_t4q9bc2)[(_yt91qe3[j]+ *_vfua9bq)][(*_t4q9bc2)[(_yt91qe3[j]+ *_vfua9bq)][0]] = (_v4zwl3z + *_vfua9bq);
                            }
                            _2h2p39y[_jcmrsez[k]][_7qqyzqi] = _hx6zdcu[1][1];
                } } }
                else
                { for (k = 0u; k < _fftouhy; k++)
                    { if ( (_ldalmje[_jcmrsez[k]][9] >= _h61iaoi[1]) && (_ldalmje[_jcmrsez[k]][10] >= _h61iaoi[4]) )
                        { if (_hx6zdcu[0][0] == 0u)
                            { _v4zwl3z += 1u;
                                _hx6zdcu[0][0] = 1u; _hx6zdcu[1][0] = _v4zwl3z;
                                _ldalmje[_50b8y5d][4] = _h61iaoi[1]; _ldalmje[_50b8y5d][5] = _h61iaoi[2];
                                _ldalmje[_50b8y5d][6] = _h61iaoi[4]; _ldalmje[_50b8y5d][7] = _h61iaoi[5];
                                _pj8t22c[_50b8y5d] = _v4zwl3z; _50b8y5d += 1u;
                                (*_t4q9bc2)[(_yt91qe3[j]+ *_vfua9bq)][0] +=1u;
                                (*_t4q9bc2)[(_yt91qe3[j]+ *_vfua9bq)][(*_t4q9bc2)[(_yt91qe3[j]+ *_vfua9bq)][0]] = (_v4zwl3z + *_vfua9bq);
                            }
                            _2h2p39y[_jcmrsez[k]][_7qqyzqi] = _hx6zdcu[1][0];
                        }
                        else if ( (_ldalmje[_jcmrsez[k]][9] < _h61iaoi[1]) && (_ldalmje[_jcmrsez[k]][10] >= _h61iaoi[4]) )
                        { if (_hx6zdcu[0][1] == 0u)
                            { _v4zwl3z += 1u;
                                _hx6zdcu[0][1] = 1u; _hx6zdcu[1][1] = _v4zwl3z;
                                _ldalmje[_50b8y5d][4] = _h61iaoi[0]; _ldalmje[_50b8y5d][5] = _h61iaoi[1];
                                _ldalmje[_50b8y5d][6] = _h61iaoi[4]; _ldalmje[_50b8y5d][7] = _h61iaoi[5];
                                _pj8t22c[_50b8y5d] = _v4zwl3z; _50b8y5d += 1u;
                                (*_t4q9bc2)[(_yt91qe3[j]+ *_vfua9bq)][0] +=1u;
                                (*_t4q9bc2)[(_yt91qe3[j]+ *_vfua9bq)][(*_t4q9bc2)[(_yt91qe3[j]+ *_vfua9bq)][0]] = (_v4zwl3z + *_vfua9bq);
                            }
                            _2h2p39y[_jcmrsez[k]][_7qqyzqi] = _hx6zdcu[1][1];
                        }
                        else if ( (_ldalmje[_jcmrsez[k]][9] < _h61iaoi[1]) && (_ldalmje[_jcmrsez[k]][10] < _h61iaoi[4]) )
                        { if (_hx6zdcu[0][2] == 0u)
                            { _v4zwl3z += 1u;
                                _hx6zdcu[0][2] = 1u; _hx6zdcu[1][2] = _v4zwl3z;
                                _ldalmje[_50b8y5d][4] = _h61iaoi[0]; _ldalmje[_50b8y5d][5] = _h61iaoi[1];
                                _ldalmje[_50b8y5d][6] = _h61iaoi[3]; _ldalmje[_50b8y5d][7] = _h61iaoi[4];
                                _pj8t22c[_50b8y5d] = _v4zwl3z; _50b8y5d += 1u;
                                (*_t4q9bc2)[(_yt91qe3[j]+ *_vfua9bq)][0] +=1u;
                                (*_t4q9bc2)[(_yt91qe3[j]+ *_vfua9bq)][(*_t4q9bc2)[(_yt91qe3[j]+ *_vfua9bq)][0]] = (_v4zwl3z + *_vfua9bq);
                            }
                            _2h2p39y[_jcmrsez[k]][_7qqyzqi] = _hx6zdcu[1][2];
                        }
                        else if ( (_ldalmje[_jcmrsez[k]][9] >= _h61iaoi[1]) && (_ldalmje[_jcmrsez[k]][10] < _h61iaoi[4]) )
                        { if (_hx6zdcu[0][3] == 0u)
                            { _v4zwl3z += 1u;
                                _hx6zdcu[0][3] = 1u; _hx6zdcu[1][3] = _v4zwl3z;
                                _ldalmje[_50b8y5d][4] = _h61iaoi[1]; _ldalmje[_50b8y5d][5] = _h61iaoi[2];
                                _ldalmje[_50b8y5d][6] = _h61iaoi[3]; _ldalmje[_50b8y5d][7] = _h61iaoi[4];
                                _pj8t22c[_50b8y5d] = _v4zwl3z; _50b8y5d += 1u;
                                (*_t4q9bc2)[(_yt91qe3[j]+ *_vfua9bq)][0] +=1u;
                                (*_t4q9bc2)[(_yt91qe3[j]+ *_vfua9bq)][(*_t4q9bc2)[(_yt91qe3[j]+ *_vfua9bq)][0]] = (_v4zwl3z + *_vfua9bq);
                            }
                            _2h2p39y[_jcmrsez[k]][_7qqyzqi] = _hx6zdcu[1][3];
                } } }
            }
            __builtin___memcpy_chk (_yt91qe3, _pj8t22c, _50b8y5d* sizeof(unsigned int), __builtin_object_size (_yt91qe3, 0));
            for (j = 0u; j < _50b8y5d; j++) { __builtin___memcpy_chk (&_ldalmje[j][0], &_ldalmje[j][4], 4* sizeof(float), __builtin_object_size (&_ldalmje[j][0], 0)); }
            _i4w6ozx = _50b8y5d;
        }
        _v4zwl3z += 1u;
        _7qqyzqi += 2u;
        for (j = 0; j < _lu6k7p9[i]; j++ )
        { _2h2p39y[j][_7qqyzqi -1] = (_v4zwl3z +j);
            (*_l9ue421)[_xd356wa[i][j]][0] = _7qqyzqi;
            for (k = 0; k < _7qqyzqi; k++)
            { (*_l9ue421)[_xd356wa[i][j]][(16 -_7qqyzqi +k)] = (_2h2p39y[j][k] + *_vfua9bq);
        } }
        *_vfua9bq += (_v4zwl3z + _lu6k7p9[i]);
        for (j = 0; j < _lu6k7p9[i]; j++) { free(_2h2p39y[j]); free(_ldalmje[j]); }
        free(_yt91qe3); free(_pj8t22c); free(_jcmrsez); free(_2h2p39y); free(_ldalmje);
    }
    for (i = 0u; i < _o31yhxs; i++) { free(_xd356wa[i]); }
    free(_lu6k7p9); free(_xd356wa);
}
void _h75lxo6(unsigned int _mda9oek, const unsigned int *restrict _qw4c6lb, const float *restrict _nbhnal3, const float *restrict _r1dwwo8, unsigned int **restrict _2h2p39y, unsigned int _vfua9bq, float ***restrict _vlnknd6)
{
    long int _050ui22;
    _050ui22 = GetByteSize(64, _vfua9bq, sizeof(float *));
    *_vlnknd6 = aligned_alloc(64, _050ui22);
    _050ui22 = GetByteSize(64, 21, sizeof(float));
    for (int i = 0; i < _vfua9bq; i++)
    { (*_vlnknd6)[i] = aligned_alloc(64, _050ui22); __builtin___memset_chk ((*_vlnknd6)[i], 0, _050ui22, __builtin_object_size ((*_vlnknd6)[i], 0));
    }
    int i, j;
    unsigned int _jg2y8rg, _ybxkpm2;
    float _q52wspo, _fv177ud, _9p3jpxg, _h30153s, _flse45b, _1vsmado;
    float _e3nwxuj[3], _3ac6pe6[3], _3lbqgu3[3], _obyghgk[3], _dan085r[3], _ljefy3z[3], _wrb7f6f[3], _jssuuqo[3], _t26ys44[3], _5gdej31[3], _2ljk4q5[3];
    for (i = 0u; i < _mda9oek; i++)
    { _jg2y8rg= _2h2p39y[i][0];
        for (j = 0u; j < _jg2y8rg; j++)
        { _ybxkpm2 = _2h2p39y[i][16 -_jg2y8rg +j];
            (*_vlnknd6)[_ybxkpm2][0] += _nbhnal3[i*16 +0]; (*_vlnknd6)[_ybxkpm2][1] += _nbhnal3[i*16 +1]; (*_vlnknd6)[_ybxkpm2][2] += _nbhnal3[i*16 +2];
            (*_vlnknd6)[_ybxkpm2][3] += _nbhnal3[i*16 +15]; (*_vlnknd6)[_ybxkpm2][4] += 1.0f;
            (*_vlnknd6)[_ybxkpm2][6] += _nbhnal3[i*16 +6]; (*_vlnknd6)[_ybxkpm2][7] += _nbhnal3[i*16 +7]; (*_vlnknd6)[_ybxkpm2][8] += _nbhnal3[i*16 +8];
    } }
    for (i = 0u; i < _vfua9bq; i++)
    { (*_vlnknd6)[i][0] /= (*_vlnknd6)[i][4]; (*_vlnknd6)[i][1] /= (*_vlnknd6)[i][4]; (*_vlnknd6)[i][2] /= (*_vlnknd6)[i][4];
        NrmlzeVect( &(*_vlnknd6)[i][6] );
        (*_vlnknd6)[i][5] = -3.40282347e+38F;
        (*_vlnknd6)[i][9] = 3.40282347e+38F; (*_vlnknd6)[i][10] = -3.40282347e+38F;
        (*_vlnknd6)[i][11] = 3.40282347e+38F; (*_vlnknd6)[i][12] = -3.40282347e+38F;
    }
    for (i = 0u; i < _mda9oek; i++)
    { __builtin___memcpy_chk (_e3nwxuj, &_r1dwwo8[_qw4c6lb[i*8 +0]*4 +0], 3*sizeof(float), __builtin_object_size (_e3nwxuj, 0));
        __builtin___memcpy_chk (_3ac6pe6, &_r1dwwo8[_qw4c6lb[i*8 +1]*4 +0], 3*sizeof(float), __builtin_object_size (_3ac6pe6, 0));
        __builtin___memcpy_chk (_3lbqgu3, &_r1dwwo8[_qw4c6lb[i*8 +2]*4 +0], 3*sizeof(float), __builtin_object_size (_3lbqgu3, 0));
        _jg2y8rg = _2h2p39y[i][0];
        for (j = 0u; j < _jg2y8rg; j++)
        { _ybxkpm2 = _2h2p39y[i][16 -_jg2y8rg +j];
            __builtin___memcpy_chk (_dan085r, &(*_vlnknd6)[_ybxkpm2][0], 3*sizeof(float), __builtin_object_size (_dan085r, 0));
            __builtin___memcpy_chk (_ljefy3z, &(*_vlnknd6)[_ybxkpm2][6], 3*sizeof(float), __builtin_object_size (_ljefy3z, 0));
            GetStkAndDipVect(_ljefy3z, _wrb7f6f, _jssuuqo);
            SubtractVect3(_t26ys44, _e3nwxuj, _dan085r);
            _5gdej31[0] = _wrb7f6f[0]*_t26ys44[0] + _wrb7f6f[1]*_t26ys44[1] + _wrb7f6f[2]*_t26ys44[2];
            _5gdej31[1] = _jssuuqo[0]*_t26ys44[0] + _jssuuqo[1]*_t26ys44[1] + _jssuuqo[2]*_t26ys44[2];
            (*_vlnknd6)[_ybxkpm2][9] = ( (((*_vlnknd6)[_ybxkpm2][9]) < (_5gdej31[0])) *((*_vlnknd6)[_ybxkpm2][9]) + (((*_vlnknd6)[_ybxkpm2][9]) >= (_5gdej31[0])) *(_5gdej31[0]) );
            (*_vlnknd6)[_ybxkpm2][10] = ( (((*_vlnknd6)[_ybxkpm2][10]) > (_5gdej31[0])) *((*_vlnknd6)[_ybxkpm2][10]) + (((*_vlnknd6)[_ybxkpm2][10]) <= (_5gdej31[0])) *(_5gdej31[0]) );
            (*_vlnknd6)[_ybxkpm2][11] = ( (((*_vlnknd6)[_ybxkpm2][11]) < (_5gdej31[1])) *((*_vlnknd6)[_ybxkpm2][11]) + (((*_vlnknd6)[_ybxkpm2][11]) >= (_5gdej31[1])) *(_5gdej31[1]) );
            (*_vlnknd6)[_ybxkpm2][12] = ( (((*_vlnknd6)[_ybxkpm2][12]) > (_5gdej31[1])) *((*_vlnknd6)[_ybxkpm2][12]) + (((*_vlnknd6)[_ybxkpm2][12]) <= (_5gdej31[1])) *(_5gdej31[1]) );
            SubtractVect3(_t26ys44, _3ac6pe6, _dan085r);
            _5gdej31[0] = _wrb7f6f[0]*_t26ys44[0] + _wrb7f6f[1]*_t26ys44[1] + _wrb7f6f[2]*_t26ys44[2];
            _5gdej31[1] = _jssuuqo[0]*_t26ys44[0] + _jssuuqo[1]*_t26ys44[1] + _jssuuqo[2]*_t26ys44[2];
            (*_vlnknd6)[_ybxkpm2][9] = ( (((*_vlnknd6)[_ybxkpm2][9]) < (_5gdej31[0])) *((*_vlnknd6)[_ybxkpm2][9]) + (((*_vlnknd6)[_ybxkpm2][9]) >= (_5gdej31[0])) *(_5gdej31[0]) );
            (*_vlnknd6)[_ybxkpm2][10] = ( (((*_vlnknd6)[_ybxkpm2][10]) > (_5gdej31[0])) *((*_vlnknd6)[_ybxkpm2][10]) + (((*_vlnknd6)[_ybxkpm2][10]) <= (_5gdej31[0])) *(_5gdej31[0]) );
            (*_vlnknd6)[_ybxkpm2][11] = ( (((*_vlnknd6)[_ybxkpm2][11]) < (_5gdej31[1])) *((*_vlnknd6)[_ybxkpm2][11]) + (((*_vlnknd6)[_ybxkpm2][11]) >= (_5gdej31[1])) *(_5gdej31[1]) );
            (*_vlnknd6)[_ybxkpm2][12] = ( (((*_vlnknd6)[_ybxkpm2][12]) > (_5gdej31[1])) *((*_vlnknd6)[_ybxkpm2][12]) + (((*_vlnknd6)[_ybxkpm2][12]) <= (_5gdej31[1])) *(_5gdej31[1]) );
            SubtractVect3(_t26ys44, _3lbqgu3, _dan085r);
            _5gdej31[0] = _wrb7f6f[0]*_t26ys44[0] + _wrb7f6f[1]*_t26ys44[1] + _wrb7f6f[2]*_t26ys44[2];
            _5gdej31[1] = _jssuuqo[0]*_t26ys44[0] + _jssuuqo[1]*_t26ys44[1] + _jssuuqo[2]*_t26ys44[2];
            (*_vlnknd6)[_ybxkpm2][9] = ( (((*_vlnknd6)[_ybxkpm2][9]) < (_5gdej31[0])) *((*_vlnknd6)[_ybxkpm2][9]) + (((*_vlnknd6)[_ybxkpm2][9]) >= (_5gdej31[0])) *(_5gdej31[0]) );
            (*_vlnknd6)[_ybxkpm2][10] = ( (((*_vlnknd6)[_ybxkpm2][10]) > (_5gdej31[0])) *((*_vlnknd6)[_ybxkpm2][10]) + (((*_vlnknd6)[_ybxkpm2][10]) <= (_5gdej31[0])) *(_5gdej31[0]) );
            (*_vlnknd6)[_ybxkpm2][11] = ( (((*_vlnknd6)[_ybxkpm2][11]) < (_5gdej31[1])) *((*_vlnknd6)[_ybxkpm2][11]) + (((*_vlnknd6)[_ybxkpm2][11]) >= (_5gdej31[1])) *(_5gdej31[1]) );
            (*_vlnknd6)[_ybxkpm2][12] = ( (((*_vlnknd6)[_ybxkpm2][12]) > (_5gdej31[1])) *((*_vlnknd6)[_ybxkpm2][12]) + (((*_vlnknd6)[_ybxkpm2][12]) <= (_5gdej31[1])) *(_5gdej31[1]) );
    } }
    for (i = 0u; i < _vfua9bq; i++)
    { _e3nwxuj[0] = (*_vlnknd6)[i][9]; _e3nwxuj[1] = (*_vlnknd6)[i][11]; _e3nwxuj[2] = 0.0f;
        _3ac6pe6[0] = (*_vlnknd6)[i][10]; _3ac6pe6[1] = (*_vlnknd6)[i][11]; _3ac6pe6[2] = 0.0f;
        _3lbqgu3[0] = (*_vlnknd6)[i][10]; _3lbqgu3[1] = (*_vlnknd6)[i][12]; _3lbqgu3[2] = 0.0f;
        _obyghgk[0] = (*_vlnknd6)[i][9]; _obyghgk[1] = (*_vlnknd6)[i][12]; _obyghgk[2] = 0.0f;
        _q52wspo = fabsf(_3ac6pe6[0]-_e3nwxuj[0])/fabsf(_3lbqgu3[1]-_3ac6pe6[1]);
        _fv177ud = sqrtf((*_vlnknd6)[i][3]/_q52wspo);
        _9p3jpxg = (*_vlnknd6)[i][3]/_fv177ud;
        _h30153s = fabsf(_3ac6pe6[0]-_e3nwxuj[0])/_9p3jpxg;
        _flse45b = fabsf(_3lbqgu3[1]-_3ac6pe6[1])/_fv177ud;
        _e3nwxuj[0] /= _h30153s; _e3nwxuj[1] /= _flse45b;
        _3ac6pe6[0] /= _h30153s; _3ac6pe6[1] /= _flse45b;
        _3lbqgu3[0] /= _h30153s; _3lbqgu3[1] /= _flse45b;
        _obyghgk[0] /= _h30153s; _obyghgk[1] /= _flse45b;
        __builtin___memcpy_chk (_ljefy3z, &(*_vlnknd6)[i][6], 3*sizeof(float), __builtin_object_size (_ljefy3z, 0));
        GetStkAndDipVect(_ljefy3z, _wrb7f6f, _jssuuqo);
        _t26ys44[0] = _wrb7f6f[0]*_e3nwxuj[0] + _jssuuqo[0]*_e3nwxuj[1] + _ljefy3z[0]*_e3nwxuj[2] + (*_vlnknd6)[i][0];
        _t26ys44[1] = _wrb7f6f[1]*_e3nwxuj[0] + _jssuuqo[1]*_e3nwxuj[1] + _ljefy3z[1]*_e3nwxuj[2] + (*_vlnknd6)[i][1];
        _t26ys44[2] = _wrb7f6f[2]*_e3nwxuj[0] + _jssuuqo[2]*_e3nwxuj[1] + _ljefy3z[2]*_e3nwxuj[2] + (*_vlnknd6)[i][2];
        _t26ys44[2] = ( ((_t26ys44[2]) < (-0.0f)) *(_t26ys44[2]) + ((_t26ys44[2]) >= (-0.0f)) *(-0.0f) );
        __builtin___memcpy_chk (&(*_vlnknd6)[i][9], _t26ys44, 3*sizeof(float), __builtin_object_size (&(*_vlnknd6)[i][9], 0));
        _t26ys44[0] = _wrb7f6f[0]*_3ac6pe6[0] + _jssuuqo[0]*_3ac6pe6[1] + _ljefy3z[0]*_3ac6pe6[2] + (*_vlnknd6)[i][0];
        _t26ys44[1] = _wrb7f6f[1]*_3ac6pe6[0] + _jssuuqo[1]*_3ac6pe6[1] + _ljefy3z[1]*_3ac6pe6[2] + (*_vlnknd6)[i][1];
        _t26ys44[2] = _wrb7f6f[2]*_3ac6pe6[0] + _jssuuqo[2]*_3ac6pe6[1] + _ljefy3z[2]*_3ac6pe6[2] + (*_vlnknd6)[i][2];
        _t26ys44[2] = ( ((_t26ys44[2]) < (-0.0f)) *(_t26ys44[2]) + ((_t26ys44[2]) >= (-0.0f)) *(-0.0f) );
        __builtin___memcpy_chk (&(*_vlnknd6)[i][12], _t26ys44, 3*sizeof(float), __builtin_object_size (&(*_vlnknd6)[i][12], 0));
        _t26ys44[0] = _wrb7f6f[0]*_3lbqgu3[0] + _jssuuqo[0]*_3lbqgu3[1] + _ljefy3z[0]*_3lbqgu3[2] + (*_vlnknd6)[i][0];
        _t26ys44[1] = _wrb7f6f[1]*_3lbqgu3[0] + _jssuuqo[1]*_3lbqgu3[1] + _ljefy3z[1]*_3lbqgu3[2] + (*_vlnknd6)[i][1];
        _t26ys44[2] = _wrb7f6f[2]*_3lbqgu3[0] + _jssuuqo[2]*_3lbqgu3[1] + _ljefy3z[2]*_3lbqgu3[2] + (*_vlnknd6)[i][2];
        _t26ys44[2] = ( ((_t26ys44[2]) < (-0.0f)) *(_t26ys44[2]) + ((_t26ys44[2]) >= (-0.0f)) *(-0.0f) );
        __builtin___memcpy_chk (&(*_vlnknd6)[i][15], _t26ys44, 3*sizeof(float), __builtin_object_size (&(*_vlnknd6)[i][15], 0));
        _t26ys44[0] = _wrb7f6f[0]*_obyghgk[0] + _jssuuqo[0]*_obyghgk[1] + _ljefy3z[0]*_obyghgk[2] + (*_vlnknd6)[i][0];
        _t26ys44[1] = _wrb7f6f[1]*_obyghgk[0] + _jssuuqo[1]*_obyghgk[1] + _ljefy3z[1]*_obyghgk[2] + (*_vlnknd6)[i][1];
        _t26ys44[2] = _wrb7f6f[2]*_obyghgk[0] + _jssuuqo[2]*_obyghgk[1] + _ljefy3z[2]*_obyghgk[2] + (*_vlnknd6)[i][2];
        _t26ys44[2] = ( ((_t26ys44[2]) < (-0.0f)) *(_t26ys44[2]) + ((_t26ys44[2]) >= (-0.0f)) *(-0.0f) );
        __builtin___memcpy_chk (&(*_vlnknd6)[i][18], _t26ys44, 3*sizeof(float), __builtin_object_size (&(*_vlnknd6)[i][18], 0));
        SubtractVect3(_2ljk4q5, &(*_vlnknd6)[i][0], &(*_vlnknd6)[i][9]); _1vsmado = VectorLength(_2ljk4q5);
        (*_vlnknd6)[i][5] = ( ((_1vsmado) > ((*_vlnknd6)[i][5])) *(_1vsmado) + ((_1vsmado) <= ((*_vlnknd6)[i][5])) *((*_vlnknd6)[i][5]) );
        SubtractVect3(_2ljk4q5, &(*_vlnknd6)[i][0], &(*_vlnknd6)[i][12]); _1vsmado = VectorLength(_2ljk4q5);
        (*_vlnknd6)[i][5] = ( ((_1vsmado) > ((*_vlnknd6)[i][5])) *(_1vsmado) + ((_1vsmado) <= ((*_vlnknd6)[i][5])) *((*_vlnknd6)[i][5]) );
        SubtractVect3(_2ljk4q5, &(*_vlnknd6)[i][0], &(*_vlnknd6)[i][15]); _1vsmado = VectorLength(_2ljk4q5);
        (*_vlnknd6)[i][5] = ( ((_1vsmado) > ((*_vlnknd6)[i][5])) *(_1vsmado) + ((_1vsmado) <= ((*_vlnknd6)[i][5])) *((*_vlnknd6)[i][5]) );
        SubtractVect3(_2ljk4q5, &(*_vlnknd6)[i][0], &(*_vlnknd6)[i][18]); _1vsmado = VectorLength(_2ljk4q5);
        (*_vlnknd6)[i][5] = ( ((_1vsmado) > ((*_vlnknd6)[i][5])) *(_1vsmado) + ((_1vsmado) <= ((*_vlnknd6)[i][5])) *((*_vlnknd6)[i][5]) );
    }
    return;
}
void _j4t0sek(const char *restrict _slwf90k, unsigned int _2vpgstq, const float *restrict _q01tee1, unsigned int **restrict _qc6dr08, unsigned int _qsq0p82, float **restrict _zh13cd5, unsigned int _04k0lth, const float *restrict _o6l1f90, unsigned int **restrict _1zlbori, unsigned int _6parsg7, float **restrict _s3jonou)
{
    int i, _gcsotdh = 16;
    char _z049enu[512];
    float *_zkwdeek = ((void *)0), *_p2u4dtb = ((void *)0);
    __builtin___strcpy_chk (_z049enu, _slwf90k, __builtin_object_size (_z049enu, 2 > 1 ? 1 : 0)); __builtin___strcat_chk (_z049enu, ".brch", __builtin_object_size (_z049enu, 2 > 1 ? 1 : 0));
    FILE *_e7eobjl;
    _zkwdeek = malloc(3*_2vpgstq*sizeof(float));
    for (i = 0; i < _2vpgstq; i++) { __builtin___memcpy_chk (&_zkwdeek[i*3 +0], &_q01tee1[i*16 +0],3*sizeof(float), __builtin_object_size (&_zkwdeek[i*3 +0], 0)); }
    if (_04k0lth > 0u)
    { _p2u4dtb = malloc(3*_04k0lth*sizeof(float));
        for (i = 0; i < _04k0lth; i++) { __builtin___memcpy_chk (&_p2u4dtb[i*3 +0], &_o6l1f90[i*16 +0],3*sizeof(float), __builtin_object_size (&_p2u4dtb[i*3 +0], 0)); }
    }
    if ((_e7eobjl = fopen(_z049enu,"wb")) == ((void *)0)) { printf("Error -cant open BranchInfo file for writing   %s\n", _z049enu); exit(10); }
    int _c54wsqt = (int)202511;
    fwrite( &_c54wsqt, sizeof(int), 1, _e7eobjl);
    fwrite( &_2vpgstq, sizeof(unsigned int), 1, _e7eobjl);
    fwrite( &_qsq0p82, sizeof(unsigned int), 1, _e7eobjl);
    for (i = 0; i < _2vpgstq; i++) { fwrite(_qc6dr08[i], sizeof(unsigned int), _gcsotdh, _e7eobjl); }
    for (i = 0; i < _qsq0p82; i++) { fwrite(_zh13cd5[i], sizeof(float), 21, _e7eobjl); }
    fwrite(_zkwdeek, sizeof(float), 3*_2vpgstq, _e7eobjl);
    if (_04k0lth > 0u)
    { fwrite( &_04k0lth, sizeof(unsigned int), 1, _e7eobjl);
        fwrite( &_6parsg7, sizeof(unsigned int), 1, _e7eobjl);
        for (i = 0; i < _04k0lth; i++) { fwrite(_1zlbori[i], sizeof(unsigned int), _gcsotdh, _e7eobjl); }
        for (i = 0; i < _6parsg7; i++) { fwrite(_s3jonou[i], sizeof(float), 21, _e7eobjl); }
        fwrite(_p2u4dtb, sizeof(float), 3*_04k0lth, _e7eobjl);
    }
    fclose(_e7eobjl);
    free(_zkwdeek); free(_p2u4dtb);
    return;
}
