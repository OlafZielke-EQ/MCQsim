
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
typedef long int ptrdiff_t;
typedef long double max_align_t;
typedef long BLASLONG;
typedef unsigned long BLASULONG;
typedef uint16_t bfloat16;
typedef uint16_t hfloat16;
typedef int blasint;
extern float _Complex cacosf(float _Complex);
extern double _Complex cacos(double _Complex);
extern long double _Complex cacosl(long double _Complex);
extern float _Complex casinf(float _Complex);
extern double _Complex casin(double _Complex);
extern long double _Complex casinl(long double _Complex);
extern float _Complex catanf(float _Complex);
extern double _Complex catan(double _Complex);
extern long double _Complex catanl(long double _Complex);
extern float _Complex ccosf(float _Complex);
extern double _Complex ccos(double _Complex);
extern long double _Complex ccosl(long double _Complex);
extern float _Complex csinf(float _Complex);
extern double _Complex csin(double _Complex);
extern long double _Complex csinl(long double _Complex);
extern float _Complex ctanf(float _Complex);
extern double _Complex ctan(double _Complex);
extern long double _Complex ctanl(long double _Complex);
extern float _Complex cacoshf(float _Complex);
extern double _Complex cacosh(double _Complex);
extern long double _Complex cacoshl(long double _Complex);
extern float _Complex casinhf(float _Complex);
extern double _Complex casinh(double _Complex);
extern long double _Complex casinhl(long double _Complex);
extern float _Complex catanhf(float _Complex);
extern double _Complex catanh(double _Complex);
extern long double _Complex catanhl(long double _Complex);
extern float _Complex ccoshf(float _Complex);
extern double _Complex ccosh(double _Complex);
extern long double _Complex ccoshl(long double _Complex);
extern float _Complex csinhf(float _Complex);
extern double _Complex csinh(double _Complex);
extern long double _Complex csinhl(long double _Complex);
extern float _Complex ctanhf(float _Complex);
extern double _Complex ctanh(double _Complex);
extern long double _Complex ctanhl(long double _Complex);
extern float _Complex cexpf(float _Complex);
extern double _Complex cexp(double _Complex);
extern long double _Complex cexpl(long double _Complex);
extern float _Complex clogf(float _Complex);
extern double _Complex clog(double _Complex);
extern long double _Complex clogl(long double _Complex);
extern float cabsf(float _Complex);
extern double cabs(double _Complex);
extern long double cabsl(long double _Complex);
extern float _Complex cpowf(float _Complex, float _Complex);
extern double _Complex cpow(double _Complex, double _Complex);
extern long double _Complex cpowl(long double _Complex, long double _Complex);
extern float _Complex csqrtf(float _Complex);
extern double _Complex csqrt(double _Complex);
extern long double _Complex csqrtl(long double _Complex);
extern float cargf(float _Complex);
extern double carg(double _Complex);
extern long double cargl(long double _Complex);
extern float cimagf(float _Complex);
extern double cimag(double _Complex);
extern long double cimagl(long double _Complex);
extern float _Complex conjf(float _Complex);
extern double _Complex conj(double _Complex);
extern long double _Complex conjl(long double _Complex);
extern float _Complex cprojf(float _Complex);
extern double _Complex cproj(double _Complex);
extern long double _Complex cprojl(long double _Complex);
extern float crealf(float _Complex);
extern double creal(double _Complex);
extern long double creall(long double _Complex);

  typedef float _Complex openblas_complex_float;
  typedef double _Complex openblas_complex_double;
  typedef double _Complex openblas_complex_xdouble;
void openblas_set_num_threads(int num_threads);
void goto_set_num_threads(int num_threads);
int openblas_set_num_threads_local(int num_threads);
int openblas_get_num_threads(void);
int openblas_get_num_procs(void);
char* openblas_get_config(void);
char* openblas_get_corename(void);
typedef void (*openblas_dojob_callback)(int thread_num, void *jobdata, int dojob_data);
typedef void (*openblas_threads_callback)(int sync, openblas_dojob_callback dojob, int numjobs, size_t jobdata_elsize, void *jobdata, int dojob_data);
void openblas_set_threads_callback_function(openblas_threads_callback callback);
int openblas_get_parallel(void);
typedef enum CBLAS_ORDER {CblasRowMajor=101, CblasColMajor=102} CBLAS_ORDER;
typedef enum CBLAS_TRANSPOSE {CblasNoTrans=111, CblasTrans=112, CblasConjTrans=113, CblasConjNoTrans=114} CBLAS_TRANSPOSE;
typedef enum CBLAS_UPLO {CblasUpper=121, CblasLower=122} CBLAS_UPLO;
typedef enum CBLAS_DIAG {CblasNonUnit=131, CblasUnit=132} CBLAS_DIAG;
typedef enum CBLAS_SIDE {CblasLeft=141, CblasRight=142} CBLAS_SIDE;
typedef CBLAS_ORDER CBLAS_LAYOUT;
float cblas_sdsdot(const blasint n, const float alpha, const float *x, const blasint incx, const float *y, const blasint incy);
double cblas_dsdot (const blasint n, const float *x, const blasint incx, const float *y, const blasint incy);
float cblas_sdot(const blasint n, const float *x, const blasint incx, const float *y, const blasint incy);
double cblas_ddot(const blasint n, const double *x, const blasint incx, const double *y, const blasint incy);
openblas_complex_float cblas_cdotu(const blasint n, const void *x, const blasint incx, const void *y, const blasint incy);
openblas_complex_float cblas_cdotc(const blasint n, const void *x, const blasint incx, const void *y, const blasint incy);
openblas_complex_double cblas_zdotu(const blasint n, const void *x, const blasint incx, const void *y, const blasint incy);
openblas_complex_double cblas_zdotc(const blasint n, const void *x, const blasint incx, const void *y, const blasint incy);
void cblas_cdotu_sub(const blasint n, const void *x, const blasint incx, const void *y, const blasint incy, void *ret);
void cblas_cdotc_sub(const blasint n, const void *x, const blasint incx, const void *y, const blasint incy, void *ret);
void cblas_zdotu_sub(const blasint n, const void *x, const blasint incx, const void *y, const blasint incy, void *ret);
void cblas_zdotc_sub(const blasint n, const void *x, const blasint incx, const void *y, const blasint incy, void *ret);
float cblas_sasum (const blasint n, const float *x, const blasint incx);
double cblas_dasum (const blasint n, const double *x, const blasint incx);
float cblas_scasum(const blasint n, const void *x, const blasint incx);
double cblas_dzasum(const blasint n, const void *x, const blasint incx);
float cblas_ssum (const blasint n, const float *x, const blasint incx);
double cblas_dsum (const blasint n, const double *x, const blasint incx);
float cblas_scsum(const blasint n, const void *x, const blasint incx);
double cblas_dzsum(const blasint n, const void *x, const blasint incx);
float cblas_snrm2 (const blasint N, const float *X, const blasint incX);
double cblas_dnrm2 (const blasint N, const double *X, const blasint incX);
float cblas_scnrm2(const blasint N, const void *X, const blasint incX);
double cblas_dznrm2(const blasint N, const void *X, const blasint incX);
size_t cblas_isamax(const blasint n, const float *x, const blasint incx);
size_t cblas_idamax(const blasint n, const double *x, const blasint incx);
size_t cblas_icamax(const blasint n, const void *x, const blasint incx);
size_t cblas_izamax(const blasint n, const void *x, const blasint incx);
size_t cblas_isamin(const blasint n, const float *x, const blasint incx);
size_t cblas_idamin(const blasint n, const double *x, const blasint incx);
size_t cblas_icamin(const blasint n, const void *x, const blasint incx);
size_t cblas_izamin(const blasint n, const void *x, const blasint incx);
float cblas_samax(const blasint n, const float *x, const blasint incx);
double cblas_damax(const blasint n, const double *x, const blasint incx);
float cblas_scamax(const blasint n, const void *x, const blasint incx);
double cblas_dzamax(const blasint n, const void *x, const blasint incx);
float cblas_samin(const blasint n, const float *x, const blasint incx);
double cblas_damin(const blasint n, const double *x, const blasint incx);
float cblas_scamin(const blasint n, const void *x, const blasint incx);
double cblas_dzamin(const blasint n, const void *x, const blasint incx);
size_t cblas_ismax(const blasint n, const float *x, const blasint incx);
size_t cblas_idmax(const blasint n, const double *x, const blasint incx);
size_t cblas_icmax(const blasint n, const void *x, const blasint incx);
size_t cblas_izmax(const blasint n, const void *x, const blasint incx);
size_t cblas_ismin(const blasint n, const float *x, const blasint incx);
size_t cblas_idmin(const blasint n, const double *x, const blasint incx);
size_t cblas_icmin(const blasint n, const void *x, const blasint incx);
size_t cblas_izmin(const blasint n, const void *x, const blasint incx);
void cblas_saxpy(const blasint n, const float alpha, const float *x, const blasint incx, float *y, const blasint incy);
void cblas_daxpy(const blasint n, const double alpha, const double *x, const blasint incx, double *y, const blasint incy);
void cblas_caxpy(const blasint n, const void *alpha, const void *x, const blasint incx, void *y, const blasint incy);
void cblas_zaxpy(const blasint n, const void *alpha, const void *x, const blasint incx, void *y, const blasint incy);
void cblas_caxpyc(const blasint n, const void *alpha, const void *x, const blasint incx, void *y, const blasint incy);
void cblas_zaxpyc(const blasint n, const void *alpha, const void *x, const blasint incx, void *y, const blasint incy);
void cblas_scopy(const blasint n, const float *x, const blasint incx, float *y, const blasint incy);
void cblas_dcopy(const blasint n, const double *x, const blasint incx, double *y, const blasint incy);
void cblas_ccopy(const blasint n, const void *x, const blasint incx, void *y, const blasint incy);
void cblas_zcopy(const blasint n, const void *x, const blasint incx, void *y, const blasint incy);
void cblas_sswap(const blasint n, float *x, const blasint incx, float *y, const blasint incy);
void cblas_dswap(const blasint n, double *x, const blasint incx, double *y, const blasint incy);
void cblas_cswap(const blasint n, void *x, const blasint incx, void *y, const blasint incy);
void cblas_zswap(const blasint n, void *x, const blasint incx, void *y, const blasint incy);
void cblas_srot(const blasint N, float *X, const blasint incX, float *Y, const blasint incY, const float c, const float s);
void cblas_drot(const blasint N, double *X, const blasint incX, double *Y, const blasint incY, const double c, const double s);
void cblas_csrot(const blasint n, const void *x, const blasint incx, void *y, const blasint incY, const float c, const float s);
void cblas_zdrot(const blasint n, const void *x, const blasint incx, void *y, const blasint incY, const double c, const double s);
void cblas_srotg(float *a, float *b, float *c, float *s);
void cblas_drotg(double *a, double *b, double *c, double *s);
void cblas_crotg(void *a, void *b, float *c, void *s);
void cblas_zrotg(void *a, void *b, double *c, void *s);
void cblas_srotm(const blasint N, float *X, const blasint incX, float *Y, const blasint incY, const float *P);
void cblas_drotm(const blasint N, double *X, const blasint incX, double *Y, const blasint incY, const double *P);
void cblas_srotmg(float *d1, float *d2, float *b1, const float b2, float *P);
void cblas_drotmg(double *d1, double *d2, double *b1, const double b2, double *P);
void cblas_sscal(const blasint N, const float alpha, float *X, const blasint incX);
void cblas_dscal(const blasint N, const double alpha, double *X, const blasint incX);
void cblas_cscal(const blasint N, const void *alpha, void *X, const blasint incX);
void cblas_zscal(const blasint N, const void *alpha, void *X, const blasint incX);
void cblas_csscal(const blasint N, const float alpha, void *X, const blasint incX);
void cblas_zdscal(const blasint N, const double alpha, void *X, const blasint incX);
void cblas_sgemv(const enum CBLAS_ORDER order, const enum CBLAS_TRANSPOSE trans, const blasint m, const blasint n,
   const float alpha, const float *a, const blasint lda, const float *x, const blasint incx, const float beta, float *y, const blasint incy);
void cblas_dgemv(const enum CBLAS_ORDER order, const enum CBLAS_TRANSPOSE trans, const blasint m, const blasint n,
   const double alpha, const double *a, const blasint lda, const double *x, const blasint incx, const double beta, double *y, const blasint incy);
void cblas_cgemv(const enum CBLAS_ORDER order, const enum CBLAS_TRANSPOSE trans, const blasint m, const blasint n,
   const void *alpha, const void *a, const blasint lda, const void *x, const blasint incx, const void *beta, void *y, const blasint incy);
void cblas_zgemv(const enum CBLAS_ORDER order, const enum CBLAS_TRANSPOSE trans, const blasint m, const blasint n,
   const void *alpha, const void *a, const blasint lda, const void *x, const blasint incx, const void *beta, void *y, const blasint incy);
void cblas_sger (const enum CBLAS_ORDER order, const blasint M, const blasint N, const float alpha, const float *X, const blasint incX, const float *Y, const blasint incY, float *A, const blasint lda);
void cblas_dger (const enum CBLAS_ORDER order, const blasint M, const blasint N, const double alpha, const double *X, const blasint incX, const double *Y, const blasint incY, double *A, const blasint lda);
void cblas_cgeru(const enum CBLAS_ORDER order, const blasint M, const blasint N, const void *alpha, const void *X, const blasint incX, const void *Y, const blasint incY, void *A, const blasint lda);
void cblas_cgerc(const enum CBLAS_ORDER order, const blasint M, const blasint N, const void *alpha, const void *X, const blasint incX, const void *Y, const blasint incY, void *A, const blasint lda);
void cblas_zgeru(const enum CBLAS_ORDER order, const blasint M, const blasint N, const void *alpha, const void *X, const blasint incX, const void *Y, const blasint incY, void *A, const blasint lda);
void cblas_zgerc(const enum CBLAS_ORDER order, const blasint M, const blasint N, const void *alpha, const void *X, const blasint incX, const void *Y, const blasint incY, void *A, const blasint lda);
void cblas_strsv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag, const blasint N, const float *A, const blasint lda, float *X, const blasint incX);
void cblas_dtrsv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag, const blasint N, const double *A, const blasint lda, double *X, const blasint incX);
void cblas_ctrsv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag, const blasint N, const void *A, const blasint lda, void *X, const blasint incX);
void cblas_ztrsv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag, const blasint N, const void *A, const blasint lda, void *X, const blasint incX);
void cblas_strmv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag, const blasint N, const float *A, const blasint lda, float *X, const blasint incX);
void cblas_dtrmv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag, const blasint N, const double *A, const blasint lda, double *X, const blasint incX);
void cblas_ctrmv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag, const blasint N, const void *A, const blasint lda, void *X, const blasint incX);
void cblas_ztrmv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag, const blasint N, const void *A, const blasint lda, void *X, const blasint incX);
void cblas_ssyr(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const blasint N, const float alpha, const float *X, const blasint incX, float *A, const blasint lda);
void cblas_dsyr(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const blasint N, const double alpha, const double *X, const blasint incX, double *A, const blasint lda);
void cblas_cher(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const blasint N, const float alpha, const void *X, const blasint incX, void *A, const blasint lda);
void cblas_zher(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const blasint N, const double alpha, const void *X, const blasint incX, void *A, const blasint lda);
void cblas_ssyr2(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,const blasint N, const float alpha, const float *X,
                const blasint incX, const float *Y, const blasint incY, float *A, const blasint lda);
void cblas_dsyr2(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const blasint N, const double alpha, const double *X,
                const blasint incX, const double *Y, const blasint incY, double *A, const blasint lda);
void cblas_cher2(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const blasint N, const void *alpha, const void *X, const blasint incX,
                const void *Y, const blasint incY, void *A, const blasint lda);
void cblas_zher2(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const blasint N, const void *alpha, const void *X, const blasint incX,
                const void *Y, const blasint incY, void *A, const blasint lda);
void cblas_sgbmv(const enum CBLAS_ORDER order, const enum CBLAS_TRANSPOSE TransA, const blasint M, const blasint N,
                 const blasint KL, const blasint KU, const float alpha, const float *A, const blasint lda, const float *X, const blasint incX, const float beta, float *Y, const blasint incY);
void cblas_dgbmv(const enum CBLAS_ORDER order, const enum CBLAS_TRANSPOSE TransA, const blasint M, const blasint N,
                 const blasint KL, const blasint KU, const double alpha, const double *A, const blasint lda, const double *X, const blasint incX, const double beta, double *Y, const blasint incY);
void cblas_cgbmv(const enum CBLAS_ORDER order, const enum CBLAS_TRANSPOSE TransA, const blasint M, const blasint N,
                 const blasint KL, const blasint KU, const void *alpha, const void *A, const blasint lda, const void *X, const blasint incX, const void *beta, void *Y, const blasint incY);
void cblas_zgbmv(const enum CBLAS_ORDER order, const enum CBLAS_TRANSPOSE TransA, const blasint M, const blasint N,
                 const blasint KL, const blasint KU, const void *alpha, const void *A, const blasint lda, const void *X, const blasint incX, const void *beta, void *Y, const blasint incY);
void cblas_ssbmv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const blasint N, const blasint K, const float alpha, const float *A,
                 const blasint lda, const float *X, const blasint incX, const float beta, float *Y, const blasint incY);
void cblas_dsbmv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const blasint N, const blasint K, const double alpha, const double *A,
                 const blasint lda, const double *X, const blasint incX, const double beta, double *Y, const blasint incY);
void cblas_stbmv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag,
                 const blasint N, const blasint K, const float *A, const blasint lda, float *X, const blasint incX);
void cblas_dtbmv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag,
                 const blasint N, const blasint K, const double *A, const blasint lda, double *X, const blasint incX);
void cblas_ctbmv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag,
                 const blasint N, const blasint K, const void *A, const blasint lda, void *X, const blasint incX);
void cblas_ztbmv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag,
                 const blasint N, const blasint K, const void *A, const blasint lda, void *X, const blasint incX);
void cblas_stbsv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag,
                 const blasint N, const blasint K, const float *A, const blasint lda, float *X, const blasint incX);
void cblas_dtbsv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag,
                 const blasint N, const blasint K, const double *A, const blasint lda, double *X, const blasint incX);
void cblas_ctbsv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag,
                 const blasint N, const blasint K, const void *A, const blasint lda, void *X, const blasint incX);
void cblas_ztbsv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag,
                 const blasint N, const blasint K, const void *A, const blasint lda, void *X, const blasint incX);
void cblas_stpmv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag,
                 const blasint N, const float *Ap, float *X, const blasint incX);
void cblas_dtpmv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag,
                 const blasint N, const double *Ap, double *X, const blasint incX);
void cblas_ctpmv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag,
                 const blasint N, const void *Ap, void *X, const blasint incX);
void cblas_ztpmv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag,
                 const blasint N, const void *Ap, void *X, const blasint incX);
void cblas_stpsv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag,
                 const blasint N, const float *Ap, float *X, const blasint incX);
void cblas_dtpsv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag,
                 const blasint N, const double *Ap, double *X, const blasint incX);
void cblas_ctpsv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag,
                 const blasint N, const void *Ap, void *X, const blasint incX);
void cblas_ztpsv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag,
                 const blasint N, const void *Ap, void *X, const blasint incX);
void cblas_ssymv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const blasint N, const float alpha, const float *A,
                 const blasint lda, const float *X, const blasint incX, const float beta, float *Y, const blasint incY);
void cblas_dsymv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const blasint N, const double alpha, const double *A,
                 const blasint lda, const double *X, const blasint incX, const double beta, double *Y, const blasint incY);
void cblas_chemv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const blasint N, const void *alpha, const void *A,
                 const blasint lda, const void *X, const blasint incX, const void *beta, void *Y, const blasint incY);
void cblas_zhemv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const blasint N, const void *alpha, const void *A,
                 const blasint lda, const void *X, const blasint incX, const void *beta, void *Y, const blasint incY);
void cblas_sspmv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const blasint N, const float alpha, const float *Ap,
                 const float *X, const blasint incX, const float beta, float *Y, const blasint incY);
void cblas_dspmv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const blasint N, const double alpha, const double *Ap,
                 const double *X, const blasint incX, const double beta, double *Y, const blasint incY);
void cblas_sspr(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const blasint N, const float alpha, const float *X, const blasint incX, float *Ap);
void cblas_dspr(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const blasint N, const double alpha, const double *X, const blasint incX, double *Ap);
void cblas_chpr(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const blasint N, const float alpha, const void *X, const blasint incX, void *A);
void cblas_zhpr(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const blasint N, const double alpha, const void *X,const blasint incX, void *A);
void cblas_sspr2(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const blasint N, const float alpha, const float *X, const blasint incX, const float *Y, const blasint incY, float *A);
void cblas_dspr2(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const blasint N, const double alpha, const double *X, const blasint incX, const double *Y, const blasint incY, double *A);
void cblas_chpr2(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const blasint N, const void *alpha, const void *X, const blasint incX, const void *Y, const blasint incY, void *Ap);
void cblas_zhpr2(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const blasint N, const void *alpha, const void *X, const blasint incX, const void *Y, const blasint incY, void *Ap);
void cblas_chbmv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const blasint N, const blasint K,
   const void *alpha, const void *A, const blasint lda, const void *X, const blasint incX, const void *beta, void *Y, const blasint incY);
void cblas_zhbmv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const blasint N, const blasint K,
   const void *alpha, const void *A, const blasint lda, const void *X, const blasint incX, const void *beta, void *Y, const blasint incY);
void cblas_chpmv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const blasint N,
   const void *alpha, const void *Ap, const void *X, const blasint incX, const void *beta, void *Y, const blasint incY);
void cblas_zhpmv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const blasint N,
   const void *alpha, const void *Ap, const void *X, const blasint incX, const void *beta, void *Y, const blasint incY);
void cblas_sgemm(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_TRANSPOSE TransB, const blasint M, const blasint N, const blasint K,
   const float alpha, const float *A, const blasint lda, const float *B, const blasint ldb, const float beta, float *C, const blasint ldc);
void cblas_dgemm(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_TRANSPOSE TransB, const blasint M, const blasint N, const blasint K,
   const double alpha, const double *A, const blasint lda, const double *B, const blasint ldb, const double beta, double *C, const blasint ldc);
void cblas_cgemm(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_TRANSPOSE TransB, const blasint M, const blasint N, const blasint K,
   const void *alpha, const void *A, const blasint lda, const void *B, const blasint ldb, const void *beta, void *C, const blasint ldc);
void cblas_cgemm3m(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_TRANSPOSE TransB, const blasint M, const blasint N, const blasint K,
   const void *alpha, const void *A, const blasint lda, const void *B, const blasint ldb, const void *beta, void *C, const blasint ldc);
void cblas_zgemm(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_TRANSPOSE TransB, const blasint M, const blasint N, const blasint K,
   const void *alpha, const void *A, const blasint lda, const void *B, const blasint ldb, const void *beta, void *C, const blasint ldc);
void cblas_zgemm3m(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_TRANSPOSE TransB, const blasint M, const blasint N, const blasint K,
   const void *alpha, const void *A, const blasint lda, const void *B, const blasint ldb, const void *beta, void *C, const blasint ldc);
void cblas_sgemmt(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_TRANSPOSE TransB, const blasint M, const blasint K,
   const float alpha, const float *A, const blasint lda, const float *B, const blasint ldb, const float beta, float *C, const blasint ldc);
void cblas_dgemmt(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_TRANSPOSE TransB, const blasint M, const blasint K,
   const double alpha, const double *A, const blasint lda, const double *B, const blasint ldb, const double beta, double *C, const blasint ldc);
void cblas_cgemmt(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_TRANSPOSE TransB, const blasint M, const blasint K,
   const void *alpha, const void *A, const blasint lda, const void *B, const blasint ldb, const void *beta, void *C, const blasint ldc);
void cblas_zgemmt(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_TRANSPOSE TransB, const blasint M, const blasint K,
   const void *alpha, const void *A, const blasint lda, const void *B, const blasint ldb, const void *beta, void *C, const blasint ldc);
void cblas_ssymm(const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side, const enum CBLAS_UPLO Uplo, const blasint M, const blasint N,
                 const float alpha, const float *A, const blasint lda, const float *B, const blasint ldb, const float beta, float *C, const blasint ldc);
void cblas_dsymm(const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side, const enum CBLAS_UPLO Uplo, const blasint M, const blasint N,
                 const double alpha, const double *A, const blasint lda, const double *B, const blasint ldb, const double beta, double *C, const blasint ldc);
void cblas_csymm(const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side, const enum CBLAS_UPLO Uplo, const blasint M, const blasint N,
                 const void *alpha, const void *A, const blasint lda, const void *B, const blasint ldb, const void *beta, void *C, const blasint ldc);
void cblas_zsymm(const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side, const enum CBLAS_UPLO Uplo, const blasint M, const blasint N,
                 const void *alpha, const void *A, const blasint lda, const void *B, const blasint ldb, const void *beta, void *C, const blasint ldc);
void cblas_ssyrk(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE Trans,
   const blasint N, const blasint K, const float alpha, const float *A, const blasint lda, const float beta, float *C, const blasint ldc);
void cblas_dsyrk(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE Trans,
   const blasint N, const blasint K, const double alpha, const double *A, const blasint lda, const double beta, double *C, const blasint ldc);
void cblas_csyrk(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE Trans,
   const blasint N, const blasint K, const void *alpha, const void *A, const blasint lda, const void *beta, void *C, const blasint ldc);
void cblas_zsyrk(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE Trans,
   const blasint N, const blasint K, const void *alpha, const void *A, const blasint lda, const void *beta, void *C, const blasint ldc);
void cblas_ssyr2k(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE Trans,
    const blasint N, const blasint K, const float alpha, const float *A, const blasint lda, const float *B, const blasint ldb, const float beta, float *C, const blasint ldc);
void cblas_dsyr2k(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE Trans,
    const blasint N, const blasint K, const double alpha, const double *A, const blasint lda, const double *B, const blasint ldb, const double beta, double *C, const blasint ldc);
void cblas_csyr2k(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE Trans,
    const blasint N, const blasint K, const void *alpha, const void *A, const blasint lda, const void *B, const blasint ldb, const void *beta, void *C, const blasint ldc);
void cblas_zsyr2k(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE Trans,
    const blasint N, const blasint K, const void *alpha, const void *A, const blasint lda, const void *B, const blasint ldb, const void *beta, void *C, const blasint ldc);
void cblas_strmm(const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA,
                 const enum CBLAS_DIAG Diag, const blasint M, const blasint N, const float alpha, const float *A, const blasint lda, float *B, const blasint ldb);
void cblas_dtrmm(const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA,
                 const enum CBLAS_DIAG Diag, const blasint M, const blasint N, const double alpha, const double *A, const blasint lda, double *B, const blasint ldb);
void cblas_ctrmm(const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA,
                 const enum CBLAS_DIAG Diag, const blasint M, const blasint N, const void *alpha, const void *A, const blasint lda, void *B, const blasint ldb);
void cblas_ztrmm(const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA,
                 const enum CBLAS_DIAG Diag, const blasint M, const blasint N, const void *alpha, const void *A, const blasint lda, void *B, const blasint ldb);
void cblas_strsm(const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA,
                 const enum CBLAS_DIAG Diag, const blasint M, const blasint N, const float alpha, const float *A, const blasint lda, float *B, const blasint ldb);
void cblas_dtrsm(const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA,
                 const enum CBLAS_DIAG Diag, const blasint M, const blasint N, const double alpha, const double *A, const blasint lda, double *B, const blasint ldb);
void cblas_ctrsm(const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA,
                 const enum CBLAS_DIAG Diag, const blasint M, const blasint N, const void *alpha, const void *A, const blasint lda, void *B, const blasint ldb);
void cblas_ztrsm(const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA,
                 const enum CBLAS_DIAG Diag, const blasint M, const blasint N, const void *alpha, const void *A, const blasint lda, void *B, const blasint ldb);
void cblas_chemm(const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side, const enum CBLAS_UPLO Uplo, const blasint M, const blasint N,
                 const void *alpha, const void *A, const blasint lda, const void *B, const blasint ldb, const void *beta, void *C, const blasint ldc);
void cblas_zhemm(const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side, const enum CBLAS_UPLO Uplo, const blasint M, const blasint N,
                 const void *alpha, const void *A, const blasint lda, const void *B, const blasint ldb, const void *beta, void *C, const blasint ldc);
void cblas_cherk(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE Trans, const blasint N, const blasint K,
                 const float alpha, const void *A, const blasint lda, const float beta, void *C, const blasint ldc);
void cblas_zherk(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE Trans, const blasint N, const blasint K,
                 const double alpha, const void *A, const blasint lda, const double beta, void *C, const blasint ldc);
void cblas_cher2k(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE Trans, const blasint N, const blasint K,
                  const void *alpha, const void *A, const blasint lda, const void *B, const blasint ldb, const float beta, void *C, const blasint ldc);
void cblas_zher2k(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE Trans, const blasint N, const blasint K,
                  const void *alpha, const void *A, const blasint lda, const void *B, const blasint ldb, const double beta, void *C, const blasint ldc);
void cblas_xerbla(blasint p, const char *rout, const char *form, ...);
void cblas_saxpby(const blasint n, const float alpha, const float *x, const blasint incx,const float beta, float *y, const blasint incy);
void cblas_daxpby(const blasint n, const double alpha, const double *x, const blasint incx,const double beta, double *y, const blasint incy);
void cblas_caxpby(const blasint n, const void *alpha, const void *x, const blasint incx,const void *beta, void *y, const blasint incy);
void cblas_zaxpby(const blasint n, const void *alpha, const void *x, const blasint incx,const void *beta, void *y, const blasint incy);
void cblas_somatcopy(const enum CBLAS_ORDER CORDER, const enum CBLAS_TRANSPOSE CTRANS, const blasint crows, const blasint ccols, const float calpha, const float *a,
       const blasint clda, float *b, const blasint cldb);
void cblas_domatcopy(const enum CBLAS_ORDER CORDER, const enum CBLAS_TRANSPOSE CTRANS, const blasint crows, const blasint ccols, const double calpha, const double *a,
       const blasint clda, double *b, const blasint cldb);
void cblas_comatcopy(const enum CBLAS_ORDER CORDER, const enum CBLAS_TRANSPOSE CTRANS, const blasint crows, const blasint ccols, const float* calpha, const float* a,
       const blasint clda, float*b, const blasint cldb);
void cblas_zomatcopy(const enum CBLAS_ORDER CORDER, const enum CBLAS_TRANSPOSE CTRANS, const blasint crows, const blasint ccols, const double* calpha, const double* a,
       const blasint clda, double *b, const blasint cldb);
void cblas_simatcopy(const enum CBLAS_ORDER CORDER, const enum CBLAS_TRANSPOSE CTRANS, const blasint crows, const blasint ccols, const float calpha, float *a,
       const blasint clda, const blasint cldb);
void cblas_dimatcopy(const enum CBLAS_ORDER CORDER, const enum CBLAS_TRANSPOSE CTRANS, const blasint crows, const blasint ccols, const double calpha, double *a,
       const blasint clda, const blasint cldb);
void cblas_cimatcopy(const enum CBLAS_ORDER CORDER, const enum CBLAS_TRANSPOSE CTRANS, const blasint crows, const blasint ccols, const float* calpha, float* a,
       const blasint clda, const blasint cldb);
void cblas_zimatcopy(const enum CBLAS_ORDER CORDER, const enum CBLAS_TRANSPOSE CTRANS, const blasint crows, const blasint ccols, const double* calpha, double* a,
       const blasint clda, const blasint cldb);
void cblas_sgeadd(const enum CBLAS_ORDER CORDER,const blasint crows, const blasint ccols, const float calpha, const float *a, const blasint clda, const float cbeta,
    float *c, const blasint cldc);
void cblas_dgeadd(const enum CBLAS_ORDER CORDER,const blasint crows, const blasint ccols, const double calpha, const double *a, const blasint clda, const double cbeta,
    double *c, const blasint cldc);
void cblas_cgeadd(const enum CBLAS_ORDER CORDER,const blasint crows, const blasint ccols, const float *calpha, const float *a, const blasint clda, const float *cbeta,
    float *c, const blasint cldc);
void cblas_zgeadd(const enum CBLAS_ORDER CORDER,const blasint crows, const blasint ccols, const double *calpha, const double *a, const blasint clda, const double *cbeta,
    double *c, const blasint cldc);
void cblas_sgemm_batch(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE * TransA_array, const enum CBLAS_TRANSPOSE * TransB_array, const blasint * M_array, const blasint * N_array, const blasint * K_array,
         const float * alpha_array, const float ** A_array, const blasint * lda_array, const float ** B_array, const blasint * ldb_array, const float * beta_array, float ** C_array, const blasint * ldc_array, const blasint group_count, const blasint * group_size);
void cblas_dgemm_batch(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE * TransA_array, const enum CBLAS_TRANSPOSE * TransB_array, const blasint * M_array, const blasint * N_array, const blasint * K_array,
         const double * alpha_array, const double ** A_array, const blasint * lda_array, const double ** B_array, const blasint * ldb_array, const double * beta_array, double ** C_array, const blasint * ldc_array, const blasint group_count, const blasint * group_size);
void cblas_cgemm_batch(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE * TransA_array, const enum CBLAS_TRANSPOSE * TransB_array, const blasint * M_array, const blasint * N_array, const blasint * K_array,
         const void * alpha_array, const void ** A_array, const blasint * lda_array, const void ** B_array, const blasint * ldb_array, const void * beta_array, void ** C_array, const blasint * ldc_array, const blasint group_count, const blasint * group_size);
void cblas_zgemm_batch(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE * TransA_array, const enum CBLAS_TRANSPOSE * TransB_array, const blasint * M_array, const blasint * N_array, const blasint * K_array,
         const void * alpha_array, const void ** A_array, const blasint * lda_array, const void ** B_array, const blasint * ldb_array, const void * beta_array, void ** C_array, const blasint * ldc_array, const blasint group_count, const blasint * group_size);
void cblas_sgemm_batch_strided(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_TRANSPOSE TransB, const blasint M, const blasint N, const blasint K, const float alpha, const float * A, const blasint lda, const blasint stridea, const float * B, const blasint ldb, const blasint strideb, const float beta, float * C, const blasint ldc, const blasint stridec, const blasint group_size);
void cblas_dgemm_batch_strided(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_TRANSPOSE TransB, const blasint M, const blasint N, const blasint K, const double alpha, const double * A, const blasint lda, const blasint stridea, const double * B, const blasint ldb, const blasint strideb, const double beta, double * C, const blasint ldc, const blasint stridec, const blasint group_size);
void cblas_cgemm_batch_strided(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_TRANSPOSE TransB, const blasint M, const blasint N, const blasint K, const void * alpha, const void * A, const blasint lda, const blasint stridea, const void * B, const blasint ldb, const blasint strideb, const void * beta, void * C, const blasint ldc, const blasint stridec, const blasint group_size);
void cblas_zgemm_batch_strided(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_TRANSPOSE TransB, const blasint M, const blasint N, const blasint K, const void * alpha, const void * A, const blasint lda, const blasint stridea, const void * B, const blasint ldb, const blasint strideb, const void * beta, void * C, const blasint ldc, const blasint stridec, const blasint group_size);
void cblas_sbstobf16(const blasint n, const float *in, const blasint incin, bfloat16 *out, const blasint incout);
void cblas_sbdtobf16(const blasint n, const double *in, const blasint incin, bfloat16 *out, const blasint incout);
void cblas_sbf16tos(const blasint n, const bfloat16 *in, const blasint incin, float *out, const blasint incout);
void cblas_dbf16tod(const blasint n, const bfloat16 *in, const blasint incin, double *out, const blasint incout);
void cblas_bgemv(const enum CBLAS_ORDER order, const enum CBLAS_TRANSPOSE trans, const blasint m, const blasint n, const bfloat16 alpha, const bfloat16 *a, const blasint lda, const bfloat16 *x, const blasint incx, const bfloat16 beta, bfloat16 *y, const blasint incy);
float cblas_sbdot(const blasint n, const bfloat16 *x, const blasint incx, const bfloat16 *y, const blasint incy);
void cblas_sbgemv(const enum CBLAS_ORDER order, const enum CBLAS_TRANSPOSE trans, const blasint m, const blasint n, const float alpha, const bfloat16 *a, const blasint lda, const bfloat16 *x, const blasint incx, const float beta, float *y, const blasint incy);
void cblas_bgemm(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_TRANSPOSE TransB, const blasint M, const blasint N, const blasint K,
      const bfloat16 alpha, const bfloat16 *A, const blasint lda, const bfloat16 *B, const blasint ldb, const bfloat16 beta, bfloat16 *C, const blasint ldc);
void cblas_sbgemm(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_TRANSPOSE TransB, const blasint M, const blasint N, const blasint K,
      const float alpha, const bfloat16 *A, const blasint lda, const bfloat16 *B, const blasint ldb, const float beta, float *C, const blasint ldc);
void cblas_sbgemm_batch(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE * TransA_array, const enum CBLAS_TRANSPOSE * TransB_array, const blasint * M_array, const blasint * N_array, const blasint * K_array,
         const float * alpha_array, const bfloat16 ** A_array, const blasint * lda_array, const bfloat16 ** B_array, const blasint * ldb_array, const float * beta_array, float ** C_array, const blasint * ldc_array, const blasint group_count, const blasint * group_size);
void cblas_sbgemm_batch_strided(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_TRANSPOSE TransB, const blasint M, const blasint N, const blasint K, const float alpha, const bfloat16 * A, const blasint lda, const blasint stridea, const bfloat16 * B, const blasint ldb, const blasint strideb, const float beta, float * C, const blasint ldc, const blasint stridec, const blasint group_size);
void cblas_shgemm(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_TRANSPOSE TransB, const blasint M, const blasint N, const blasint K,
      const float alpha, const hfloat16 *A, const blasint lda, const hfloat16 *B, const blasint ldb, const float beta, float *C, const blasint ldc);
struct _r1yxr2d
{ long _k5ogz2u;
    double _q30o30u;
    float _xpbwhvn;
    float _dntcagr;
    float _pkxvhag;
    float _uug95rl;
    float _9dde3cu;
    float _4yfds0d;
    float _gfu8iog;
    float _at9mx2x;
    float _rtj5pht;
    float _ktakpux;
    unsigned int _148p0tv;
};
struct _q0oxfq8
{ float _yj9boq0;
    float _p2u4dtb;
    float _x5qjicj;
    float _e9ftaxp;
    float _6snglh6;
    float _3ac6pe6;
    unsigned int _jp61lfc;
    unsigned int _urvqt3g;
    unsigned int _o31yhxs;
    unsigned int _l9ue421;
    unsigned int _2h2p39y;
    unsigned short int _nwngucz;
    unsigned short int _zk94yxf;
    unsigned short int _qsq0p82;
};
struct _wrdmly4
{ float _n8la100;
    float _lw2lm34;
    float _6fokzix;
    float _jssuuqo;
    float _q01tee1;
    float _8o1bj9r;
    float _rdihsa7;
    float _ihair4c;
    float _owgzrix;
    float _a11dpg8;
    unsigned int _3ru7myo;
    unsigned int _5izolnn;
    unsigned short int _hx6zdcu;
    unsigned short int _fftouhy;
    unsigned short int _aq15o2n;
    unsigned short int _cs4kek4;
    unsigned short int _2qgv3dr;
    unsigned short int _26q3xlr;
};
extern void LoadInputParaFile(int _20n8tte, const char *restrict _sbelrqv, char *restrict _slwf90k, char *restrict _8dalfyo, struct _r1yxr2d *restrict _5zja8dz , struct _q0oxfq8 *restrict _z3nx39r, struct _wrdmly4 *restrict _lyk83tq);
extern void LoadRghStghBnd(int _20n8tte, const char *restrict _8dalfyo, const struct _q0oxfq8 *restrict _z3nx39r, struct _wrdmly4 *restrict _lyk83tq, const int *restrict _b6gs8g8, const unsigned int *restrict _xo100vh, unsigned int *restrict _t4q9bc2, unsigned int *restrict _mgvu4is, float *restrict _qw5dsk7, float *restrict _rcpprcj, float *restrict _z39pecn, float *restrict _8e7flqh, float *restrict _mzwffhh);
extern void AssignRefValues(int _20n8tte, struct _r1yxr2d *restrict _5zja8dz, struct _q0oxfq8 *restrict _z3nx39r, struct _wrdmly4 *restrict _lyk83tq, const int *restrict _b6gs8g8, const int *restrict _zfflb7w, const unsigned int *restrict _xo100vh, const unsigned int *restrict _g1ved32, float *restrict _cmktlzk, float *restrict _ldalmje, unsigned int *restrict _t4q9bc2, unsigned int *restrict _mgvu4is, float *restrict _qw5dsk7, float *restrict _rcpprcj, const float *restrict _z39pecn, const float *restrict _8e7flqh, float *restrict _3n3bzt9, float *restrict _mzwffhh, float *restrict _3p91i3v, float *restrict _m7pyxdw);
extern void MakeNew_Khmatrix(int _avt8k06, int _20n8tte, const char *restrict _8dalfyo, struct _r1yxr2d *_5zja8dz, struct _q0oxfq8 *restrict _z3nx39r, const struct _wrdmly4 *restrict _lyk83tq, const int *restrict _b6gs8g8, const int *restrict _zfflb7w, const unsigned int *restrict _xo100vh, const unsigned int *restrict _g1ved32, const float *restrict _cmktlzk, const float *restrict _ldalmje, const unsigned int *restrict _t4q9bc2, const unsigned int *restrict _mgvu4is, const float *restrict _qw5dsk7, const float *restrict _rcpprcj, const float *restrict _z39pecn, const float *restrict _8e7flqh, unsigned int *restrict _6zjq7w7, float *restrict _3n3bzt9, float *restrict _mzwffhh, float *restrict _m7pyxdw, unsigned short int *restrict _iob0nlw, unsigned short int *restrict _iymyrsv, unsigned short int *restrict _btkf6xd, unsigned short int *restrict _9f77lu8, unsigned short int ***restrict _3b9so0u, unsigned short int ***restrict _1ajo091, unsigned short int ***restrict _86lthet, unsigned short int ***restrict _op767j2, float ***restrict _0dtirc9, float ***restrict _qk4ltah, float ***restrict _h61iaoi, float ***restrict _nf74t60, float ***restrict _mrunjdu, float ***restrict _jbc2sws, float ***restrict _4xkss0k, float ***restrict _t997e4u, float ***restrict _ge2j51b, float ***restrict _sh3tow4, float ***restrict _pqvl00e);
extern void Save_Khmatrix(int _20n8tte, const char *restrict _slwf90k, const struct _q0oxfq8 *restrict _z3nx39r, const struct _wrdmly4 *restrict _lyk83tq, const int *restrict _dtz91x7, const int *restrict _gcsotdh, const int *restrict _b6gs8g8, const int *restrict _zfflb7w, const unsigned int *restrict _xo100vh, const unsigned int *restrict _g1ved32, const float *restrict _3n3bzt9, const float *restrict _m7pyxdw, unsigned short int *restrict _iob0nlw, unsigned short int *restrict _iymyrsv, unsigned short int *restrict _btkf6xd, unsigned short int *restrict _9f77lu8, unsigned short int **restrict _3b9so0u, unsigned short int **restrict _1ajo091, unsigned short int **restrict _86lthet, unsigned short int **restrict _op767j2, float **restrict _0dtirc9, float **restrict _qk4ltah, float **restrict _h61iaoi, float **restrict _nf74t60, float **restrict _mrunjdu, float **restrict _jbc2sws, float **restrict _4xkss0k, float **restrict _t997e4u, float **restrict _ge2j51b, float **restrict _sh3tow4, float **restrict _pqvl00e);
extern void Load_Khmatrix(int _avt8k06, int _20n8tte, const char *restrict _slwf90k, struct _r1yxr2d *restrict _5zja8dz, struct _q0oxfq8 *restrict _z3nx39r, const struct _wrdmly4 *restrict _lyk83tq, const int *restrict _dtz91x7, const int *restrict _gcsotdh, const int *restrict _b6gs8g8, const int *restrict _zfflb7w, const unsigned int *restrict _xo100vh, unsigned int *restrict _6zjq7w7, float *restrict _3n3bzt9, float *restrict _mzwffhh, const unsigned int *restrict _g1ved32, float *restrict _m7pyxdw, unsigned short int *restrict _iob0nlw, unsigned short int *restrict _iymyrsv, unsigned short int *restrict _btkf6xd, unsigned short int *restrict _9f77lu8, unsigned short int ***restrict _3b9so0u, unsigned short int ***restrict _1ajo091, unsigned short int ***restrict _86lthet, unsigned short int ***restrict _op767j2, float ***restrict _0dtirc9, float ***restrict _qk4ltah, float ***restrict _h61iaoi, float ***restrict _nf74t60, float ***restrict _mrunjdu, float ***restrict _jbc2sws, float ***restrict _4xkss0k, float ***restrict _t997e4u, float ***restrict _ge2j51b, float ***restrict _sh3tow4, float ***restrict _pqvl00e);
extern void TectonicLoading(int _avt8k06, int _20n8tte, const struct _r1yxr2d *restrict _5zja8dz, const struct _q0oxfq8 *restrict _z3nx39r, const int *restrict _b6gs8g8, const int *restrict _zfflb7w, const unsigned int *restrict _xo100vh, const unsigned int *restrict _g1ved32, float *restrict _3n3bzt9, float *restrict _3p91i3v, const float *restrict _m7pyxdw, float *restrict _rcpprcj, const unsigned short int *restrict _iob0nlw, const unsigned short int *restrict _iymyrsv, const unsigned short int *restrict _btkf6xd, const unsigned short int *restrict _9f77lu8, const unsigned short int **restrict _3b9so0u, const unsigned short int **restrict _1ajo091, const unsigned short int **restrict _86lthet, const unsigned short int **restrict _op767j2, const float **restrict _0dtirc9, const float **restrict _qk4ltah, const float **restrict _h61iaoi, const float **restrict _nf74t60, const float **restrict _mrunjdu, const float **restrict _jbc2sws, const float **restrict _4xkss0k, const float **restrict _t997e4u, const float **restrict _ge2j51b, const float **restrict _sh3tow4, const float **restrict _pqvl00e);
extern long int GetByteSize(int _6zx1tv1, int _tk63n2q, int _qtkcd0r);
extern void SubtractVect3(float *restrict _zh13cd5, const float *restrict _3lbqgu3, const float *restrict _li076a3);
extern void CrossProduct(float *restrict _bk2l0d1, const float *restrict _3dzeedh, const float *restrict _hjshzkb);
extern float VectorLength(const float *restrict _nvg3kk0);
extern void NrmlzeVect(float *restrict _nvg3kk0);
extern void GetStkAndDipVect(const float *restrict _5bwikot, float *restrict _pt5r5rc, float *restrict _glqj5sx);
extern void RotStressG2L(float *restrict _2292mtr, const float *restrict _088o98b, const float *restrict _kkrzd5e);
extern void ScaleStress6(float *restrict _e1g983e, float _wrb7f6f);
extern float ran0(long *_e05qmkt);
extern float ran1(long *_e05qmkt);
extern void ReAssignElemVals(struct _r1yxr2d *restrict _5zja8dz, const struct _q0oxfq8 *restrict _z3nx39r, unsigned int i, unsigned int *restrict _6zjq7w7, float *restrict _3n3bzt9, float *restrict _mzwffhh, unsigned int _7qqyzqi);
extern void ResetAfterEvent(int _20n8tte, const int *restrict _b6gs8g8, struct _r1yxr2d *restrict _5zja8dz, unsigned int *restrict _6zjq7w7, float *restrict _3n3bzt9, float *restrict _mzwffhh, float *restrict _3p91i3v);
extern void WriteFltComp(int _20n8tte, int _8kgyggs, int _npuevjl, const char *restrict _m1zxyvg, unsigned int _qy8yb7w, const int *restrict _b6gs8g8, const unsigned int *restrict _xo100vh, const unsigned int *restrict _t4q9bc2, const float *restrict _qw5dsk7, const float *restrict _z39pecn, const float *restrict _mzwffhh, const float *restrict _3n3bzt9, const float *restrict _3p91i3v);
extern void WritePreRunData(int _20n8tte, const char *restrict _m1zxyvg, const struct _q0oxfq8 *restrict _z3nx39r, const struct _wrdmly4 *restrict _lyk83tq, const int *restrict _b6gs8g8, const unsigned int *restrict _xo100vh, const unsigned int *restrict _t4q9bc2, const float *restrict _qw5dsk7, const float *restrict _z39pecn, const unsigned int *restrict _6zjq7w7, float *restrict _3n3bzt9, const float *restrict _mzwffhh);
extern void ContinueCatalog(int _20n8tte, const char *restrict _8dalfyo, struct _r1yxr2d *restrict _5zja8dz, const struct _q0oxfq8 *restrict _z3nx39r, const struct _wrdmly4 *restrict _lyk83tq, unsigned int *restrict _jzr25md, double *restrict _txhi15c, const int *restrict _b6gs8g8, const int *restrict _zfflb7w, const unsigned int *restrict _xo100vh, const unsigned int *restrict _g1ved32, const unsigned int *restrict _t4q9bc2, const float *restrict _qw5dsk7, const float *restrict _z39pecn, unsigned int *restrict _6zjq7w7, float *restrict _3n3bzt9, float *restrict _mzwffhh, float *restrict _3p91i3v, float *restrict _m7pyxdw);
int main(int _m9j98cw, char **_n2y8fv2)
{ if ( (_m9j98cw > 3 ) || (_m9j98cw < 2) ) { fprintf(__stdoutp,"Input Error\n Please start the code in the following way:\n\n mpirun -np 4 ./MCQsim25 RunParaFile.txt"); exit(1); }
    int _20n8tte, _avt8k06;
    MPI_Init(&_m9j98cw, &_n2y8fv2);
    MPI_Comm_rank(((MPI_Comm)0x44000000), &_20n8tte);
    MPI_Comm_size(((MPI_Comm)0x44000000), &_avt8k06);
    struct _r1yxr2d __attribute__((aligned(64))) _5zja8dz; __builtin___memset_chk (&_5zja8dz, 0, sizeof(_5zja8dz), __builtin_object_size (&_5zja8dz, 0));
    struct _q0oxfq8 __attribute__((aligned(64))) _z3nx39r; __builtin___memset_chk (&_z3nx39r, 0, sizeof(_z3nx39r), __builtin_object_size (&_z3nx39r, 0));
    struct _wrdmly4 __attribute__((aligned(64))) _lyk83tq; __builtin___memset_chk (&_lyk83tq, 0, sizeof(_lyk83tq), __builtin_object_size (&_lyk83tq, 0));
    int *_0fqldb2 = ((void*)0), *_hlmh1t9 = ((void*)0), *_tacy15p = ((void*)0), *_xrum2lh = ((void*)0);
    unsigned int *_mygty3a = ((void*)0), *_lx2toow = ((void*)0);
    const int *_dtz91x7 = ((void*)0), *_gcsotdh = ((void*)0), *_b6gs8g8 = ((void*)0), *_zfflb7w = ((void*)0);
    const unsigned int *_xo100vh = ((void*)0), *_g1ved32 = ((void*)0);
    unsigned int *_t4q9bc2 = ((void*)0), *_mgvu4is = ((void*)0), *_6zjq7w7 = ((void*)0);
    float *_qw5dsk7 = ((void*)0), *_z39pecn = ((void*)0), *_rcpprcj = ((void*)0), *_8e7flqh = ((void*)0), *_cmktlzk = ((void*)0), *_ldalmje = ((void*)0);
    float *_mzwffhh = ((void*)0), *_3n3bzt9 = ((void*)0), *_m7pyxdw = ((void*)0), *_3p91i3v = ((void*)0);
    unsigned short int *_0vs52sq = ((void*)0), *_zn4regb = ((void*)0), *_n8wwiif = ((void*)0), *_el5015o = ((void*)0);
    unsigned short int **_swbkok3 = ((void*)0), **_c795rlq = ((void*)0), **_25lzse9 = ((void*)0), **_dg1wjui = ((void*)0);
    float **_zd76szn = ((void*)0), **_plozjcn = ((void*)0), **_ozctwjl = ((void*)0), **_3jczucq = ((void*)0);
    float **_ukwqvr1 = ((void*)0), **_0isb5iq = ((void*)0), **_b68rerv = ((void*)0), **_ycz04cd = ((void*)0);
    float **_m458ax4 = ((void*)0), **_ny181ac = ((void*)0), **_zzf8s2z = ((void*)0);
    const unsigned short int *_iob0nlw = ((void*)0), *_iymyrsv = ((void*)0), *_btkf6xd = ((void*)0), *_9f77lu8 = ((void*)0);
    const unsigned short int **_3b9so0u = ((void*)0), **_1ajo091 = ((void*)0), **_86lthet = ((void*)0), **_op767j2 = ((void*)0);
    const float **_0dtirc9 = ((void*)0), **_nf74t60 = ((void*)0), **_jbc2sws = ((void*)0), **_ge2j51b = ((void*)0);
    const float **_qk4ltah = ((void*)0), **_mrunjdu = ((void*)0), **_4xkss0k = ((void*)0), **_sh3tow4 = ((void*)0);
    const float **_h61iaoi = ((void*)0), **_t997e4u = ((void*)0), **_pqvl00e = ((void*)0);
    float **_ntdr5xy = ((void*)0);
    unsigned int *_k00ilu6 = ((void*)0), *_vfua9bq = ((void*)0), *_7c6qc6k = ((void*)0), *_8djzr5s = ((void*)0), *_ddqhkga = ((void*)0);
    float *_28xg8y7= ((void*)0), *_7odkef8 = ((void*)0), *_v71448n = ((void*)0), *_2ivglye = ((void*)0), *_j3swebn= ((void*)0), *_vg5vmf6 = ((void*)0), *_0yknfu8 = ((void*)0), *_4xyjoa2 = ((void*)0), *_795e5iz = ((void*)0), *_vlnknd6 = ((void*)0);
    float *_4pa8hye = ((void*)0), *_hj3muur = ((void*)0), *_qxmcthp = ((void*)0), *_wd1xmxb = ((void*)0), *_z2c9ros = ((void*)0), *_32orwix = ((void*)0);
    const float *_mybd2qe = ((void*)0), *_7yt2m6t = ((void*)0), *_xf7jxjn = ((void*)0), *_9yo18pa = ((void*)0);
    long int _39ayq6w;
    char _sbelrqv[512], _slwf90k[512], _8dalfyo[512];
    __builtin___strcpy_chk (_sbelrqv, _n2y8fv2[1], __builtin_object_size (_sbelrqv, 2 > 1 ? 1 : 0));
    LoadInputParaFile(_20n8tte, _sbelrqv, _slwf90k, _8dalfyo, &_5zja8dz, &_z3nx39r, &_lyk83tq);
    if (_z3nx39r._jp61lfc > 50000)
    { if (_20n8tte == 0) { fprintf(__stdoutp,"\n\nThis version of MCQsim restricts the maximum number of fault elements to %d.\nThe imported model however contains more elements (%d).\nReduce element number to run simulation.\n\n",50000,_z3nx39r._jp61lfc); }
        return 0;
    }
    { unsigned int i;
        int _d7x48wz = (int)(_z3nx39r._jp61lfc/_avt8k06);
        int _8p74mjd = (int)(_z3nx39r._jp61lfc%_avt8k06);
        int _t9nm3v2 = (int)(_z3nx39r._urvqt3g/_avt8k06);
        int _25qbrmq = (int)(_z3nx39r._urvqt3g%_avt8k06);
        _39ayq6w = GetByteSize(64, _avt8k06, sizeof(int));
        _0fqldb2 = aligned_alloc(64, _39ayq6w);
        _hlmh1t9 = aligned_alloc(64, _39ayq6w);
        _tacy15p = aligned_alloc(64, _39ayq6w);
        _xrum2lh = aligned_alloc(64, _39ayq6w);
        _0fqldb2[0] = 0;
        _hlmh1t9[0] = 0;
        for (i = 0u; i < _avt8k06; i++) { _tacy15p[i] = _d7x48wz; }
        for (i = 0u; i < _8p74mjd; i++) { _tacy15p[i] += 1; }
        for (i = 1u; i < _avt8k06; i++) { _0fqldb2[i] = _0fqldb2[i-1] + _tacy15p[i-1]; }
        for (i = 0u; i < _avt8k06; i++) { _xrum2lh[i] = _t9nm3v2; }
        for (i = 0u; i < _25qbrmq; i++) { _xrum2lh[i] += 1; }
        for (i = 1u; i < _avt8k06; i++) { _hlmh1t9[i] = _hlmh1t9[i-1] + _xrum2lh[i-1]; }
        _39ayq6w = GetByteSize(64, _tacy15p[_20n8tte], sizeof(unsigned int));
        _mygty3a = aligned_alloc(64, _39ayq6w);
        for (i = 0u; i < _tacy15p[_20n8tte]; i++) { _mygty3a[i] = i*_avt8k06 +_20n8tte; }
        _dtz91x7 = (const int *)_0fqldb2; _b6gs8g8 = (const int *)_tacy15p; _xo100vh = (const unsigned int *)_mygty3a;
        if (_z3nx39r._urvqt3g > 0u)
        { _39ayq6w = GetByteSize(64, _xrum2lh[_20n8tte], sizeof(unsigned int));
            _lx2toow = aligned_alloc(64, _39ayq6w);
            for (i = 0u; i < _xrum2lh[_20n8tte]; i++) { _lx2toow[i] = i*_avt8k06 +_20n8tte; }
        }
        _gcsotdh = (const int *)_hlmh1t9; _zfflb7w = (const int *)_xrum2lh; _g1ved32 = (const unsigned int *)_lx2toow;
        _39ayq6w = GetByteSize(64, 8*_z3nx39r._jp61lfc, sizeof(unsigned int)); _t4q9bc2 = aligned_alloc(64, _39ayq6w); __builtin___memset_chk (_t4q9bc2, 0, _39ayq6w, __builtin_object_size (_t4q9bc2, 0));
        _39ayq6w = GetByteSize(64, 16*_z3nx39r._jp61lfc, sizeof(float)); _qw5dsk7 = aligned_alloc(64, _39ayq6w); __builtin___memset_chk (_qw5dsk7, 0, _39ayq6w, __builtin_object_size (_qw5dsk7, 0));
        _39ayq6w = GetByteSize(64, 4*_lyk83tq._3ru7myo, sizeof(float)); _z39pecn = aligned_alloc(64, _39ayq6w); __builtin___memset_chk (_z39pecn, 0, _39ayq6w, __builtin_object_size (_z39pecn, 0));
        _39ayq6w = GetByteSize(64, 16*_b6gs8g8[_20n8tte], sizeof(float)); _mzwffhh = aligned_alloc(64, _39ayq6w); __builtin___memset_chk (_mzwffhh, 0, _39ayq6w, __builtin_object_size (_mzwffhh, 0));
        _39ayq6w = GetByteSize(64, 8*_b6gs8g8[_20n8tte], sizeof(float)); _3p91i3v = aligned_alloc(64, _39ayq6w); __builtin___memset_chk (_3p91i3v, 0, _39ayq6w, __builtin_object_size (_3p91i3v, 0));
        _39ayq6w = GetByteSize(64, 4*_b6gs8g8[_20n8tte], sizeof(unsigned int)); _6zjq7w7 = aligned_alloc(64, _39ayq6w); __builtin___memset_chk (_6zjq7w7, 0, _39ayq6w, __builtin_object_size (_6zjq7w7, 0));
        _39ayq6w = GetByteSize(64, 16*_b6gs8g8[_20n8tte], sizeof(float)); _3n3bzt9 = aligned_alloc(64, _39ayq6w); __builtin___memset_chk (_3n3bzt9, 0, _39ayq6w, __builtin_object_size (_3n3bzt9, 0));
        _39ayq6w = GetByteSize(64, _b6gs8g8[_20n8tte], sizeof(unsigned short int)); _0vs52sq = aligned_alloc(64, _39ayq6w); __builtin___memset_chk (_0vs52sq, 0, _39ayq6w, __builtin_object_size (_0vs52sq, 0));
        if (_z3nx39r._urvqt3g > 0u)
        {
            _39ayq6w = GetByteSize(64, 8*_z3nx39r._urvqt3g, sizeof(unsigned int)); _mgvu4is = aligned_alloc(64, _39ayq6w); __builtin___memset_chk (_mgvu4is, 0, _39ayq6w, __builtin_object_size (_mgvu4is, 0));
            _39ayq6w = GetByteSize(64, 16*_z3nx39r._urvqt3g, sizeof(float)); _rcpprcj = aligned_alloc(64, _39ayq6w); __builtin___memset_chk (_rcpprcj, 0, _39ayq6w, __builtin_object_size (_rcpprcj, 0));
            _39ayq6w = GetByteSize(64, 4*_lyk83tq._5izolnn, sizeof(float)); _8e7flqh = aligned_alloc(64, _39ayq6w); __builtin___memset_chk (_8e7flqh, 0, _39ayq6w, __builtin_object_size (_8e7flqh, 0));
            _39ayq6w = GetByteSize(64, 16*_zfflb7w[_20n8tte], sizeof(float)); _m7pyxdw = aligned_alloc(64, _39ayq6w); __builtin___memset_chk (_m7pyxdw, 0, _39ayq6w, __builtin_object_size (_m7pyxdw, 0));
            _39ayq6w = GetByteSize(64, _b6gs8g8[_20n8tte], sizeof(unsigned short int)); _zn4regb = aligned_alloc(64, _39ayq6w); __builtin___memset_chk (_zn4regb, 0, _39ayq6w, __builtin_object_size (_zn4regb, 0));
            _39ayq6w = GetByteSize(64, _zfflb7w[_20n8tte], sizeof(unsigned short int)); _n8wwiif = aligned_alloc(64, _39ayq6w); __builtin___memset_chk (_n8wwiif, 0, _39ayq6w, __builtin_object_size (_n8wwiif, 0));
                                                                                                _el5015o = aligned_alloc(64, _39ayq6w); __builtin___memset_chk (_el5015o, 0, _39ayq6w, __builtin_object_size (_el5015o, 0));
        }
    }
    LoadRghStghBnd(_20n8tte, _8dalfyo, &_z3nx39r, &_lyk83tq, _b6gs8g8, _xo100vh, _t4q9bc2, _mgvu4is, _qw5dsk7, _rcpprcj, _z39pecn, _8e7flqh, _mzwffhh);
    { _39ayq6w = GetByteSize(64, 8*_lyk83tq._hx6zdcu, sizeof(float)); _cmktlzk = (float *) aligned_alloc(64, _39ayq6w); __builtin___memset_chk (_cmktlzk, 0, _39ayq6w, __builtin_object_size (_cmktlzk, 0));
        if (_z3nx39r._urvqt3g > 0u)
        { _39ayq6w = GetByteSize(64, 8*_lyk83tq._fftouhy, sizeof(float)); _ldalmje = (float *) aligned_alloc(64, _39ayq6w); __builtin___memset_chk (_ldalmje, 0, _39ayq6w, __builtin_object_size (_ldalmje, 0));
        }
    }
    AssignRefValues(_20n8tte, &_5zja8dz, &_z3nx39r, &_lyk83tq, _b6gs8g8, _zfflb7w, _xo100vh, _g1ved32, _cmktlzk, _ldalmje, _t4q9bc2, _mgvu4is, _qw5dsk7, _rcpprcj, _z39pecn, _8e7flqh, _3n3bzt9, _mzwffhh, _3p91i3v, _m7pyxdw);
    if (_20n8tte == 0)
    { float _9p3jpxg;
        fprintf(__stdoutp,"\n\n\n------------------------------------------------------------------\n");
        fprintf(__stdoutp,"MCQsim version number:  %u\n",202511);
        fprintf(__stdoutp,"----------------------\n");
        fprintf(__stdoutp,"Number of RANKS: %u\n",_avt8k06);
        fprintf(__stdoutp,"----------------------\n");
        fprintf(__stdoutp,"System info: Byte Size for FLOAT: %lu     unsigned INT: %lu    long INT: %lu\n", sizeof(float), sizeof(unsigned int), sizeof(long int));
        fprintf(__stdoutp,"----------------------\n");
        fprintf(__stdoutp,"FileName:               %s\n",_8dalfyo); fprintf(__stdoutp,"RunNumber:              %u\n",_lyk83tq._aq15o2n);
        fprintf(__stdoutp,"PlotCat2Screen:         %u\n",_z3nx39r._zk94yxf);
        fprintf(__stdoutp,"MinElem4Catalog:        %u\n",_z3nx39r._o31yhxs); fprintf(__stdoutp,"ContPreviousRun:        %u\n",_lyk83tq._cs4kek4);
        fprintf(__stdoutp,"LoadPrev_Kh_mat:        %u\n",_lyk83tq._2qgv3dr); fprintf(__stdoutp,"Kh_mat file name:       %s\n",_slwf90k);
        fprintf(__stdoutp,"----------------------\n");
        fprintf(__stdoutp,"SaveSTF4LargeEQ:        %u\n",_z3nx39r._nwngucz); fprintf(__stdoutp,"MinMag2SaveSTF:         %f\n",_z3nx39r._yj9boq0);
        fprintf(__stdoutp,"UseHalfSpace:           %u\n",_lyk83tq._26q3xlr);
        fprintf(__stdoutp,"----------------------\n");
        fprintf(__stdoutp,"ViscAftSlipTime:        %fyears\n",_z3nx39r._p2u4dtb);
        fprintf(__stdoutp,"----------------------\n");
        fprintf(__stdoutp,"CoSeisHealFract:        %f\n",_z3nx39r._e9ftaxp); fprintf(__stdoutp,"OvershootFract:         %f\n",_z3nx39r._6snglh6);
        fprintf(__stdoutp,"PreStressFract:         %f\n",_lyk83tq._rdihsa7); fprintf(__stdoutp,"MinCoSeisSlipRate(m/s): %f\n",(_z3nx39r._3ac6pe6/_5zja8dz._gfu8iog));
        fprintf(__stdoutp,"----------------------\n");
        fprintf(__stdoutp,"RecLength:              %lf\n",_5zja8dz._q30o30u);
        fprintf(__stdoutp,"----------------------\n");
        fprintf(__stdoutp,"Elastic properties\n");
        fprintf(__stdoutp,"Medium density (kg/m^3):   %e\n",_5zja8dz._xpbwhvn); fprintf(__stdoutp,"AddedNormalStress (MPa):   %e\n",_5zja8dz._dntcagr);
        fprintf(__stdoutp,"ShearModulus (Pa):         %e\n",_5zja8dz._pkxvhag); fprintf(__stdoutp,"PoissonRatio:              %e\n",_5zja8dz._uug95rl);
        fprintf(__stdoutp,"ChangeFricBtwEQs:          %u\n",_z3nx39r._qsq0p82); fprintf(__stdoutp,"Lambda (Pa):               %e\n",_5zja8dz._9dde3cu);
        fprintf(__stdoutp,"S-waveVelocity (m/s):      %e\n",_5zja8dz._4yfds0d);
        fprintf(__stdoutp,"unit slip (m, fault):      %e\n",_lyk83tq._6fokzix);
        fprintf(__stdoutp,"min. slip 2 start EQ (m):  %e\n",_z3nx39r._3ac6pe6); fprintf(__stdoutp,"coseis. timestep (s):      %e\n",_5zja8dz._gfu8iog);
        fprintf(__stdoutp,"------------------------------------------------------------------\n\n");
        _9p3jpxg = 100.0f - expf(-1.0f/_z3nx39r._p2u4dtb)*100.0f;
        fprintf(__stdoutp,"\nFractional post-seismic change during first year of after-slip:      %2.3f%% released \n",_9p3jpxg);
        fprintf(__stdoutp,"FaultPatchNumber %d     FaultVertexNumber %d\n", _z3nx39r._jp61lfc,_lyk83tq._3ru7myo);
       fprintf(__stdoutp,"FaultCent2EdgeLgth and ElementArea: %5.1fm  and %2.4fkm^2\n\n", _lyk83tq._n8la100, _lyk83tq._q01tee1*1.0E-6f);
    }
    int _qwlvunl = 0;
    if (_qwlvunl == 1)
    { char _m1zxyvg[] = "ArrayEntry2Check.txt";
        int _8kgyggs = 4;
        int _npuevjl = 8;
        WriteFltComp(_20n8tte, _8kgyggs, _npuevjl, _m1zxyvg, _z3nx39r._jp61lfc, _b6gs8g8, _xo100vh, _t4q9bc2, _qw5dsk7, _z39pecn, _mzwffhh, _3n3bzt9, _3p91i3v);
    }
    if (_lyk83tq._2qgv3dr == 0)
    {
        MakeNew_Khmatrix(_avt8k06, _20n8tte, _8dalfyo, &_5zja8dz, &_z3nx39r, &_lyk83tq, _b6gs8g8, _zfflb7w, _xo100vh, _g1ved32, _cmktlzk, _ldalmje, _t4q9bc2, _mgvu4is, _qw5dsk7, _rcpprcj, _z39pecn, _8e7flqh, _6zjq7w7, _3n3bzt9, _mzwffhh, _m7pyxdw, _0vs52sq, _zn4regb, _n8wwiif, _el5015o, &_swbkok3, &_c795rlq, &_25lzse9, &_dg1wjui, &_zd76szn, &_ukwqvr1, &_m458ax4, &_plozjcn, &_0isb5iq, &_ozctwjl, &_b68rerv, &_ny181ac, &_3jczucq, &_ycz04cd, &_zzf8s2z);
    }
    else if (_lyk83tq._2qgv3dr == 1)
    {
        MakeNew_Khmatrix(_avt8k06, _20n8tte, _8dalfyo, &_5zja8dz, &_z3nx39r, &_lyk83tq, _b6gs8g8, _zfflb7w, _xo100vh, _g1ved32, _cmktlzk, _ldalmje, _t4q9bc2, _mgvu4is, _qw5dsk7, _rcpprcj, _z39pecn, _8e7flqh, _6zjq7w7, _3n3bzt9, _mzwffhh, _m7pyxdw, _0vs52sq, _zn4regb, _n8wwiif, _el5015o, &_swbkok3, &_c795rlq, &_25lzse9, &_dg1wjui, &_zd76szn, &_ukwqvr1, &_m458ax4, &_plozjcn, &_0isb5iq, &_ozctwjl, &_b68rerv, &_ny181ac, &_3jczucq, &_ycz04cd, &_zzf8s2z);
        Save_Khmatrix(_20n8tte, _slwf90k, &_z3nx39r, &_lyk83tq, _dtz91x7, _gcsotdh, _b6gs8g8, _zfflb7w, _xo100vh, _g1ved32, _3n3bzt9, _m7pyxdw, _0vs52sq, _zn4regb, _n8wwiif, _el5015o, _swbkok3, _c795rlq, _25lzse9, _dg1wjui, _zd76szn, _ukwqvr1, _m458ax4, _plozjcn, _0isb5iq, _ozctwjl, _b68rerv, _ny181ac, _3jczucq, _ycz04cd, _zzf8s2z);
    }
    else if (_lyk83tq._2qgv3dr == 2)
    { Load_Khmatrix(_avt8k06, _20n8tte, _slwf90k, &_5zja8dz, &_z3nx39r, &_lyk83tq, _dtz91x7, _gcsotdh, _b6gs8g8, _zfflb7w, _xo100vh, _6zjq7w7, _3n3bzt9, _mzwffhh, _g1ved32, _m7pyxdw, _0vs52sq, _zn4regb, _n8wwiif, _el5015o, &_swbkok3, &_c795rlq, &_25lzse9, &_dg1wjui, &_zd76szn, &_ukwqvr1, &_m458ax4, &_plozjcn, &_0isb5iq, &_ozctwjl, &_b68rerv, &_ny181ac, &_3jczucq, &_ycz04cd, &_zzf8s2z);
    }
    else if (_lyk83tq._2qgv3dr == 3)
    {
    }
    _iob0nlw = (const unsigned short int *)_0vs52sq; _iymyrsv = (const unsigned short int *)_zn4regb; _btkf6xd = (const unsigned short int *)_n8wwiif; _9f77lu8 = (const unsigned short int *)_el5015o;
    _3b9so0u = (const unsigned short int **)_swbkok3; _1ajo091 = (const unsigned short int **)_c795rlq; _86lthet = (const unsigned short int **)_25lzse9; _op767j2 = (const unsigned short int **)_dg1wjui;
    _0dtirc9 = (const float **)_zd76szn; _qk4ltah = (const float **)_ukwqvr1; _h61iaoi = (const float **)_m458ax4;
    _nf74t60 = (const float **)_plozjcn; _mrunjdu = (const float **)_0isb5iq;
    _jbc2sws = (const float **)_ozctwjl; _4xkss0k = (const float **)_b68rerv; _t997e4u = (const float **)_ny181ac;
    _ge2j51b = (const float **)_3jczucq; _sh3tow4 = (const float **)_ycz04cd; _pqvl00e = (const float **)_zzf8s2z;
    TectonicLoading(_avt8k06, _20n8tte, &_5zja8dz, &_z3nx39r, _b6gs8g8, _zfflb7w, _xo100vh, _g1ved32, _3n3bzt9, _3p91i3v, _m7pyxdw, _rcpprcj, _iob0nlw, _iymyrsv, _btkf6xd, _9f77lu8, _3b9so0u, _1ajo091, _86lthet, _op767j2, _0dtirc9, _qk4ltah, _h61iaoi, _nf74t60, _mrunjdu, _jbc2sws, _4xkss0k, _t997e4u, _ge2j51b, _sh3tow4, _pqvl00e);
    WritePreRunData(_20n8tte,_8dalfyo, &_z3nx39r, &_lyk83tq, _b6gs8g8, _xo100vh, _t4q9bc2, _qw5dsk7, _z39pecn, _6zjq7w7, _3n3bzt9, _mzwffhh);
    float _qijc1yd = 3.40282347e+38F;
    float _53j2d54 = 1.0E-6f;
    { unsigned int i;
        float _0i9g5dz, _obyghgk, _9p3jpxg[3] = {0.0f, 0.0f, 0.0f};
        for (i = 0u; i < _b6gs8g8[_20n8tte]; i++)
        { _0i9g5dz = (-1.0f*_3n3bzt9[i*16 +2]*_3n3bzt9[i*16 +3]) *_53j2d54;
            _obyghgk = _0i9g5dz /(sqrtf(_3n3bzt9[i*16 +9]*_3n3bzt9[i*16 +9] + _3n3bzt9[i*16 +10]*_3n3bzt9[i*16 +10])) *365.25f;
            _qijc1yd = ( ((_obyghgk) < (_qijc1yd)) *(_obyghgk) + ((_obyghgk) >= (_qijc1yd)) *(_qijc1yd) );
            _9p3jpxg[0] = (_6zjq7w7[i*4 +0] == 1u)*(_9p3jpxg[0] +1.0f) + (_6zjq7w7[i*4 +0] != 1u)*_9p3jpxg[0];
            _9p3jpxg[1] = (_6zjq7w7[i*4 +0] == 2u)*(_9p3jpxg[1] +1.0f) + (_6zjq7w7[i*4 +0] != 2u)*_9p3jpxg[1];
            _9p3jpxg[2] = (_6zjq7w7[i*4 +0] == 3u)*(_9p3jpxg[2] +1.0f) + (_6zjq7w7[i*4 +0] != 3u)*_9p3jpxg[2];
        }
        MPI_Allreduce((void *) -1, _9p3jpxg, 3, ((MPI_Datatype)0x4c00040a), (MPI_Op)(0x58000003), ((MPI_Comm)0x44000000));
        MPI_Allreduce((void *) -1, &_qijc1yd, 1, ((MPI_Datatype)0x4c00040a), (MPI_Op)(0x58000002), ((MPI_Comm)0x44000000));
        _9p3jpxg[0] /= (float)_z3nx39r._jp61lfc; _9p3jpxg[1] /= (float)_z3nx39r._jp61lfc; _9p3jpxg[2] /= (float)_z3nx39r._jp61lfc;
        if (_20n8tte == 0) { fprintf(__stdoutp,"\nStability at catalog start:  unstable %2.2f%%     cond. stable %2.2f%%     stable %2.2f%%\n\n",_9p3jpxg[0]*100.0f, _9p3jpxg[1]*100.0f, _9p3jpxg[2]*100.0f); }
        _qijc1yd *= 5.0f;
    }
    FILE *_szbffzz;
    char _y9v20ya[512], _m3tzk3n[512];
    __builtin___strcpy_chk (_y9v20ya, _8dalfyo, __builtin_object_size (_y9v20ya, 2 > 1 ? 1 : 0)); __builtin___snprintf_chk (_m3tzk3n, sizeof(_m3tzk3n), 0, __builtin_object_size (_m3tzk3n, 2 > 1 ? 1 : 0), "_%u.cat",_lyk83tq._aq15o2n); __builtin___strcat_chk (_y9v20ya, _m3tzk3n, __builtin_object_size (_y9v20ya, 2 > 1 ? 1 : 0));
    unsigned int _jzr25md = 0u;
    double _txhi15c = 0.0;
    ContinueCatalog(_20n8tte, _8dalfyo, &_5zja8dz, &_z3nx39r, &_lyk83tq, &_jzr25md, &_txhi15c, _b6gs8g8, _zfflb7w, _xo100vh, _g1ved32, _t4q9bc2, _qw5dsk7, _z39pecn, _6zjq7w7, _3n3bzt9, _mzwffhh, _3p91i3v, _m7pyxdw);
    _5zja8dz._q30o30u += _txhi15c;
    MPI_Barrier( ((MPI_Comm)0x44000000) );
    if (_20n8tte == 0)
    { if ((_szbffzz = fopen(_y9v20ya,"rb+")) == ((void*)0)) { perror("Error -cant open/write RawCatalogFile...\n"); exit(1); }
    }
    free(_t4q9bc2); free(_z39pecn); free(_qw5dsk7); free(_cmktlzk);
    free(_mgvu4is); free(_8e7flqh); free(_rcpprcj); free(_ldalmje);
    MPI_Status _py6b2kn;
    MPI_File _qx8ertb;
    unsigned int i, j, k, _ufwbucj, _g8vuazx, _fu7bfxq, _7bzegbl, _ybxkpm2, _jg2y8rg, _rgb93wk, _kb700uc, _yt91qe3, _v4zwl3z, _f7fr91i, _blmntfk;
    float _ckeh5sw, _wyx9trf, _d1vn8lf, _8sryo4r, _0ljch2v;
    double _qzcucts = _txhi15c;
    double _uflsehu = _txhi15c;
    _0ljch2v = -3.40282347e+38F;
    double _0s1f537 = 0.0, _q8hke22 = 0.0, _1zlbori = 0.0, _ds3qsqc = 0.0, _bo119xp;
    struct timespec _7dknncp, _6nouhfj, _y79kvls, _rg2gzr7, _tkezjss, _7upfln4, _zwd0w9r, _rxyj38x;
    int __attribute__((aligned(64))) _gg24479[_avt8k06]; __builtin___memset_chk (_gg24479, 0, _avt8k06*sizeof(int), __builtin_object_size (_gg24479, 0));
    int __attribute__((aligned(64))) _o2t8fy6[_avt8k06]; __builtin___memset_chk (_o2t8fy6, 0, _avt8k06*sizeof(int), __builtin_object_size (_o2t8fy6, 0));
    int __attribute__((aligned(64))) _moxsfh4[_avt8k06]; __builtin___memset_chk (_moxsfh4, 0, _avt8k06*sizeof(int), __builtin_object_size (_moxsfh4, 0));
    int __attribute__((aligned(64))) _1poxrte[_avt8k06]; __builtin___memset_chk (_1poxrte, 0, _avt8k06*sizeof(int), __builtin_object_size (_1poxrte, 0));
    long long int __attribute__((aligned(64))) _ti22nv0[_avt8k06]; __builtin___memset_chk (_ti22nv0, 0, _avt8k06*sizeof(long long int), __builtin_object_size (_ti22nv0, 0));
    long long int __attribute__((aligned(64))) _1zderi4[_avt8k06]; __builtin___memset_chk (_1zderi4, 0, _avt8k06*sizeof(long long int), __builtin_object_size (_1zderi4, 0));
    unsigned int __attribute__((aligned(64))) _h9ier7u[16];
    float __attribute__((aligned(64))) _9p3jpxg[16];
    float __attribute__((aligned(64))) _0wfoy4d[3];
    float __attribute__((aligned(64))) _lqbmae3[6];
    { _39ayq6w = GetByteSize(64, (_b6gs8g8[_20n8tte]+1), sizeof(unsigned int)); _ddqhkga = aligned_alloc(64, _39ayq6w); __builtin___memset_chk (_ddqhkga, 0, _39ayq6w, __builtin_object_size (_ddqhkga, 0));
        _39ayq6w = GetByteSize(64, 4*_b6gs8g8[_20n8tte], sizeof(float)); _795e5iz = aligned_alloc(64,_39ayq6w); __builtin___memset_chk (_795e5iz, 0, _39ayq6w, __builtin_object_size (_795e5iz, 0));
        _39ayq6w = GetByteSize(64, 3*_b6gs8g8[_20n8tte], sizeof(float)); _4pa8hye = aligned_alloc(64, _39ayq6w); __builtin___memset_chk (_4pa8hye, 0, _39ayq6w, __builtin_object_size (_4pa8hye, 0));
        _39ayq6w = GetByteSize(64, 3*_z3nx39r._jp61lfc, sizeof(float)); _hj3muur = aligned_alloc(64, _39ayq6w); __builtin___memset_chk (_hj3muur, 0, _39ayq6w, __builtin_object_size (_hj3muur, 0));
        _39ayq6w = GetByteSize(64, 2*_z3nx39r._l9ue421, sizeof(float)); _qxmcthp = aligned_alloc(64, _39ayq6w); __builtin___memset_chk (_qxmcthp, 0, _39ayq6w, __builtin_object_size (_qxmcthp, 0));
        _39ayq6w = GetByteSize(64, _b6gs8g8[_20n8tte], sizeof(unsigned int)); _k00ilu6 = aligned_alloc(64, _39ayq6w); __builtin___memset_chk (_k00ilu6, 0, _39ayq6w, __builtin_object_size (_k00ilu6, 0));
                                                                                            _vfua9bq = aligned_alloc(64, _39ayq6w); __builtin___memset_chk (_vfua9bq, 0, _39ayq6w, __builtin_object_size (_vfua9bq, 0));
        _39ayq6w = GetByteSize(64, _b6gs8g8[_20n8tte], sizeof(float)); _28xg8y7 = aligned_alloc(64, _39ayq6w); __builtin___memset_chk (_28xg8y7, 0, _39ayq6w, __builtin_object_size (_28xg8y7, 0));
                                                                                            _7odkef8 = aligned_alloc(64, _39ayq6w); __builtin___memset_chk (_7odkef8, 0, _39ayq6w, __builtin_object_size (_7odkef8, 0));
                                                                                            _v71448n = aligned_alloc(64, _39ayq6w); __builtin___memset_chk (_v71448n, 0, _39ayq6w, __builtin_object_size (_v71448n, 0));
                                                                                            _2ivglye = aligned_alloc(64, _39ayq6w); __builtin___memset_chk (_2ivglye, 0, _39ayq6w, __builtin_object_size (_2ivglye, 0));
        _39ayq6w = GetByteSize(64, _z3nx39r._jp61lfc, sizeof(unsigned int)); _7c6qc6k= aligned_alloc(64, _39ayq6w); __builtin___memset_chk (_7c6qc6k, 0, _39ayq6w, __builtin_object_size (_7c6qc6k, 0));
                                                                                            _8djzr5s= aligned_alloc(64, _39ayq6w); __builtin___memset_chk (_8djzr5s, 0, _39ayq6w, __builtin_object_size (_8djzr5s, 0));
        _39ayq6w = GetByteSize(64, _z3nx39r._jp61lfc, sizeof(float)); _j3swebn= aligned_alloc(64, _39ayq6w); __builtin___memset_chk (_j3swebn, 0, _39ayq6w, __builtin_object_size (_j3swebn, 0));
                                                                                            _vg5vmf6= aligned_alloc(64, _39ayq6w); __builtin___memset_chk (_vg5vmf6, 0, _39ayq6w, __builtin_object_size (_vg5vmf6, 0));
                                                                                            _0yknfu8= aligned_alloc(64, _39ayq6w); __builtin___memset_chk (_0yknfu8, 0, _39ayq6w, __builtin_object_size (_0yknfu8, 0));
                                                                                            _4xyjoa2= aligned_alloc(64, _39ayq6w); __builtin___memset_chk (_4xyjoa2, 0, _39ayq6w, __builtin_object_size (_4xyjoa2, 0));
        _39ayq6w = GetByteSize(64, _b6gs8g8[_20n8tte], sizeof(float *)); _ntdr5xy = aligned_alloc(64, _39ayq6w);
        _39ayq6w = GetByteSize(64, 2*_5zja8dz._148p0tv, sizeof(float));
        for (i = 0u; i < _b6gs8g8[_20n8tte]; i++) { _ntdr5xy[i] = aligned_alloc(64, _39ayq6w); __builtin___memset_chk (_ntdr5xy[i], 0, _39ayq6w, __builtin_object_size (_ntdr5xy[i], 0)); }
        if (_z3nx39r._urvqt3g > 0u)
        { _39ayq6w = GetByteSize(64, 3*_zfflb7w[_20n8tte], sizeof(float)); _vlnknd6 = aligned_alloc(64, _39ayq6w); __builtin___memset_chk (_vlnknd6, 0, _39ayq6w, __builtin_object_size (_vlnknd6, 0));
            _39ayq6w = GetByteSize(64, 4*_zfflb7w[_20n8tte], sizeof(float)); _wd1xmxb = aligned_alloc(64, _39ayq6w); __builtin___memset_chk (_wd1xmxb, 0, _39ayq6w, __builtin_object_size (_wd1xmxb, 0));
            _39ayq6w = GetByteSize(64, 4*_z3nx39r._urvqt3g, sizeof(float)); _z2c9ros = aligned_alloc(64, _39ayq6w); __builtin___memset_chk (_z2c9ros, 0, _39ayq6w, __builtin_object_size (_z2c9ros, 0));
            _39ayq6w = GetByteSize(64, 3*_z3nx39r._2h2p39y, sizeof(float)); _32orwix = aligned_alloc(64, _39ayq6w); __builtin___memset_chk (_32orwix, 0, _39ayq6w, __builtin_object_size (_32orwix, 0));
        }
        _mybd2qe = (const float *)_hj3muur; _7yt2m6t = (const float *)_z2c9ros;
        _xf7jxjn = (const float *)_qxmcthp; _9yo18pa = (const float *)_32orwix;
    }
    if (_20n8tte == 0) { fprintf(__stdoutp,"Starting the catalog....\n"); }
    for (i = 0u; i < _zfflb7w[_20n8tte]; i++) { _m7pyxdw[i*16 +0] = 0.0f; _m7pyxdw[i*16 +1] = 0.0f; _m7pyxdw[i*16 +2] = 0.0f; }
    clock_gettime(_CLOCK_REALTIME, &_7dknncp);
    unsigned int _w2pcn9t = 0u;
    unsigned int _zfk8ut0 = 0u;
    unsigned int _dswm31y = 1u;
    while (_qzcucts <= _5zja8dz._q30o30u)
    {
        clock_gettime(_CLOCK_REALTIME, &_y79kvls);
        for (i = 0u; i < _b6gs8g8[_20n8tte]; i++) { __builtin___memcpy_chk (&_795e5iz[i*4 +0], &_3n3bzt9[i*16 +0], 4*sizeof(float), __builtin_object_size (&_795e5iz[i*4 +0], 0)); }
        for (i = 0u; i < _zfflb7w[_20n8tte]; i++) { __builtin___memcpy_chk (&_vlnknd6[i*3 +0], &_m7pyxdw[i*16 +0], 3*sizeof(float), __builtin_object_size (&_vlnknd6[i*3 +0], 0)); }
        _ckeh5sw = _qijc1yd;
        k = 1;
        _g8vuazx = 0;
        while (k == 1)
        { _g8vuazx++;
            _9p3jpxg[12] = ( ((0.0f) > (expf(-1.0f*_ckeh5sw/_z3nx39r._p2u4dtb))) *(0.0f) + ((0.0f) <= (expf(-1.0f*_ckeh5sw/_z3nx39r._p2u4dtb))) *(expf(-1.0f*_ckeh5sw/_z3nx39r._p2u4dtb)) );
            _9p3jpxg[13] = ( ((0.0f) > (expf(-1.0f*_ckeh5sw/_z3nx39r._x5qjicj))) *(0.0f) + ((0.0f) <= (expf(-1.0f*_ckeh5sw/_z3nx39r._x5qjicj))) *(expf(-1.0f*_ckeh5sw/_z3nx39r._x5qjicj)) );
            _9p3jpxg[14] = (_5zja8dz._pkxvhag/(1.0f*_5zja8dz._4yfds0d)) *(_z3nx39r._3ac6pe6/_5zja8dz._gfu8iog);
            _9p3jpxg[15] = 3.40282347e+38F;
            for (i = 0u; i < _b6gs8g8[_20n8tte]; i++)
            { _bo119xp = (_qzcucts + (double)_ckeh5sw)*(double)_3p91i3v[i*8 +5];
                _3n3bzt9[i*16 +0] = (_3p91i3v[i*8 +3] <= _bo119xp)*(_795e5iz[i*4 +0] +_ckeh5sw*_3n3bzt9[i*16 +9]) + (_3p91i3v[i*8 +3] > _bo119xp)*_795e5iz[i*4 +0];
                _3n3bzt9[i*16 +1] = (_3p91i3v[i*8 +3] <= _bo119xp)*(_795e5iz[i*4 +1] +_ckeh5sw*_3n3bzt9[i*16 +10]) + (_3p91i3v[i*8 +3] > _bo119xp)*_795e5iz[i*4 +1];
                _3n3bzt9[i*16 +2] = _795e5iz[i*4 +2];
                _3n3bzt9[i*16 +3] = _795e5iz[i*4 +3] - (1.0f - _9p3jpxg[12])*_mzwffhh[i*16 +7];
                _3n3bzt9[i*16 +15] = 0.0f;
            }
            for (i = 0u; i < _zfflb7w[_20n8tte]; i++)
            { __builtin___memcpy_chk (&_m7pyxdw[i*16 +0], &_vlnknd6[i*3 +0], 3*sizeof(float), __builtin_object_size (&_m7pyxdw[i*16 +0], 0));
            }
            for (_ufwbucj = 0; _ufwbucj < _dswm31y; _ufwbucj++)
            {
                _ybxkpm2 = 0u;
                _rgb93wk = 0u;
                for (i = 0u; i < _b6gs8g8[_20n8tte]; i++)
                {
                    _9p3jpxg[0] = sqrtf( _3n3bzt9[i*16 +0]*_3n3bzt9[i*16 +0] + _3n3bzt9[i*16 +1]*_3n3bzt9[i*16 +1] );
                    _9p3jpxg[0] = (_6zjq7w7[i*4 +0] == 1u)*0.0f + (_6zjq7w7[i*4 +0] != 1u)*_9p3jpxg[0];
                    _9p3jpxg[1] = _3n3bzt9[i*16 +3] *-1.0f*_3n3bzt9[i*16 +2] ;
                    _9p3jpxg[3] = ( ((0.0f) > ((_9p3jpxg[0] -_9p3jpxg[1]))) *(0.0f) + ((0.0f) <= ((_9p3jpxg[0] -_9p3jpxg[1]))) *((_9p3jpxg[0] -_9p3jpxg[1])) );
                    _9p3jpxg[0] = (_9p3jpxg[0] <= 0.0f)*-1.0f + (_9p3jpxg[0] > 0.0f)*_9p3jpxg[0];
                    _9p3jpxg[1] = -1.0f*(((_9p3jpxg[3]/_9p3jpxg[0])*_3n3bzt9[i*16 +0])/_3n3bzt9[i*16 +5]);
                    _9p3jpxg[2] = -1.0f*(((_9p3jpxg[3]/_9p3jpxg[0])*_3n3bzt9[i*16 +1])/_3n3bzt9[i*16 +6]);
                    _3n3bzt9[i*16 +15] = (_9p3jpxg[3] <= _5zja8dz._at9mx2x)*_3n3bzt9[i*16 +15] + (_9p3jpxg[3] > _5zja8dz._at9mx2x)*(_3n3bzt9[i*16 +15] + sqrtf(_9p3jpxg[1]*_9p3jpxg[1] + _9p3jpxg[2]*_9p3jpxg[2]));
                    _4pa8hye[_ybxkpm2*3 +0] = (_9p3jpxg[3] <= _5zja8dz._at9mx2x)*_4pa8hye[_ybxkpm2*3 +0] + (_9p3jpxg[3] > _5zja8dz._at9mx2x)*(float)(_xo100vh[i]);
                    _4pa8hye[_ybxkpm2*3 +1] = (_9p3jpxg[3] <= _5zja8dz._at9mx2x)*_4pa8hye[_ybxkpm2*3 +1] + (_9p3jpxg[3] > _5zja8dz._at9mx2x)*_9p3jpxg[1]*_3n3bzt9[i*16 +8];
                    _4pa8hye[_ybxkpm2*3 +2] = (_9p3jpxg[3] <= _5zja8dz._at9mx2x)*_4pa8hye[_ybxkpm2*3 +2] + (_9p3jpxg[3] > _5zja8dz._at9mx2x)*_9p3jpxg[2]*_3n3bzt9[i*16 +8];
                    _ybxkpm2 = (_9p3jpxg[3] <= _5zja8dz._at9mx2x)*_ybxkpm2 + (_9p3jpxg[3] > _5zja8dz._at9mx2x)*(_ybxkpm2 + 1u);
                }
                _jg2y8rg = 0u;
                _kb700uc = 0u;
                for (i = 0u; i < _zfflb7w[_20n8tte]; i++)
                { _9p3jpxg[0] = sqrtf( _m7pyxdw[i*16 +0]*_m7pyxdw[i*16 +0] + _m7pyxdw[i*16 +1]*_m7pyxdw[i*16 +1] );
                    _9p3jpxg[1] = fabsf(_m7pyxdw[i*16 +2]);
                    _9p3jpxg[2] = ( ((0.0f) > ((1.0f - _9p3jpxg[13])*_9p3jpxg[0])) *(0.0f) + ((0.0f) <= ((1.0f - _9p3jpxg[13])*_9p3jpxg[0])) *((1.0f - _9p3jpxg[13])*_9p3jpxg[0]) );
                    _9p3jpxg[3] = ( ((0.0f) > ((1.0f - _9p3jpxg[13])*_9p3jpxg[1])) *(0.0f) + ((0.0f) <= ((1.0f - _9p3jpxg[13])*_9p3jpxg[1])) *((1.0f - _9p3jpxg[13])*_9p3jpxg[1]) );
                    _9p3jpxg[4] = sqrtf(_9p3jpxg[2]*_9p3jpxg[2] + _9p3jpxg[3]*_9p3jpxg[3]);
                    _9p3jpxg[0] = (_9p3jpxg[0] <= 0.0f)*-1.0f + (_9p3jpxg[0] > 0.0f)*_9p3jpxg[0];
                    _9p3jpxg[1] = (_9p3jpxg[1] <= 0.0f)*-1.0f + (_9p3jpxg[1] > 0.0f)*_9p3jpxg[1];
                    _9p3jpxg[5] = (_9p3jpxg[0] <= 0.0f)*0.0f + (_9p3jpxg[0] > 0.0f)*(-1.0f*((_9p3jpxg[2]/_9p3jpxg[0])*_m7pyxdw[i*16 +0])/_m7pyxdw[i*16 +5]);
                    _9p3jpxg[6] = (_9p3jpxg[0] <= 0.0f)*0.0f + (_9p3jpxg[0] > 0.0f)*(-1.0f*((_9p3jpxg[2]/_9p3jpxg[0])*_m7pyxdw[i*16 +1])/_m7pyxdw[i*16 +6]);
                    _9p3jpxg[7] = (_9p3jpxg[1] <= 0.0f)*0.0f + (_9p3jpxg[1] > 0.0f)*(-1.0f*((_9p3jpxg[3]/_9p3jpxg[1])*_m7pyxdw[i*16 +2])/_m7pyxdw[i*16 +7]);
                    _wd1xmxb[_jg2y8rg*4 +0] = (_9p3jpxg[4] <= _5zja8dz._at9mx2x)*_wd1xmxb[_jg2y8rg*4 +0] + (_9p3jpxg[4] > _5zja8dz._at9mx2x)*(float)(_g1ved32[i]);
                    _wd1xmxb[_jg2y8rg*4 +1] = (_9p3jpxg[4] <= _5zja8dz._at9mx2x)*_wd1xmxb[_jg2y8rg*4 +1] + (_9p3jpxg[4] > _5zja8dz._at9mx2x)*_9p3jpxg[5]*_m7pyxdw[i*16 +8];
                    _wd1xmxb[_jg2y8rg*4 +2] = (_9p3jpxg[4] <= _5zja8dz._at9mx2x)*_wd1xmxb[_jg2y8rg*4 +2] + (_9p3jpxg[4] > _5zja8dz._at9mx2x)*_9p3jpxg[6]*_m7pyxdw[i*16 +8];
                    _wd1xmxb[_jg2y8rg*4 +3] = (_9p3jpxg[4] <= _5zja8dz._at9mx2x)*_wd1xmxb[_jg2y8rg*4 +3] + (_9p3jpxg[4] > _5zja8dz._at9mx2x)*_9p3jpxg[7]*_m7pyxdw[i*16 +8];
                    _jg2y8rg = (_9p3jpxg[4] <= _5zja8dz._at9mx2x)*_jg2y8rg + (_9p3jpxg[4] > _5zja8dz._at9mx2x)*(_jg2y8rg + 1u);
                }
                MPI_Allgather(&_ybxkpm2, 1, ((MPI_Datatype)0x4c000406), _o2t8fy6, 1, ((MPI_Datatype)0x4c000405), ((MPI_Comm)0x44000000));
                _gg24479[0] = 0;
                _rgb93wk = _o2t8fy6[0];
                for (i = 1; i < _avt8k06; i++) { _rgb93wk += _o2t8fy6[i]; _o2t8fy6[i-1] *= 3; _gg24479[i] = _gg24479[i-1] + _o2t8fy6[i-1]; }
                _o2t8fy6[_avt8k06-1] *= 3;
                MPI_Allgatherv(_4pa8hye, _o2t8fy6[_20n8tte], ((MPI_Datatype)0x4c00040a), _hj3muur, _o2t8fy6, _gg24479, ((MPI_Datatype)0x4c00040a), ((MPI_Comm)0x44000000));
                if (_z3nx39r._urvqt3g > 0u)
                { MPI_Allgather(&_jg2y8rg, 1, ((MPI_Datatype)0x4c000406), _1poxrte, 1, ((MPI_Datatype)0x4c000405), ((MPI_Comm)0x44000000));
                    _moxsfh4[0] = 0u;
                    _kb700uc = _1poxrte[0];
                    for (i = 1; i < _avt8k06; i++) { _kb700uc += _1poxrte[i]; _1poxrte[i-1] *= 4; _moxsfh4[i] = _moxsfh4[i-1] + _1poxrte[i-1]; }
                    _1poxrte[_avt8k06-1] *= 4;
                    MPI_Allgatherv(_wd1xmxb, _1poxrte[_20n8tte], ((MPI_Datatype)0x4c00040a), _z2c9ros, _1poxrte, _moxsfh4, ((MPI_Datatype)0x4c00040a), ((MPI_Comm)0x44000000));
                }
                for (i = 0u; i < _b6gs8g8[_20n8tte]; i++)
                {
                    if (_z3nx39r._urvqt3g > 0u)
                    { _h9ier7u[0] = (_kb700uc <= _iymyrsv[i])*1u + (_kb700uc > _iymyrsv[i])*2u;
                        _h9ier7u[0] = (_kb700uc <= 0u)*0u + (_kb700uc > 0u)*_h9ier7u[0];
                        if (_h9ier7u[0] == 1u)
                        { _9p3jpxg[0] = 0.0f; _9p3jpxg[1] = 0.0f;
                            for (j = 0u; j < _kb700uc; j++)
                            { _h9ier7u[2] = (unsigned int)_7yt2m6t[j*4 +0];
                                _h9ier7u[1] = _1ajo091[i][_h9ier7u[2]] *3;
                                _h9ier7u[2] = j*4;
                                _9p3jpxg[0] += (_7yt2m6t[_h9ier7u[2] +1] *_nf74t60[i][_h9ier7u[1] +0] + _7yt2m6t[_h9ier7u[2] +2] *_nf74t60[i][_h9ier7u[1] +1] + _7yt2m6t[_h9ier7u[2]+3] *_nf74t60[i][_h9ier7u[1] +2]);
                                _9p3jpxg[1] += (_7yt2m6t[_h9ier7u[2] +1] *_mrunjdu[i][_h9ier7u[1] +0] + _7yt2m6t[_h9ier7u[2] +2] *_mrunjdu[i][_h9ier7u[1] +1] + _7yt2m6t[_h9ier7u[2]+3] *_mrunjdu[i][_h9ier7u[1] +2]);
                            }
                            _3n3bzt9[i*16 +0] += _9p3jpxg[0];
                            _3n3bzt9[i*16 +1] += _9p3jpxg[1];
                        }
                        else if (_h9ier7u[0] == 2u)
                        { __builtin___memset_chk (_32orwix, 0, 3*_iymyrsv[i]*sizeof(float), __builtin_object_size (_32orwix, 0));
                            for (j = 0u; j < _kb700uc; j++)
                            { _h9ier7u[2] = (unsigned int)_7yt2m6t[j*4 +0];
                                _h9ier7u[1] = _1ajo091[i][_h9ier7u[2]] *3;
                                _h9ier7u[2] = j*4;
                                _32orwix[_h9ier7u[1] +0] += _7yt2m6t[_h9ier7u[2] +1];
                                _32orwix[_h9ier7u[1] +1] += _7yt2m6t[_h9ier7u[2] +2];
                                _32orwix[_h9ier7u[1] +2] += _7yt2m6t[_h9ier7u[2] +3];
                            }
                            _3n3bzt9[i*16 +0] += cblas_sdot(3*_iymyrsv[i], _9yo18pa, 1, _nf74t60[i], 1);
                            _3n3bzt9[i*16 +1] += cblas_sdot(3*_iymyrsv[i], _9yo18pa, 1, _mrunjdu[i], 1);
                        }
                    }
                    _h9ier7u[0] = (_rgb93wk <= _iob0nlw[i])*1u + (_rgb93wk > _iob0nlw[i])*2u;
                    _h9ier7u[0] = (_rgb93wk <= 0u)*0u + (_rgb93wk > 0u)*_h9ier7u[0];
                    if (_h9ier7u[0] == 1u)
                    { _9p3jpxg[0] = 0.0f; _9p3jpxg[1] = 0.0f;
                        for (j = 0u; j < _rgb93wk; j++)
                        { _h9ier7u[2] = (unsigned int)_mybd2qe[j*3 +0];
                            _h9ier7u[1] = _3b9so0u[i][_h9ier7u[2]] *2;
                            _h9ier7u[2] = j*3;
                            _9p3jpxg[0] += (_mybd2qe[_h9ier7u[2] +1]*_0dtirc9[i][_h9ier7u[1] +0] + _mybd2qe[_h9ier7u[2] +2]*_0dtirc9[i][_h9ier7u[1] +1] );
                            _9p3jpxg[1] += (_mybd2qe[_h9ier7u[2] +1]*_qk4ltah[i][_h9ier7u[1] +0] + _mybd2qe[_h9ier7u[2] +2]*_qk4ltah[i][_h9ier7u[1] +1] );
                        }
                        _3n3bzt9[i*16 +0] += _9p3jpxg[0];
                        _3n3bzt9[i*16 +1] += _9p3jpxg[1];
                    }
                    else if (_h9ier7u[0] == 2u)
                    { __builtin___memset_chk (_qxmcthp, 0, 2*_iob0nlw[i]*sizeof(float), __builtin_object_size (_qxmcthp, 0));
                        for (j = 0u; j < _rgb93wk; j++)
                        { _h9ier7u[2] = (unsigned int)_mybd2qe[j*3 +0];
                            _h9ier7u[1] = _3b9so0u[i][_h9ier7u[2]] *2;
                            _h9ier7u[2] = j*3;
                            _qxmcthp[_h9ier7u[1] +0] += _mybd2qe[_h9ier7u[2] +1];
                            _qxmcthp[_h9ier7u[1] +1] += _mybd2qe[_h9ier7u[2] +2];
                        }
                        _3n3bzt9[i*16 +0] += cblas_sdot(2*_iob0nlw[i], _xf7jxjn, 1, _0dtirc9[i], 1);
                        _3n3bzt9[i*16 +1] += cblas_sdot(2*_iob0nlw[i], _xf7jxjn, 1, _qk4ltah[i], 1);
                    }
                }
                for (i = 0u; i < _zfflb7w[_20n8tte]; i++)
                { _h9ier7u[0] = (_kb700uc <= _btkf6xd[i])*1u + (_kb700uc > _btkf6xd[i])*2u;
                    _h9ier7u[0] = (_kb700uc <= 0u)*0u + (_kb700uc > 0u)*_h9ier7u[0];
                    if ( _h9ier7u[0] == 1)
                    { _9p3jpxg[0] = 0.0f; _9p3jpxg[1] = 0.0f; _9p3jpxg[2] = 0.0f;
                        for (j = 0u; j < _kb700uc; j++)
                        { _h9ier7u[2] = (unsigned int)_7yt2m6t[j*4 +0];
                            _h9ier7u[1] = _86lthet[i][_h9ier7u[2]] *3;
                            _h9ier7u[2] = j*4;
                            _9p3jpxg[0] += (_7yt2m6t[_h9ier7u[2] +1]*_jbc2sws[i][_h9ier7u[1] +0] + _7yt2m6t[_h9ier7u[2] +2]*_jbc2sws[i][_h9ier7u[1] +1] + _7yt2m6t[_h9ier7u[2] +3]*_jbc2sws[i][_h9ier7u[1] +2]);
                            _9p3jpxg[1] += (_7yt2m6t[_h9ier7u[2] +1]*_4xkss0k[i][_h9ier7u[1] +0] + _7yt2m6t[_h9ier7u[2] +2]*_4xkss0k[i][_h9ier7u[1] +1] + _7yt2m6t[_h9ier7u[2] +3]*_4xkss0k[i][_h9ier7u[1] +2]);
                            _9p3jpxg[2] += (_7yt2m6t[_h9ier7u[2] +1]*_t997e4u[i][_h9ier7u[1] +0] + _7yt2m6t[_h9ier7u[2] +2]*_t997e4u[i][_h9ier7u[1] +1] + _7yt2m6t[_h9ier7u[2] +3]*_t997e4u[i][_h9ier7u[1] +2]);
                        }
                        _m7pyxdw[i*16 +0] += _9p3jpxg[0];
                        _m7pyxdw[i*16 +1] += _9p3jpxg[1];
                        _m7pyxdw[i*16 +2] += _9p3jpxg[2];
                    }
                    else if ( _h9ier7u[0] == 2)
                    { __builtin___memset_chk (_32orwix, 0, 3*_btkf6xd[i]*sizeof(float), __builtin_object_size (_32orwix, 0));
                        for (j = 0u; j < _kb700uc; j++)
                        { _h9ier7u[2] = (unsigned int)_7yt2m6t[j*4 +0];
                            _h9ier7u[1] = _86lthet[i][_h9ier7u[2]] *3;
                            _h9ier7u[2] = j*4;
                            _32orwix[_h9ier7u[1] +0] += _7yt2m6t[_h9ier7u[2] +1];
                            _32orwix[_h9ier7u[1] +1] += _7yt2m6t[_h9ier7u[2] +2];
                            _32orwix[_h9ier7u[1] +2] += _7yt2m6t[_h9ier7u[2] +3];
                        }
                        _m7pyxdw[i*16 +0] += cblas_sdot(3*_btkf6xd[i], _9yo18pa, 1, _jbc2sws[i], 1);
                        _m7pyxdw[i*16 +1] += cblas_sdot(3*_btkf6xd[i], _9yo18pa, 1, _4xkss0k[i], 1);
                        _m7pyxdw[i*16 +2] += cblas_sdot(3*_btkf6xd[i], _9yo18pa, 1, _t997e4u[i], 1);
                    }
                    _h9ier7u[0] = (_rgb93wk <= _9f77lu8[i])*1u + (_rgb93wk > _9f77lu8[i])*2u;
                    _h9ier7u[0] = (_rgb93wk <= 0u)*0u + (_rgb93wk > 0u)*_h9ier7u[0];
                    if ( _h9ier7u[0] == 1)
                    { _9p3jpxg[0] = 0.0f; _9p3jpxg[1] = 0.0f; _9p3jpxg[2] = 0.0f;
                        for (j = 0u; j < _rgb93wk; j++)
                        { _h9ier7u[2] = (unsigned int)_mybd2qe[j*3 +0];
                            _h9ier7u[1] = _op767j2[i][_h9ier7u[2]] *2;
                            _h9ier7u[2] = j*3;
                            _9p3jpxg[0] += (_mybd2qe[_h9ier7u[2] +1]*_ge2j51b[i][_h9ier7u[1] +0] + _mybd2qe[_h9ier7u[2] +2]*_ge2j51b[i][_h9ier7u[1] +1]);
                            _9p3jpxg[1] += (_mybd2qe[_h9ier7u[2] +1]*_sh3tow4[i][_h9ier7u[1] +0] + _mybd2qe[_h9ier7u[2] +2]*_sh3tow4[i][_h9ier7u[1] +1]);
                            _9p3jpxg[2] += (_mybd2qe[_h9ier7u[2] +1]*_pqvl00e[i][_h9ier7u[1] +0] + _mybd2qe[_h9ier7u[2] +2]*_pqvl00e[i][_h9ier7u[1] +1]);
                        }
                        _m7pyxdw[i*16 +0] += _9p3jpxg[0];
                        _m7pyxdw[i*16 +1] += _9p3jpxg[1];
                        _m7pyxdw[i*16 +2] += _9p3jpxg[2];
                    }
                    else if ( _h9ier7u[0] == 2)
                    { __builtin___memset_chk (_qxmcthp, 0, 2*_9f77lu8[i]*sizeof(float), __builtin_object_size (_qxmcthp, 0));
                        for (j = 0u; j < _rgb93wk; j++)
                        { _h9ier7u[2] = (unsigned int)_mybd2qe[j*3 +0];
                            _h9ier7u[1] = _op767j2[i][_h9ier7u[2]] *2;
                            _h9ier7u[2] = j*3;
                            _qxmcthp[_h9ier7u[1] +0] += _mybd2qe[_h9ier7u[2] +1];
                            _qxmcthp[_h9ier7u[1] +1] += _mybd2qe[_h9ier7u[2] +2];
                        }
                        _m7pyxdw[i*16 +0] += cblas_sdot(2*_9f77lu8[i], _xf7jxjn, 1, _ge2j51b[i], 1);
                        _m7pyxdw[i*16 +1] += cblas_sdot(2*_9f77lu8[i], _xf7jxjn, 1, _sh3tow4[i], 1);
                        _m7pyxdw[i*16 +2] += cblas_sdot(2*_9f77lu8[i], _xf7jxjn, 1, _pqvl00e[i], 1);
                    }
                }
            }
            for (i = 0u; i < _b6gs8g8[_20n8tte]; i++)
            {
                if (_6zjq7w7[i*4 +0] == 1u)
                { _9p3jpxg[0] = sqrtf(_795e5iz[i*4 +0]*_795e5iz[i*4 +0] + _795e5iz[i*4 +1]*_795e5iz[i*4 +1]);
                    _9p3jpxg[1] = _795e5iz[i*4 +3]* -1.0f*_795e5iz[i*4 +2] + 1.05f*_9p3jpxg[14];
                    _9p3jpxg[2] = sqrtf(_3n3bzt9[i*16 + 0]*_3n3bzt9[i*16 + 0] + _3n3bzt9[i*16 + 1]*_3n3bzt9[i*16 + 1]);
                    _9p3jpxg[3] = _3n3bzt9[i*16 +3]* -1.0f*_3n3bzt9[i*16 +2] + 1.05f*_9p3jpxg[14];
                    _9p3jpxg[4] = _9p3jpxg[0] - _9p3jpxg[1];
                    _9p3jpxg[5] = _9p3jpxg[2] - _9p3jpxg[3];
                    _9p3jpxg[8] = sqrtf(_3n3bzt9[i*16 +9]*_3n3bzt9[i*16 +9] + _3n3bzt9[i*16 +10]*_3n3bzt9[i*16 +10]);
                    _9p3jpxg[8] = (_9p3jpxg[8] > 0.0f)*_9p3jpxg[8] + (_9p3jpxg[8] > 0.0f)*-1.0f;
                    _9p3jpxg[6] = (( ((_9p3jpxg[0]) > (_9p3jpxg[1])) *(_9p3jpxg[0]) + ((_9p3jpxg[0]) <= (_9p3jpxg[1])) *(_9p3jpxg[1]) )*_53j2d54)/_9p3jpxg[8];
                    _9p3jpxg[6] = (_9p3jpxg[8] > 0.0f)*_9p3jpxg[6] + (_9p3jpxg[8] <= 0.0f)*-3.40282347e+38F;
                    _9p3jpxg[8] = (_9p3jpxg[5] > _9p3jpxg[4])*(_9p3jpxg[5] - _9p3jpxg[4]) + (_9p3jpxg[5] <= _9p3jpxg[4])*-1.0f;
                    _9p3jpxg[7] = (_9p3jpxg[5] > _9p3jpxg[4])*((-1.0f*_9p3jpxg[4])/_9p3jpxg[8] *_ckeh5sw) + (_9p3jpxg[5] <= _9p3jpxg[4])*3.40282347e+38F;
                    _9p3jpxg[7] = ( ((_9p3jpxg[7]) > (0.0f)) *(_9p3jpxg[7]) + ((_9p3jpxg[7]) <= (0.0f)) *(0.0f) );
                    _h9ier7u[0] = (_9p3jpxg[7] > _ckeh5sw)*1u + (_9p3jpxg[7] <= _ckeh5sw)*0u;
                    _h9ier7u[1] = (fabsf(_9p3jpxg[7] - _ckeh5sw) <= _9p3jpxg[6])*1u + (fabsf(_9p3jpxg[7] - _ckeh5sw) > _9p3jpxg[6])*0u;
                    _9p3jpxg[7] = (_h9ier7u[0] == 1u)*( ((_ckeh5sw+_9p3jpxg[6]) > ((_9p3jpxg[7]))) *(_ckeh5sw+_9p3jpxg[6]) + ((_ckeh5sw+_9p3jpxg[6]) <= ((_9p3jpxg[7]))) *((_9p3jpxg[7])) ) + (_h9ier7u[0] != 1u)*((_h9ier7u[1] == 1u)*_ckeh5sw + (_h9ier7u[1] != 1u)*_9p3jpxg[7]);
                    _9p3jpxg[7] = (_9p3jpxg[0] >= _9p3jpxg[1])*0.0f + (_9p3jpxg[0] < _9p3jpxg[1])*_9p3jpxg[7];
                    _9p3jpxg[15]= ( ((_9p3jpxg[7]) < (_9p3jpxg[15])) *(_9p3jpxg[7]) + ((_9p3jpxg[7]) >= (_9p3jpxg[15])) *(_9p3jpxg[15]) );
                 }
            }
            MPI_Allreduce((void *) -1, &_9p3jpxg[15], 1, ((MPI_Datatype)0x4c00040a), (MPI_Op)(0x58000002), ((MPI_Comm)0x44000000));
            if (_9p3jpxg[15] == 0.0f)
            { k = 0;
                _ckeh5sw = 0.0f;
                for (i = 0u; i < _b6gs8g8[_20n8tte]; i++) { __builtin___memcpy_chk (&_3n3bzt9[i*16 +0], &_795e5iz[i*4 +0], 4*sizeof(float), __builtin_object_size (&_3n3bzt9[i*16 +0], 0)); _3n3bzt9[i*16 +15]= 0.0f; }
                for (i = 0u; i < _zfflb7w[_20n8tte]; i++) { __builtin___memcpy_chk (&_m7pyxdw[i*16 +0], &_vlnknd6[i*3 +0], 3*sizeof(float), __builtin_object_size (&_m7pyxdw[i*16 +0], 0)); }
            }
            else
            {
                if (_9p3jpxg[15] > _ckeh5sw)
                { _qzcucts += (double)_ckeh5sw;
                    _ckeh5sw = (_9p3jpxg[15] -_ckeh5sw);
                    for (i = 0u; i < _b6gs8g8[_20n8tte]; i++) { __builtin___memcpy_chk (&_795e5iz[i*4 +0], &_3n3bzt9[i*16 +0], 4*sizeof(float), __builtin_object_size (&_795e5iz[i*4 +0], 0)); _mzwffhh[i*16 +7] = _3n3bzt9[i*16 +3] - _mzwffhh[i*16+2]; _3p91i3v[i*8 +3] += _3n3bzt9[i*16 +15]; }
                    for (i = 0u; i < _zfflb7w[_20n8tte]; i++) { __builtin___memcpy_chk (&_vlnknd6[i*3 +0], &_m7pyxdw[i*16 +0], 3*sizeof(float), __builtin_object_size (&_vlnknd6[i*3 +0], 0)); }
                    _w2pcn9t += 1u;
                }
                else if (_9p3jpxg[15] == _ckeh5sw) { k = 0; }
                else if (_9p3jpxg[15] < _ckeh5sw) { _ckeh5sw = _9p3jpxg[15]; }
            }
        }
        _qzcucts += (double)_ckeh5sw;
        for (i = 0u; i < _b6gs8g8[_20n8tte]; i++) { _3p91i3v[i*8 +3] += _3n3bzt9[i*16 +15]; }
        _7bzegbl = 0u;
        _fu7bfxq = 0u;
        _d1vn8lf = 0.0f;
        __builtin___memset_chk (_0wfoy4d, 0, 3*sizeof(float), __builtin_object_size (_0wfoy4d, 0));
        for (i = 0u; i < _b6gs8g8[_20n8tte]; i++)
        { _9p3jpxg[0] = sqrtf( _3n3bzt9[i*16 +0]*_3n3bzt9[i*16 +0] + _3n3bzt9[i*16 +1]*_3n3bzt9[i*16 +1] );
            _9p3jpxg[0] = (_6zjq7w7[i*4 +0] == 1u)*_9p3jpxg[0] + (_6zjq7w7[i*4 +0] != 1u)*0.0f;
            _9p3jpxg[1] = _3n3bzt9[i*16 +3] *-1.0f*_3n3bzt9[i*16 +2];
            _9p3jpxg[2] = ((1.0f*_5zja8dz._4yfds0d)/_5zja8dz._pkxvhag) *( ((0.0f) > ((_9p3jpxg[0] -_9p3jpxg[1]))) *(0.0f) + ((0.0f) <= ((_9p3jpxg[0] -_9p3jpxg[1]))) *((_9p3jpxg[0] -_9p3jpxg[1])) ) *_5zja8dz._gfu8iog;
            _h9ier7u[0] = (_9p3jpxg[2] >= _z3nx39r._3ac6pe6)*1u + (_9p3jpxg[2] < _z3nx39r._3ac6pe6)*0u;
            _9p3jpxg[2] = (_h9ier7u[0] == 1u)*_9p3jpxg[2] + (_h9ier7u[0] != 1u)*0.0f;
            _0wfoy4d[0] = (_9p3jpxg[2] > _d1vn8lf)*_3p91i3v[i*8 +0] + (_9p3jpxg[2] <= _d1vn8lf)*_0wfoy4d[0];
            _0wfoy4d[1] = (_9p3jpxg[2] > _d1vn8lf)*_3p91i3v[i*8 +1] + (_9p3jpxg[2] <= _d1vn8lf)*_0wfoy4d[1];
            _0wfoy4d[2] = (_9p3jpxg[2] > _d1vn8lf)*_3p91i3v[i*8 +2] + (_9p3jpxg[2] <= _d1vn8lf)*_0wfoy4d[2];
            _fu7bfxq = (_9p3jpxg[2] > _d1vn8lf)*i + (_9p3jpxg[2] <= _d1vn8lf)*_fu7bfxq;
            _d1vn8lf = (_9p3jpxg[2] > _d1vn8lf)*_9p3jpxg[2] + (_9p3jpxg[2] <= _d1vn8lf)*_d1vn8lf;
        }
        { struct { float _s1gd1t5; int _ubods1e; } _ov6tt9q[1], _163yk03[1];
            _ov6tt9q[0]._s1gd1t5 = _d1vn8lf; _ov6tt9q[0]._ubods1e = _20n8tte;
            MPI_Allreduce(_ov6tt9q, _163yk03, 1, ((MPI_Datatype)0x8c000000), (MPI_Op)(0x5800000c), ((MPI_Comm)0x44000000));
            MPI_Bcast(&_163yk03[0]._ubods1e, 1, ((MPI_Datatype)0x4c000405), _163yk03[0]._ubods1e, ((MPI_Comm)0x44000000));
            MPI_Bcast(_0wfoy4d, 3, ((MPI_Datatype)0x4c00040a), _163yk03[0]._ubods1e, ((MPI_Comm)0x44000000));
            _7bzegbl = _163yk03[0]._ubods1e;
        }
        _f7fr91i = (_d1vn8lf > 0.0f)*1u + (_d1vn8lf <= 0.0f)*0u;
        MPI_Allreduce((void *) -1, &_f7fr91i, 1, ((MPI_Datatype)0x4c000406), (MPI_Op)(0x58000003), ((MPI_Comm)0x44000000));
        if (_f7fr91i == 0u) { perror("something went wrong, left Load2NextEQ without finding next event"); exit(1); }
        clock_gettime(_CLOCK_REALTIME, &_rg2gzr7);
        _q8hke22 += (_rg2gzr7.tv_sec - _y79kvls.tv_sec) + (_rg2gzr7.tv_nsec - _y79kvls.tv_nsec)/1.0E+9;
        for (i = 0u; i < _b6gs8g8[_20n8tte]; i++)
        { __builtin___memcpy_chk (&_795e5iz[i*4 +0], &_3n3bzt9[i*16 +0], 2*sizeof(float), __builtin_object_size (&_795e5iz[i*4 +0], 0));
            _9p3jpxg[0] = (_mzwffhh[i*16 +3] *-1.0f*_3n3bzt9[i*16 +2]) + _mzwffhh[i*16 +1]*_3n3bzt9[i*16 +4];
            _9p3jpxg[0] /= (-1.0f *_3n3bzt9[i*16 +2]);
            _3n3bzt9[i*16 +3] = ( ((_3n3bzt9[i*16 +3]) > (_9p3jpxg[0])) *(_3n3bzt9[i*16 +3]) + ((_3n3bzt9[i*16 +3]) <= (_9p3jpxg[0])) *(_9p3jpxg[0]) );
            _mzwffhh[i*16 +5] = _3n3bzt9[i*16 +3];
            _mzwffhh[i*16 +6] = _3n3bzt9[i*16 +3];
            _3n3bzt9[i*16 +2] = (_20n8tte != _7bzegbl)*((i != _fu7bfxq)*(_3n3bzt9[i*16 +2] + -50000.0f) + (i == _fu7bfxq)*_3n3bzt9[i*16 +2]) + (_20n8tte == _7bzegbl)*_3n3bzt9[i*16 +2];
        }
        clock_gettime(_CLOCK_REALTIME, &_tkezjss);
        clock_gettime(_CLOCK_REALTIME, &_zwd0w9r);
        _rgb93wk = 1u;
        _yt91qe3 = 0u;
        _blmntfk = 0u;
        _8sryo4r = 0.0f;
        while (_rgb93wk > 0u)
        {
            _ybxkpm2 = 0u;
            _rgb93wk = 0u;
            for (i = 0u; i < _b6gs8g8[_20n8tte]; i++)
            { _3n3bzt9[i*16 +12] = (_3n3bzt9[i*16 +11] <= 0.0f)*0.0f + (_3n3bzt9[i*16 +11] > 0.0f)*_3n3bzt9[i*16 +12];
                _mzwffhh[i*16 +5] = (_3n3bzt9[i*16 +11] <= 0.0f)*( ((_mzwffhh[i*16 +3]) > (_mzwffhh[i*16 +5])) *(_mzwffhh[i*16 +3]) + ((_mzwffhh[i*16 +3]) <= (_mzwffhh[i*16 +5])) *(_mzwffhh[i*16 +5]) ) + (_3n3bzt9[i*16 +11] > 0.0f)*_mzwffhh[i*16 +5];
                _3n3bzt9[i*16 +3] = (_3n3bzt9[i*16 +11] <= 0.0f)*( ((_mzwffhh[i*16 +3]) > (_3n3bzt9[i*16 +3])) *(_mzwffhh[i*16 +3]) + ((_mzwffhh[i*16 +3]) <= (_3n3bzt9[i*16 +3])) *(_3n3bzt9[i*16 +3]) ) + (_3n3bzt9[i*16 +11] > 0.0f)*_3n3bzt9[i*16 +3];
                _9p3jpxg[1] = (_mzwffhh[i*16 +3] - _mzwffhh[i*16 +5])/(1.0f*_3n3bzt9[i*16 +4]);
                _9p3jpxg[2] = (_mzwffhh[i*16 +4] - _mzwffhh[i*16 +3])/(0.5f*_3n3bzt9[i*16 +4]);
                _9p3jpxg[3] = _mzwffhh[i*16 +3] - _9p3jpxg[2]*_3n3bzt9[i*16 +4];
                _9p3jpxg[4] = (_3n3bzt9[i*16 +12] <= _3n3bzt9[i*16 +4]) *(_mzwffhh[i*16 +5] + _9p3jpxg[1]*_3n3bzt9[i*16 +12]) + (_3n3bzt9[i*16 +12] > _3n3bzt9[i*16 +4])*(_9p3jpxg[3] + _9p3jpxg[2]*_3n3bzt9[i*16 +12]);
                _9p3jpxg[4] = ( ((_9p3jpxg[4]) > (_mzwffhh[i*16 +4])) *(_9p3jpxg[4]) + ((_9p3jpxg[4]) <= (_mzwffhh[i*16 +4])) *(_mzwffhh[i*16 +4]) );
                _9p3jpxg[5] = _3n3bzt9[i*16 +3] + _z3nx39r._e9ftaxp*(_mzwffhh[i*16 +6] - _3n3bzt9[i*16 +3]);
                _3n3bzt9[i*16 +3] = (_6zjq7w7[i*4 +1] == 0u)*_3n3bzt9[i*16 +3] + (_6zjq7w7[i*4 +1] != 0u)*( (_3n3bzt9[i*16 +11] <= 0.0f)*_9p3jpxg[5] + (_3n3bzt9[i*16 +11] > 0.0f)*_9p3jpxg[4] );
                _mzwffhh[i*16 +5] = (_3n3bzt9[i*16 +11] <= 0.0f)*_9p3jpxg[5] + (_3n3bzt9[i*16 +11] > 0.0f)*_mzwffhh[i*16 +5];
                _3n3bzt9[i*16 +11] = 0.0f;
                _9p3jpxg[0] = sqrtf( _3n3bzt9[i*16 +0]*_3n3bzt9[i*16 +0] + _3n3bzt9[i*16 +1]*_3n3bzt9[i*16 +1] );
                _3n3bzt9[i*16 +2] = ( ((0.0f) < (_3n3bzt9[i*16 +2])) *(0.0f) + ((0.0f) >= (_3n3bzt9[i*16 +2])) *(_3n3bzt9[i*16 +2]) );
                _9p3jpxg[15] = _5zja8dz._rtj5pht + _5zja8dz._ktakpux*-1.0f*_3n3bzt9[i*16 +2];
                _9p3jpxg[1] = _3n3bzt9[i*16 +3] *-1.0f*_3n3bzt9[i*16 +2];
                _9p3jpxg[3] = _9p3jpxg[0] - _9p3jpxg[1];
                _9p3jpxg[5] = ((1.0f*_5zja8dz._4yfds0d)/_5zja8dz._pkxvhag) *_9p3jpxg[3] *_5zja8dz._gfu8iog;
                _h9ier7u[0] = (_9p3jpxg[5] < _z3nx39r._3ac6pe6)*0u + (_9p3jpxg[5] >= _z3nx39r._3ac6pe6)*1u;
                _h9ier7u[1] = (_6zjq7w7[i*4 +1] == 2u)*2u + (_6zjq7w7[i*4 +1] != 2u)*1u;
                _h9ier7u[0] *= _h9ier7u[1];
                if (_h9ier7u[0] == 1u)
                {
                    _ddqhkga[_yt91qe3] = (_6zjq7w7[i*4 +1] == 0u)*i + (_6zjq7w7[i*4 +1] != 0u)*_ddqhkga[_yt91qe3];
                    _yt91qe3 = (_6zjq7w7[i*4 +1] == 0u)*(_yt91qe3 +1u) + (_6zjq7w7[i*4 +1] != 0u)*_yt91qe3;
                    _6zjq7w7[i*4 +2] = (_6zjq7w7[i*4 +1] == 0u)*_blmntfk + (_6zjq7w7[i*4 +1] != 0u)*_6zjq7w7[i*4 +2];
                    _6zjq7w7[i*4 +1] = (_6zjq7w7[i*4 +1] == 0u)*1u + (_6zjq7w7[i*4 +1] != 0u)*_6zjq7w7[i*4 +1];
                    _6zjq7w7[i*4 +1] = (_9p3jpxg[0] > _9p3jpxg[15])*2u + (_9p3jpxg[0] <= _9p3jpxg[15])*_6zjq7w7[i*4 +1];
                    _6zjq7w7[i*4 +3] = _blmntfk +1u - _6zjq7w7[i*4 +2];
                    _3n3bzt9[i*16 +15] = fabsf(_9p3jpxg[5]);
                    _9p3jpxg[5] = -1.0f*((_9p3jpxg[3]/_9p3jpxg[0]) *_3n3bzt9[i*16 +0]) / _3n3bzt9[i*16 +5];
                    _9p3jpxg[6] = -1.0f*((_9p3jpxg[3]/_9p3jpxg[0]) *_3n3bzt9[i*16 +1]) / _3n3bzt9[i*16 +6];
                    _9p3jpxg[7] = _3n3bzt9[i*16 +15]/sqrtf(_9p3jpxg[5]*_9p3jpxg[5] +_9p3jpxg[6]*_9p3jpxg[6]);
                    _9p3jpxg[7] = ( ((_9p3jpxg[7]) < (1.0f)) *(_9p3jpxg[7]) + ((_9p3jpxg[7]) >= (1.0f)) *(1.0f) );
                    _9p3jpxg[5] *= _9p3jpxg[7];
                    _9p3jpxg[6] *= _9p3jpxg[7];
                    _3n3bzt9[i*16 +11] = sqrtf(_9p3jpxg[5]*_9p3jpxg[5] +_9p3jpxg[6]*_9p3jpxg[6]);
                    _3n3bzt9[i*16 +12] += _3n3bzt9[i*16 +11];
                    _3n3bzt9[i*16 +13] += _9p3jpxg[5];
                    _3n3bzt9[i*16 +14] += _9p3jpxg[6];
                    _8sryo4r += (_3n3bzt9[i*16 +11]*_3n3bzt9[i*16 +8]);
                    _ntdr5xy[i][ (2u*(_blmntfk - _6zjq7w7[i*4 +2])) +0] = _9p3jpxg[5];
                    _ntdr5xy[i][ (2u*(_blmntfk - _6zjq7w7[i*4 +2])) +1] = _9p3jpxg[6];
                    _4pa8hye[_ybxkpm2*3 +0] = (float)(_xo100vh[i]);
                    _4pa8hye[_ybxkpm2*3 +1] = _9p3jpxg[5]*_3n3bzt9[i*16 +8];
                    _4pa8hye[_ybxkpm2*3 +2] = _9p3jpxg[6]*_3n3bzt9[i*16 +8];
                    _ybxkpm2 += 1u;
                }
                else if (_h9ier7u[0] == 2u)
                { _6zjq7w7[i*4 +3] = _blmntfk +1u - _6zjq7w7[i*4 +2];
                    _3n3bzt9[i*16 +15] = fabsf(_9p3jpxg[5]);
                    _9p3jpxg[5] = -1.0f*((_9p3jpxg[3]/_9p3jpxg[0]) *_3n3bzt9[i*16 +0]) / _3n3bzt9[i*16 +5];
                    _9p3jpxg[6] = -1.0f*((_9p3jpxg[3]/_9p3jpxg[0]) *_3n3bzt9[i*16 +1]) / _3n3bzt9[i*16 +6];
                    _9p3jpxg[7] = _3n3bzt9[i*16 +15]/sqrtf(_9p3jpxg[5]*_9p3jpxg[5] +_9p3jpxg[6]*_9p3jpxg[6]);
                    _9p3jpxg[7] = ( ((_9p3jpxg[7]) < (1.0f)) *(_9p3jpxg[7]) + ((_9p3jpxg[7]) >= (1.0f)) *(1.0f) );
                    _9p3jpxg[5] *= _9p3jpxg[7];
                    _9p3jpxg[6] *= _9p3jpxg[7];
                    _3n3bzt9[i*16 +11] = sqrtf(_9p3jpxg[5]*_9p3jpxg[5] +_9p3jpxg[6]*_9p3jpxg[6]);
                    _3n3bzt9[i*16 +12] += _3n3bzt9[i*16 +11];
                    _3n3bzt9[i*16 +13] += _9p3jpxg[5];
                    _3n3bzt9[i*16 +14] += _9p3jpxg[6];
                    _8sryo4r += (_3n3bzt9[i*16 +11]*_3n3bzt9[i*16 +8]);
                    _ntdr5xy[i][ (2u*(_blmntfk - _6zjq7w7[i*4 +2])) +0] = _9p3jpxg[5];
                    _ntdr5xy[i][ (2u*(_blmntfk - _6zjq7w7[i*4 +2])) +1] = _9p3jpxg[6];
                    _3n3bzt9[i*16 +0] += (_9p3jpxg[5] *_3n3bzt9[i*16 +5]);
                    _3n3bzt9[i*16 +1] += (_9p3jpxg[6] *_3n3bzt9[i*16 +6]);
            } }
            MPI_Allgather(&_ybxkpm2, 1, ((MPI_Datatype)0x4c000406), _o2t8fy6, 1, ((MPI_Datatype)0x4c000405), ((MPI_Comm)0x44000000));
            _gg24479[0] = 0u;
            _rgb93wk = _o2t8fy6[0];
            for (i = 1; i < _avt8k06; i++) { _rgb93wk += _o2t8fy6[i]; _o2t8fy6[i-1] *= 3; _gg24479[i] = _gg24479[i-1] + _o2t8fy6[i-1]; }
            _o2t8fy6[_avt8k06-1] *= 3;
            MPI_Allgatherv(_4pa8hye, _o2t8fy6[_20n8tte], ((MPI_Datatype)0x4c00040a), _hj3muur, _o2t8fy6, _gg24479, ((MPI_Datatype)0x4c00040a), ((MPI_Comm)0x44000000));
            for (i = 0u; i < _b6gs8g8[_20n8tte]; i++)
            {
                _h9ier7u[0] = (_rgb93wk <= _iob0nlw[i])*1u + (_rgb93wk > _iob0nlw[i])*2u;
                _h9ier7u[0] = (_rgb93wk <= 0u)*0u + (_rgb93wk > 0u)*_h9ier7u[0];
                if (_h9ier7u[0] == 1u)
                { _9p3jpxg[0] = 0.0f; _9p3jpxg[1] = 0.0f; _9p3jpxg[2] = 0.0f;
                    for (j = 0u; j < _rgb93wk; j++)
                    { _h9ier7u[2] = (unsigned int)_mybd2qe[j*3 +0];
                        _h9ier7u[1] = _3b9so0u[i][_h9ier7u[2]] *2;
                        _h9ier7u[2] = j*3;
                        _9p3jpxg[0] += (_mybd2qe[_h9ier7u[2] +1]*_0dtirc9[i][_h9ier7u[1] +0] + _mybd2qe[_h9ier7u[2] +2]*_0dtirc9[i][_h9ier7u[1] +1] );
                        _9p3jpxg[1] += (_mybd2qe[_h9ier7u[2] +1]*_qk4ltah[i][_h9ier7u[1] +0] + _mybd2qe[_h9ier7u[2] +2]*_qk4ltah[i][_h9ier7u[1] +1] );
                        _9p3jpxg[2] += (_mybd2qe[_h9ier7u[2] +1]*_h61iaoi[i][_h9ier7u[1] +0] + _mybd2qe[_h9ier7u[2] +2]*_h61iaoi[i][_h9ier7u[1] +1] );
                    }
                    _3n3bzt9[i*16 +0] += _9p3jpxg[0];
                    _3n3bzt9[i*16 +1] += _9p3jpxg[1];
                    _3n3bzt9[i*16 +2] += _9p3jpxg[2];
                }
                else if (_h9ier7u[0] == 2u)
                { __builtin___memset_chk (_qxmcthp, 0, 2*_iob0nlw[i]*sizeof(float), __builtin_object_size (_qxmcthp, 0));
                    for (j = 0u; j < _rgb93wk; j++)
                    { _h9ier7u[2] = (unsigned int)_mybd2qe[j*3 +0];
                        _h9ier7u[1] = _3b9so0u[i][_h9ier7u[2]] *2;
                        _h9ier7u[2] = j*3;
                        _qxmcthp[_h9ier7u[1] +0] += _mybd2qe[_h9ier7u[2] +1];
                        _qxmcthp[_h9ier7u[1] +1] += _mybd2qe[_h9ier7u[2] +2];
                    }
                    _3n3bzt9[i*16 +0] += cblas_sdot(2*_iob0nlw[i], _xf7jxjn, 1, _0dtirc9[i], 1);
                    _3n3bzt9[i*16 +1] += cblas_sdot(2*_iob0nlw[i], _xf7jxjn, 1, _qk4ltah[i], 1);
                    _3n3bzt9[i*16 +2] += cblas_sdot(2*_iob0nlw[i], _xf7jxjn, 1, _h61iaoi[i], 1);
            } }
            _rgb93wk = (_blmntfk < _5zja8dz._148p0tv-1)*_rgb93wk + (_blmntfk >= _5zja8dz._148p0tv-1)*0;
            _blmntfk += 1u;
        }
        if (_z3nx39r._urvqt3g > 0u)
        {
            _ybxkpm2 = 0u;
            for (i = 0u; i < _yt91qe3; i++)
            { _h9ier7u[0] =((fabsf(_3n3bzt9[_ddqhkga[i]*16 +13]) + fabsf(_3n3bzt9[_ddqhkga[i]*16 +14])) > 0.0f)*1u + ((fabsf(_3n3bzt9[_ddqhkga[i]*16 +13]) + fabsf(_3n3bzt9[_ddqhkga[i]*16 +14])) <= 0.0f)*0u;
                _4pa8hye[_ybxkpm2*3 +0] = (_h9ier7u[0] == 0u)*_4pa8hye[_ybxkpm2*3 +0] + (_h9ier7u[0] != 0u)*(float)(_xo100vh[_ddqhkga[i]]);
                _4pa8hye[_ybxkpm2*3 +1] = (_h9ier7u[0] == 0u)*_4pa8hye[_ybxkpm2*3 +1] + (_h9ier7u[0] != 0u)*_3n3bzt9[_ddqhkga[i]*16 +13]*_3n3bzt9[_ddqhkga[i]*16 +8];
                _4pa8hye[_ybxkpm2*3 +2] = (_h9ier7u[0] == 0u)*_4pa8hye[_ybxkpm2*3 +2] + (_h9ier7u[0] != 0u)*_3n3bzt9[_ddqhkga[i]*16 +14]*_3n3bzt9[_ddqhkga[i]*16 +8];
                _ybxkpm2 = (_h9ier7u[0] == 0u)*_ybxkpm2 + (_h9ier7u[0] != 0u)*(_ybxkpm2 + 1u);
            }
            MPI_Allgather(&_ybxkpm2, 1, ((MPI_Datatype)0x4c000406), _o2t8fy6, 1, ((MPI_Datatype)0x4c000405), ((MPI_Comm)0x44000000));
            _gg24479[0] = 0u;
            _rgb93wk = _o2t8fy6[0];
            for (i = 1; i < _avt8k06; i++) { _rgb93wk += _o2t8fy6[i]; _o2t8fy6[i-1] *= 3; _gg24479[i] = _gg24479[i-1] + _o2t8fy6[i-1]; }
            _o2t8fy6[_avt8k06-1] *= 3;
            MPI_Allgatherv(_4pa8hye, _o2t8fy6[_20n8tte], ((MPI_Datatype)0x4c00040a), _hj3muur, _o2t8fy6, _gg24479, ((MPI_Datatype)0x4c00040a), ((MPI_Comm)0x44000000));
            for (i = 0u; i < _zfflb7w[_20n8tte]; i++)
            { _h9ier7u[0] = (_rgb93wk <= _9f77lu8[i])*1u + (_rgb93wk > _9f77lu8[i])*2u;
                _h9ier7u[0] = (_rgb93wk <= 0u)*0u + (_rgb93wk > 0u)*_h9ier7u[0];
                if (_h9ier7u[0] == 1)
                { _9p3jpxg[0] = 0.0f; _9p3jpxg[1] = 0.0f; _9p3jpxg[2] = 0.0f;
                    for (j = 0u; j < _rgb93wk; j++)
                    { _h9ier7u[2] = (unsigned int)_mybd2qe[j*3 +0];
                        _h9ier7u[1] = _op767j2[i][_h9ier7u[2]] *2;
                        _h9ier7u[2] = j*3;
                        _9p3jpxg[0] += (_mybd2qe[_h9ier7u[2] +1]*_ge2j51b[i][_h9ier7u[1] +0] + _mybd2qe[_h9ier7u[2] +2]*_ge2j51b[i][_h9ier7u[1] +1]);
                        _9p3jpxg[1] += (_mybd2qe[_h9ier7u[2]+ 1]*_sh3tow4[i][_h9ier7u[1] +0] + _mybd2qe[_h9ier7u[2] +2]*_sh3tow4[i][_h9ier7u[1] +1]);
                        _9p3jpxg[2] += (_mybd2qe[_h9ier7u[2] +1]*_pqvl00e[i][_h9ier7u[1] +0] + _mybd2qe[_h9ier7u[2] +2]*_pqvl00e[i][_h9ier7u[1] +1]);
                    }
                    _m7pyxdw[i*16 +0] += _9p3jpxg[0];
                    _m7pyxdw[i*16 +1] += _9p3jpxg[1];
                    _m7pyxdw[i*16 +2] += _9p3jpxg[2];
                }
                else if (_h9ier7u[0] == 2)
                { __builtin___memset_chk (_qxmcthp, 0, 2*_9f77lu8[i]*sizeof(float), __builtin_object_size (_qxmcthp, 0));
                    for (j = 0u; j < _rgb93wk; j++)
                    { _h9ier7u[2] = (unsigned int)_mybd2qe[j*3 +0];
                        _h9ier7u[1] = _op767j2[i][_h9ier7u[2]];
                        _h9ier7u[2] = j*3;
                        _qxmcthp[_h9ier7u[1] +0] += _mybd2qe[_h9ier7u[2] +1];
                        _qxmcthp[_h9ier7u[1] +1] += _mybd2qe[_h9ier7u[2] +2];
                    }
                    _m7pyxdw[i*16 +0] += cblas_sdot(2*_9f77lu8[i], _xf7jxjn, 1, _ge2j51b[i], 1);
                    _m7pyxdw[i*16 +1] += cblas_sdot(2*_9f77lu8[i], _xf7jxjn, 1, _sh3tow4[i], 1);
                    _m7pyxdw[i*16 +2] += cblas_sdot(2*_9f77lu8[i], _xf7jxjn, 1, _pqvl00e[i], 1);
                }
            }
        }
        _v4zwl3z = _yt91qe3;
        MPI_Allreduce((void *) -1, &_v4zwl3z, 1, ((MPI_Datatype)0x4c000406), (MPI_Op)(0x58000003), ((MPI_Comm)0x44000000));
        MPI_Allreduce((void *) -1, &_8sryo4r, 1, ((MPI_Datatype)0x4c00040a), (MPI_Op)(0x58000003), ((MPI_Comm)0x44000000));
        _8sryo4r *= _5zja8dz._pkxvhag;
        _wyx9trf = (log10f(_8sryo4r) -9.1f)/1.5f;
        if (_v4zwl3z >= _z3nx39r._o31yhxs)
        { int _w8yl68n = 0;
            _jzr25md += 1u;
            _zfk8ut0 += 1u;
            __builtin___memset_chk (_lqbmae3, 0, 6*sizeof(float), __builtin_object_size (_lqbmae3, 0));
            for (i = 0u; i < _yt91qe3; i++)
            { _h9ier7u[1] = _ddqhkga[i];
                _9p3jpxg[0] = sqrtf(_3n3bzt9[_ddqhkga[i]*16 +13]*_3n3bzt9[_ddqhkga[i]*16 +13] +_3n3bzt9[_ddqhkga[i]*16 +14]*_3n3bzt9[_ddqhkga[i]*16 +14]);
                _lqbmae3[0] += (_3p91i3v[_ddqhkga[i]*8 +0]*_9p3jpxg[0]*_3n3bzt9[_ddqhkga[i]*16 +8]);
                _lqbmae3[1] += (_3p91i3v[_ddqhkga[i]*8 +1]*_9p3jpxg[0]*_3n3bzt9[_ddqhkga[i]*16 +8]);
                _lqbmae3[2] += (_3p91i3v[_ddqhkga[i]*8 +2]*_9p3jpxg[0]*_3n3bzt9[_ddqhkga[i]*16 +8]);
                _lqbmae3[3] += _3n3bzt9[_ddqhkga[i]*16 +8];
                _lqbmae3[4] += _9p3jpxg[0];
                _lqbmae3[5] += sqrtf(_795e5iz[_ddqhkga[i]*4 +0]*_795e5iz[_ddqhkga[i]*4 +0] + _795e5iz[_ddqhkga[i]*4 +1]*_795e5iz[_ddqhkga[i]*4 +1]) - sqrtf(_3n3bzt9[_ddqhkga[i]*16 +0]*_3n3bzt9[_ddqhkga[i]*16 +0] + _3n3bzt9[_ddqhkga[i]*16 +1]*_3n3bzt9[_ddqhkga[i]*16 +1]);
                _k00ilu6[i] = _xo100vh[_ddqhkga[i]];
                _vfua9bq[i] = _6zjq7w7[_ddqhkga[i]*4 +2];
                _28xg8y7[i] = _3n3bzt9[_ddqhkga[i]*16 +13];
                _7odkef8[i] = _3n3bzt9[_ddqhkga[i]*16 +14];
                _v71448n[i] = _3p91i3v[_ddqhkga[i]*8 +3] + sqrtf(_28xg8y7[i]*_28xg8y7[i] + _7odkef8[i]*_7odkef8[i]);
                _2ivglye[i] = sqrtf(_795e5iz[_ddqhkga[i]*4 +0]*_795e5iz[_ddqhkga[i]*4 +0] + _795e5iz[_ddqhkga[i]*4 +1]*_795e5iz[_ddqhkga[i]*4 +1]) - sqrtf(_3n3bzt9[_ddqhkga[i]*16 +0]*_3n3bzt9[_ddqhkga[i]*16 +0] + _3n3bzt9[_ddqhkga[i]*16 +1]*_3n3bzt9[_ddqhkga[i]*16 +1]);
                _w8yl68n += (_6zjq7w7[_ddqhkga[i]*4 +1] == 2)*1 + (_6zjq7w7[_ddqhkga[i]*4 +1] != 2)*0;
            }
            MPI_Allreduce((void *) -1, _lqbmae3, 6, ((MPI_Datatype)0x4c00040a), (MPI_Op)(0x58000003), ((MPI_Comm)0x44000000));
            MPI_Allreduce((void *) -1, &_w8yl68n, 1, ((MPI_Datatype)0x4c000405), (MPI_Op)(0x58000003), ((MPI_Comm)0x44000000));
            _lqbmae3[0] /= (_8sryo4r/_5zja8dz._pkxvhag); _lqbmae3[1] /= (_8sryo4r/_5zja8dz._pkxvhag); _lqbmae3[2] /= (_8sryo4r/_5zja8dz._pkxvhag);
            _lqbmae3[4] /= (float)_v4zwl3z; _lqbmae3[5] /= (float)_v4zwl3z;
            MPI_Allgather(&_yt91qe3, 1, ((MPI_Datatype)0x4c000406), _o2t8fy6, 1, ((MPI_Datatype)0x4c000405), ((MPI_Comm)0x44000000));
            _gg24479[0] = 0;
            for (i = 1; i < _avt8k06; i++) { _gg24479[i] = _gg24479[i-1] + _o2t8fy6[i-1]; }
            MPI_Gatherv(_k00ilu6, _o2t8fy6[_20n8tte], ((MPI_Datatype)0x4c000406), _7c6qc6k, _o2t8fy6, _gg24479, ((MPI_Datatype)0x4c000406), 0, ((MPI_Comm)0x44000000));
            MPI_Gatherv(_vfua9bq, _o2t8fy6[_20n8tte], ((MPI_Datatype)0x4c000406), _8djzr5s, _o2t8fy6, _gg24479, ((MPI_Datatype)0x4c000406), 0, ((MPI_Comm)0x44000000));
            MPI_Gatherv(_28xg8y7, _o2t8fy6[_20n8tte], ((MPI_Datatype)0x4c00040a), _j3swebn, _o2t8fy6, _gg24479, ((MPI_Datatype)0x4c00040a), 0, ((MPI_Comm)0x44000000));
            MPI_Gatherv(_7odkef8, _o2t8fy6[_20n8tte], ((MPI_Datatype)0x4c00040a), _vg5vmf6, _o2t8fy6, _gg24479, ((MPI_Datatype)0x4c00040a), 0, ((MPI_Comm)0x44000000));
            MPI_Gatherv(_v71448n, _o2t8fy6[_20n8tte], ((MPI_Datatype)0x4c00040a), _0yknfu8, _o2t8fy6, _gg24479, ((MPI_Datatype)0x4c00040a), 0, ((MPI_Comm)0x44000000));
            MPI_Gatherv(_2ivglye, _o2t8fy6[_20n8tte], ((MPI_Datatype)0x4c00040a), _4xyjoa2, _o2t8fy6, _gg24479, ((MPI_Datatype)0x4c00040a), 0, ((MPI_Comm)0x44000000));
            if (_20n8tte == 0)
            { _0ljch2v = (_wyx9trf > _0ljch2v)*_wyx9trf + (_wyx9trf <= _0ljch2v)*_0ljch2v;
                if (_z3nx39r._zk94yxf == 1u) { fprintf(__stdoutp,"%7u  Time: %9.4lf  (%7.2lf days since last)   intSeisLoops: %2u   Element#: %6u (%3u broken)    RA: %8.2f   mSlip: %7.3f      MRF %5u     M: %4.2f      M_max: %4.2f\n", _jzr25md, _qzcucts, (_qzcucts-_uflsehu)*365.25, _g8vuazx, _v4zwl3z, _w8yl68n, (_lqbmae3[3]*1.0E-6f), _lqbmae3[4], (_blmntfk-1u), _wyx9trf, _0ljch2v); }
                fseek(_szbffzz, 0, 2);
                fwrite(&_qzcucts, sizeof(double), 1, _szbffzz);
                fwrite(&_8sryo4r, sizeof(float), 1, _szbffzz);
                fwrite(&_wyx9trf, sizeof(float), 1, _szbffzz);
                fwrite(&_blmntfk, sizeof(unsigned int), 1, _szbffzz);
                fwrite(_0wfoy4d, sizeof(float), 3, _szbffzz);
                fwrite(_lqbmae3, sizeof(float), 6, _szbffzz);
                fwrite(&_v4zwl3z, sizeof(unsigned), 1, _szbffzz);
                fwrite(_7c6qc6k, sizeof(unsigned int), _v4zwl3z, _szbffzz);
                fwrite(_8djzr5s, sizeof(unsigned int), _v4zwl3z, _szbffzz);
                fwrite(_j3swebn, sizeof(float), _v4zwl3z, _szbffzz);
                fwrite(_vg5vmf6, sizeof(float), _v4zwl3z, _szbffzz);
                fwrite(_0yknfu8, sizeof(float), _v4zwl3z, _szbffzz);
                fwrite(_4xyjoa2, sizeof(float), _v4zwl3z, _szbffzz);
                fseek(_szbffzz, 0, 0);
                fwrite(&_jzr25md, sizeof(unsigned int), 1, _szbffzz);
            }
            _ds3qsqc += (_jzr25md > 1u)*(_qzcucts - _uflsehu) + (_jzr25md <= 1u)*0.0;
            _uflsehu = _qzcucts;
        }
        MPI_Barrier( ((MPI_Comm)0x44000000) );
        if ((_z3nx39r._nwngucz == 1u) && (_wyx9trf >= _z3nx39r._yj9boq0))
        { float *_fc4kl8c = (float *) malloc( _5zja8dz._148p0tv *sizeof(float));
            float *_gn0t6dv = (float *) malloc( _5zja8dz._148p0tv *sizeof(float));
            long long int _eaxylnf = 0;
            long long int _361xsnb = 0;
            char _1750ceh[512];
            for (i = 0; i < _b6gs8g8[_20n8tte]; i++)
            { _361xsnb = (_6zjq7w7[i*4 +1] == 0u)*_361xsnb + (_6zjq7w7[i*4 +1] != 0u)*(_361xsnb + 2*sizeof(float) +4*sizeof(unsigned int) +2*_6zjq7w7[i*4 +3]*sizeof(float) );
            }
            MPI_Allgather(&_361xsnb, 1, ((MPI_Datatype)0x4c000809), _ti22nv0, 1, ((MPI_Datatype)0x4c000809), ((MPI_Comm)0x44000000));
            _1zderi4[0] = 0;
            for (i = 1; i < _avt8k06; i++) { _1zderi4[i] = _1zderi4[i-1] +_ti22nv0[i-1]; }
            __builtin___strcpy_chk (_1750ceh, _8dalfyo, __builtin_object_size (_1750ceh, 2 > 1 ? 1 : 0)); __builtin___strcat_chk (_1750ceh, "_t", __builtin_object_size (_1750ceh, 2 > 1 ? 1 : 0)); __builtin___snprintf_chk (_m3tzk3n, sizeof(_m3tzk3n), 0, __builtin_object_size (_m3tzk3n, 2 > 1 ? 1 : 0), "%lf",_qzcucts); __builtin___strcat_chk (_1750ceh, _m3tzk3n, __builtin_object_size (_1750ceh, 2 > 1 ? 1 : 0)); __builtin___strcat_chk (_1750ceh, "_M", __builtin_object_size (_1750ceh, 2 > 1 ? 1 : 0)); __builtin___snprintf_chk (_m3tzk3n, sizeof(_m3tzk3n), 0, __builtin_object_size (_m3tzk3n, 2 > 1 ? 1 : 0), "%f",_wyx9trf); __builtin___strcat_chk (_1750ceh, _m3tzk3n, __builtin_object_size (_1750ceh, 2 > 1 ? 1 : 0)); __builtin___strcat_chk (_1750ceh, ".srfb", __builtin_object_size (_1750ceh, 2 > 1 ? 1 : 0));
            MPI_File_delete(_1750ceh, ((MPI_Info)0x1c000000));
            MPI_File_open(((MPI_Comm)0x44000000), _1750ceh, 1|4, ((MPI_Info)0x1c000000), &_qx8ertb);
            if (_20n8tte == 0)
            { int _dyk0x65 = (int)202511;
                MPI_File_write_at(_qx8ertb, _eaxylnf, &_dyk0x65, 1, ((MPI_Datatype)0x4c000405), &_py6b2kn); _eaxylnf += 1*sizeof(int);
                MPI_File_write_at(_qx8ertb, _eaxylnf, &_wyx9trf, 1, ((MPI_Datatype)0x4c00040a), &_py6b2kn); _eaxylnf += 1*sizeof(float);
                _9p3jpxg[0] = (float)_qzcucts;
                MPI_File_write_at(_qx8ertb, _eaxylnf, &_9p3jpxg[0], 1, ((MPI_Datatype)0x4c00040a), &_py6b2kn); _eaxylnf += 1*sizeof(float);
                MPI_File_write_at(_qx8ertb, _eaxylnf, &_v4zwl3z, 1, ((MPI_Datatype)0x4c000406), &_py6b2kn); _eaxylnf += 1*sizeof(unsigned int);
                MPI_File_write_at(_qx8ertb, _eaxylnf, &_5zja8dz._gfu8iog, 1, ((MPI_Datatype)0x4c00040a), &_py6b2kn); _eaxylnf += 1*sizeof(float);
                MPI_File_write_at(_qx8ertb, _eaxylnf, &_5zja8dz._4yfds0d, 1, ((MPI_Datatype)0x4c00040a), &_py6b2kn); _eaxylnf += 1*sizeof(float);
                MPI_File_write_at(_qx8ertb, _eaxylnf, &_5zja8dz._xpbwhvn, 1, ((MPI_Datatype)0x4c00040a), &_py6b2kn); _eaxylnf += 1*sizeof(float);
            }
            _eaxylnf = 1*sizeof(int) +1*sizeof(unsigned int) +5*sizeof(float);
            for (i = 0; i < _b6gs8g8[_20n8tte]; i++)
            { if (_6zjq7w7[i*4 +1] != 0u)
                { _h9ier7u[0] = _xo100vh[i];
                    MPI_File_write_at(_qx8ertb, (_1zderi4[_20n8tte]+_eaxylnf), &_h9ier7u[0], 1, ((MPI_Datatype)0x4c000406), &_py6b2kn); _eaxylnf += 1*sizeof(unsigned int);
                    MPI_File_write_at(_qx8ertb, (_1zderi4[_20n8tte]+_eaxylnf), &_6zjq7w7[i*4 +2], 1, ((MPI_Datatype)0x4c000406), &_py6b2kn); _eaxylnf += 1*sizeof(unsigned int);
                    MPI_File_write_at(_qx8ertb, (_1zderi4[_20n8tte]+_eaxylnf), &_3n3bzt9[i*16 +13], 1, ((MPI_Datatype)0x4c00040a), &_py6b2kn); _eaxylnf += 1*sizeof(float);
                    MPI_File_write_at(_qx8ertb, (_1zderi4[_20n8tte]+_eaxylnf), &_6zjq7w7[i*4 +3], 1, ((MPI_Datatype)0x4c000406), &_py6b2kn); _eaxylnf += 1*sizeof(unsigned int);
                    MPI_File_write_at(_qx8ertb, (_1zderi4[_20n8tte]+_eaxylnf), &_3n3bzt9[i*16 +14], 1, ((MPI_Datatype)0x4c00040a), &_py6b2kn); _eaxylnf += 1*sizeof(float);
                    MPI_File_write_at(_qx8ertb, (_1zderi4[_20n8tte]+_eaxylnf), &_6zjq7w7[i*4 +3], 1, ((MPI_Datatype)0x4c000406), &_py6b2kn); _eaxylnf += 1*sizeof(unsigned int);
                    for (j = 0; j < _6zjq7w7[i*4 +3]; j++) { _fc4kl8c[j] = _ntdr5xy[i][2u*j +0]; _gn0t6dv[j] = _ntdr5xy[i][2u*j +1]; }
                    MPI_File_write_at(_qx8ertb, (_1zderi4[_20n8tte]+_eaxylnf), _fc4kl8c, _6zjq7w7[i*4 +3], ((MPI_Datatype)0x4c00040a), &_py6b2kn); _eaxylnf += _6zjq7w7[i*4 +3]*sizeof(float);
                    MPI_File_write_at(_qx8ertb, (_1zderi4[_20n8tte]+_eaxylnf), _gn0t6dv, _6zjq7w7[i*4 +3], ((MPI_Datatype)0x4c00040a), &_py6b2kn); _eaxylnf += _6zjq7w7[i*4 +3]*sizeof(float);
            } }
            free(_fc4kl8c); free(_gn0t6dv);
            MPI_Barrier( ((MPI_Comm)0x44000000) );
            MPI_File_close(&_qx8ertb);
            double _eimvihc;
            clock_gettime(_CLOCK_REALTIME, &_rxyj38x);
            _eimvihc = (double)_avt8k06 *(_rxyj38x.tv_sec - _zwd0w9r.tv_sec) + (_rxyj38x.tv_nsec - _zwd0w9r.tv_nsec)/1.0E+9;
            if (_20n8tte == 0) { fprintf(__stdoutp,"CPU time (hours) for STF event %s:  %lf\n",_1750ceh, _eimvihc/3600.0); }
        }
        for (i = 0u; i < _yt91qe3; i++) { __builtin___memset_chk (_ntdr5xy[_ddqhkga[i]], 0, 2*_5zja8dz._148p0tv*sizeof(float), __builtin_object_size (_ntdr5xy[_ddqhkga[i]], 0)); }
        if (_z3nx39r._qsq0p82 == 1u)
        {
            for (i = 0u; i < _yt91qe3; i++)
            {
                ReAssignElemVals(&_5zja8dz, &_z3nx39r, _ddqhkga[i], _6zjq7w7, _3n3bzt9, _mzwffhh, 0);
        } }
        ResetAfterEvent(_20n8tte, _b6gs8g8, &_5zja8dz, _6zjq7w7, _3n3bzt9, _mzwffhh, _3p91i3v);
        clock_gettime(_CLOCK_REALTIME, &_7upfln4);
        _1zlbori += (_7upfln4.tv_sec - _tkezjss.tv_sec) + (_7upfln4.tv_nsec - _tkezjss.tv_nsec)/1.0E+9;
    }
    clock_gettime(_CLOCK_REALTIME, &_6nouhfj);
    _0s1f537 = (_6nouhfj.tv_sec - _7dknncp.tv_sec) + (_6nouhfj.tv_nsec - _7dknncp.tv_nsec)/1.0E+9;
    if (_20n8tte == 0)
    {
        _ds3qsqc /= (double)_jzr25md;
        fprintf(__stdoutp,"\nTotal RunTime for EQcycle (minutes): %6.4f       and mean inter-event time (days): %3.2lf\n",_0s1f537/60.0f, _ds3qsqc*365.25 );
        fprintf(__stdoutp,"Time spent in interseismic: %4.2lf   Time spent in coseismic: %4.2lf\n",_q8hke22/60.0f, _1zlbori/60.0f);
        fprintf(__stdoutp,"\n average iteration number: %f\n",(float)_w2pcn9t*_dswm31y/_zfk8ut0);
        fclose(_szbffzz);
    }
    MPI_Barrier( ((MPI_Comm)0x44000000) );
    free(_mygty3a); free(_lx2toow); free(_6zjq7w7); free(_3n3bzt9); free(_mzwffhh); free(_3p91i3v); free(_m7pyxdw); free(_qxmcthp); free(_32orwix); free(_7c6qc6k); free(_8djzr5s);
    free(_4pa8hye); free(_hj3muur); free(_wd1xmxb); free(_z2c9ros); free(_795e5iz); free(_vlnknd6); free(_k00ilu6); free(_28xg8y7); free(_7odkef8); free(_v71448n);
    free(_j3swebn); free(_vg5vmf6); free(_0yknfu8); free(_4xyjoa2); free(_ddqhkga); free(_vfua9bq); free(_2ivglye);
    for (i = 0; i < _b6gs8g8[_20n8tte]; i++)
    { free(_ntdr5xy[i]);
        free(_swbkok3[i]); free(_zd76szn[i]); free(_ukwqvr1[i]); free(_m458ax4[i]);
    }
    if (_z3nx39r._urvqt3g > 0u)
    { for (i = 0; i < _b6gs8g8[_20n8tte]; i++)
        { free(_c795rlq[i]); free(_plozjcn[i]); free(_0isb5iq[i]);
        }
        for (i = 0; i < _zfflb7w[_20n8tte]; i++)
        { free(_25lzse9[i]); free(_dg1wjui[i]); free(_ozctwjl[i]); free(_3jczucq[i]); free(_b68rerv[i]); free(_ycz04cd[i]); free(_ny181ac[i]); free(_zzf8s2z[i]);
        }
    }
    free(_0fqldb2); free(_hlmh1t9); free(_tacy15p); free(_xrum2lh); free(_ntdr5xy);
    free(_0vs52sq); free(_swbkok3); free(_zd76szn); free(_ukwqvr1); free(_m458ax4);
    free(_zn4regb); free(_c795rlq); free(_plozjcn); free(_0isb5iq);
    free(_n8wwiif); free(_25lzse9); free(_ozctwjl); free(_b68rerv); free(_ny181ac);
    free(_el5015o); free(_dg1wjui); free(_3jczucq); free(_ycz04cd); free(_zzf8s2z);
    MPI_Barrier( ((MPI_Comm)0x44000000) );
    MPI_Finalize( );
    return 0;
}
