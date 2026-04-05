
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

void _25sq8cr(double _5lf6a4f[6], double _e8f8409[6], double X, double Y, double Z, double _atj9tze[3], double _b7r3ocw[3], double _v2isqp4[3], double _wyzx4w6, double _frhocn9, double _yrybova, double _r965r96, double _xj4pnjy);
void _oyr6ccv(double _4qxxva2[6],double _0get8t6[6], double X, double Y, double Z, double _atj9tze[3] ,double _b7r3ocw[3], double _v2isqp4[3], double _wyzx4w6, double _frhocn9, double _yrybova, double _r965r96, double _xj4pnjy);
void _tb94874(double _yr06nkw[3], double _kvl9ptk, double _112hxrj, double _mdmmmxi, double _s3kkv9r[3][3]);
void _hpmnz1u(int _rpj1epk[1], double x,double y,double z,double _jz0ifrm,double _l29q2l4, double _33e1i1p,double _0bq1s9r,double _jwh55he,double _wlgje2a);
void _z8vmuhx(double _ihx74d0[6], double _z36u8ue[6], double B[3][3]);
void _r2sds8b(double x,double y,double z,double _tj1gm25,double _wk8r2yk,double _ccprhvl,double _lqjidkn,double _gdc5tao, double _dbfwyjx[3],double _vxpj4hn[3],double e[6]);
void _nnhj7xv(double x, double y, double z, double _tj1gm25, double _wk8r2yk, double _ccprhvl, double _lqjidkn, double _gdc5tao, double e[6]);
void _g0fn9xu(double _p0fodro[6],double _yyv37di[6], double X,double Y,double Z,double _z7w2seh,double _e2atn2t,double _wdl9zs1,double _ll8os9x[3], double _togn1w1[3], double _r965r96,double _xj4pnjy);
void _t3w5uzb(double _4gc4ryg, double _ovd79c9, double _mv9e40k, double _52xg8p6, double _jjq9k6r, double _e4q6qdy, double _4lmct64, double _gdc5tao, double a, double _82n2atq[6]);
void StrainHS_Nikkhoo(float _d5eqy08[6], float _82n2atq[6], float _ji9kqdg, float _b4xgmpw, float _5d3go8j, float _2c9t199[3], float _dxj37c9[3], float _h42m44g[3], float _dtg5nqj, float _gyspeng, float _0vgnxb5, const float _8bixoiq, const float _yn4u9bx)
{
    double X = (double)_ji9kqdg; double Y = (double)_b4xgmpw; double Z = (double)_5d3go8j;
    double _wyzx4w6 = (double)_dtg5nqj; double _frhocn9= (double)_gyspeng; double _yrybova= (double)_0vgnxb5;
    double _r965r96 = (double)_8bixoiq; double _xj4pnjy= (double)_yn4u9bx;
    double _atj9tze[3], _b7r3ocw[3], _v2isqp4[3];
    _atj9tze[0] = (double)_2c9t199[0]; _atj9tze[1] = (double)_2c9t199[1]; _atj9tze[2] = (double)_2c9t199[2];
    _b7r3ocw[0] = (double)_dxj37c9[0]; _b7r3ocw[1] = (double)_dxj37c9[1]; _b7r3ocw[2] = (double)_dxj37c9[2];
    _v2isqp4[0] = (double)_h42m44g[0]; _v2isqp4[1] = (double)_h42m44g[1]; _v2isqp4[2] = (double)_h42m44g[2];
    int i;
    double _5lf6a4f[6], _4qxxva2[6], _ieuc1hh[6];
    double _e8f8409[6], _0get8t6[6], _3j3r30m[6];
    double _5r5uqpg[3], _bfuk9e1[3], _9kwv08h[3];
    _5r5uqpg[0] = _atj9tze[0]; _bfuk9e1[0] = _b7r3ocw[0]; _9kwv08h[0] = _v2isqp4[0];
    _5r5uqpg[1] = _atj9tze[1]; _bfuk9e1[1] = _b7r3ocw[1]; _9kwv08h[1] = _v2isqp4[1];
    _5r5uqpg[2] = -1.0*_atj9tze[2]; _bfuk9e1[2] = -1.0*_b7r3ocw[2]; _9kwv08h[2] = -1.0*_v2isqp4[2];
    for (i = 0; i < 6; i++)
    { _5lf6a4f[i] = 0.0; _e8f8409[i] = 0.0;
        _4qxxva2[i] = 0.0; _0get8t6[i] = 0.0;
        _ieuc1hh[i] = 0.0; _3j3r30m[i] = 0.0;
    }
    if ((Z > 0.0) || (_atj9tze[2] > 0.0) || (_b7r3ocw[2] > 0.0) || (_v2isqp4[2] > 0.0))
    { fprintf(__stdoutp,"A triangle vertex or the center of testing location have positive z-position (above half-space) => abort\n");
        fprintf(__stdoutp,"Z: %f    P1[2]: %f    P2[2]: %f      P3[2]: %f   \n",Z, _atj9tze[2], _b7r3ocw[2],_v2isqp4[2]);
        exit(10);
    }
    if ((_atj9tze[2] == 0.0) && (_b7r3ocw[2] == 0.0) && (_v2isqp4[2] == 0.0))
    { for (i = 0; i < 6; i++)
        { _d5eqy08[i] = 0.0;
            _82n2atq[i] = 0.0;
    } }
    else
    {
        _25sq8cr( _5lf6a4f, _e8f8409, X, Y, Z, _atj9tze, _b7r3ocw, _v2isqp4, _wyzx4w6, _frhocn9, _yrybova, _r965r96, _xj4pnjy);
        _oyr6ccv(_4qxxva2, _0get8t6, X, Y, Z, _atj9tze, _b7r3ocw, _v2isqp4, _wyzx4w6, _frhocn9, _yrybova, _r965r96, _xj4pnjy);
        _25sq8cr( _ieuc1hh, _3j3r30m, X, Y, Z, _5r5uqpg, _bfuk9e1, _9kwv08h, _wyzx4w6, _frhocn9, _yrybova, _r965r96, _xj4pnjy);
        if ((_5r5uqpg[2] == 0.0) && (_bfuk9e1[2] == 0.0) && (_9kwv08h[2] == 0.0))
        { _ieuc1hh[4] = -1.0*_ieuc1hh[4]; _ieuc1hh[5] = -1.0*_ieuc1hh[5];
            _3j3r30m[4] = -1.0*_3j3r30m[4]; _3j3r30m[5] = -1.0*_3j3r30m[5];
        }
        for (i = 0; i < 6; i++)
        { _d5eqy08[i] = (_5lf6a4f[i] +_ieuc1hh[i] +_4qxxva2[i]);
            _82n2atq[i] = (_e8f8409[i] +_3j3r30m[i] +_0get8t6[i]);
    } }
    return;
}
void _25sq8cr(double _5lf6a4f[6], double _e8f8409[6], double X, double Y, double Z, double _atj9tze[3], double _b7r3ocw[3], double _v2isqp4[3], double _ccprhvl, double _lqjidkn, double _wk8r2yk, double _r965r96, double _xj4pnjy)
{
    int _t7zewvi, _kzmpjrh;
    double A, B, C, x, y, z, _gdc5tao, _3dd03uv;
    int _3jgzlma[1];
    double _m1enf31[3], _gippsfn[3], _c4d145b[3];
    double _laym0wb[3], _glmy62x[3];
    double _wgabtgd[3], _8c1cqg4[3], _yllwb7s[3];
    double _h5csk30[3], _lvylkfx[3], _hkc856w[3];
    double _pr97y2l[3], _ofn1r4k[6];
    double _gjyv5g1[6], _0ur5l8a[6], _78e9riv[6];
    double _4v30t9n[3][3];
    _gdc5tao = _xj4pnjy/(2.0*(_xj4pnjy+_r965r96));
    _wgabtgd[0] = 0.0; _wgabtgd[1] = 0.0; _wgabtgd[2] = 0.0;
    _8c1cqg4[0] = 0.0; _8c1cqg4[1] = 0.0; _8c1cqg4[2] = 0.0;
    _yllwb7s[0] = 0.0; _yllwb7s[1] = 0.0; _yllwb7s[2] = 0.0;
    _pr97y2l[0] = 0.0; _pr97y2l[1] = 0.0; _pr97y2l[2] = 1.0;
    _laym0wb[0] = _b7r3ocw[0] -_atj9tze[0]; _laym0wb[1] = _b7r3ocw[1] -_atj9tze[1]; _laym0wb[2] = _b7r3ocw[2] -_atj9tze[2];
    _glmy62x[0] = _v2isqp4[0] -_atj9tze[0]; _glmy62x[1] = _v2isqp4[1] -_atj9tze[1]; _glmy62x[2] = _v2isqp4[2] -_atj9tze[2];
    _c4d145b[0] = _laym0wb[1]*_glmy62x[2] - _laym0wb[2]*_glmy62x[1];
    _c4d145b[1] = _laym0wb[2]*_glmy62x[0] - _laym0wb[0]*_glmy62x[2];
    _c4d145b[2] = _laym0wb[0]*_glmy62x[1] - _laym0wb[1]*_glmy62x[0];
    _3dd03uv = sqrt( _c4d145b[0]*_c4d145b[0] +_c4d145b[1]*_c4d145b[1] + _c4d145b[2]*_c4d145b[2]);
    _c4d145b[0] = _c4d145b[0]/_3dd03uv; _c4d145b[1] = _c4d145b[1]/_3dd03uv; _c4d145b[2] = _c4d145b[2]/_3dd03uv;
    _m1enf31[0] = _pr97y2l[1]*_c4d145b[2] - _pr97y2l[2]*_c4d145b[1];
    _m1enf31[1] = _pr97y2l[2]*_c4d145b[0] - _pr97y2l[0]*_c4d145b[2];
    _m1enf31[2] = _pr97y2l[0]*_c4d145b[1] - _pr97y2l[1]*_c4d145b[0];
    _3dd03uv = sqrt( _m1enf31[0]*_m1enf31[0] +_m1enf31[1]*_m1enf31[1] + _m1enf31[2]*_m1enf31[2]);
    if (_3dd03uv < 1.19209290e-7F)
    { _m1enf31[0] = 0.0 ; _m1enf31[1] = 1.0; _m1enf31[2] = 0.0;
    }
    else
    { _m1enf31[0]/= _3dd03uv; _m1enf31[1]/= _3dd03uv; _m1enf31[2]/= _3dd03uv;
    }
    _gippsfn[0] = _c4d145b[1]*_m1enf31[2] - _c4d145b[2]*_m1enf31[1];
    _gippsfn[1] = _c4d145b[2]*_m1enf31[0] - _c4d145b[0]*_m1enf31[2];
    _gippsfn[2] = _c4d145b[0]*_m1enf31[1] - _c4d145b[1]*_m1enf31[0];
    _3dd03uv = sqrt( _gippsfn[0]*_gippsfn[0] +_gippsfn[1]*_gippsfn[1] + _gippsfn[2]*_gippsfn[2]);
    _gippsfn[0] = _gippsfn[0]/_3dd03uv; _gippsfn[1] = _gippsfn[1]/_3dd03uv; _gippsfn[2] = _gippsfn[2]/_3dd03uv;
    _4v30t9n[0][0] = _c4d145b[0]; _4v30t9n[0][1] = _c4d145b[1]; _4v30t9n[0][2] = _c4d145b[2];
    _4v30t9n[1][0] = _m1enf31[0]; _4v30t9n[1][1] = _m1enf31[1]; _4v30t9n[1][2] = _m1enf31[2];
    _4v30t9n[2][0] = _gippsfn[0]; _4v30t9n[2][1] = _gippsfn[1]; _4v30t9n[2][2] = _gippsfn[2];
    _tb94874(_laym0wb, (X-_b7r3ocw[0]), (Y-_b7r3ocw[1]), (Z-_b7r3ocw[2]), _4v30t9n);
    x = _laym0wb[0]; y = _laym0wb[1]; z = _laym0wb[2];
    _tb94874(_laym0wb, (_atj9tze[0]-_b7r3ocw[0]), (_atj9tze[1]-_b7r3ocw[1]), (_atj9tze[2]-_b7r3ocw[2]), _4v30t9n);
    _wgabtgd[0] = _laym0wb[0]; _wgabtgd[1] = _laym0wb[1]; _wgabtgd[2] = _laym0wb[2];
    _tb94874(_glmy62x, (_v2isqp4[0]-_b7r3ocw[0]), (_v2isqp4[1]-_b7r3ocw[1]), (_v2isqp4[2]-_b7r3ocw[2]), _4v30t9n);
    _yllwb7s[0] = _glmy62x[0]; _yllwb7s[1] = _glmy62x[1]; _yllwb7s[2] = _glmy62x[2];
    _3dd03uv = sqrt( (_8c1cqg4[0]-_wgabtgd[0])*(_8c1cqg4[0]-_wgabtgd[0]) +(_8c1cqg4[1]-_wgabtgd[1])*(_8c1cqg4[1]-_wgabtgd[1]) + (_8c1cqg4[2]-_wgabtgd[2])*(_8c1cqg4[2]-_wgabtgd[2]));
    _h5csk30[0] = (_8c1cqg4[0]-_wgabtgd[0])/_3dd03uv; _h5csk30[1] = (_8c1cqg4[1]-_wgabtgd[1])/_3dd03uv; _h5csk30[2] = (_8c1cqg4[2]-_wgabtgd[2])/_3dd03uv;
    _3dd03uv = sqrt( (_yllwb7s[0]-_wgabtgd[0])*(_yllwb7s[0]-_wgabtgd[0]) +(_yllwb7s[1]-_wgabtgd[1])*(_yllwb7s[1]-_wgabtgd[1]) + (_yllwb7s[2]-_wgabtgd[2])*(_yllwb7s[2]-_wgabtgd[2]));
    _lvylkfx[0] = (_yllwb7s[0]-_wgabtgd[0])/_3dd03uv; _lvylkfx[1] = (_yllwb7s[1]-_wgabtgd[1])/_3dd03uv; _lvylkfx[2] = (_yllwb7s[2]-_wgabtgd[2])/_3dd03uv;
    _3dd03uv = sqrt( (_yllwb7s[0]-_8c1cqg4[0])*(_yllwb7s[0]-_8c1cqg4[0]) +(_yllwb7s[1]-_8c1cqg4[1])*(_yllwb7s[1]-_8c1cqg4[1]) + (_yllwb7s[2]-_8c1cqg4[2])*(_yllwb7s[2]-_8c1cqg4[2]));
    _hkc856w[0] = (_yllwb7s[0]-_8c1cqg4[0])/_3dd03uv; _hkc856w[1] = (_yllwb7s[1]-_8c1cqg4[1])/_3dd03uv; _hkc856w[2] = (_yllwb7s[2]-_8c1cqg4[2])/_3dd03uv;
    _3dd03uv = _h5csk30[0]*_lvylkfx[0] + _h5csk30[1]*_lvylkfx[1] + _h5csk30[2]*_lvylkfx[2];
    A = acos(_3dd03uv) ;
    _3dd03uv = -1.0*_h5csk30[0]*_hkc856w[0] + -1.0*_h5csk30[1]*_hkc856w[1] + -1.0*_h5csk30[2]*_hkc856w[2];
    B = acos(_3dd03uv) ;
    _3dd03uv = _hkc856w[0]*_lvylkfx[0] + _hkc856w[1]*_lvylkfx[1] + _hkc856w[2]*_lvylkfx[2];
    C = acos(_3dd03uv) ;
    _hpmnz1u(_3jgzlma,y,z,x,_wgabtgd[1],_wgabtgd[2], _8c1cqg4[1], _8c1cqg4[2], _yllwb7s[1], _yllwb7s[2]);
    if (_3jgzlma[0] == 1) { _t7zewvi = 1; _kzmpjrh = 0; }
    if (_3jgzlma[0] ==-1) { _t7zewvi = 0; _kzmpjrh = 1; }
    if (_3jgzlma[0] == 0) { _t7zewvi = 0; _kzmpjrh = 0; }
    if (_t7zewvi == 1)
    {
        _laym0wb[0] = -1.0*_lvylkfx[0]; _laym0wb[1] = -1.0*_lvylkfx[1]; _laym0wb[2] = -1.0*_lvylkfx[2];
        _r2sds8b(x, y, z, A, _wk8r2yk, _ccprhvl, _lqjidkn, _gdc5tao, _wgabtgd, _laym0wb, _gjyv5g1);
        _r2sds8b(x, y, z, B, _wk8r2yk, _ccprhvl, _lqjidkn, _gdc5tao, _8c1cqg4, _h5csk30, _0ur5l8a);
        _r2sds8b(x, y, z, C, _wk8r2yk, _ccprhvl, _lqjidkn, _gdc5tao, _yllwb7s, _hkc856w, _78e9riv);
    }
    if (_kzmpjrh == 1)
    {
        _r2sds8b(x, y, z, A, _wk8r2yk, _ccprhvl, _lqjidkn, _gdc5tao, _wgabtgd, _lvylkfx, _gjyv5g1);
        _laym0wb[0] = -1.0*_h5csk30[0]; _laym0wb[1] = -1.0*_h5csk30[1]; _laym0wb[2] = -1.0*_h5csk30[2];
        _r2sds8b(x, y, z, B, _wk8r2yk, _ccprhvl, _lqjidkn, _gdc5tao, _8c1cqg4, _laym0wb, _0ur5l8a);
        _laym0wb[0] = -1.0*_hkc856w[0]; _laym0wb[1] = -1.0*_hkc856w[1]; _laym0wb[2] = -1.0*_hkc856w[2];
        _r2sds8b(x, y, z, C, _wk8r2yk, _ccprhvl, _lqjidkn, _gdc5tao, _yllwb7s, _laym0wb, _78e9riv);
    }
    if ((_kzmpjrh == 1) || (_t7zewvi == 1))
    {
        _ofn1r4k[0] = _gjyv5g1[0]+_0ur5l8a[0]+_78e9riv[0];
        _ofn1r4k[1] = _gjyv5g1[1]+_0ur5l8a[1]+_78e9riv[1];
        _ofn1r4k[2] = _gjyv5g1[2]+_0ur5l8a[2]+_78e9riv[2];
        _ofn1r4k[3] = _gjyv5g1[3]+_0ur5l8a[3]+_78e9riv[3];
        _ofn1r4k[4] = _gjyv5g1[4]+_0ur5l8a[4]+_78e9riv[4];
        _ofn1r4k[5] = _gjyv5g1[5]+_0ur5l8a[5]+_78e9riv[5];
    }
    else
    {
        _ofn1r4k[0] = (__builtin_nanf(""));
        _ofn1r4k[1] = (__builtin_nanf(""));
        _ofn1r4k[2] = (__builtin_nanf(""));
        _ofn1r4k[3] = (__builtin_nanf(""));
        _ofn1r4k[4] = (__builtin_nanf(""));
        _ofn1r4k[5] = (__builtin_nanf(""));
    }
    _4v30t9n[0][0] = _c4d145b[0]; _4v30t9n[0][1] = _m1enf31[0]; _4v30t9n[0][2] = _gippsfn[0];
    _4v30t9n[1][0] = _c4d145b[1]; _4v30t9n[1][1] = _m1enf31[1]; _4v30t9n[1][2] = _gippsfn[1];
    _4v30t9n[2][0] = _c4d145b[2]; _4v30t9n[2][1] = _m1enf31[2]; _4v30t9n[2][2] = _gippsfn[2];
    _z8vmuhx(_e8f8409, _ofn1r4k, _4v30t9n);
    _5lf6a4f[0] = 2.0*_r965r96*_e8f8409[0]+_xj4pnjy*(_e8f8409[0]+_e8f8409[3]+_e8f8409[5]);
    _5lf6a4f[3] = 2.0*_r965r96*_e8f8409[3]+_xj4pnjy*(_e8f8409[0]+_e8f8409[3]+_e8f8409[5]);
    _5lf6a4f[5] = 2.0*_r965r96*_e8f8409[5]+_xj4pnjy*(_e8f8409[0]+_e8f8409[3]+_e8f8409[5]);
    _5lf6a4f[1] = 2.0*_r965r96*_e8f8409[1];
    _5lf6a4f[2] = 2.0*_r965r96*_e8f8409[2];
    _5lf6a4f[4] = 2.0*_r965r96*_e8f8409[4];
    return;
}
void _oyr6ccv(double _4qxxva2[6], double _0get8t6[6], double X, double Y, double Z, double _atj9tze[3], double _b7r3ocw[3], double _v2isqp4[3], double _ccprhvl, double _lqjidkn, double _wk8r2yk, double _r965r96, double _xj4pnjy)
{
    int i;
    double _laym0wb[3], _glmy62x[3];
    double _p0fodro[6], _8bnp84d[6], _yt9moge[6];
    double _yyv37di[6], _a5u6hc3[6], _6kof8kr[6];
    double _4v30t9n[3][3];
    double _3dd03uv, _pr97y2l[3];
    double _c4d145b[3], _m1enf31[3], _gippsfn[3];
    _pr97y2l[0] = 0.0; _pr97y2l[1] = 0.0; _pr97y2l[2] = 1.0;
    _laym0wb[0] = _b7r3ocw[0] -_atj9tze[0]; _laym0wb[1] = _b7r3ocw[1] -_atj9tze[1]; _laym0wb[2] = _b7r3ocw[2] -_atj9tze[2];
    _glmy62x[0] = _v2isqp4[0] -_atj9tze[0]; _glmy62x[1] = _v2isqp4[1] -_atj9tze[1]; _glmy62x[2] = _v2isqp4[2] -_atj9tze[2];
    _c4d145b[0] = _laym0wb[1]*_glmy62x[2] - _laym0wb[2]*_glmy62x[1];
    _c4d145b[1] = _laym0wb[2]*_glmy62x[0] - _laym0wb[0]*_glmy62x[2];
    _c4d145b[2] = _laym0wb[0]*_glmy62x[1] - _laym0wb[1]*_glmy62x[0];
    _3dd03uv = sqrt( _c4d145b[0]*_c4d145b[0] +_c4d145b[1]*_c4d145b[1] + _c4d145b[2]*_c4d145b[2]);
    _c4d145b[0] = _c4d145b[0]/_3dd03uv; _c4d145b[1] = _c4d145b[1]/_3dd03uv; _c4d145b[2] = _c4d145b[2]/_3dd03uv;
    _m1enf31[0] = _pr97y2l[1]*_c4d145b[2] - _pr97y2l[2]*_c4d145b[1];
    _m1enf31[1] = _pr97y2l[2]*_c4d145b[0] - _pr97y2l[0]*_c4d145b[2];
    _m1enf31[2] = _pr97y2l[0]*_c4d145b[1] - _pr97y2l[1]*_c4d145b[0];
    _3dd03uv = sqrt( _m1enf31[0]*_m1enf31[0] +_m1enf31[1]*_m1enf31[1] + _m1enf31[2]*_m1enf31[2]);
    if (_3dd03uv < 1.19209290e-7F)
    { _m1enf31[0] = 0.0 ; _m1enf31[1] = 1.0; _m1enf31[2] = 0.0;
    }
    else
    { _m1enf31[0]/= _3dd03uv; _m1enf31[1]/= _3dd03uv; _m1enf31[2]/= _3dd03uv;
    }
    _gippsfn[0] = _c4d145b[1]*_m1enf31[2] - _c4d145b[2]*_m1enf31[1];
    _gippsfn[1] = _c4d145b[2]*_m1enf31[0] - _c4d145b[0]*_m1enf31[2];
    _gippsfn[2] = _c4d145b[0]*_m1enf31[1] - _c4d145b[1]*_m1enf31[0];
    _3dd03uv = sqrt( _gippsfn[0]*_gippsfn[0] +_gippsfn[1]*_gippsfn[1] + _gippsfn[2]*_gippsfn[2]);
    _gippsfn[0] = _gippsfn[0]/_3dd03uv; _gippsfn[1] = _gippsfn[1]/_3dd03uv; _gippsfn[2] = _gippsfn[2]/_3dd03uv;
    _4v30t9n[0][0] = _c4d145b[0]; _4v30t9n[0][1] = _m1enf31[0]; _4v30t9n[0][2] = _gippsfn[0];
    _4v30t9n[1][0] = _c4d145b[1]; _4v30t9n[1][1] = _m1enf31[1]; _4v30t9n[1][2] = _gippsfn[1];
    _4v30t9n[2][0] = _c4d145b[2]; _4v30t9n[2][1] = _m1enf31[2]; _4v30t9n[2][2] = _gippsfn[2];
    _tb94874(_laym0wb, _wk8r2yk, _ccprhvl, _lqjidkn, _4v30t9n);
    _g0fn9xu(_p0fodro,_yyv37di, X,Y,Z,_laym0wb[0],_laym0wb[1],_laym0wb[2],_atj9tze,_b7r3ocw,_r965r96,_xj4pnjy);
    _g0fn9xu(_8bnp84d,_a5u6hc3, X,Y,Z,_laym0wb[0],_laym0wb[1],_laym0wb[2],_b7r3ocw,_v2isqp4,_r965r96,_xj4pnjy);
    _g0fn9xu(_yt9moge,_6kof8kr, X,Y,Z,_laym0wb[0],_laym0wb[1],_laym0wb[2],_v2isqp4,_atj9tze,_r965r96,_xj4pnjy);
    for (i = 0; i < 6; i++)
    { _4qxxva2[i] = _p0fodro[i] + _8bnp84d[i] + _yt9moge[i];
        _0get8t6[i] = _yyv37di[i] + _a5u6hc3[i] + _6kof8kr[i];
    }
    return;
}
void _z8vmuhx(double _ihx74d0[6], double _z36u8ue[6], double B[3][3])
{
    double _khejla0, _5ba75uu, _37d4ivs, _s6m033i;
    double _r11q7uw, _jcslqc3, _gnzs43c, _zpeal2w;
    double _olmlgse, _fhv1xg1, _w5wfcpk, _c54wsqt;
    double A[9];
    _khejla0 = _z36u8ue[0]; _5ba75uu = _z36u8ue[1]; _37d4ivs = _z36u8ue[2];
    _s6m033i = _z36u8ue[3]; _r11q7uw = _z36u8ue[4]; _jcslqc3 = _z36u8ue[5];
    A[0] = B[0][0]; A[1] = B[1][0]; A[2] = B[2][0];
    A[3] = B[0][1]; A[4] = B[1][1]; A[5] = B[2][1];
    A[6] = B[0][2]; A[7] = B[1][2]; A[8] = B[2][2];
    _gnzs43c = A[0]*A[0]*_khejla0 + 2.0*A[0]*A[3] *_5ba75uu + 2.0*A[0]*A[6] *_37d4ivs + 2.0*A[3]*A[6] *_r11q7uw + A[3]*A[3]*_s6m033i + A[6]*A[6]*_jcslqc3;
    _fhv1xg1 = A[1]*A[1]*_khejla0 + 2.0*A[1]*A[4] *_5ba75uu + 2.0*A[1]*A[7] *_37d4ivs + 2.0*A[4]*A[7] *_r11q7uw + A[4]*A[4]*_s6m033i + A[7]*A[7]*_jcslqc3;
    _c54wsqt = A[2]*A[2]*_khejla0 + 2.0*A[2]*A[5] *_5ba75uu + 2.0*A[2]*A[8] *_37d4ivs + 2.0*A[5]*A[8] *_r11q7uw + A[5]*A[5]*_s6m033i + A[8]*A[8]*_jcslqc3;
    _zpeal2w = A[0]*A[1]*_khejla0 + (A[0]*A[4]+ A[1]*A[3])*_5ba75uu + (A[0]*A[7] + A[1]*A[6])*_37d4ivs + (A[7]*A[3] + A[6]*A[4])*_r11q7uw + A[4]*A[3]*_s6m033i + A[6]*A[7]*_jcslqc3;
    _olmlgse = A[0]*A[2]*_khejla0 + (A[0]*A[5]+ A[2]*A[3])*_5ba75uu + (A[0]*A[8] + A[2]*A[6])*_37d4ivs + (A[8]*A[3] + A[6]*A[5])*_r11q7uw + A[5]*A[3]*_s6m033i + A[6]*A[8]*_jcslqc3;
    _w5wfcpk = A[1]*A[2]*_khejla0 + (A[2]*A[4]+ A[1]*A[5])*_5ba75uu + (A[2]*A[7] + A[1]*A[8])*_37d4ivs + (A[7]*A[5] + A[8]*A[4])*_r11q7uw + A[4]*A[5]*_s6m033i + A[7]*A[8]*_jcslqc3;
    _ihx74d0[0] = _gnzs43c; _ihx74d0[1] = _zpeal2w; _ihx74d0[2] = _olmlgse;
    _ihx74d0[3] = _fhv1xg1; _ihx74d0[4] = _w5wfcpk; _ihx74d0[5] = _c54wsqt;
    return;
}
void _tb94874(double _yr06nkw[3], double _znvp8gm, double _nhie9ed, double _e1zuble, double A[3][3])
{
    _yr06nkw[0] = A[0][0]*_znvp8gm + A[0][1]*_nhie9ed + A[0][2]*_e1zuble; _yr06nkw[1] = A[1][0]*_znvp8gm + A[1][1]*_nhie9ed + A[1][2]*_e1zuble; _yr06nkw[2] = A[2][0]*_znvp8gm + A[2][1]*_nhie9ed + A[2][2]*_e1zuble;
    return;
}
void _hpmnz1u(int _3jgzlma[1], double x,double y,double z,double _jz0ifrm,double _l29q2l4, double _33e1i1p,double _0bq1s9r,double _jwh55he,double _wlgje2a)
{
    double a, b, c;
    a = ((_0bq1s9r-_wlgje2a)*(x-_jwh55he) +(_jwh55he-_33e1i1p)*(y-_wlgje2a)) / ((_0bq1s9r-_wlgje2a)*(_jz0ifrm-_jwh55he) +(_jwh55he-_33e1i1p)*(_l29q2l4-_wlgje2a));
    b = ((_wlgje2a-_l29q2l4)*(x-_jwh55he) +(_jz0ifrm-_jwh55he)*(y-_wlgje2a)) / ((_0bq1s9r-_wlgje2a)*(_jz0ifrm-_jwh55he) +(_jwh55he-_33e1i1p)*(_l29q2l4-_wlgje2a));
    c = 1.0 -a -b;
    _3jgzlma[0] = 1;
    if ((a <= 0.0) && (b > c) && (c > a)) { _3jgzlma[0] = -1; }
    if ((a > b) && (b <= 0.0) && (c > a)) { _3jgzlma[0] = -1; }
    if ((a > b) && (b > c) && (c <= 0.0)) { _3jgzlma[0] = -1; }
    if ((a == 0.0) && (b >= 0.0) && (c >= 0.0)) { _3jgzlma[0] = 0; }
    if ((a >= 0.0) && (b == 0.0) && (c >= 0.0)) { _3jgzlma[0] = 0; }
    if ((a >= 0.0) && (b >= 0.0) && (c == 0.0)) { _3jgzlma[0] = 0; }
    if ((_3jgzlma[0] == 0) && (z != 0.0)) { _3jgzlma[0] = 1; }
    return;
}
void _r2sds8b(double x,double y,double z,double _tj1gm25,double _wk8r2yk,double _ccprhvl,double _lqjidkn,double _gdc5tao, double _dbfwyjx[3],double _vxpj4hn[3],double _ihx74d0[6])
{
    double A[2][2];
    double B[3][3];
    double _ln1p0vu[2];
    double _h17j8qo[2];
    double _4gc4ryg;
    double _kzrezrd;
    double _z049enu;
    double _wi5flul;
    double _laym0wb[2];
    double e[6];
    A[0][0] = _vxpj4hn[2]; A[0][1] = -1.0*_vxpj4hn[1];
    A[1][0] = _vxpj4hn[1]; A[1][1] = _vxpj4hn[2];
    _laym0wb[0] = y -_dbfwyjx[1]; _laym0wb[1] = z -_dbfwyjx[2];
    _ln1p0vu[0] = A[0][0]*_laym0wb[0] +A[0][1]*_laym0wb[1];
    _ln1p0vu[1] = A[1][0]*_laym0wb[0] +A[1][1]*_laym0wb[1];
    _4gc4ryg = _ln1p0vu[0]; _kzrezrd = _ln1p0vu[1];
    _laym0wb[0] = _ccprhvl; _laym0wb[1] = _lqjidkn;
    _h17j8qo[0] = A[0][0]*_laym0wb[0] +A[0][1]*_laym0wb[1];
    _h17j8qo[1] = A[1][0]*_laym0wb[0] +A[1][1]*_laym0wb[1];
    _z049enu = _h17j8qo[0]; _wi5flul = _h17j8qo[1];
    _nnhj7xv(x, _4gc4ryg, _kzrezrd, (-1.0*3.14159265358979323846264338327950288 +_tj1gm25), _wk8r2yk, _z049enu, _wi5flul, _gdc5tao, e);
    B[0][0] = 1.0; B[0][1] = 0.0; B[0][2] = 0.0;
    B[1][0] = 0.0; B[1][1] = A[0][0]; B[1][2] = A[1][0];
    B[2][0] = 0.0; B[2][1] = A[0][1]; B[2][2] = A[1][1];
    _z8vmuhx(_ihx74d0, e, B);
    return;
}
void _g0fn9xu(double _d5eqy08[6],double _82n2atq[6], double X,double Y,double Z,double _z7w2seh,double _e2atn2t,double _wdl9zs1,double _hxfgcaq[3], double _vkizzqy[3],double _r965r96,double _xj4pnjy)
{
    int i, I;
    double _gdc5tao, _52xg8p6, _cjwtl9t, _qvdz4x1;
    double _vxpj4hn[3], _pr97y2l[3], _y36ddfg[3], _s3jonou[3], _9iy4x41[3];
    double _zhqimqo[3], _dpxnf7c[3], _rwedbcy[3], _8nyfwdc[3];
    double _uqaaj84[6], _rx6pryi[6], _jp6mkt7[6];
    double A[3][3], _xaji0y6[3][3];
    _gdc5tao = _xj4pnjy/(_r965r96+_xj4pnjy)/2.0;
    _vxpj4hn[0] = _vkizzqy[0]-_hxfgcaq[0]; _vxpj4hn[1] = _vkizzqy[1]-_hxfgcaq[1]; _vxpj4hn[2] = _vkizzqy[2]-_hxfgcaq[2];
    _pr97y2l[0] = 0.0; _pr97y2l[1] = 0.0; _pr97y2l[2] = 1.0;
    _cjwtl9t = sqrt( _vxpj4hn[0]*_vxpj4hn[0] +_vxpj4hn[1]*_vxpj4hn[1] +_vxpj4hn[2]*_vxpj4hn[2]);
    _qvdz4x1 = -1.0*(_vxpj4hn[0]*_pr97y2l[0] +_vxpj4hn[1]*_pr97y2l[1] +_vxpj4hn[2]*_pr97y2l[2]);
    _52xg8p6 = acos(_qvdz4x1/_cjwtl9t);
    if (fabs(cos(_52xg8p6)/sin(_52xg8p6)) > 5.0e+5*(3.14159265358979323846264338327950288/360.0))
    { for (i = 0; i < 6; i++) { _d5eqy08[i] = 0.0; _82n2atq[i] = 0.0; }
    }
    else
    {
        _y36ddfg[0] = _vxpj4hn[0]; _y36ddfg[1] = _vxpj4hn[1]; _y36ddfg[2] = 0.0;
        _cjwtl9t = sqrt( _y36ddfg[0]*_y36ddfg[0] +_y36ddfg[1]*_y36ddfg[1] +_y36ddfg[2]*_y36ddfg[2]);
        _y36ddfg[0] /= _cjwtl9t; _y36ddfg[1] /= _cjwtl9t; _y36ddfg[2] /= _cjwtl9t;
        _9iy4x41[0] = -1.0*_pr97y2l[0]; _9iy4x41[1] = -1.0*_pr97y2l[1]; _9iy4x41[2] = -1.0*_pr97y2l[2];
        _s3jonou[0] = _9iy4x41[1]*_y36ddfg[2] -_9iy4x41[2]*_y36ddfg[1];
        _s3jonou[1] = _9iy4x41[2]*_y36ddfg[0] -_9iy4x41[0]*_y36ddfg[2];
        _s3jonou[2] = _9iy4x41[0]*_y36ddfg[1] -_9iy4x41[1]*_y36ddfg[0];
        _cjwtl9t = sqrt( _s3jonou[0]*_s3jonou[0] +_s3jonou[1]*_s3jonou[1] +_s3jonou[2]*_s3jonou[2]);
        _s3jonou[0] /= _cjwtl9t; _s3jonou[1] /= _cjwtl9t; _s3jonou[2] /= _cjwtl9t;
        A[0][0] = _y36ddfg[0]; A[0][1] = _s3jonou[0]; A[0][2] = _9iy4x41[0];
        A[1][0] = _y36ddfg[1]; A[1][1] = _s3jonou[1]; A[1][2] = _9iy4x41[1];
        A[2][0] = _y36ddfg[2]; A[2][1] = _s3jonou[2]; A[2][2] = _9iy4x41[2];
        _xaji0y6[0][0] = _y36ddfg[0]; _xaji0y6[0][1] = _y36ddfg[1]; _xaji0y6[0][2] = _y36ddfg[2];
        _xaji0y6[1][0] = _s3jonou[0]; _xaji0y6[1][1] = _s3jonou[1]; _xaji0y6[1][2] = _s3jonou[2];
        _xaji0y6[2][0] = _9iy4x41[0]; _xaji0y6[2][1] = _9iy4x41[1]; _xaji0y6[2][2] = _9iy4x41[2];
        _tb94874(_rwedbcy, (X-_hxfgcaq[0]), (Y-_hxfgcaq[1]), (Z-_hxfgcaq[2]), A);
        _tb94874(_zhqimqo, _vxpj4hn[0], _vxpj4hn[1], _vxpj4hn[2], A);
        _8nyfwdc[0] = _rwedbcy[0] - _zhqimqo[0];
        _8nyfwdc[1] = _rwedbcy[1] - _zhqimqo[1];
        _8nyfwdc[2] = _rwedbcy[2] - _zhqimqo[2];
        _tb94874(_dpxnf7c, _z7w2seh, _e2atn2t, _wdl9zs1, A);
        I = (_52xg8p6*_rwedbcy[0]) >=0.0 ? 1 : 0;
        for (i = 0; i < 6; i++) { _uqaaj84[i] = 0.0; _rx6pryi[i] = 0.0; }
        if (I == 1)
        { _t3w5uzb(_rwedbcy[0],_rwedbcy[1], _rwedbcy[2], (-3.14159265358979323846264338327950288 +_52xg8p6),_dpxnf7c[0],_dpxnf7c[1], _dpxnf7c[2], _gdc5tao,(-1.0*_hxfgcaq[2]), _uqaaj84);
            _t3w5uzb(_8nyfwdc[0],_8nyfwdc[1], _8nyfwdc[2], (-3.14159265358979323846264338327950288 +_52xg8p6),_dpxnf7c[0],_dpxnf7c[1], _dpxnf7c[2], _gdc5tao,(-1.0*_vkizzqy[2]), _rx6pryi);
        }
        else if (I == 0)
        { _t3w5uzb(_rwedbcy[0],_rwedbcy[1], _rwedbcy[2], _52xg8p6,_dpxnf7c[0],_dpxnf7c[1], _dpxnf7c[2], _gdc5tao,(-1.0*_hxfgcaq[2]), _uqaaj84);
            _t3w5uzb(_8nyfwdc[0],_8nyfwdc[1], _8nyfwdc[2], _52xg8p6,_dpxnf7c[0],_dpxnf7c[1], _dpxnf7c[2], _gdc5tao,(-1.0*_vkizzqy[2]), _rx6pryi);
        }
         for (i = 0; i < 6; i++) { _jp6mkt7[i] = _rx6pryi[i] - _uqaaj84[i]; }
         _z8vmuhx(_82n2atq, _jp6mkt7, _xaji0y6);
         _d5eqy08[0] = 2.0*_r965r96*_82n2atq[0] +_xj4pnjy*(_82n2atq[0] +_82n2atq[3] +_82n2atq[5]);
              _d5eqy08[3] = 2.0*_r965r96*_82n2atq[3] +_xj4pnjy*(_82n2atq[0] +_82n2atq[3] +_82n2atq[5]);
         _d5eqy08[5] = 2.0*_r965r96*_82n2atq[5] +_xj4pnjy*(_82n2atq[0] +_82n2atq[3] +_82n2atq[5]);
         _d5eqy08[1] = 2.0*_r965r96*_82n2atq[1];
         _d5eqy08[2] = 2.0*_r965r96*_82n2atq[2];
          _d5eqy08[4] = 2.0*_r965r96*_82n2atq[4];
    }
    return;
}
void _nnhj7xv(double x, double y, double z, double _tj1gm25, double _wk8r2yk, double _ccprhvl, double _lqjidkn, double _gdc5tao, double e[6])
{ double _qnsn2s8, _u1g463s, _nwywrra, _lq5u6jz;
    double _nhie9ed, _ovd79c9, _or8wm9l, _h17j8qo;
    double r, _h4a3yjt, _j4x5a5x, _09vs9kk;
    double _u7fbmen, W, _k0qtl3q, _ioxd9cs;
    double _4n89qkn, _s9toqw9, _s832ejx, C;
    double S, _ui3do70, _02g3yu6, _vk0ucc1;
    double _zhwjr64, _dpja1wj, _0tpac56, _t4ufel6;
    double _246hid2, _5owf759;
    double _lz15bcw, _mv9e40k, _gv44ftj, _xlul6lg;
    double _8drc11e, _dygflij;
    _qnsn2s8 = sin(_tj1gm25); _u1g463s = cos(_tj1gm25);
    _nwywrra = y*_u1g463s - z*_qnsn2s8; _lq5u6jz = y*_qnsn2s8 + z*_u1g463s;
    _nhie9ed = x*x; _ovd79c9 = y*y;
    _or8wm9l = z*z; _h17j8qo = _nhie9ed+_ovd79c9+_or8wm9l;
    r = sqrt(_h17j8qo); _h4a3yjt = r*r*r;
    _j4x5a5x = r*(r-z); _09vs9kk = _h17j8qo*(r-z)*(r-z);
    _u7fbmen = _h4a3yjt*(r-z);
    W = _lq5u6jz-r; _k0qtl3q = W*W;
    _ioxd9cs = W*r; _4n89qkn = _k0qtl3q*r;
    _s9toqw9 = W*_h4a3yjt; _s832ejx = _k0qtl3q*_h17j8qo;
    C = (r*_u1g463s -z)/_ioxd9cs; S = (r*_qnsn2s8 -y)/_ioxd9cs;
    _lz15bcw = S*S; _mv9e40k = y*y*y;
    _8drc11e = C*C;
    _gv44ftj = (1.0-_gdc5tao); _xlul6lg = (2.0*_gdc5tao+1.0);
    _dygflij= _u1g463s*_u1g463s;
    _ui3do70 = (_nwywrra/r/(r-_lq5u6jz)-y/r/(r-z))/4.0/3.14159265358979323846264338327950288;
    _02g3yu6 = (x/r/(r-z)-_u1g463s*x/r/(r-_lq5u6jz))/4.0/3.14159265358979323846264338327950288;
    _vk0ucc1 = (_qnsn2s8*x/r/(r-_lq5u6jz))/4.0/3.14159265358979323846264338327950288;
    _zhwjr64 = _wk8r2yk*(_ui3do70) + _wk8r2yk/8.0/3.14159265358979323846264338327950288/_gv44ftj*(_nwywrra/_ioxd9cs+_nwywrra*_nhie9ed/_s832ejx-_nwywrra*_nhie9ed/_s9toqw9+y/_j4x5a5x-_nhie9ed*y/_09vs9kk-_nhie9ed*y/_u7fbmen) - _ccprhvl*x/8.0/3.14159265358979323846264338327950288/_gv44ftj*((_xlul6lg/_ioxd9cs+_nhie9ed/_s832ejx-_nhie9ed/_s9toqw9)*_u1g463s+_xlul6lg/_j4x5a5x-_nhie9ed/_09vs9kk-_nhie9ed/_u7fbmen) + _lqjidkn*x*_qnsn2s8/8.0/3.14159265358979323846264338327950288/_gv44ftj*(_xlul6lg/_ioxd9cs+_nhie9ed/_s832ejx-_nhie9ed/_s9toqw9);
    _t4ufel6 = _ccprhvl*(_02g3yu6) + _wk8r2yk/8.0/3.14159265358979323846264338327950288/_gv44ftj*((1.0/_ioxd9cs+_lz15bcw-_ovd79c9/_s9toqw9)*_nwywrra+_xlul6lg*y/_j4x5a5x-_mv9e40k/_09vs9kk-_mv9e40k/_u7fbmen-2.0*_gdc5tao*_u1g463s*S) - _ccprhvl*x/8.0/3.14159265358979323846264338327950288/_gv44ftj*(1.0/_j4x5a5x-_ovd79c9/_09vs9kk-_ovd79c9/_u7fbmen+(1.0/_ioxd9cs+_lz15bcw-_ovd79c9/_s9toqw9)*_u1g463s) + _lqjidkn*x*_qnsn2s8/8.0/3.14159265358979323846264338327950288/_gv44ftj*(1.0/_ioxd9cs+_lz15bcw-_ovd79c9/_s9toqw9);
    _5owf759 = _lqjidkn*(_vk0ucc1) + _wk8r2yk/8.0/3.14159265358979323846264338327950288/_gv44ftj*(_nwywrra/W/r+_nwywrra*_8drc11e-_nwywrra*_or8wm9l/_s9toqw9+y*z/_h4a3yjt+2.0*_gdc5tao*_qnsn2s8*C) - _ccprhvl*x/8.0/3.14159265358979323846264338327950288/_gv44ftj*((1.0/_ioxd9cs+_8drc11e-_or8wm9l/_s9toqw9)*_u1g463s+z/_h4a3yjt) + _lqjidkn*x*_qnsn2s8/8.0/3.14159265358979323846264338327950288/_gv44ftj*(1.0/_ioxd9cs+_8drc11e-_or8wm9l/_s9toqw9);
    _dpja1wj = _wk8r2yk*(_02g3yu6)/2.0+_ccprhvl*(_ui3do70)/2.0 - _wk8r2yk/8.0/3.14159265358979323846264338327950288/_gv44ftj*(x*_ovd79c9/_09vs9kk-_gdc5tao*x/_j4x5a5x+x*_ovd79c9/_u7fbmen-_gdc5tao*x*_u1g463s/_ioxd9cs+_nwywrra*x*S/_ioxd9cs+_nwywrra*x*y/_s9toqw9) + _ccprhvl/8.0/3.14159265358979323846264338327950288/_gv44ftj*(_nhie9ed*y/_09vs9kk-_gdc5tao*y/_j4x5a5x+_nhie9ed*y/_u7fbmen+_gdc5tao*_u1g463s*S+_nhie9ed*y*_u1g463s/_s9toqw9+_nhie9ed*_u1g463s*S/_ioxd9cs) - _lqjidkn*_qnsn2s8/8.0/3.14159265358979323846264338327950288/_gv44ftj*(_gdc5tao*S+_nhie9ed*S/_ioxd9cs+_nhie9ed*y/_s9toqw9);
    _0tpac56 = _wk8r2yk*(_vk0ucc1)/2.0+_lqjidkn*(_ui3do70)/2.0 - _wk8r2yk/8.0/3.14159265358979323846264338327950288/_gv44ftj*(-x*y/_h4a3yjt+_gdc5tao*x*_qnsn2s8/_ioxd9cs+_nwywrra*x*C/_ioxd9cs+_nwywrra*x*z/_s9toqw9) + _ccprhvl/8.0/3.14159265358979323846264338327950288/_gv44ftj*(-_nhie9ed/_h4a3yjt+_gdc5tao/r+_gdc5tao*_u1g463s*C+_nhie9ed*z*_u1g463s/_s9toqw9+_nhie9ed*_u1g463s*C/_ioxd9cs) - _lqjidkn*_qnsn2s8/8.0/3.14159265358979323846264338327950288/_gv44ftj*(_gdc5tao*C+_nhie9ed*C/_ioxd9cs+_nhie9ed*z/_s9toqw9);
    _246hid2 = _ccprhvl*(_vk0ucc1)/2.0+_lqjidkn*(_02g3yu6)/2.0 + _wk8r2yk/8.0/3.14159265358979323846264338327950288/_gv44ftj*(_ovd79c9/_h4a3yjt-_gdc5tao/r-_gdc5tao*_u1g463s*C+_gdc5tao*_qnsn2s8*S+_nwywrra*_qnsn2s8*_u1g463s/_k0qtl3q-_nwywrra*(y*_u1g463s+z*_qnsn2s8)/_4n89qkn+_nwywrra*y*z/_s832ejx-_nwywrra*y*z/_s9toqw9) - _ccprhvl*x/8.0/3.14159265358979323846264338327950288/_gv44ftj*(y/_h4a3yjt+_qnsn2s8*_dygflij/_k0qtl3q-_u1g463s*(y*_u1g463s+z*_qnsn2s8)/_4n89qkn+y*z*_u1g463s/_s832ejx-y*z*_u1g463s/_s9toqw9) - _lqjidkn*x*_qnsn2s8/8.0/3.14159265358979323846264338327950288/_gv44ftj*(y*z/_s9toqw9-_qnsn2s8*_u1g463s/_k0qtl3q+(y*_u1g463s+z*_qnsn2s8)/_4n89qkn-y*z/_s832ejx);
    e[0] = _zhwjr64; e[1] = _dpja1wj; e[2] = _0tpac56;
    e[3] = _t4ufel6; e[4] = _246hid2; e[5] = _5owf759;
    return;
}
void _t3w5uzb(double _4gc4ryg, double _ovd79c9, double _mv9e40k, double _52xg8p6, double _jjq9k6r, double _e4q6qdy, double _4lmct64, double _gdc5tao, double a, double _82n2atq[6])
{
    double _nn16uxf, _i18k0ax, _ef0votq, _zuj2xor, _yxyur4s;
    double _bb4m0a6, _pmlvak9, _ic1l7dy, _3fpoyip, _k0qtl3q;
    double _5td44k1, _j6fp8hq, _mabsi9n, _3zezqac, _j4t0sek;
    double _kbp2qd6, _vasifpw, _hylurtp, _va1nkbt, _m3o55ix;
    double _fshb0oq, _r59c81d, _o6ouz1e, _gntb34k, _ej9c0uu;
    double _5wryesu, _4fcl5ef;
    double _rdu4xsz, _nsrtwlx, _0j36tcf, _6c2cu7o, _kjnawng;
    double _88d0tgs, _dsc4e5h, _nz3jbck, _jnvc7hr, _iwg0odw;
    double _a41izxt, _7yvdqbs, _cv3tlcx, _9rdpxdf, _w9pfb3h;
    double _xw1tbrw, _8zgvwwu, _b74x2ek, _vork1dj, _sxlmt2p;
    _nn16uxf = sin(_52xg8p6); _i18k0ax = cos(_52xg8p6); _ef0votq = _i18k0ax/_nn16uxf;
    _zuj2xor = _mv9e40k+a+a; _yxyur4s = _4gc4ryg*_i18k0ax+_zuj2xor*_nn16uxf; _bb4m0a6 = -_4gc4ryg*_nn16uxf+_zuj2xor*_i18k0ax;
    _pmlvak9 = _4gc4ryg*_4gc4ryg + _ovd79c9*_ovd79c9 + _zuj2xor*_zuj2xor; _ic1l7dy = sqrt(_pmlvak9); _hylurtp = 1.0-_gdc5tao-_gdc5tao;
    _3fpoyip = _ic1l7dy*_i18k0ax+_zuj2xor; _k0qtl3q = _i18k0ax +a/_ic1l7dy; _5td44k1 = _i18k0ax +_zuj2xor/_ic1l7dy;
    _j6fp8hq = _gdc5tao +a/_ic1l7dy; _mabsi9n = _gdc5tao+_gdc5tao +a/_ic1l7dy; _3zezqac = _ic1l7dy +_zuj2xor;
    _j4t0sek = _ic1l7dy +_bb4m0a6; _kbp2qd6 = _mv9e40k +a; _vasifpw = 1.0 +a/_ic1l7dy/_i18k0ax;
    _va1nkbt = _yxyur4s/_ic1l7dy/(_ic1l7dy+_bb4m0a6) -_4gc4ryg/_ic1l7dy/(_ic1l7dy+_zuj2xor);
    _m3o55ix = _ovd79c9/_ic1l7dy/(_ic1l7dy+_zuj2xor) -_i18k0ax*_ovd79c9/_ic1l7dy/(_ic1l7dy+_bb4m0a6);
    _fshb0oq = -_nn16uxf*_ovd79c9/_ic1l7dy/(_ic1l7dy+_bb4m0a6);
    _rdu4xsz = _ic1l7dy*_ic1l7dy*_ic1l7dy; _nsrtwlx = _ef0votq*_ef0votq; _0j36tcf = _3zezqac*_3zezqac;
    _6c2cu7o = _4gc4ryg*_4gc4ryg; _kjnawng = _j4t0sek*_j4t0sek; _88d0tgs = _ovd79c9*_ovd79c9;
    _dsc4e5h = _pmlvak9*_pmlvak9; _nz3jbck = _4gc4ryg*_4gc4ryg*_4gc4ryg; _jnvc7hr = _ic1l7dy*_ic1l7dy*_ic1l7dy*_ic1l7dy*_ic1l7dy;
    _iwg0odw = _ovd79c9*_ovd79c9*_ovd79c9; _a41izxt = _i18k0ax*_i18k0ax; _9rdpxdf = (_4gc4ryg/_ic1l7dy-_nn16uxf);
    _7yvdqbs = (2.0-_gdc5tao-_gdc5tao); _cv3tlcx = (a/_pmlvak9+1.0/_3zezqac); _w9pfb3h = (_ic1l7dy*_nn16uxf-_4gc4ryg);
    _xw1tbrw = (-2.0+_gdc5tao+_gdc5tao); _8zgvwwu = (1.0-_mabsi9n); _b74x2ek = (_mabsi9n-1.0);
    _vork1dj = (_zuj2xor/_ic1l7dy+1.0); _sxlmt2p = (1.0-_gdc5tao);
    _r59c81d = _jjq9k6r*(0.25*(_xw1tbrw*_hylurtp*_m3o55ix*_nsrtwlx-_hylurtp*_ovd79c9/_0j36tcf*(_8zgvwwu*_ef0votq-_4gc4ryg/_3zezqac*_j6fp8hq)/_ic1l7dy*_4gc4ryg+_hylurtp*_ovd79c9/_3zezqac*(a/_rdu4xsz*_4gc4ryg*_ef0votq-1.0/_3zezqac*_j6fp8hq+_6c2cu7o/_0j36tcf*_j6fp8hq/_ic1l7dy+_6c2cu7o/_3zezqac*a/_rdu4xsz)-_hylurtp*_ovd79c9*_i18k0ax*_ef0votq/_kjnawng*_k0qtl3q*_9rdpxdf-_hylurtp*_ovd79c9*_i18k0ax*_ef0votq/_j4t0sek*a/_rdu4xsz*_4gc4ryg-3.0*a*_ovd79c9*_kbp2qd6*_ef0votq/_jnvc7hr*_4gc4ryg-_ovd79c9*_kbp2qd6/_rdu4xsz/_3zezqac*(-_hylurtp*_ef0votq+_4gc4ryg/_3zezqac*_mabsi9n+a*_4gc4ryg/_pmlvak9)*_4gc4ryg-_ovd79c9*_kbp2qd6/_pmlvak9/_0j36tcf*(-_hylurtp*_ef0votq+_4gc4ryg/_3zezqac*_mabsi9n+a*_4gc4ryg/_pmlvak9)*_4gc4ryg+_ovd79c9*_kbp2qd6/_ic1l7dy/_3zezqac*(1.0/_3zezqac*_mabsi9n-_6c2cu7o/_0j36tcf*_mabsi9n/_ic1l7dy-_6c2cu7o/_3zezqac*a/_rdu4xsz+a/_pmlvak9-2.0*a*_6c2cu7o/_dsc4e5h)-_ovd79c9*_kbp2qd6/_rdu4xsz/_j4t0sek*(_i18k0ax/_j4t0sek*(_3fpoyip*(_hylurtp*_i18k0ax-a/_ic1l7dy)*_ef0votq+_7yvdqbs*_w9pfb3h*_i18k0ax)-a*_zuj2xor*_i18k0ax*_ef0votq/_pmlvak9)*_4gc4ryg-_ovd79c9*_kbp2qd6/_ic1l7dy/_kjnawng*(_i18k0ax/_j4t0sek*(_3fpoyip*(_hylurtp*_i18k0ax-a/_ic1l7dy)*_ef0votq+_7yvdqbs*_w9pfb3h*_i18k0ax)-a*_zuj2xor*_i18k0ax*_ef0votq/_pmlvak9)*_9rdpxdf+_ovd79c9*_kbp2qd6/_ic1l7dy/_j4t0sek*(-_i18k0ax/_kjnawng*(_3fpoyip*(_hylurtp*_i18k0ax-a/_ic1l7dy)*_ef0votq+_7yvdqbs*_w9pfb3h*_i18k0ax)*_9rdpxdf+_i18k0ax/_j4t0sek*(1.0/_ic1l7dy*_i18k0ax*_4gc4ryg*(_hylurtp*_i18k0ax-a/_ic1l7dy)*_ef0votq+_3fpoyip*a/_rdu4xsz*_4gc4ryg*_ef0votq+_7yvdqbs*(1.0/_ic1l7dy*_nn16uxf*_4gc4ryg-1.0)*_i18k0ax)+2.0*a*_zuj2xor*_i18k0ax*_ef0votq/_dsc4e5h*_4gc4ryg))/3.14159265358979323846264338327950288/_sxlmt2p)
        + _e4q6qdy*(0.25*(_hylurtp*((_7yvdqbs*_nsrtwlx+_gdc5tao)/_ic1l7dy*_4gc4ryg/_3zezqac-(_7yvdqbs*_nsrtwlx+1.0)*_i18k0ax*_9rdpxdf/_j4t0sek)-_hylurtp/_0j36tcf*(-_hylurtp*_4gc4ryg*_ef0votq+_gdc5tao*_zuj2xor-a+a*_4gc4ryg*_ef0votq/_ic1l7dy+_6c2cu7o/_3zezqac*_j6fp8hq)/_ic1l7dy*_4gc4ryg+_hylurtp/_3zezqac*(-_hylurtp*_ef0votq+a*_ef0votq/_ic1l7dy-a*_6c2cu7o*_ef0votq/_rdu4xsz+2.0*_4gc4ryg/_3zezqac*_j6fp8hq-_nz3jbck/_0j36tcf*_j6fp8hq/_ic1l7dy-_nz3jbck/_3zezqac*a/_rdu4xsz)+_hylurtp*_ef0votq/_kjnawng*(_yxyur4s*_i18k0ax-a*_w9pfb3h/_ic1l7dy/_i18k0ax)*_9rdpxdf-_hylurtp*_ef0votq/_j4t0sek*(_a41izxt-a*(1.0/_ic1l7dy*_nn16uxf*_4gc4ryg-1.0)/_ic1l7dy/_i18k0ax+a*_w9pfb3h/_rdu4xsz/_i18k0ax*_4gc4ryg)-a*_kbp2qd6*_ef0votq/_rdu4xsz+3.0*a*_6c2cu7o*_kbp2qd6*_ef0votq/_jnvc7hr-_kbp2qd6/_0j36tcf*(2.0*_gdc5tao+1.0/_ic1l7dy*(_hylurtp*_4gc4ryg*_ef0votq+a)-_6c2cu7o/_ic1l7dy/_3zezqac*_mabsi9n-a*_6c2cu7o/_rdu4xsz)/_ic1l7dy*_4gc4ryg+_kbp2qd6/_3zezqac*(-1.0/_rdu4xsz*(_hylurtp*_4gc4ryg*_ef0votq+a)*_4gc4ryg+1.0/_ic1l7dy*_hylurtp*_ef0votq-2.0*_4gc4ryg/_ic1l7dy/_3zezqac*_mabsi9n+_nz3jbck/_rdu4xsz/_3zezqac*_mabsi9n+_nz3jbck/_pmlvak9/_0j36tcf*_mabsi9n+_nz3jbck/_dsc4e5h/_3zezqac*a-2.0*a/_rdu4xsz*_4gc4ryg+3.0*a*_nz3jbck/_jnvc7hr)-_kbp2qd6*_ef0votq/_kjnawng*(-_i18k0ax*_nn16uxf+a*_4gc4ryg*_zuj2xor/_rdu4xsz/_i18k0ax+_w9pfb3h/_ic1l7dy*(_7yvdqbs*_i18k0ax-_3fpoyip/_j4t0sek*_vasifpw))*_9rdpxdf+_kbp2qd6*_ef0votq/_j4t0sek*(a*_zuj2xor/_rdu4xsz/_i18k0ax-3.0*a*_6c2cu7o*_zuj2xor/_jnvc7hr/_i18k0ax+(1.0/_ic1l7dy*_nn16uxf*_4gc4ryg-1.0)/_ic1l7dy*(_7yvdqbs*_i18k0ax-_3fpoyip/_j4t0sek*_vasifpw)-_w9pfb3h/_rdu4xsz*(_7yvdqbs*_i18k0ax-_3fpoyip/_j4t0sek*_vasifpw)*_4gc4ryg+_w9pfb3h/_ic1l7dy*(-1.0/_ic1l7dy*_i18k0ax*_4gc4ryg/_j4t0sek*_vasifpw+_3fpoyip/_kjnawng*_vasifpw*_9rdpxdf+_3fpoyip/_j4t0sek*a/_rdu4xsz/_i18k0ax*_4gc4ryg)))/3.14159265358979323846264338327950288/_sxlmt2p)
        + _4lmct64*(0.25*(_hylurtp*(-_ovd79c9/_0j36tcf*(1.0+a/_ic1l7dy)/_ic1l7dy*_4gc4ryg-_ovd79c9/_3zezqac*a/_rdu4xsz*_4gc4ryg+_ovd79c9*_i18k0ax/_kjnawng*_k0qtl3q*_9rdpxdf+_ovd79c9*_i18k0ax/_j4t0sek*a/_rdu4xsz*_4gc4ryg)+_ovd79c9*_kbp2qd6/_rdu4xsz*_cv3tlcx*_4gc4ryg-_ovd79c9*_kbp2qd6/_ic1l7dy*(-2.0*a/_dsc4e5h*_4gc4ryg-1.0/_0j36tcf/_ic1l7dy*_4gc4ryg)-_ovd79c9*_kbp2qd6*_i18k0ax/_rdu4xsz/_j4t0sek*(_3fpoyip/_j4t0sek*_k0qtl3q+a*_zuj2xor/_pmlvak9)*_4gc4ryg-_ovd79c9*_kbp2qd6*_i18k0ax/_ic1l7dy/_kjnawng*(_3fpoyip/_j4t0sek*_k0qtl3q+a*_zuj2xor/_pmlvak9)*_9rdpxdf+_ovd79c9*_kbp2qd6*_i18k0ax/_ic1l7dy/_j4t0sek*(1.0/_ic1l7dy*_i18k0ax*_4gc4ryg/_j4t0sek*_k0qtl3q-_3fpoyip/_kjnawng*_k0qtl3q*_9rdpxdf-_3fpoyip/_j4t0sek*a/_rdu4xsz*_4gc4ryg-2.0*a*_zuj2xor/_dsc4e5h*_4gc4ryg))/3.14159265358979323846264338327950288/_sxlmt2p);
    _ej9c0uu = _jjq9k6r*(0.25*(_hylurtp*((_7yvdqbs*_nsrtwlx-_gdc5tao)/_ic1l7dy*_ovd79c9/_3zezqac-(_7yvdqbs*_nsrtwlx+1.0-2.0*_gdc5tao)*_i18k0ax/_ic1l7dy*_ovd79c9/_j4t0sek)+_hylurtp/_0j36tcf*(_4gc4ryg*_ef0votq*_8zgvwwu+_gdc5tao*_zuj2xor-a+_88d0tgs/_3zezqac*_j6fp8hq)/_ic1l7dy*_ovd79c9-_hylurtp/_3zezqac*(a*_4gc4ryg*_ef0votq/_rdu4xsz*_ovd79c9+2.0*_ovd79c9/_3zezqac*_j6fp8hq-_iwg0odw/_0j36tcf*_j6fp8hq/_ic1l7dy-_iwg0odw/_3zezqac*a/_rdu4xsz)+_hylurtp*_yxyur4s*_ef0votq/_kjnawng*_k0qtl3q/_ic1l7dy*_ovd79c9+_hylurtp*_yxyur4s*_ef0votq/_j4t0sek*a/_rdu4xsz*_ovd79c9+3.0*a*_ovd79c9*_kbp2qd6*_ef0votq/_jnvc7hr*_4gc4ryg-_kbp2qd6/_0j36tcf*(-2.0*_gdc5tao+1.0/_ic1l7dy*(_hylurtp*_4gc4ryg*_ef0votq-a)+_88d0tgs/_ic1l7dy/_3zezqac*_mabsi9n+a*_88d0tgs/_rdu4xsz)/_ic1l7dy*_ovd79c9+_kbp2qd6/_3zezqac*(-1.0/_rdu4xsz*(_hylurtp*_4gc4ryg*_ef0votq-a)*_ovd79c9+2.0*_ovd79c9/_ic1l7dy/_3zezqac*_mabsi9n-_iwg0odw/_rdu4xsz/_3zezqac*_mabsi9n-_iwg0odw/_pmlvak9/_0j36tcf*_mabsi9n-_iwg0odw/_dsc4e5h/_3zezqac*a+2.0*a/_rdu4xsz*_ovd79c9-3.0*a*_iwg0odw/_jnvc7hr)-_kbp2qd6/_kjnawng*(_a41izxt-1.0/_ic1l7dy*(_hylurtp*_yxyur4s*_ef0votq+a*_i18k0ax)+a*_zuj2xor*_yxyur4s*_ef0votq/_rdu4xsz-1.0/_ic1l7dy/_j4t0sek*(_88d0tgs*_a41izxt-a*_yxyur4s*_ef0votq/_ic1l7dy*_3fpoyip))/_ic1l7dy*_ovd79c9+_kbp2qd6/_j4t0sek*(1.0/_rdu4xsz*(_hylurtp*_yxyur4s*_ef0votq+a*_i18k0ax)*_ovd79c9-3.0*a*_zuj2xor*_yxyur4s*_ef0votq/_jnvc7hr*_ovd79c9+1.0/_rdu4xsz/_j4t0sek*(_88d0tgs*_a41izxt-a*_yxyur4s*_ef0votq/_ic1l7dy*_3fpoyip)*_ovd79c9+1.0/_pmlvak9/_kjnawng*(_88d0tgs*_a41izxt-a*_yxyur4s*_ef0votq/_ic1l7dy*_3fpoyip)*_ovd79c9-1.0/_ic1l7dy/_j4t0sek*(2.0*_ovd79c9*_a41izxt+a*_yxyur4s*_ef0votq/_rdu4xsz*_3fpoyip*_ovd79c9-a*_yxyur4s*_ef0votq/_pmlvak9*_i18k0ax*_ovd79c9)))/3.14159265358979323846264338327950288/_sxlmt2p)
        + _e4q6qdy*(0.25*(_7yvdqbs*_hylurtp*_va1nkbt*_nsrtwlx+_hylurtp/_3zezqac*(_b74x2ek*_ef0votq+_4gc4ryg/_3zezqac*_j6fp8hq)-_hylurtp*_88d0tgs/_0j36tcf*(_b74x2ek*_ef0votq+_4gc4ryg/_3zezqac*_j6fp8hq)/_ic1l7dy+_hylurtp*_ovd79c9/_3zezqac*(-a/_rdu4xsz*_ovd79c9*_ef0votq-_4gc4ryg/_0j36tcf*_j6fp8hq/_ic1l7dy*_ovd79c9-_ovd79c9/_3zezqac*a/_rdu4xsz*_4gc4ryg)-_hylurtp*_ef0votq/_j4t0sek*_vasifpw+_hylurtp*_88d0tgs*_ef0votq/_kjnawng*_vasifpw/_ic1l7dy+_hylurtp*_88d0tgs*_ef0votq/_j4t0sek*a/_rdu4xsz/_i18k0ax-a*_kbp2qd6*_ef0votq/_rdu4xsz+3.0*a*_88d0tgs*_kbp2qd6*_ef0votq/_jnvc7hr+_kbp2qd6/_ic1l7dy/_3zezqac*(_hylurtp*_ef0votq-2.0*_gdc5tao*_4gc4ryg/_3zezqac-a*_4gc4ryg/_ic1l7dy*(1.0/_ic1l7dy+1.0/_3zezqac))-_88d0tgs*_kbp2qd6/_rdu4xsz/_3zezqac*(_hylurtp*_ef0votq-2.0*_gdc5tao*_4gc4ryg/_3zezqac-a*_4gc4ryg/_ic1l7dy*(1.0/_ic1l7dy+1.0/_3zezqac))-_88d0tgs*_kbp2qd6/_pmlvak9/_0j36tcf*(_hylurtp*_ef0votq-2.0*_gdc5tao*_4gc4ryg/_3zezqac-a*_4gc4ryg/_ic1l7dy*(1.0/_ic1l7dy+1.0/_3zezqac))+_ovd79c9*_kbp2qd6/_ic1l7dy/_3zezqac*(2.0*_gdc5tao*_4gc4ryg/_0j36tcf/_ic1l7dy*_ovd79c9+a*_4gc4ryg/_rdu4xsz*(1.0/_ic1l7dy+1.0/_3zezqac)*_ovd79c9-a*_4gc4ryg/_ic1l7dy*(-1.0/_rdu4xsz*_ovd79c9-1.0/_0j36tcf/_ic1l7dy*_ovd79c9))+_kbp2qd6*_ef0votq/_ic1l7dy/_j4t0sek*(_xw1tbrw*_i18k0ax+_3fpoyip/_j4t0sek*_vasifpw+a*_zuj2xor/_pmlvak9/_i18k0ax)-_88d0tgs*_kbp2qd6*_ef0votq/_rdu4xsz/_j4t0sek*(_xw1tbrw*_i18k0ax+_3fpoyip/_j4t0sek*_vasifpw+a*_zuj2xor/_pmlvak9/_i18k0ax)-_88d0tgs*_kbp2qd6*_ef0votq/_pmlvak9/_kjnawng*(_xw1tbrw*_i18k0ax+_3fpoyip/_j4t0sek*_vasifpw+a*_zuj2xor/_pmlvak9/_i18k0ax)+_ovd79c9*_kbp2qd6*_ef0votq/_ic1l7dy/_j4t0sek*(1.0/_ic1l7dy*_i18k0ax*_ovd79c9/_j4t0sek*_vasifpw-_3fpoyip/_kjnawng*_vasifpw/_ic1l7dy*_ovd79c9-_3fpoyip/_j4t0sek*a/_rdu4xsz/_i18k0ax*_ovd79c9-2.0*a*_zuj2xor/_dsc4e5h/_i18k0ax*_ovd79c9))/3.14159265358979323846264338327950288/_sxlmt2p)
        + _4lmct64*(0.25*(_hylurtp*(-_nn16uxf/_ic1l7dy*_ovd79c9/_j4t0sek+_ovd79c9/_0j36tcf*(1.0+a/_ic1l7dy)/_ic1l7dy*_4gc4ryg+_ovd79c9/_3zezqac*a/_rdu4xsz*_4gc4ryg-_yxyur4s/_kjnawng*_k0qtl3q/_ic1l7dy*_ovd79c9-_yxyur4s/_j4t0sek*a/_rdu4xsz*_ovd79c9)-_ovd79c9*_kbp2qd6/_rdu4xsz*_cv3tlcx*_4gc4ryg+_4gc4ryg*_kbp2qd6/_ic1l7dy*(-2.0*a/_dsc4e5h*_ovd79c9-1.0/_0j36tcf/_ic1l7dy*_ovd79c9)+_kbp2qd6/_kjnawng*(_nn16uxf*(_i18k0ax-a/_ic1l7dy)+_yxyur4s/_ic1l7dy*(1.0+a*_zuj2xor/_pmlvak9)-1.0/_ic1l7dy/_j4t0sek*(_88d0tgs*_i18k0ax*_nn16uxf-a*_yxyur4s/_ic1l7dy*_3fpoyip))/_ic1l7dy*_ovd79c9-_kbp2qd6/_j4t0sek*(_nn16uxf*a/_rdu4xsz*_ovd79c9-_yxyur4s/_rdu4xsz*(1.0+a*_zuj2xor/_pmlvak9)*_ovd79c9-2.0*_yxyur4s/_jnvc7hr*a*_zuj2xor*_ovd79c9+1.0/_rdu4xsz/_j4t0sek*(_88d0tgs*_i18k0ax*_nn16uxf-a*_yxyur4s/_ic1l7dy*_3fpoyip)*_ovd79c9+1.0/_pmlvak9/_kjnawng*(_88d0tgs*_i18k0ax*_nn16uxf-a*_yxyur4s/_ic1l7dy*_3fpoyip)*_ovd79c9-1/_ic1l7dy/_j4t0sek*(2*_ovd79c9*_i18k0ax*_nn16uxf+a*_yxyur4s/_rdu4xsz*_3fpoyip*_ovd79c9-a*_yxyur4s/_pmlvak9*_i18k0ax*_ovd79c9)))/3.14159265358979323846264338327950288/_sxlmt2p);
    _4fcl5ef = _jjq9k6r*(0.25*(_7yvdqbs*(_hylurtp*_fshb0oq*_ef0votq-_ovd79c9/_0j36tcf*_mabsi9n*_vork1dj-0.5*_ovd79c9/_3zezqac*a/_rdu4xsz*2.0*_zuj2xor+_ovd79c9*_i18k0ax/_kjnawng*_k0qtl3q*_5td44k1+0.5*_ovd79c9*_i18k0ax/_j4t0sek*a/_rdu4xsz*2.0*_zuj2xor)+_ovd79c9/_ic1l7dy*(2.0*_gdc5tao/_3zezqac+a/_pmlvak9)-0.5*_ovd79c9*_kbp2qd6/_rdu4xsz*(2.0*_gdc5tao/_3zezqac+a/_pmlvak9)*2.0*_zuj2xor+_ovd79c9*_kbp2qd6/_ic1l7dy*(-2.0*_gdc5tao/_0j36tcf*_vork1dj-a/_dsc4e5h*2.0*_zuj2xor)+_ovd79c9*_i18k0ax/_ic1l7dy/_j4t0sek*(1.0-2.0*_gdc5tao-_3fpoyip/_j4t0sek*_k0qtl3q-a*_zuj2xor/_pmlvak9)-0.5*_ovd79c9*_kbp2qd6*_i18k0ax/_rdu4xsz/_j4t0sek*(1.0-2.0*_gdc5tao-_3fpoyip/_j4t0sek*_k0qtl3q-a*_zuj2xor/_pmlvak9)*2.0*_zuj2xor-_ovd79c9*_kbp2qd6*_i18k0ax/_ic1l7dy/_kjnawng*(1.0-2.0*_gdc5tao-_3fpoyip/_j4t0sek*_k0qtl3q-a*_zuj2xor/_pmlvak9)*_5td44k1+_ovd79c9*_kbp2qd6*_i18k0ax/_ic1l7dy/_j4t0sek*(-(_i18k0ax*_zuj2xor/_ic1l7dy+1.0)/_j4t0sek*_k0qtl3q+_3fpoyip/_kjnawng*_k0qtl3q*_5td44k1+0.5*_3fpoyip/_j4t0sek*a/_rdu4xsz*2.0*_zuj2xor-a/_pmlvak9+a*_zuj2xor/_dsc4e5h*2.0*_zuj2xor))/3.14159265358979323846264338327950288/_sxlmt2p)
        + _e4q6qdy*(0.25*(_xw1tbrw*_hylurtp*_ef0votq*(_vork1dj/_3zezqac-_i18k0ax*_5td44k1/_j4t0sek)+_7yvdqbs*_4gc4ryg/_0j36tcf*_mabsi9n*_vork1dj+0.5*_7yvdqbs*_4gc4ryg/_3zezqac*a/_rdu4xsz*2.0*_zuj2xor+_7yvdqbs*_nn16uxf/_j4t0sek*_k0qtl3q-_7yvdqbs*_yxyur4s/_kjnawng*_k0qtl3q*_5td44k1-0.5*_7yvdqbs*_yxyur4s/_j4t0sek*a/_rdu4xsz*2.0*_zuj2xor+1.0/_ic1l7dy*(_hylurtp*_ef0votq-2.0*_gdc5tao*_4gc4ryg/_3zezqac-a*_4gc4ryg/_pmlvak9)-0.5*_kbp2qd6/_rdu4xsz*(_hylurtp*_ef0votq-2.0*_gdc5tao*_4gc4ryg/_3zezqac-a*_4gc4ryg/_pmlvak9)*2.0*_zuj2xor+_kbp2qd6/_ic1l7dy*(2.0*_gdc5tao*_4gc4ryg/_0j36tcf*_vork1dj+a*_4gc4ryg/_dsc4e5h*2.0*_zuj2xor)-1.0/_j4t0sek*(_i18k0ax*_nn16uxf+_3fpoyip*_ef0votq/_ic1l7dy*(_7yvdqbs*_i18k0ax-_3fpoyip/_j4t0sek)+a/_ic1l7dy*(_nn16uxf-_zuj2xor*_yxyur4s/_pmlvak9-_yxyur4s*_3fpoyip/_ic1l7dy/_j4t0sek))+_kbp2qd6/_kjnawng*(_i18k0ax*_nn16uxf+_3fpoyip*_ef0votq/_ic1l7dy*(_7yvdqbs*_i18k0ax-_3fpoyip/_j4t0sek)+a/_ic1l7dy*(_nn16uxf-_zuj2xor*_yxyur4s/_pmlvak9-_yxyur4s*_3fpoyip/_ic1l7dy/_j4t0sek))*_5td44k1-_kbp2qd6/_j4t0sek*((_i18k0ax*_zuj2xor/_ic1l7dy+1.0)*_ef0votq/_ic1l7dy*(_7yvdqbs*_i18k0ax-_3fpoyip/_j4t0sek)-0.5*_3fpoyip*_ef0votq/_rdu4xsz*(_7yvdqbs*_i18k0ax-_3fpoyip/_j4t0sek)*2.0*_zuj2xor+_3fpoyip*_ef0votq/_ic1l7dy*(-(_i18k0ax*_zuj2xor/_ic1l7dy+1.0)/_j4t0sek+_3fpoyip/_kjnawng*_5td44k1)-0.5*a/_rdu4xsz*(_nn16uxf-_zuj2xor*_yxyur4s/_pmlvak9-_yxyur4s*_3fpoyip/_ic1l7dy/_j4t0sek)*2.0*_zuj2xor+a/_ic1l7dy*(-_yxyur4s/_pmlvak9-_zuj2xor*_nn16uxf/_pmlvak9+_zuj2xor*_yxyur4s/_dsc4e5h*2.0*_zuj2xor-_nn16uxf*_3fpoyip/_ic1l7dy/_j4t0sek-_yxyur4s*(_i18k0ax*_zuj2xor/_ic1l7dy+1.0)/_ic1l7dy/_j4t0sek+0.5*_yxyur4s*_3fpoyip/_rdu4xsz/_j4t0sek*2.0*_zuj2xor+_yxyur4s*_3fpoyip/_ic1l7dy/_kjnawng*_5td44k1)))/3.14159265358979323846264338327950288/_sxlmt2p)
        + _4lmct64*(0.25*(_7yvdqbs*_fshb0oq-_7yvdqbs*_ovd79c9*_nn16uxf/_kjnawng*_k0qtl3q*_5td44k1-0.5*_7yvdqbs*_ovd79c9*_nn16uxf/_j4t0sek*a/_rdu4xsz*2.0*_zuj2xor+_ovd79c9*_nn16uxf/_ic1l7dy/_j4t0sek*(1.0+_3fpoyip/_j4t0sek*_k0qtl3q+a*_zuj2xor/_pmlvak9)-0.5*_ovd79c9*_kbp2qd6*_nn16uxf/_rdu4xsz/_j4t0sek*(1.0+_3fpoyip/_j4t0sek*_k0qtl3q+a*_zuj2xor/_pmlvak9)*2.0*_zuj2xor-_ovd79c9*_kbp2qd6*_nn16uxf/_ic1l7dy/_kjnawng*(1.0+_3fpoyip/_j4t0sek*_k0qtl3q+a*_zuj2xor/_pmlvak9)*_5td44k1+_ovd79c9*_kbp2qd6*_nn16uxf/_ic1l7dy/_j4t0sek*((_i18k0ax*_zuj2xor/_ic1l7dy+1.0)/_j4t0sek*_k0qtl3q-_3fpoyip/_kjnawng*_k0qtl3q*_5td44k1-0.5*_3fpoyip/_j4t0sek*a/_rdu4xsz*2.0*_zuj2xor+a/_pmlvak9-a*_zuj2xor/_dsc4e5h*2.0*_zuj2xor))/3.14159265358979323846264338327950288/_sxlmt2p);
    _o6ouz1e = _jjq9k6r/2.0*(0.25*(_xw1tbrw*_hylurtp*_va1nkbt*_nsrtwlx+_hylurtp/_3zezqac*(_8zgvwwu*_ef0votq-_4gc4ryg/_3zezqac*_j6fp8hq)-_hylurtp*_88d0tgs/_0j36tcf*(_8zgvwwu*_ef0votq-_4gc4ryg/_3zezqac*_j6fp8hq)/_ic1l7dy+_hylurtp*_ovd79c9/_3zezqac*(a/_rdu4xsz*_ovd79c9*_ef0votq+_4gc4ryg/_0j36tcf*_j6fp8hq/_ic1l7dy*_ovd79c9+_ovd79c9/_3zezqac*a/_rdu4xsz*_4gc4ryg)+_hylurtp*_i18k0ax*_ef0votq/_j4t0sek*_k0qtl3q-_hylurtp*_88d0tgs*_i18k0ax*_ef0votq/_kjnawng*_k0qtl3q/_ic1l7dy-_hylurtp*_88d0tgs*_i18k0ax*_ef0votq/_j4t0sek*a/_rdu4xsz+a*_kbp2qd6*_ef0votq/_rdu4xsz-3.0*a*_88d0tgs*_kbp2qd6*_ef0votq/_jnvc7hr+_kbp2qd6/_ic1l7dy/_3zezqac*(-_hylurtp*_ef0votq+_4gc4ryg/_3zezqac*_mabsi9n+a*_4gc4ryg/_pmlvak9)-_88d0tgs*_kbp2qd6/_rdu4xsz/_3zezqac*(-_hylurtp*_ef0votq+_4gc4ryg/_3zezqac*_mabsi9n+a*_4gc4ryg/_pmlvak9)-_88d0tgs*_kbp2qd6/_pmlvak9/_0j36tcf*(-_hylurtp*_ef0votq+_4gc4ryg/_3zezqac*_mabsi9n+a*_4gc4ryg/_pmlvak9)+_ovd79c9*_kbp2qd6/_ic1l7dy/_3zezqac*(-_4gc4ryg/_0j36tcf*_mabsi9n/_ic1l7dy*_ovd79c9-_ovd79c9/_3zezqac*a/_rdu4xsz*_4gc4ryg-2.0*a*_4gc4ryg/_dsc4e5h*_ovd79c9)+_kbp2qd6/_ic1l7dy/_j4t0sek*(_i18k0ax/_j4t0sek*(_3fpoyip*(_hylurtp*_i18k0ax-a/_ic1l7dy)*_ef0votq+_7yvdqbs*_w9pfb3h*_i18k0ax)-a*_zuj2xor*_i18k0ax*_ef0votq/_pmlvak9)-_88d0tgs*_kbp2qd6/_rdu4xsz/_j4t0sek*(_i18k0ax/_j4t0sek*(_3fpoyip*(_hylurtp*_i18k0ax-a/_ic1l7dy)*_ef0votq+_7yvdqbs*_w9pfb3h*_i18k0ax)-a*_zuj2xor*_i18k0ax*_ef0votq/_pmlvak9)-_88d0tgs*_kbp2qd6/_pmlvak9/_kjnawng*(_i18k0ax/_j4t0sek*(_3fpoyip*(_hylurtp*_i18k0ax-a/_ic1l7dy)*_ef0votq+_7yvdqbs*_w9pfb3h*_i18k0ax)-a*_zuj2xor*_i18k0ax*_ef0votq/_pmlvak9)+_ovd79c9*_kbp2qd6/_ic1l7dy/_j4t0sek*(-_i18k0ax/_kjnawng*(_3fpoyip*(_hylurtp*_i18k0ax-a/_ic1l7dy)*_ef0votq+_7yvdqbs*_w9pfb3h*_i18k0ax)/_ic1l7dy*_ovd79c9+_i18k0ax/_j4t0sek*(1.0/_ic1l7dy*_i18k0ax*_ovd79c9*(_hylurtp*_i18k0ax-a/_ic1l7dy)*_ef0votq+_3fpoyip*a/_rdu4xsz*_ovd79c9*_ef0votq+_7yvdqbs/_ic1l7dy*_nn16uxf*_ovd79c9*_i18k0ax)+2.0*a*_zuj2xor*_i18k0ax*_ef0votq/_dsc4e5h*_ovd79c9))/3.14159265358979323846264338327950288/_sxlmt2p)
        + _e4q6qdy/2.0*(0.25*(_hylurtp*((_7yvdqbs*_nsrtwlx+_gdc5tao)/_ic1l7dy*_ovd79c9/_3zezqac-(_7yvdqbs*_nsrtwlx+1.0)*_i18k0ax/_ic1l7dy*_ovd79c9/_j4t0sek)-_hylurtp/_0j36tcf*(-_hylurtp*_4gc4ryg*_ef0votq+_gdc5tao*_zuj2xor-a+a*_4gc4ryg*_ef0votq/_ic1l7dy+_6c2cu7o/_3zezqac*_j6fp8hq)/_ic1l7dy*_ovd79c9+_hylurtp/_3zezqac*(-a*_4gc4ryg*_ef0votq/_rdu4xsz*_ovd79c9-_6c2cu7o/_0j36tcf*_j6fp8hq/_ic1l7dy*_ovd79c9-_6c2cu7o/_3zezqac*a/_rdu4xsz*_ovd79c9)+_hylurtp*_ef0votq/_kjnawng*(_yxyur4s*_i18k0ax-a*_w9pfb3h/_ic1l7dy/_i18k0ax)/_ic1l7dy*_ovd79c9-_hylurtp*_ef0votq/_j4t0sek*(-a/_pmlvak9*_nn16uxf*_ovd79c9/_i18k0ax+a*_w9pfb3h/_rdu4xsz/_i18k0ax*_ovd79c9)+3.0*a*_ovd79c9*_kbp2qd6*_ef0votq/_jnvc7hr*_4gc4ryg-_kbp2qd6/_0j36tcf*(2.0*_gdc5tao+1.0/_ic1l7dy*(_hylurtp*_4gc4ryg*_ef0votq+a)-_6c2cu7o/_ic1l7dy/_3zezqac*_mabsi9n-a*_6c2cu7o/_rdu4xsz)/_ic1l7dy*_ovd79c9+_kbp2qd6/_3zezqac*(-1.0/_rdu4xsz*(_hylurtp*_4gc4ryg*_ef0votq+a)*_ovd79c9+_6c2cu7o/_rdu4xsz/_3zezqac*_mabsi9n*_ovd79c9+_6c2cu7o/_pmlvak9/_0j36tcf*_mabsi9n*_ovd79c9+_6c2cu7o/_dsc4e5h/_3zezqac*a*_ovd79c9+3.0*a*_6c2cu7o/_jnvc7hr*_ovd79c9)-_kbp2qd6*_ef0votq/_kjnawng*(-_i18k0ax*_nn16uxf+a*_4gc4ryg*_zuj2xor/_rdu4xsz/_i18k0ax+_w9pfb3h/_ic1l7dy*(_7yvdqbs*_i18k0ax-_3fpoyip/_j4t0sek*_vasifpw))/_ic1l7dy*_ovd79c9+_kbp2qd6*_ef0votq/_j4t0sek*(-3.0*a*_4gc4ryg*_zuj2xor/_jnvc7hr/_i18k0ax*_ovd79c9+1.0/_pmlvak9*_nn16uxf*_ovd79c9*(_7yvdqbs*_i18k0ax-_3fpoyip/_j4t0sek*_vasifpw)-_w9pfb3h/_rdu4xsz*(_7yvdqbs*_i18k0ax-_3fpoyip/_j4t0sek*_vasifpw)*_ovd79c9+_w9pfb3h/_ic1l7dy*(-1.0/_ic1l7dy*_i18k0ax*_ovd79c9/_j4t0sek*_vasifpw+_3fpoyip/_kjnawng*_vasifpw/_ic1l7dy*_ovd79c9+_3fpoyip/_j4t0sek*a/_rdu4xsz/_i18k0ax*_ovd79c9)))/3.14159265358979323846264338327950288/_sxlmt2p)
        + _4lmct64/2.0*(0.25*(_hylurtp*(1.0/_3zezqac*(1.0+a/_ic1l7dy)-_88d0tgs/_0j36tcf*(1.0+a/_ic1l7dy)/_ic1l7dy-_88d0tgs/_3zezqac*a/_rdu4xsz-_i18k0ax/_j4t0sek*_k0qtl3q+_88d0tgs*_i18k0ax/_kjnawng*_k0qtl3q/_ic1l7dy+_88d0tgs*_i18k0ax/_j4t0sek*a/_rdu4xsz)-_kbp2qd6/_ic1l7dy*_cv3tlcx+_88d0tgs*_kbp2qd6/_rdu4xsz*_cv3tlcx-_ovd79c9*_kbp2qd6/_ic1l7dy*(-2.0*a/_dsc4e5h*_ovd79c9-1.0/_0j36tcf/_ic1l7dy*_ovd79c9)+_kbp2qd6*_i18k0ax/_ic1l7dy/_j4t0sek*(_3fpoyip/_j4t0sek*_k0qtl3q+a*_zuj2xor/_pmlvak9)-_88d0tgs*_kbp2qd6*_i18k0ax/_rdu4xsz/_j4t0sek*(_3fpoyip/_j4t0sek*_k0qtl3q+a*_zuj2xor/_pmlvak9)-_88d0tgs*_kbp2qd6*_i18k0ax/_pmlvak9/_kjnawng*(_3fpoyip/_j4t0sek*_k0qtl3q+a*_zuj2xor/_pmlvak9)+_ovd79c9*_kbp2qd6*_i18k0ax/_ic1l7dy/_j4t0sek*(1.0/_ic1l7dy*_i18k0ax*_ovd79c9/_j4t0sek*_k0qtl3q-_3fpoyip/_kjnawng*_k0qtl3q/_ic1l7dy*_ovd79c9-_3fpoyip/_j4t0sek*a/_rdu4xsz*_ovd79c9-2.0*a*_zuj2xor/_dsc4e5h*_ovd79c9))/3.14159265358979323846264338327950288/_sxlmt2p)
        + _jjq9k6r/2.0*(0.25*(_hylurtp*((_7yvdqbs*_nsrtwlx-_gdc5tao)/_ic1l7dy*_4gc4ryg/_3zezqac-(_7yvdqbs*_nsrtwlx+1.0-2.0*_gdc5tao)*_i18k0ax*_9rdpxdf/_j4t0sek)+_hylurtp/_0j36tcf*(_4gc4ryg*_ef0votq*_8zgvwwu+_gdc5tao*_zuj2xor-a+_88d0tgs/_3zezqac*_j6fp8hq)/_ic1l7dy*_4gc4ryg-_hylurtp/_3zezqac*(_8zgvwwu*_ef0votq+a*_6c2cu7o*_ef0votq/_rdu4xsz-_88d0tgs/_0j36tcf*_j6fp8hq/_ic1l7dy*_4gc4ryg-_88d0tgs/_3zezqac*a/_rdu4xsz*_4gc4ryg)-_hylurtp*_i18k0ax*_ef0votq/_j4t0sek*_k0qtl3q+_hylurtp*_yxyur4s*_ef0votq/_kjnawng*_k0qtl3q*_9rdpxdf+_hylurtp*_yxyur4s*_ef0votq/_j4t0sek*a/_rdu4xsz*_4gc4ryg-a*_kbp2qd6*_ef0votq/_rdu4xsz+3.0*a*_6c2cu7o*_kbp2qd6*_ef0votq/_jnvc7hr-_kbp2qd6/_0j36tcf*(-2.0*_gdc5tao+1.0/_ic1l7dy*(_hylurtp*_4gc4ryg*_ef0votq-a)+_88d0tgs/_ic1l7dy/_3zezqac*_mabsi9n+a*_88d0tgs/_rdu4xsz)/_ic1l7dy*_4gc4ryg+_kbp2qd6/_3zezqac*(-1.0/_rdu4xsz*(_hylurtp*_4gc4ryg*_ef0votq-a)*_4gc4ryg+1.0/_ic1l7dy*_hylurtp*_ef0votq-_88d0tgs/_rdu4xsz/_3zezqac*_mabsi9n*_4gc4ryg-_88d0tgs/_pmlvak9/_0j36tcf*_mabsi9n*_4gc4ryg-_88d0tgs/_dsc4e5h/_3zezqac*a*_4gc4ryg-3.0*a*_88d0tgs/_jnvc7hr*_4gc4ryg)-_kbp2qd6/_kjnawng*(_a41izxt-1.0/_ic1l7dy*(_hylurtp*_yxyur4s*_ef0votq+a*_i18k0ax)+a*_zuj2xor*_yxyur4s*_ef0votq/_rdu4xsz-1.0/_ic1l7dy/_j4t0sek*(_88d0tgs*_a41izxt-a*_yxyur4s*_ef0votq/_ic1l7dy*_3fpoyip))*_9rdpxdf+_kbp2qd6/_j4t0sek*(1.0/_rdu4xsz*(_hylurtp*_yxyur4s*_ef0votq+a*_i18k0ax)*_4gc4ryg-1.0/_ic1l7dy*_hylurtp*_i18k0ax*_ef0votq+a*_zuj2xor*_i18k0ax*_ef0votq/_rdu4xsz-3.0*a*_zuj2xor*_yxyur4s*_ef0votq/_jnvc7hr*_4gc4ryg+1.0/_rdu4xsz/_j4t0sek*(_88d0tgs*_a41izxt-a*_yxyur4s*_ef0votq/_ic1l7dy*_3fpoyip)*_4gc4ryg+1.0/_ic1l7dy/_kjnawng*(_88d0tgs*_a41izxt-a*_yxyur4s*_ef0votq/_ic1l7dy*_3fpoyip)*_9rdpxdf-1.0/_ic1l7dy/_j4t0sek*(-a*_i18k0ax*_ef0votq/_ic1l7dy*_3fpoyip+a*_yxyur4s*_ef0votq/_rdu4xsz*_3fpoyip*_4gc4ryg-a*_yxyur4s*_ef0votq/_pmlvak9*_i18k0ax*_4gc4ryg)))/3.14159265358979323846264338327950288/_sxlmt2p)
        + _e4q6qdy/2.0*(0.25*(_7yvdqbs*_hylurtp*_m3o55ix*_nsrtwlx-_hylurtp*_ovd79c9/_0j36tcf*(_b74x2ek*_ef0votq+_4gc4ryg/_3zezqac*_j6fp8hq)/_ic1l7dy*_4gc4ryg+_hylurtp*_ovd79c9/_3zezqac*(-a/_rdu4xsz*_4gc4ryg*_ef0votq+1.0/_3zezqac*_j6fp8hq-_6c2cu7o/_0j36tcf*_j6fp8hq/_ic1l7dy-_6c2cu7o/_3zezqac*a/_rdu4xsz)+_hylurtp*_ovd79c9*_ef0votq/_kjnawng*_vasifpw*_9rdpxdf+_hylurtp*_ovd79c9*_ef0votq/_j4t0sek*a/_rdu4xsz/_i18k0ax*_4gc4ryg+3.0*a*_ovd79c9*_kbp2qd6*_ef0votq/_jnvc7hr*_4gc4ryg-_ovd79c9*_kbp2qd6/_rdu4xsz/_3zezqac*(_hylurtp*_ef0votq-2.0*_gdc5tao*_4gc4ryg/_3zezqac-a*_4gc4ryg/_ic1l7dy*(1.0/_ic1l7dy+1.0/_3zezqac))*_4gc4ryg-_ovd79c9*_kbp2qd6/_pmlvak9/_0j36tcf*(_hylurtp*_ef0votq-2.0*_gdc5tao*_4gc4ryg/_3zezqac-a*_4gc4ryg/_ic1l7dy*(1.0/_ic1l7dy+1.0/_3zezqac))*_4gc4ryg+_ovd79c9*_kbp2qd6/_ic1l7dy/_3zezqac*(-2.0*_gdc5tao/_3zezqac+2.0*_gdc5tao*_6c2cu7o/_0j36tcf/_ic1l7dy-a/_ic1l7dy*(1.0/_ic1l7dy+1.0/_3zezqac)+a*_6c2cu7o/_rdu4xsz*(1.0/_ic1l7dy+1.0/_3zezqac)-a*_4gc4ryg/_ic1l7dy*(-1.0/_rdu4xsz*_4gc4ryg-1.0/_0j36tcf/_ic1l7dy*_4gc4ryg))-_ovd79c9*_kbp2qd6*_ef0votq/_rdu4xsz/_j4t0sek*(_xw1tbrw*_i18k0ax+_3fpoyip/_j4t0sek*_vasifpw+a*_zuj2xor/_pmlvak9/_i18k0ax)*_4gc4ryg-_ovd79c9*_kbp2qd6*_ef0votq/_ic1l7dy/_kjnawng*(_xw1tbrw*_i18k0ax+_3fpoyip/_j4t0sek*_vasifpw+a*_zuj2xor/_pmlvak9/_i18k0ax)*_9rdpxdf+_ovd79c9*_kbp2qd6*_ef0votq/_ic1l7dy/_j4t0sek*(1.0/_ic1l7dy*_i18k0ax*_4gc4ryg/_j4t0sek*_vasifpw-_3fpoyip/_kjnawng*_vasifpw*_9rdpxdf-_3fpoyip/_j4t0sek*a/_rdu4xsz/_i18k0ax*_4gc4ryg-2.0*a*_zuj2xor/_dsc4e5h/_i18k0ax*_4gc4ryg))/3.14159265358979323846264338327950288/_sxlmt2p)
        + _4lmct64/2.0*(0.25*(_hylurtp*(-_nn16uxf*_9rdpxdf/_j4t0sek-1.0/_3zezqac*(1.0+a/_ic1l7dy)+_6c2cu7o/_0j36tcf*(1.0+a/_ic1l7dy)/_ic1l7dy+_6c2cu7o/_3zezqac*a/_rdu4xsz+_i18k0ax/_j4t0sek*_k0qtl3q-_yxyur4s/_kjnawng*_k0qtl3q*_9rdpxdf-_yxyur4s/_j4t0sek*a/_rdu4xsz*_4gc4ryg)+_kbp2qd6/_ic1l7dy*_cv3tlcx-_6c2cu7o*_kbp2qd6/_rdu4xsz*_cv3tlcx+_4gc4ryg*_kbp2qd6/_ic1l7dy*(-2.0*a/_dsc4e5h*_4gc4ryg-1.0/_0j36tcf/_ic1l7dy*_4gc4ryg)+_kbp2qd6/_kjnawng*(_nn16uxf*(_i18k0ax-a/_ic1l7dy)+_yxyur4s/_ic1l7dy*(1.0+a*_zuj2xor/_pmlvak9)-1.0/_ic1l7dy/_j4t0sek*(_88d0tgs*_i18k0ax*_nn16uxf-a*_yxyur4s/_ic1l7dy*_3fpoyip))*_9rdpxdf-_kbp2qd6/_j4t0sek*(_nn16uxf*a/_rdu4xsz*_4gc4ryg+_i18k0ax/_ic1l7dy*(1.0+a*_zuj2xor/_pmlvak9)-_yxyur4s/_rdu4xsz*(1.0+a*_zuj2xor/_pmlvak9)*_4gc4ryg-2.0*_yxyur4s/_jnvc7hr*a*_zuj2xor*_4gc4ryg+1.0/_rdu4xsz/_j4t0sek*(_88d0tgs*_i18k0ax*_nn16uxf-a*_yxyur4s/_ic1l7dy*_3fpoyip)*_4gc4ryg+1.0/_ic1l7dy/_kjnawng*(_88d0tgs*_i18k0ax*_nn16uxf-a*_yxyur4s/_ic1l7dy*_3fpoyip)*_9rdpxdf-1.0/_ic1l7dy/_j4t0sek*(-a*_i18k0ax/_ic1l7dy*_3fpoyip+a*_yxyur4s/_rdu4xsz*_3fpoyip*_4gc4ryg-a*_yxyur4s/_pmlvak9*_i18k0ax*_4gc4ryg)))/3.14159265358979323846264338327950288/_sxlmt2p);
    _gntb34k = _jjq9k6r/2.0*(0.25*(_xw1tbrw*_hylurtp*_fshb0oq*_nsrtwlx-_hylurtp*_ovd79c9/_0j36tcf*(_8zgvwwu*_ef0votq-_4gc4ryg/_3zezqac*_j6fp8hq)*_vork1dj+_hylurtp*_ovd79c9/_3zezqac*(0.5*a/_rdu4xsz*2.0*_zuj2xor*_ef0votq+_4gc4ryg/_0j36tcf*_j6fp8hq*_vork1dj+0.5*_4gc4ryg/_3zezqac*a/_rdu4xsz*2.0*_zuj2xor)-_hylurtp*_ovd79c9*_i18k0ax*_ef0votq/_kjnawng*_k0qtl3q*_5td44k1-0.5*_hylurtp*_ovd79c9*_i18k0ax*_ef0votq/_j4t0sek*a/_rdu4xsz*2.0*_zuj2xor+a/_rdu4xsz*_ovd79c9*_ef0votq-1.5*a*_ovd79c9*_kbp2qd6*_ef0votq/_jnvc7hr*2.0*_zuj2xor+_ovd79c9/_ic1l7dy/_3zezqac*(-_hylurtp*_ef0votq+_4gc4ryg/_3zezqac*_mabsi9n+a*_4gc4ryg/_pmlvak9)-0.5*_ovd79c9*_kbp2qd6/_rdu4xsz/_3zezqac*(-_hylurtp*_ef0votq+_4gc4ryg/_3zezqac*_mabsi9n+a*_4gc4ryg/_pmlvak9)*2.0*_zuj2xor-_ovd79c9*_kbp2qd6/_ic1l7dy/_0j36tcf*(-_hylurtp*_ef0votq+_4gc4ryg/_3zezqac*_mabsi9n+a*_4gc4ryg/_pmlvak9)*_vork1dj+_ovd79c9*_kbp2qd6/_ic1l7dy/_3zezqac*(-_4gc4ryg/_0j36tcf*_mabsi9n*_vork1dj-0.5*_4gc4ryg/_3zezqac*a/_rdu4xsz*2.0*_zuj2xor-a*_4gc4ryg/_dsc4e5h*2.0*_zuj2xor)+_ovd79c9/_ic1l7dy/_j4t0sek*(_i18k0ax/_j4t0sek*(_3fpoyip*(_hylurtp*_i18k0ax-a/_ic1l7dy)*_ef0votq+_7yvdqbs*_w9pfb3h*_i18k0ax)-a*_zuj2xor*_i18k0ax*_ef0votq/_pmlvak9)-0.5*_ovd79c9*_kbp2qd6/_rdu4xsz/_j4t0sek*(_i18k0ax/_j4t0sek*(_3fpoyip*(_hylurtp*_i18k0ax-a/_ic1l7dy)*_ef0votq+_7yvdqbs*_w9pfb3h*_i18k0ax)-a*_zuj2xor*_i18k0ax*_ef0votq/_pmlvak9)*2.0*_zuj2xor-_ovd79c9*_kbp2qd6/_ic1l7dy/_kjnawng*(_i18k0ax/_j4t0sek*(_3fpoyip*(_hylurtp*_i18k0ax-a/_ic1l7dy)*_ef0votq+_7yvdqbs*_w9pfb3h*_i18k0ax)-a*_zuj2xor*_i18k0ax*_ef0votq/_pmlvak9)*_5td44k1+_ovd79c9*_kbp2qd6/_ic1l7dy/_j4t0sek*(-_i18k0ax/_kjnawng*(_3fpoyip*(_hylurtp*_i18k0ax-a/_ic1l7dy)*_ef0votq+_7yvdqbs*_w9pfb3h*_i18k0ax)*_5td44k1+_i18k0ax/_j4t0sek*((_i18k0ax*_zuj2xor/_ic1l7dy+1.0)*(_hylurtp*_i18k0ax-a/_ic1l7dy)*_ef0votq+0.5*_3fpoyip*a/_rdu4xsz*2.0*_zuj2xor*_ef0votq+0.5*_7yvdqbs/_ic1l7dy*_nn16uxf*2.0*_zuj2xor*_i18k0ax)-a*_i18k0ax*_ef0votq/_pmlvak9+a*_zuj2xor*_i18k0ax*_ef0votq/_dsc4e5h*2.0*_zuj2xor))/3.14159265358979323846264338327950288/_sxlmt2p)
        + _e4q6qdy/2.0*(0.25*(_hylurtp*((_7yvdqbs*_nsrtwlx+_gdc5tao)*_vork1dj/_3zezqac-(_7yvdqbs*_nsrtwlx+1.0)*_i18k0ax*_5td44k1/_j4t0sek)-_hylurtp/_0j36tcf*(-_hylurtp*_4gc4ryg*_ef0votq+_gdc5tao*_zuj2xor-a+a*_4gc4ryg*_ef0votq/_ic1l7dy+_6c2cu7o/_3zezqac*_j6fp8hq)*_vork1dj+_hylurtp/_3zezqac*(_gdc5tao-0.5*a*_4gc4ryg*_ef0votq/_rdu4xsz*2.0*_zuj2xor-_6c2cu7o/_0j36tcf*_j6fp8hq*_vork1dj-0.5*_6c2cu7o/_3zezqac*a/_rdu4xsz*2.0*_zuj2xor)+_hylurtp*_ef0votq/_kjnawng*(_yxyur4s*_i18k0ax-a*_w9pfb3h/_ic1l7dy/_i18k0ax)*_5td44k1-_hylurtp*_ef0votq/_j4t0sek*(_i18k0ax*_nn16uxf-0.5*a/_pmlvak9*_nn16uxf*2.0*_zuj2xor/_i18k0ax+0.5*a*_w9pfb3h/_rdu4xsz/_i18k0ax*2.0*_zuj2xor)-a/_rdu4xsz*_4gc4ryg*_ef0votq+1.5*a*_4gc4ryg*_kbp2qd6*_ef0votq/_jnvc7hr*2.0*_zuj2xor+1.0/_3zezqac*(2.0*_gdc5tao+1.0/_ic1l7dy*(_hylurtp*_4gc4ryg*_ef0votq+a)-_6c2cu7o/_ic1l7dy/_3zezqac*_mabsi9n-a*_6c2cu7o/_rdu4xsz)-_kbp2qd6/_0j36tcf*(2.0*_gdc5tao+1.0/_ic1l7dy*(_hylurtp*_4gc4ryg*_ef0votq+a)-_6c2cu7o/_ic1l7dy/_3zezqac*_mabsi9n-a*_6c2cu7o/_rdu4xsz)*_vork1dj+_kbp2qd6/_3zezqac*(-0.5/_rdu4xsz*(_hylurtp*_4gc4ryg*_ef0votq+a)*2.0*_zuj2xor+0.5*_6c2cu7o/_rdu4xsz/_3zezqac*_mabsi9n*2.0*_zuj2xor+_6c2cu7o/_ic1l7dy/_0j36tcf*_mabsi9n*_vork1dj+0.5*_6c2cu7o/_dsc4e5h/_3zezqac*a*2.0*_zuj2xor+1.5*a*_6c2cu7o/_jnvc7hr*2.0*_zuj2xor)+_ef0votq/_j4t0sek*(-_i18k0ax*_nn16uxf+a*_4gc4ryg*_zuj2xor/_rdu4xsz/_i18k0ax+_w9pfb3h/_ic1l7dy*(_7yvdqbs*_i18k0ax-_3fpoyip/_j4t0sek*_vasifpw))-_kbp2qd6*_ef0votq/_kjnawng*(-_i18k0ax*_nn16uxf+a*_4gc4ryg*_zuj2xor/_rdu4xsz/_i18k0ax+_w9pfb3h/_ic1l7dy*(_7yvdqbs*_i18k0ax-_3fpoyip/_j4t0sek*_vasifpw))*_5td44k1+_kbp2qd6*_ef0votq/_j4t0sek*(a/_rdu4xsz/_i18k0ax*_4gc4ryg-1.5*a*_4gc4ryg*_zuj2xor/_jnvc7hr/_i18k0ax*2.0*_zuj2xor+0.5/_pmlvak9*_nn16uxf*2.0*_zuj2xor*(_7yvdqbs*_i18k0ax-_3fpoyip/_j4t0sek*_vasifpw)-0.5*_w9pfb3h/_rdu4xsz*(_7yvdqbs*_i18k0ax-_3fpoyip/_j4t0sek*_vasifpw)*2.0*_zuj2xor+_w9pfb3h/_ic1l7dy*(-(_i18k0ax*_zuj2xor/_ic1l7dy+1.0)/_j4t0sek*_vasifpw+_3fpoyip/_kjnawng*_vasifpw*_5td44k1+0.5*_3fpoyip/_j4t0sek*a/_rdu4xsz/_i18k0ax*2.0*_zuj2xor)))/3.14159265358979323846264338327950288/_sxlmt2p)
        + _4lmct64/2.0*(0.25*(_hylurtp*(-_ovd79c9/_0j36tcf*(1.0+a/_ic1l7dy)*_vork1dj-0.5*_ovd79c9/_3zezqac*a/_rdu4xsz*2.0*_zuj2xor+_ovd79c9*_i18k0ax/_kjnawng*_k0qtl3q*_5td44k1+0.5*_ovd79c9*_i18k0ax/_j4t0sek*a/_rdu4xsz*2.0*_zuj2xor)-_ovd79c9/_ic1l7dy*_cv3tlcx+0.5*_ovd79c9*_kbp2qd6/_rdu4xsz*_cv3tlcx*2.0*_zuj2xor-_ovd79c9*_kbp2qd6/_ic1l7dy*(-a/_dsc4e5h*2.0*_zuj2xor-1.0/_0j36tcf*_vork1dj)+_ovd79c9*_i18k0ax/_ic1l7dy/_j4t0sek*(_3fpoyip/_j4t0sek*_k0qtl3q+a*_zuj2xor/_pmlvak9)-0.5*_ovd79c9*_kbp2qd6*_i18k0ax/_rdu4xsz/_j4t0sek*(_3fpoyip/_j4t0sek*_k0qtl3q+a*_zuj2xor/_pmlvak9)*2.0*_zuj2xor-_ovd79c9*_kbp2qd6*_i18k0ax/_ic1l7dy/_kjnawng*(_3fpoyip/_j4t0sek*_k0qtl3q+a*_zuj2xor/_pmlvak9)*_5td44k1+_ovd79c9*_kbp2qd6*_i18k0ax/_ic1l7dy/_j4t0sek*((_i18k0ax*_zuj2xor/_ic1l7dy+1.0)/_j4t0sek*_k0qtl3q-_3fpoyip/_kjnawng*_k0qtl3q*_5td44k1-0.5*_3fpoyip/_j4t0sek*a/_rdu4xsz*2.0*_zuj2xor+a/_pmlvak9-a*_zuj2xor/_dsc4e5h*2.0*_zuj2xor))/3.14159265358979323846264338327950288/_sxlmt2p)
        + _jjq9k6r/2.0*(0.25*(_7yvdqbs*(_hylurtp*_m3o55ix*_ef0votq-_4gc4ryg/_0j36tcf*_mabsi9n/_ic1l7dy*_ovd79c9-_ovd79c9/_3zezqac*a/_rdu4xsz*_4gc4ryg+_ovd79c9*_i18k0ax/_kjnawng*_k0qtl3q*_9rdpxdf+_ovd79c9*_i18k0ax/_j4t0sek*a/_rdu4xsz*_4gc4ryg)-_ovd79c9*_kbp2qd6/_rdu4xsz*(2.0*_gdc5tao/_3zezqac+a/_pmlvak9)*_4gc4ryg+_ovd79c9*_kbp2qd6/_ic1l7dy*(-2.0*_gdc5tao/_0j36tcf/_ic1l7dy*_4gc4ryg-2.0*a/_dsc4e5h*_4gc4ryg)-_ovd79c9*_kbp2qd6*_i18k0ax/_rdu4xsz/_j4t0sek*(1.0-2.0*_gdc5tao-_3fpoyip/_j4t0sek*_k0qtl3q-a*_zuj2xor/_pmlvak9)*_4gc4ryg-_ovd79c9*_kbp2qd6*_i18k0ax/_ic1l7dy/_kjnawng*(1.0-2.0*_gdc5tao-_3fpoyip/_j4t0sek*_k0qtl3q-a*_zuj2xor/_pmlvak9)*_9rdpxdf+_ovd79c9*_kbp2qd6*_i18k0ax/_ic1l7dy/_j4t0sek*(-1.0/_ic1l7dy*_i18k0ax*_4gc4ryg/_j4t0sek*_k0qtl3q+_3fpoyip/_kjnawng*_k0qtl3q*_9rdpxdf+_3fpoyip/_j4t0sek*a/_rdu4xsz*_4gc4ryg+2.0*a*_zuj2xor/_dsc4e5h*_4gc4ryg))/3.14159265358979323846264338327950288/_sxlmt2p)
        + _e4q6qdy/2.0*(0.25*(_xw1tbrw*_hylurtp*_ef0votq*(1.0/_ic1l7dy*_4gc4ryg/_3zezqac-_i18k0ax*_9rdpxdf/_j4t0sek)-_7yvdqbs/_3zezqac*_mabsi9n+_7yvdqbs*_6c2cu7o/_0j36tcf*_mabsi9n/_ic1l7dy+_7yvdqbs*_6c2cu7o/_3zezqac*a/_rdu4xsz+_7yvdqbs*_i18k0ax/_j4t0sek*_k0qtl3q-_7yvdqbs*_yxyur4s/_kjnawng*_k0qtl3q*_9rdpxdf-_7yvdqbs*_yxyur4s/_j4t0sek*a/_rdu4xsz*_4gc4ryg-_kbp2qd6/_rdu4xsz*(_hylurtp*_ef0votq-2.0*_gdc5tao*_4gc4ryg/_3zezqac-a*_4gc4ryg/_pmlvak9)*_4gc4ryg+_kbp2qd6/_ic1l7dy*(-2.0*_gdc5tao/_3zezqac+2.0*_gdc5tao*_6c2cu7o/_0j36tcf/_ic1l7dy-a/_pmlvak9+2.0*a*_6c2cu7o/_dsc4e5h)+_kbp2qd6/_kjnawng*(_i18k0ax*_nn16uxf+_3fpoyip*_ef0votq/_ic1l7dy*(_7yvdqbs*_i18k0ax-_3fpoyip/_j4t0sek)+a/_ic1l7dy*(_nn16uxf-_zuj2xor*_yxyur4s/_pmlvak9-_yxyur4s*_3fpoyip/_ic1l7dy/_j4t0sek))*_9rdpxdf-_kbp2qd6/_j4t0sek*(1.0/_pmlvak9*_i18k0ax*_4gc4ryg*_ef0votq*(_7yvdqbs*_i18k0ax-_3fpoyip/_j4t0sek)-_3fpoyip*_ef0votq/_rdu4xsz*(_7yvdqbs*_i18k0ax-_3fpoyip/_j4t0sek)*_4gc4ryg+_3fpoyip*_ef0votq/_ic1l7dy*(-1.0/_ic1l7dy*_i18k0ax*_4gc4ryg/_j4t0sek+_3fpoyip/_kjnawng*_9rdpxdf)-a/_rdu4xsz*(_nn16uxf-_zuj2xor*_yxyur4s/_pmlvak9-_yxyur4s*_3fpoyip/_ic1l7dy/_j4t0sek)*_4gc4ryg+a/_ic1l7dy*(-_zuj2xor*_i18k0ax/_pmlvak9+2.0*_zuj2xor*_yxyur4s/_dsc4e5h*_4gc4ryg-_i18k0ax*_3fpoyip/_ic1l7dy/_j4t0sek-_yxyur4s/_pmlvak9*_i18k0ax*_4gc4ryg/_j4t0sek+_yxyur4s*_3fpoyip/_rdu4xsz/_j4t0sek*_4gc4ryg+_yxyur4s*_3fpoyip/_ic1l7dy/_kjnawng*_9rdpxdf)))/3.14159265358979323846264338327950288/_sxlmt2p)
        + _4lmct64/2.0*(0.25*(_7yvdqbs*_m3o55ix-_7yvdqbs*_ovd79c9*_nn16uxf/_kjnawng*_k0qtl3q*_9rdpxdf-_7yvdqbs*_ovd79c9*_nn16uxf/_j4t0sek*a/_rdu4xsz*_4gc4ryg-_ovd79c9*_kbp2qd6*_nn16uxf/_rdu4xsz/_j4t0sek*(1.0+_3fpoyip/_j4t0sek*_k0qtl3q+a*_zuj2xor/_pmlvak9)*_4gc4ryg-_ovd79c9*_kbp2qd6*_nn16uxf/_ic1l7dy/_kjnawng*(1.0+_3fpoyip/_j4t0sek*_k0qtl3q+a*_zuj2xor/_pmlvak9)*_9rdpxdf+_ovd79c9*_kbp2qd6*_nn16uxf/_ic1l7dy/_j4t0sek*(1.0/_ic1l7dy*_i18k0ax*_4gc4ryg/_j4t0sek*_k0qtl3q-_3fpoyip/_kjnawng*_k0qtl3q*_9rdpxdf-_3fpoyip/_j4t0sek*a/_rdu4xsz*_4gc4ryg-2.0*a*_zuj2xor/_dsc4e5h*_4gc4ryg))/3.14159265358979323846264338327950288/_sxlmt2p);
    _5wryesu = _jjq9k6r/2.0*(0.25*(_hylurtp*((_7yvdqbs*_nsrtwlx-_gdc5tao)*_vork1dj/_3zezqac-(_7yvdqbs*_nsrtwlx+1.0-2.0*_gdc5tao)*_i18k0ax*_5td44k1/_j4t0sek)+_hylurtp/_0j36tcf*(_4gc4ryg*_ef0votq*_8zgvwwu+_gdc5tao*_zuj2xor-a+_88d0tgs/_3zezqac*_j6fp8hq)*_vork1dj-_hylurtp/_3zezqac*(0.5*a*_4gc4ryg*_ef0votq/_rdu4xsz*2.0*_zuj2xor+_gdc5tao-_88d0tgs/_0j36tcf*_j6fp8hq*_vork1dj-0.5*_88d0tgs/_3zezqac*a/_rdu4xsz*2.0*_zuj2xor)-_hylurtp*_nn16uxf*_ef0votq/_j4t0sek*_k0qtl3q+_hylurtp*_yxyur4s*_ef0votq/_kjnawng*_k0qtl3q*_5td44k1+0.5*_hylurtp*_yxyur4s*_ef0votq/_j4t0sek*a/_rdu4xsz*2.0*_zuj2xor-a/_rdu4xsz*_4gc4ryg*_ef0votq+1.5*a*_4gc4ryg*_kbp2qd6*_ef0votq/_jnvc7hr*2.0*_zuj2xor+1.0/_3zezqac*(-2.0*_gdc5tao+1.0/_ic1l7dy*(_hylurtp*_4gc4ryg*_ef0votq-a)+_88d0tgs/_ic1l7dy/_3zezqac*_mabsi9n+a*_88d0tgs/_rdu4xsz)-_kbp2qd6/_0j36tcf*(-2.0*_gdc5tao+1.0/_ic1l7dy*(_hylurtp*_4gc4ryg*_ef0votq-a)+_88d0tgs/_ic1l7dy/_3zezqac*_mabsi9n+a*_88d0tgs/_rdu4xsz)*_vork1dj+_kbp2qd6/_3zezqac*(-0.5/_rdu4xsz*(_hylurtp*_4gc4ryg*_ef0votq-a)*2.0*_zuj2xor-0.5*_88d0tgs/_rdu4xsz/_3zezqac*_mabsi9n*2.0*_zuj2xor-_88d0tgs/_ic1l7dy/_0j36tcf*_mabsi9n*_vork1dj-0.5*_88d0tgs/_dsc4e5h/_3zezqac*a*2.0*_zuj2xor-1.5*a*_88d0tgs/_jnvc7hr*2.0*_zuj2xor)+1.0/_j4t0sek*(_a41izxt-1.0/_ic1l7dy*(_hylurtp*_yxyur4s*_ef0votq+a*_i18k0ax)+a*_zuj2xor*_yxyur4s*_ef0votq/_rdu4xsz-1.0/_ic1l7dy/_j4t0sek*(_88d0tgs*_a41izxt-a*_yxyur4s*_ef0votq/_ic1l7dy*_3fpoyip))-_kbp2qd6/_kjnawng*(_a41izxt-1.0/_ic1l7dy*(_hylurtp*_yxyur4s*_ef0votq+a*_i18k0ax)+a*_zuj2xor*_yxyur4s*_ef0votq/_rdu4xsz-1.0/_ic1l7dy/_j4t0sek*(_88d0tgs*_a41izxt-a*_yxyur4s*_ef0votq/_ic1l7dy*_3fpoyip))*_5td44k1+_kbp2qd6/_j4t0sek*(0.5/_rdu4xsz*(_hylurtp*_yxyur4s*_ef0votq+a*_i18k0ax)*2.0*_zuj2xor-1.0/_ic1l7dy*_hylurtp*_nn16uxf*_ef0votq+a*_yxyur4s*_ef0votq/_rdu4xsz+a*_zuj2xor*_nn16uxf*_ef0votq/_rdu4xsz-1.5*a*_zuj2xor*_yxyur4s*_ef0votq/_jnvc7hr*2.0*_zuj2xor+0.5/_rdu4xsz/_j4t0sek*(_88d0tgs*_a41izxt-a*_yxyur4s*_ef0votq/_ic1l7dy*_3fpoyip)*2.0*_zuj2xor+1.0/_ic1l7dy/_kjnawng*(_88d0tgs*_a41izxt-a*_yxyur4s*_ef0votq/_ic1l7dy*_3fpoyip)*_5td44k1-1.0/_ic1l7dy/_j4t0sek*(-a*_nn16uxf*_ef0votq/_ic1l7dy*_3fpoyip+0.5*a*_yxyur4s*_ef0votq/_rdu4xsz*_3fpoyip*2.0*_zuj2xor-a*_yxyur4s*_ef0votq/_ic1l7dy*(_i18k0ax*_zuj2xor/_ic1l7dy+1.0))))/3.14159265358979323846264338327950288/_sxlmt2p)
        + _e4q6qdy/2.0*(0.25*(_7yvdqbs*_hylurtp*_fshb0oq*_nsrtwlx-_hylurtp*_ovd79c9/_0j36tcf*(_b74x2ek*_ef0votq+_4gc4ryg/_3zezqac*_j6fp8hq)*_vork1dj+_hylurtp*_ovd79c9/_3zezqac*(-0.5*a/_rdu4xsz*2.0*_zuj2xor*_ef0votq-_4gc4ryg/_0j36tcf*_j6fp8hq*_vork1dj-0.5*_4gc4ryg/_3zezqac*a/_rdu4xsz*2.0*_zuj2xor)+_hylurtp*_ovd79c9*_ef0votq/_kjnawng*_vasifpw*_5td44k1+0.5*_hylurtp*_ovd79c9*_ef0votq/_j4t0sek*a/_rdu4xsz/_i18k0ax*2.0*_zuj2xor-a/_rdu4xsz*_ovd79c9*_ef0votq+1.5*a*_ovd79c9*_kbp2qd6*_ef0votq/_jnvc7hr*2.0*_zuj2xor+_ovd79c9/_ic1l7dy/_3zezqac*(_hylurtp*_ef0votq-2.0*_gdc5tao*_4gc4ryg/_3zezqac-a*_4gc4ryg/_ic1l7dy*(1.0/_ic1l7dy+1.0/_3zezqac))-0.5*_ovd79c9*_kbp2qd6/_rdu4xsz/_3zezqac*(_hylurtp*_ef0votq-2.0*_gdc5tao*_4gc4ryg/_3zezqac-a*_4gc4ryg/_ic1l7dy*(1.0/_ic1l7dy+1.0/_3zezqac))*2.0*_zuj2xor-_ovd79c9*_kbp2qd6/_ic1l7dy/_0j36tcf*(_hylurtp*_ef0votq-2.0*_gdc5tao*_4gc4ryg/_3zezqac-a*_4gc4ryg/_ic1l7dy*(1.0/_ic1l7dy+1.0/_3zezqac))*_vork1dj+_ovd79c9*_kbp2qd6/_ic1l7dy/_3zezqac*(2.0*_gdc5tao*_4gc4ryg/_0j36tcf*_vork1dj+0.5*a*_4gc4ryg/_rdu4xsz*(1.0/_ic1l7dy+1.0/_3zezqac)*2.0*_zuj2xor-a*_4gc4ryg/_ic1l7dy*(-0.5/_rdu4xsz*2.0*_zuj2xor-1.0/_0j36tcf*_vork1dj))+_ovd79c9*_ef0votq/_ic1l7dy/_j4t0sek*(_xw1tbrw*_i18k0ax+_3fpoyip/_j4t0sek*_vasifpw+a*_zuj2xor/_pmlvak9/_i18k0ax)-0.5*_ovd79c9*_kbp2qd6*_ef0votq/_rdu4xsz/_j4t0sek*(_xw1tbrw*_i18k0ax+_3fpoyip/_j4t0sek*_vasifpw+a*_zuj2xor/_pmlvak9/_i18k0ax)*2.0*_zuj2xor-_ovd79c9*_kbp2qd6*_ef0votq/_ic1l7dy/_kjnawng*(_xw1tbrw*_i18k0ax+_3fpoyip/_j4t0sek*_vasifpw+a*_zuj2xor/_pmlvak9/_i18k0ax)*_5td44k1+_ovd79c9*_kbp2qd6*_ef0votq/_ic1l7dy/_j4t0sek*((_i18k0ax*_zuj2xor/_ic1l7dy+1.0)/_j4t0sek*_vasifpw-_3fpoyip/_kjnawng*_vasifpw*_5td44k1-0.5*_3fpoyip/_j4t0sek*a/_rdu4xsz/_i18k0ax*2.0*_zuj2xor+a/_pmlvak9/_i18k0ax-a*_zuj2xor/_dsc4e5h/_i18k0ax*2.0*_zuj2xor))/3.14159265358979323846264338327950288/_sxlmt2p)
        + _4lmct64/2.0*(0.25*(_hylurtp*(-_nn16uxf*_5td44k1/_j4t0sek+_4gc4ryg/_0j36tcf*(1.0+a/_ic1l7dy)*_vork1dj+0.5*_4gc4ryg/_3zezqac*a/_rdu4xsz*2.0*_zuj2xor+_nn16uxf/_j4t0sek*_k0qtl3q-_yxyur4s/_kjnawng*_k0qtl3q*_5td44k1-0.5*_yxyur4s/_j4t0sek*a/_rdu4xsz*2.0*_zuj2xor)+_4gc4ryg/_ic1l7dy*_cv3tlcx-0.5*_4gc4ryg*_kbp2qd6/_rdu4xsz*_cv3tlcx*2.0*_zuj2xor+_4gc4ryg*_kbp2qd6/_ic1l7dy*(-a/_dsc4e5h*2.0*_zuj2xor-1.0/_0j36tcf*_vork1dj)-1.0/_j4t0sek*(_nn16uxf*(_i18k0ax-a/_ic1l7dy)+_yxyur4s/_ic1l7dy*(1.0+a*_zuj2xor/_pmlvak9)-1.0/_ic1l7dy/_j4t0sek*(_88d0tgs*_i18k0ax*_nn16uxf-a*_yxyur4s/_ic1l7dy*_3fpoyip))+_kbp2qd6/_kjnawng*(_nn16uxf*(_i18k0ax-a/_ic1l7dy)+_yxyur4s/_ic1l7dy*(1.0+a*_zuj2xor/_pmlvak9)-1.0/_ic1l7dy/_j4t0sek*(_88d0tgs*_i18k0ax*_nn16uxf-a*_yxyur4s/_ic1l7dy*_3fpoyip))*_5td44k1-_kbp2qd6/_j4t0sek*(0.5*_nn16uxf*a/_rdu4xsz*2.0*_zuj2xor+_nn16uxf/_ic1l7dy*(1.0+a*_zuj2xor/_pmlvak9)-0.5*_yxyur4s/_rdu4xsz*(1.0+a*_zuj2xor/_pmlvak9)*2.0*_zuj2xor+_yxyur4s/_ic1l7dy*(a/_pmlvak9-a*_zuj2xor/_dsc4e5h*2.0*_zuj2xor)+0.5/_rdu4xsz/_j4t0sek*(_88d0tgs*_i18k0ax*_nn16uxf-a*_yxyur4s/_ic1l7dy*_3fpoyip)*2.0*_zuj2xor+1.0/_ic1l7dy/_kjnawng*(_88d0tgs*_i18k0ax*_nn16uxf-a*_yxyur4s/_ic1l7dy*_3fpoyip)*_5td44k1-1.0/_ic1l7dy/_j4t0sek*(-a*_nn16uxf/_ic1l7dy*_3fpoyip+0.5*a*_yxyur4s/_rdu4xsz*_3fpoyip*2.0*_zuj2xor-a*_yxyur4s/_ic1l7dy*(_i18k0ax*_zuj2xor/_ic1l7dy+1.0))))/3.14159265358979323846264338327950288/_sxlmt2p)
        + _jjq9k6r/2.0*(0.25*(_7yvdqbs*(_hylurtp*_va1nkbt*_ef0votq+1.0/_3zezqac*_mabsi9n-_88d0tgs/_0j36tcf*_mabsi9n/_ic1l7dy-_88d0tgs/_3zezqac*a/_rdu4xsz-_i18k0ax/_j4t0sek*_k0qtl3q+_88d0tgs*_i18k0ax/_kjnawng*_k0qtl3q/_ic1l7dy+_88d0tgs*_i18k0ax/_j4t0sek*a/_rdu4xsz)+_kbp2qd6/_ic1l7dy*(2.0*_gdc5tao/_3zezqac+a/_pmlvak9)-_88d0tgs*_kbp2qd6/_rdu4xsz*(2.0*_gdc5tao/_3zezqac+a/_pmlvak9)+_ovd79c9*_kbp2qd6/_ic1l7dy*(-2.0*_gdc5tao/_0j36tcf/_ic1l7dy*_ovd79c9-2.0*a/_dsc4e5h*_ovd79c9)+_kbp2qd6*_i18k0ax/_ic1l7dy/_j4t0sek*(1.0-2.0*_gdc5tao-_3fpoyip/_j4t0sek*_k0qtl3q-a*_zuj2xor/_pmlvak9)-_88d0tgs*_kbp2qd6*_i18k0ax/_rdu4xsz/_j4t0sek*(1.0-2.0*_gdc5tao-_3fpoyip/_j4t0sek*_k0qtl3q-a*_zuj2xor/_pmlvak9)-_88d0tgs*_kbp2qd6*_i18k0ax/_pmlvak9/_kjnawng*(1.0-2.0*_gdc5tao-_3fpoyip/_j4t0sek*_k0qtl3q-a*_zuj2xor/_pmlvak9)+_ovd79c9*_kbp2qd6*_i18k0ax/_ic1l7dy/_j4t0sek*(-1.0/_ic1l7dy*_i18k0ax*_ovd79c9/_j4t0sek*_k0qtl3q+_3fpoyip/_kjnawng*_k0qtl3q/_ic1l7dy*_ovd79c9+_3fpoyip/_j4t0sek*a/_rdu4xsz*_ovd79c9+2.0*a*_zuj2xor/_dsc4e5h*_ovd79c9))/3.14159265358979323846264338327950288/_sxlmt2p)
        + _e4q6qdy/2.0*(0.25*(_xw1tbrw*_hylurtp*_ef0votq*(1.0/_ic1l7dy*_ovd79c9/_3zezqac-_i18k0ax/_ic1l7dy*_ovd79c9/_j4t0sek)+_7yvdqbs*_4gc4ryg/_0j36tcf*_mabsi9n/_ic1l7dy*_ovd79c9+_7yvdqbs*_4gc4ryg/_3zezqac*a/_rdu4xsz*_ovd79c9-_7yvdqbs*_yxyur4s/_kjnawng*_k0qtl3q/_ic1l7dy*_ovd79c9-_7yvdqbs*_yxyur4s/_j4t0sek*a/_rdu4xsz*_ovd79c9-_kbp2qd6/_rdu4xsz*(_hylurtp*_ef0votq-2.0*_gdc5tao*_4gc4ryg/_3zezqac-a*_4gc4ryg/_pmlvak9)*_ovd79c9+_kbp2qd6/_ic1l7dy*(2.0*_gdc5tao*_4gc4ryg/_0j36tcf/_ic1l7dy*_ovd79c9+2.0*a*_4gc4ryg/_dsc4e5h*_ovd79c9)+_kbp2qd6/_kjnawng*(_i18k0ax*_nn16uxf+_3fpoyip*_ef0votq/_ic1l7dy*(_7yvdqbs*_i18k0ax-_3fpoyip/_j4t0sek)+a/_ic1l7dy*(_nn16uxf-_zuj2xor*_yxyur4s/_pmlvak9-_yxyur4s*_3fpoyip/_ic1l7dy/_j4t0sek))/_ic1l7dy*_ovd79c9-_kbp2qd6/_j4t0sek*(1.0/_pmlvak9*_i18k0ax*_ovd79c9*_ef0votq*(_7yvdqbs*_i18k0ax-_3fpoyip/_j4t0sek)-_3fpoyip*_ef0votq/_rdu4xsz*(_7yvdqbs*_i18k0ax-_3fpoyip/_j4t0sek)*_ovd79c9+_3fpoyip*_ef0votq/_ic1l7dy*(-_i18k0ax/_ic1l7dy*_ovd79c9/_j4t0sek+_3fpoyip/_kjnawng/_ic1l7dy*_ovd79c9)-a/_rdu4xsz*(_nn16uxf-_zuj2xor*_yxyur4s/_pmlvak9-_yxyur4s*_3fpoyip/_ic1l7dy/_j4t0sek)*_ovd79c9+a/_ic1l7dy*(2.0*_zuj2xor*_yxyur4s/_dsc4e5h*_ovd79c9-_yxyur4s/_pmlvak9*_i18k0ax*_ovd79c9/_j4t0sek+_yxyur4s*_3fpoyip/_rdu4xsz/_j4t0sek*_ovd79c9+_yxyur4s*_3fpoyip/_pmlvak9/_kjnawng*_ovd79c9)))/3.14159265358979323846264338327950288/_sxlmt2p)
        + _4lmct64/2.0*(0.25*(_7yvdqbs*_va1nkbt+_7yvdqbs*_nn16uxf/_j4t0sek*_k0qtl3q-_7yvdqbs*_88d0tgs*_nn16uxf/_kjnawng*_k0qtl3q/_ic1l7dy-_7yvdqbs*_88d0tgs*_nn16uxf/_j4t0sek*a/_rdu4xsz+_kbp2qd6*_nn16uxf/_ic1l7dy/ _j4t0sek*(1.0+_3fpoyip/_j4t0sek*_k0qtl3q+a*_zuj2xor/_pmlvak9)-_88d0tgs*_kbp2qd6*_nn16uxf/_rdu4xsz/_j4t0sek*(1.0+_3fpoyip/_j4t0sek*_k0qtl3q+a*_zuj2xor/_pmlvak9)-_88d0tgs*_kbp2qd6*_nn16uxf/_pmlvak9/_kjnawng*(1.0+_3fpoyip/_j4t0sek*_k0qtl3q+a* _zuj2xor/_pmlvak9)+_ovd79c9*_kbp2qd6*_nn16uxf/_ic1l7dy/_j4t0sek*(1.0/_ic1l7dy*_i18k0ax*_ovd79c9/_j4t0sek*_k0qtl3q-_3fpoyip/_kjnawng* _k0qtl3q/_ic1l7dy*_ovd79c9-_3fpoyip/_j4t0sek*a/_rdu4xsz*_ovd79c9-2.0*a*_zuj2xor/_dsc4e5h*_ovd79c9))/3.14159265358979323846264338327950288/_sxlmt2p);
    _82n2atq[0] = _r59c81d; _82n2atq[1] = _o6ouz1e; _82n2atq[2] = _gntb34k;
    _82n2atq[3] = _ej9c0uu; _82n2atq[4] = _5wryesu; _82n2atq[5] = _4fcl5ef;
    return;
}
