!Directives that determine compilation of the program

!Profiling setup
#define use_profiling     .true.
#define profile_code      .true. 

#define zeror             0.0_wp
#define oner              1.0_wp
#define im                (0.0_wp,1.0_wp)
#define zeroc             (0.0_wp,0.0_wp)
#define onec              (1.0_wp,0.0_wp)

#ifdef ksamp
#define ORB               complex(wp)
#define KStrans           "C"
#define KSone             (1.0_wp,0.0_wp)
#define KSzero            (0.0_wp,0.0_wp)
#define KSoner            1.0_wp
#define KSzero            0.0_wp
#else
#define ORB               real(wp)
#define KStrans           "T"
#define KSone             1.0_wp
#define KSzero            0.0_wp
#define KSoner            1.0_wp
#define KSzeror           0.0_wp
#endif
