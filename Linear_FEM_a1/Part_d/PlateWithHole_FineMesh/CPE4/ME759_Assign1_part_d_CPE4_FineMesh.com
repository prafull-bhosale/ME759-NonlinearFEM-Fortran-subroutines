from driverConstants import *
from driverStandard import StandardAnalysis
import driverUtils, sys
options = {
    'SIMExt':'.sim',
    'abaquslm_license_file':'27001@matlab.iitb.ac.in;27800@localhost',
    'academic':RESEARCH,
    'ams':OFF,
    'analysisType':STANDARD,
    'applicationName':'analysis',
    'aqua':OFF,
    'ask_delete':OFF,
    'background':None,
    'beamSectGen':OFF,
    'biorid':OFF,
    'cavityTypes':[],
    'cavitydefinition':OFF,
    'cavparallel':OFF,
    'compile_cpp':['cl', '/c', '/W0', '/MD', '/TP', '/EHsc', '/DNDEBUG', '/DWIN32', '/DTP_IP', '/D_CONSOLE', '/DNTI', '/DFLT_LIC', '/DOL_DOC', '/D__LIB__', '/DHKS_NT', '/D_WINDOWS_SOURCE', '/DFAR=', '/D_WINDOWS', '/DABQ_WIN86_64', '%P', '/I%I', '/IC:\\Program Files\\Dassault Systemes\\SimulationServices\\V6R2018x'],
    'compile_fmu':['win64CmpWrp', '-m64', '-msvc9', 'cl', '/LD', '/D_WINDOWS', '/TC', '/W0', '/I%I', '/IC:\\Program Files\\Dassault Systemes\\SimulationServices\\V6R2018x'],
    'compile_fortran':['ifort', '/Qmkl:sequential', '/free', '/c', '/DABQ_WIN86_64', '/extend-source', '/fpp', '/iface:cref', '/recursive', '/Qauto-scalar', '/QxSSE3', '/QaxAVX', '/heap-arrays:1', '/include:%I', '/include:C:\\Program Files\\Dassault Systemes\\SimulationServices\\V6R2018x', '%P', '/names:lowercase', '/names:lowercase'],
    'complexFrequency':OFF,
    'contact':OFF,
    'cosimulation':OFF,
    'coupledProcedure':OFF,
    'cpus':1,
    'cse':OFF,
    'cyclicSymmetryModel':OFF,
    'directCyclic':OFF,
    'direct_port':'53985',
    'direct_solver':DMP,
    'doc_root':'http://help.3ds.com',
    'dsa':OFF,
    'dynamic':OFF,
    'filPrt':[],
    'fils':[],
    'finitesliding':OFF,
    'foundation':OFF,
    'geostatic':OFF,
    'geotech':OFF,
    'heatTransfer':OFF,
    'hysteresis':OFF,
    'impJobExpVars':{},
    'importJobList':[],
    'importer':OFF,
    'importerParts':OFF,
    'includes':[],
    'initialConditionsFile':OFF,
    'input':'ME759_Assign1_part_d_CPE4_FineMesh',
    'inputFormat':INP,
    'job':'ME759_Assign1_part_d_CPE4_FineMesh',
    'keyword_licenses':[],
    'lanczos':OFF,
    'libs':[],
    'license_server_type':FLEXNET,
    'link_exe':['LINK', '/nologo', '/INCREMENTAL:NO', '/subsystem:console', '/machine:AMD64', '/STACK:20000000', '/NODEFAULTLIB:LIBC.LIB', '/NODEFAULTLIB:LIBCMT.LIB', '/DEFAULTLIB:OLDNAMES.LIB', '/DEFAULTLIB:LIBIFCOREMD.LIB', '/DEFAULTLIB:LIBIFPORTMD.LIB', '/DEFAULTLIB:LIBMMD.LIB', '/DEFAULTLIB:kernel32.lib', '/DEFAULTLIB:user32.lib', '/DEFAULTLIB:advapi32.lib', '/FIXED:NO', '/LARGEADDRESSAWARE', '/out:%J', '%F', '%M', '%L', '%B', '%O', 'oldnames.lib', 'user32.lib', 'ws2_32.lib', 'netapi32.lib', 'advapi32.lib', 'msvcrt.lib', 'vcruntime.lib', 'ucrt.lib'],
    'link_sl':['LINK', '/nologo', '/NOENTRY', '/INCREMENTAL:NO', '/subsystem:console', '/machine:AMD64', '/NODEFAULTLIB:LIBC.LIB', '/NODEFAULTLIB:LIBCMT.LIB', '/DEFAULTLIB:OLDNAMES.LIB', '/DEFAULTLIB:LIBIFCOREMD.LIB', '/DEFAULTLIB:LIBIFPORTMD.LIB', '/DEFAULTLIB:LIBMMD.LIB', '/DEFAULTLIB:kernel32.lib', '/DEFAULTLIB:user32.lib', '/DEFAULTLIB:advapi32.lib', '/FIXED:NO', '/dll', '/def:%E', '/out:%U', '%F', '%A', '%L', '%B', 'oldnames.lib', 'user32.lib', 'ws2_32.lib', 'netapi32.lib', 'advapi32.lib', 'msvcrt.lib', 'vcruntime.lib', 'ucrt.lib'],
    'listener_name':'LAPTOP-BRO9FT16',
    'listener_resource':'8904',
    'magnetostatic':OFF,
    'massDiffusion':OFF,
    'memory':'90%',
    'message':None,
    'messaging_mechanism':'DIRECT',
    'modifiedTet':OFF,
    'moldflowFiles':[],
    'moldflowMaterial':OFF,
    'mp_environment_export':('ABAQUSLM_LICENSE_FILE', 'ABAQUS_CCI_DEBUG', 'ABAQUS_CSE_CURRCONFIGMAPPING', 'ABAQUS_CSE_RELTIMETOLERANCE', 'ABAQUS_LANG', 'ABAQUS_MPF_DIAGNOSTIC_LEVEL', 'ABA_ADM_ALIGNMENT', 'ABA_ADM_MINIMUMDECREASE', 'ABA_ADM_MINIMUMINCREASE', 'ABA_ALL_ADB_IN_TMPDIR', 'ABA_CM_BUFFERING', 'ABA_CM_BUFFERING_LIMIT', 'ABA_CUTOFF_SLAVEFACET_ANGLE', 'ABA_DMPSOLVER_BWDPARALLELOFF', 'ABA_ELP_SURFACE_SPLIT', 'ABA_ELP_SUSPEND', 'ABA_EXT_SIMOUTPUT', 'ABA_HOME', 'ABA_ITERATIVE_SOLVER_VERBOSE', 'ABA_MEMORY_MODE', 'ABA_MPI_MESSAGE_TRACKING', 'ABA_MPI_VERBOSE_LEVEL', 'ABA_NUM_INTEGRATION_POINTS_LINE3D', 'ABA_SHARED_SAVEDIR', 'ABA_PATH', 'ABA_PRE_DECOMPOSITION', 'ABA_RESOURCE_MONITOR', 'ABA_RESOURCE_USEMALLINFO', 'ABA_SYMBOLIC_GENERALCOLLAPSE', 'ABA_SYMBOLIC_GENERAL_MAXCLIQUERANK', 'ABA_TOSCA_PROTOTYPE', 'ABA_TOSCA_SEQFILES', 'ABA_UNIT_INDEPENDENT_CONTACT', 'ABQLMHANGLIMIT', 'ABQLMIMPL', 'ABQLMQUEUE', 'ABQLMUSER', 'ABQ_ACTIVATE_PTK', 'ABQ_CRTMALLOC', 'ABQ_DATACHECK', 'ABQ_DLALLOCATOR', 'ABQ_RECOVER', 'ABQ_RESTART', 'ABQ_SPLITFILE', 'ABQ_STD_ACCUM_CSLIP', 'ABQ_STD_ACTIVATE_BEAM_ROTATION', 'ABQ_STD_ALLOW_SURFACE_TO_BEAM', 'ABQ_XFEM_POREPRESSURE', 'ABQ_XPL_PARTITIONSIZE', 'ABQ_XPL_WINDOWDUMP', 'ACML_FAST_MALLOC', 'ACML_FAST_MALLOC_CHUNK_SIZE', 'ACML_FAST_MALLOC_DEBUG', 'ACML_FAST_MALLOC_MAX_CHUNKS', 'ADB_USE_OLDSLDB', 'ADB_USE_NEWSLDB', 'CCI_RENDEZVOUS', 'DOMAIN', 'DOMAIN_CPUS', 'DOUBLE_PRECISION', 'DSLS_CONFIG', 'FLEXLM_DIAGNOSTICS', 'FOR0006', 'FOR0064', 'FOR_DISABLE_DIAGNOSTIC_DISPLAY', 'FOR_IGNORE_EXCEPTIONS', 'IPATH_NO_CPUAFFINITY', 'LD_PRELOAD', 'MALLOC_MMAP_THRESHOLD_', 'MKL_DYNAMIC', 'MKL_NUM_THREADS', 'MPCCI_CODEID', 'MPCCI_DEBUG', 'MPCCI_JOBID', 'MPCCI_NETDEVICE', 'MPCCI_SERVER', 'MPCCI_TINFO', 'MPC_GANG', 'MPIEXEC_AFFINITY_TABLE', 'MPI_FLAGS', 'MPI_FLUSH_FCACHE', 'MPI_PROPAGATE_TSTP', 'MPI_RDMA_MSGSIZE', 'MPI_RDMA_NENVELOPE', 'MPI_SOCKBUFSIZE', 'MPI_USE_MALLOPT_MMAP_MAX', 'MPI_USE_MALLOPT_MMAP_THRESHOLD', 'MPI_USE_MALLOPT_SBRK_PROTECTION', 'MPI_WORKDIR', 'MP_NUMBER_OF_THREADS', 'MPICH_ND_ZCOPY_THRESHOLD', 'NCPUS', 'OMP_DYNAMIC', 'OMP_NUM_THREADS', 'OUTDIR', 'PAIDUP', 'PARALLEL_METHOD', 'RAIDEV_NDREG_LAZYMEM', 'SMA_PARENT', 'SMA_PLATFORM', 'SMA_WS', 'SIMULIA_COSIN_PATH'),
    'mp_file_system':(DETECT, DETECT),
    'mp_mode':THREADS,
    'mp_mode_requested':MPI,
    'mp_mpi_implementation':NATIVE,
    'mp_mpi_searchpath':['Microsoft MPI', 'Microsoft HPC Pack', 'Microsoft HPC Pack 2008 R2', 'Microsoft HPC Pack 2008', 'Microsoft HPC Pack 2008 SDK', 'Microsoft HPC Pack 2012'],
    'mp_mpirun_path':{MSSDK: 'C:\\Program Files\\Microsoft MPI', NATIVE: 'C:\\Program Files\\Microsoft MPI\\bin\\mpiexec.exe'},
    'mp_num_parallel_ftps':(4, 4),
    'mp_rsh_command':'dummy %H -l %U -n %C',
    'multiphysics':OFF,
    'noDmpDirect':[],
    'noMultiHost':[],
    'noMultiHostElemLoop':[],
    'no_domain_check':1,
    'onCaeGraphicsStartup':driverUtils.decodeFunction('begin 666 -\nM8P     H    1P   $,   !SY1   \'0  &H! &H" \'T  \'0  &H! &H# \'T!\nM \'0  &H! &H$ \'T" \'0  &H! &H% \'T# \'0  &H! &H& \'T$ \'0\' \'T% \'0(\nM \'T& \'0( \'T\' \'0( \'T( \'0) \'0* \'0+ \'0, &8$ \'T) \'0- \'T* \'0. \'T+\nM \'0  &H/ &H0 \'T, \'0  &H/ &H1 \'T- \'0  &H/ &H2 \'T. \'0  &H/ &H3\nM \'T/ \'0( \'T0 \'0( \'T1 \'0( \'T2 &0! \'T3 \'0( \'T4 \'0  &H/ &H4 \'T5\nM \'0  &H/ &H5 \'T6 \'0  &H/ &H6 \'T7 \'0  &H/ &H7 \'T8 \'0  &H/ &H8\nM \'T9 \'0  &H/ &H9 \'T: \'0  &H/ &H: \'T; \'0  &H/ &H; \'T< &0" \'T=\nM &0" \'T> \'P# &0  &L" \'-0 7P# &0# &L" \'-0 7P# &0$ &L" \'+) 60%\nM \'T> \'0. \'T7 &0& \'T, &0\' \'T- \'P# &0  &L" \'(M GP$ &0" !ED  !K\nM P!R+0)\\! !D @ 99 @ :P$ <[0!? 0 9 ( &60) &L" \'+& 7P$ &0% !ED\nM!0!K  !RQ@%D!@!]# !D"@!]#0!QQ@%Q+0)N9 !\\ P!D"P!K @!RY %D# !]\nM# !D#0!]#0!N20!\\ P!D#@!K @!R+0)\\ 0!D#P!K @!R"P)D!@!]# !D$ !]\nM#0!Q+0)\\ 0!D$0 @9!( :P( <BT"9!, ?0P 9 T ?0T <2T";@  ?   9!0 \nM:P( <O4"= X ?1H ? $ 9!4 (&06 &L" \')D G0. \'T0 &07 \'T, &08 \'T-\nM \'%L#GP! &09 &L" \'*. GP> \')_ G0. \'T6 \'\'R F0: \'T, &0\' \'T- \'%L\nM#GP! &0; &L" \'*^ F0< \'T, &0= \'T- \'0. \'T0 \'0) \'0* \'0, &8# \'T)\nM \'%L#GP! &0> "!D\'P!K @!R; YD( !]# !D(0!]#0!T#@!]$ !T"0!T"@!T\nM# !F P!]"0!Q; YN=PM\\  !D(@ @9", :P( <AH(9"0 9"4 ;!T ;1X ?1\\ \nM;1\\ ?2   7P> \'(J W0. \'T7 &X& &0& \'T, \'0. \'T1 &0F \'T3 \'P" &0"\nM !ED(@!K! !R;@-T" !]$0!D) !D  !L( !](0!D)P!\\(0!J(0!D*  \\;@  \nM? ( 9 ( &60I &L$ \'.> WP" &0" !ED*0!K @!R0 5\\ @!D!0 99 4 :P4 \nM<D %9"0 9   ;"( ?2( ?"( :B, ? ( 9"D &60J (," \'TC \'0D \'PC (,!\nM &0% &L% \') !7D7 \'PB &HE \'PC &0" !F# 0!]) !7;@T  0$!9"L ?20 \nM;@$ 6\'PD &0L &L  \'(2!\'PD &0M !1]) !N  !\\) !D+@!K!0!R/05\\\'@!R\nM2 1\\ 0!D+P!K @!S/P1\\\'P!D, !\\ 0"# @!R2 1D*0!]\' !N&0!D) !D  !L\nM( !](0!D)P!\\(0!J(0!D*  \\?"0 9!, :P4 <CH%= X ?1  ?"0 9#$ :P4 \nM<J $?"0 9#( :P  <J $?!P 9"( :P0 <J $9"( ?1P ;BH ?"0 9#( :P0 \nM<LH$= @ ?1  ?"0 9 8 :P0 <LH$9 D ?1P <<H$;@  ?!X <C<%>4( 9"0 \nM9   ;"8 ?24 ?"4 :B< @P  9 ( &60I "!D,P!K @!R$05T#@!]& !T"@!T\nM"P!T# !F P!]"0!N  !7<30% 0$!= X ?1@ = H = L = P 9@, ?0D <30%\nM6\'$W!7$Z!7$]!7% !6X  \'P? &0T \'P! \'P@ (,# \')G!70. \'T: &0" \'TF\nM \'0. \'T8 &X  \'P? &0U \'P! \'P@ (,# \'.+!7P? &0V \'P! \'P@ (,# \'( \nM!W0. \'T( \'P? &0W \'P! (," \'.O!7P? &0X \'P! (," \'(\\!G0( \'T( \'0.\nM \'T: \'P? &0Y \'P! (," \'+9!60Z \'T, &0- \'T- \'%7!GP? &0[ \'P! (,"\nM \'+Q!608 \'T, \'%7!GP? &0\\ \'P! (," \'()!G0. \'T2 \'%7!GP? &0P \'P!\nM (," \'(A!G0. \'T2 \'%7!GP? &0] \'P! (," \')7!G0. \'T0 \'%7!FX; \'P?\nM &0V \'P! \'P@ (,# \')7!F0I \'T< &X  \'P? &0^ \'P! \'P@ (,# \'( !W0.\nM \'T( \'0. \'T0 \'P? &0_ \'P! \'P@ (,# \'*0!G0( \'T( \'\']!GP? &1  \'P!\nM (," \'+]!G0( \'T0 \'P? &1! \'P! (," \'+0!F0D &0  &P@ \'TA &1" \'PA\nM &HA &1# #QN  !\\\'P!D1 !\\ 0"# @!R^@9T#@!]& !T"@!T"P!T# !F P!]\nM"0!Q^@9Q_09Q  =N  !\\\'P!D10!\\ 0"# @!RU0=\\\'P!D1@!\\ 0"# @!R6@=\\\nM\'@!R2P=\\ @!DO@!K @!R5P=T#@!]& !T"P!T"@!T# !F P!]"0!Q5P=QT@=D\nM.@!]# !D#0!]#0!Q%PA\\\'P!D2 !\\ 0"# @!S> =\\\'P!D20!\\ 0"# @!RM =\\\nM\'@!RT@=\\& !T#@!K P!RGP=T"0!T"P!T"@!T# !F! !]"0!QL0=T"P!T"@!T\nM# !F P!]"0!QT@=Q%PA\\\'P!D2@!\\ 0"# @!R%PAD!@!]# !D#0!]#0!Q%PAQ\nM; Y\\\'P!D2P!\\ 0"# @!R_ =T#@!]& !T"@!T"P!T# !F P!]"0!Q; Y\\ 0!D\nM3 !K @!R; YD!P!]# !D30!]#0!Q; YN4@9\\  !D3@!K @!R/@A\\ 0!D3P!K\nM @!R; YD4 !]# !Q; YN+@9\\  !D40!K @!R;@A\\ 0!D4@!K @!R; YT#@!]\nM$ !D.@!]# !D.@!]#0!Q; YN_@5\\  !D4P!K @!RLPA\\ 0!D5 !K @!RE0AD\nM!@!]# !D50!]#0!Q; Y\\ 0!D5@!K @!R; YD5P!]# !D6 !]#0!Q; YNN05\\\nM  !D60!K @!R.@E\\ 0!D6@!K @!RV@AD6P!]# !D7 !]#0!Q; Y\\ 0!D70!K\nM @!R^PAD7@!]# !D.@!]#0!T#@!]% !Q; Y\\ 0!D7P!K @!R\' ED!@!]# !D\nM#0!]#0!T#@!]$@!Q; Y\\ 0!D8 !K @!R; YT#@!]& !T"@!]"0!Q; YN,@5\\\nM  !D80!K @!RRPET#@!]&P!T"0!T"P!T"@!T# !F! !]"0!\\ 0!D8@!K @!R\nM>0ET#@!]" !T#@!]&@!Q; Y\\ 0!D8P!K @!RH ET"0!T"@!T"P!T# !F! !]\nM"0!T#@!]$0!Q; Y\\ 0!D9  @9&4 :P( <FP.= H = L = P 9@, ?0D = X \nM?18 <6P.;J$$?   9&8 :P( <O,)? $ 9!$ (&1G &L" \')L#G0. \'T7 \'%L\nM#FYY!\'P  &1H &L" \')7#7P> \',A"F0D &0  &P@ \'TA &1I \'PA &HA &1J\nM #QN  !\\ P!D: !K @!R.0MT#@!]%P!\\ 0!D:P @9&P :P( <X,*? $ 9&L \nM(&1M &L" \'.#"GP! &05 "!D;@!K @!S@PI\\ 0!D%0 @9&\\ :P( <X,*? $ \nM9!4 (&1P &L" \'*,"F0& \'T, \'\'^"WP! &1K "!D<0!K @!RI0IT#@!]&@!Q\nM_@M\\ 0!D:P @9\'( :P( <KX*9"( ?1P <?X+? $ 9\', (&1T &L" \'+="F1U\nM \'T, &1- \'T- \'\'^"WP! &1V "!D=P!K @!R\'0MD#0!]# !D!@!]#0!\\ 0!D\nM> !K @!R-@M\\ @!DOP!K @!R-@MD>@!]\'0!Q-@MQ_@M\\ 0!D=@ @9\'L :P( \nM<OX+= X ?1H <?X+;L4 ?!X <OX+? $ 9\'P (&1] &L" \'-O"WP! &1^ "!D\nM?P!K @!S;PM\\ 0!D?@ @9(  :P( <OX+= X ?1< ? $ 9&L (&1L &L" \'*.\nM"V0- \'T, \'\'["WP! &1K "!D@0!K @!RIPMT#@!]&@!Q^PM\\ 0!D:P @9(( \nM:P( <L8+9"( ?1P 9"$ ?0P <?L+? $ 9&L (&1R &L" \'+?"V0B \'T< \'\'[\nM"WP! &1S "!D= !K @!R^PMT#@!]$ !Q^PMQ_@MN  !\\ 0!D%0 @9&\\ :P( \nM<CP,9"0 9   ;"8 ?24 ?"4 :B< @P  9 ( &62# &L" \')4#60& \'T, \'%4\nM#7%L#GP! &2$ "!DA0!K @!R50QD(@!]\' !Q; Y\\ 0!D+0 @9(8 :P( <GH,\nM= X ?1< = X ?1L 9"D ?1P <6P.? $ 9\', (&2\' &L" \'+Y#\'P> \'+>#\'0+\nM \'0* \'0, &8# \'T) &0& \'T, &0Z \'T- \'P" &3  &L" \'+V#&0& \'T- \'P!\nM &2( &L" \'+2#\'0. \'T1 \'\';#\'0. \'T( \'\'V#\'%4#60& \'T, &2) \'T- \'0.\nM \'T0 \'0. \'T: \'%L#GP! &2* &L" \'(4#60& \'T, &0- \'T- \'%L#GP! &2+\nM &L" \'(O#60Z \'T, &2, \'T- \'%L#GP! &1\\ "!DC0!K @!R; YD!@!]# !D\nM6 !]#0!T#@!]%P!Q; YN%0%\\  !D"P!K @!RM0U\\ 0!D$0 @9(X :P( <H@-\nM9 8 ?0P 9%@ ?0T = X ?1( <6P.? $ 9(\\ :P( <IT-9 8 ?0P <6P.? $ \nM9)  :P( <FP.= X ?1( <6P.;K< ?   9)$ :P( <NX-= X ?1  ? $ 9)( \nM:P( <FP.= X ?1@ = H = L = P 9@, ?0D <6P.;GX ?   9), :P( <FP.\nM? $ 9)0 :P( <Q(.? $ 9)4 :P( <FP.? ( 9,$ :P( <B<.9)< ?0P ;C, \nM9"0 9   ;"@ ?2< ?"< :BD @P  9)@ :P( <E0.= X ?0@ 9)D ?0P ;@8 \nM9)H ?0P 9)L ?0T = X ?1< <6P.;@  ? $ 9 @ (&2< &L" \')V#W0+ \'0*\nM \'0, &8# \'T) &0D &0  &P@ \'TA \'PA &HA &HJ &2= (,! \'.R#G0. \'T;\nM &X  \'P" &0" !ED!0!K @!R=@]\\ @!D!0 99 D :P  <G8/9"< ?"$ :B$ \nM9)X /\'P! &2? &L" \')S#WP  &2@ &L" \'),#WP" &0" !ED!0!K @!R;0]\\\nM @!D!0 99"D :P( <FT/? ( 9"D &61S "!DH0!K @!R;0]T"@!T"P!T# !F\nM P!]"0!D!@!]# !T#@!]&@!Q;0]Q< ]\\  !DH@!K @!R< ]D.@!]# !DHP!]\nM#0!T#@!]" !Q< ]Q<P]Q=@]N  !\\ P!D  !K @!R"!!\\ 0!D"  @9)P :P, \nM<@@0? $ 9*0 :P( <J</= X ?1L ;@  9 8 ?0X 9 < ?0\\ ? 0 9 ( &60 \nM &L# \'(4$\'P$ &0" !ED" !K 0!S\\P]\\! !D @ 99 D :P( <@40? 0 9 4 \nM&60% &L  \'(%$&0& \'T, &0* \'T- \'$%$\'$4$&X, \'P, \'T. \'P- \'T/ \'0.\nM \'P2 &L" \'(I$\'0. \'T4 &X  \'0. \'P( &L" \'(^$\'0K \'T* &X  \'0  &H/\nM &HL &2E \'P% &2F \'P& &2G \'P\' &2H \'P( &2I \'P) &2J \'P* &2K \'P+\nM &2L \'P; &2M \'P< &2N \'P, &2O \'P- &2P \'P. &2Q \'P/ &2R \'P0 &2S\nM \'P1 &2T \'P2 &2U \'P4 &2V \'P3 &2W \'P5 &2X \'P6 &2Y \'P7 &2Z \'P8\nM &2[ \'P9 &2\\ \'P: &2] \'P= (, &0%D  !3*,(   !.9_%HXXBU^-0^:0  \nM  !S\'P   $AU;6UI;F=B:7)D($-O;6UU;FEC871I;VYS($QT9"YS$    $AU\nM;6UI;F=B:7)D($QT9"YI 0   &<       #P/V<       #H/VD$    :04 \nM  !G,0BL\'%ID^S]T P   %-\'26=[%*Y\'X7JD/V<       #X/W0#    24)-\nM<PX   !\'6%0X,# @5&5X=\'5R96>%ZU&X\'H4!0&D&    = 8   !\'6%0X,#!G\nMUZ-P/0K7ZS]T!@   #-$;&%B<VD-    <PT   !7:6QD8V%T(%90.#@P9P  \nM     !9 9YJ9F9F9F>D_= D   !015)-141)03-G,S,S,S,SPS]S"P   $=,\nM24Y4(%(S(%!49YJ9F9F9F:D_9V9F9F9F9N8_:1D   !S&0   $=,24Y4(%(S\nM(%!4("L@1TQ)3E0@1V%M;6%GFIF9F9F9V3]G         $!I P   \'0#    \nM051):?____\\H @   \'0&    <V5A<F-H= H   !)1TY/4D5#05-%9RU#\'.OB\nM-@H_= $    Q=!<   !!0D%155-?04Q,3U=?051)7U1204Y33&D"    = $ \nM   @9P          9S,S,S,S,],_:0H   !G<1L-X"V0ZC]S$    $%422!&\nM:7)E1TP@5C<S,#!T!0   %8S-# P9W$]"M>C<.T_9U@YM,AVON\\_= (    S\nM,G0(    34]"24Q)5%ET!    $9I<F5T!@   %)A9&5O;G,(    1FER94=,\nM(%9T P   %!R;W0%    5C<Q,#!G        X#]T!0   %8U,3 P= 4   !6\nM,S(P,\'0%    5CDX,#!T"    $UO8FEL:71Y<Q<   !-;V)I;&ET>2!2861E\nM;VX@2$0@,S8U,\',)    1FER94=,(%8U<PP   !&:7)E1TP@5C4R,#!T 0  \nM #!T%@   $%"05%54U])4U9)4U]42%)%4TA/3$1S#    $9I<F5\'3"!6-3<P\nM,\',(    1FER92!\'3"!S"@   $9I<F4@1TP@5#)S!0   "XT-3$Y<PH   !&\nM:7)E($=,(%@Q<PH   !&:7)E($=,(%HQ<PP   !&:7)E($=,(#AX,#!S"P  \nM $9I<F4@1TPX.# P<Q0   !204=%(#$R."!0<F\\@>#@V+U-316<S,S,S,S/[\nM/W0&    0V]M<&%Q<P@   !00EA\'1"U!1&<4KD?A>A3N/W,2    1&EA;6]N\nM9"!-=6QT:6UE9&EA<P@   !&:7)E($=,,70$    14Q307,-    14Q302!%\nM4D%:3U(@6&>:F9F9F9GA/W,/    14Q302!3>6YE<F=Y($E)9[@>A>M1N,X_\nM9S,S,S,S,^L_<Q<   !(97=L971T+5!A8VMA<F0@0V]M<&%N>70\'    :\'!V\nM:7-F>&<I7(_"]2CL/V?-S,S,S,SL/W0+    ;&EB,S5A8V1A,S!G"M>C<#T*\nM[S]T"P   &QI8F1D=FES>&=L<Q4   !6:7)T=6%L($UE;6]R>2!$<FEV97)T\nM!0   $EN=&5L<PL   !);G1E;" Y-#5\'37,-    26YT96P@0V%N=&EG86D@\nM    <R    !-;V)I;&4@26YT96PH4BD@-"!397)I97,@17AP<F5S<W,6    \nM26YT97)G<F%P:"!#;W)P;W)A=&EO;G0&    =V-G9\')V<Q(   !.5DE$24$@\nM0V]R<&]R871I;VYT @   $].=!(   !!0D%155-?1$E304),15]&4$5I#@  \nM \',.    475A9\')O($98(#0V,#!S#@   %%U861R;R!&6" Q-S P<PT   !1\nM=6%D<F\\@1E@@-3@P<PT   !1=6%D<F\\@1E@@-3<P<PT   !1=6%D<F\\@1E@@\nM,S<P<PX   !1=6%D<F\\@1E@@,S4P,\',.    475A9\')O($98(#@X,$UI"P  \nM \',+    475A9\')O,B!0<F]G        !$!I#P   \',/    475A9\')O-" W\nM-3 @6$=,<Q@   !1=6%D<F\\T(#<U,"!81TPO04=0+U-313)S#P   "XQ($Y6\nM241)02 U,RXS-FDU    = <   !\'949O<F-E:0D   !S"0   %%U861R;R!&\nM6&D\'    = <   !1=6%D<F\\T= <   !1=6%D<F\\R<PX   !1=6%D<F\\@1E@@\nM,3,P,\',.    475A9\')O($98(#4W,$UT!0   #8T8FET:0P   !S#    %%U\nM861R;R Q,# P37,*    475A9\')O($Y64W,+    475A9\')O,B!-6%)S%   \nM %%U861R;S(@35A2+T%\'4"]34T4R9V9F9F9F9OX_<Q,   !\'949O<F-E(#(U\nM-B]00TDO4U-%<PX   !2259!(%1.5" H4$-)*6>%ZU&X\'H7_/W,)    4DE6\nM02!43E0R= 8   !)35!!0U1S"0   %)%4R]3+S$O,G,*    5E!23R]"+S$R\nM.\',6    4W5N($UI8W)O<WES=&5M<RP@26YC+G,-    6%92+3$R,# L(%9)\nM4W,*    0G)I86X@4&%U;\',(    365S82!8,3%S#@   $UE<V$@3V9F4V-R\nM965N<PH   !-97-A(#,N-"XR9WL4KD?A>K0_= @   !L;G@X-E\\V-&<     \nM  !)0&<      (!!0&<       #R/W0$    365S870,    04)17U%!7U!2\nM24Y4= L   !-15-!7TY/7T%337,1    365S82!\'3%@@26YD:7)E8W1S\'   \nM $UE<V$@<\')O:F5C=#H@=W=W+FUE<V$S9"YO<F=S"P   "@Q+C4@365S82 V\nM<Q8   !602!,:6YU>"!3>7-T96US+"!);F,N9U*X\'H7K4?0_<PL   !\'1$D@\nM1V5N97)I8W0.    9W)A<&AI8W-$<FEV97)T#P   &1O=6)L94)U9F9E<FEN\nM9W0/    8F%C:V9A8V5#=6QL:6YG= P   !D:7-P;&%Y3&ES=\'-T$P   &AI\nM9VAL:6=H=$UE=&AO9$AI;G1T"    &1R86=-;V1E=!(   !A=71O1FET069T\nM97)2;W1A=&5T"0   &%N=&E!;&EA<W00    =\')A;G-L=6-E;F-Y36]D9705\nM    <&]L>6=O;D]F9G-E=$-O;G-T86YT=!(   !P;VQY9V]N3V9F<V5T4VQO\nM<&5T&@   \'!R:6YT4&]L>6=O;D]F9G-E=$-O;G-T86YT=!<   !P<FEN=%!O\nM;\'EG;VY/9F9S9713;&]P970,    =F5R=&5X07)R87ES=!H   !V97)T97A!\nM<G)A>7-);D1I<W!L87E,:7-T<W0.    =&5X=\'5R94UA<\'!I;F=T$P   \'!R\nM:6YT5&5X=\'5R94UA<\'!I;F=T\'    &-O;G1O=7)286YG951E>\'1U<F50<F5C\nM:7-I;VYT#P   &1I<F5C=%)E;F1E<FEN9W04    :&%R9\'=A<F5!8V-E;&5R\nM871I;VYT$P   &%C8V5L97)A=&5/9F938W)E96YT#P   &AA<F1W87)E3W9E\nM<FQA>70/    8F%C:V=R;W5N9$-O;&]R= P   !B86-K:6YG4W1O<F5T"0  \nM %]V:7-U86Q)9"@#    :0$   !I!0   \',%    +C0U,3DH P   &D!    \nM:00   !S#P   "XQ($Y6241)02 U,RXS-B@#    :0$   !I @   $XH P  \nM &D!    :0(   !S"@   $UE<V$@,RXT+C(H+0   \'0\'    <V5S<VEO;G0,\nM    9W)A<&AI8W-);F9O= @   !G;%9E;F1O<G0*    9VQ296YD97)E<G0)\nM    9VQ697)S:6]N= \\   !G;\'A397)V97)696YD;W)T$    &=L>%-E<G9E\nM<E9E<G-I;VYT!P   $]014Y?1TQ2\'P   \'00    2$%21%=!4D5?3U9%4DQ!\nM6700    4T]&5%=!4D5?3U9%4DQ!670#    6$]2= 4   !"3$5.1\'0%    \nM05-?25-T P   $]&1G06    9&5F875L=$=R87!H:6-S3W!T:6]N<U(S    \nM4C0   !2-0   %(V    4CP   !2/0   %(^    4C\\   !20    %)!    \nM4C$   !2,@   \'0$    3F]N970"    <F52!@   %(\'    = (   !O<W0\'\nM    96YV:7)O;G0&    <W1R:6YG= 4   !S<&QI=\'0#    ;&5N= 0   !A\nM=&]F= @   !P;&%T9F]R;70,    87)C:&ET96-T=7)E= ,   !U=&ET"P  \nM &=E=%!L871F;W)M= <   !H87-?:V5Y= 0   !&05-4= D   !S971686QU\nM97,H*    %)%    4D8   !21P   %)(    4DD   !2*@   %(K    4BP \nM  !2+0   %(N    4B\\   !2,    %(S    4C0   !2-0   %(V    4C< \nM  !2.    %(Y    4CL   !2.@   %(\\    4CT   !2/@   %(_    4D  \nM  !200   %(Q    4C(   !T"    \'9I<W5A;$ED=!,   !D:7-P;&%Y:6YG\nM5&]7:6YD;W=S4@8   !2!P   %)4    4E8   !T!@   &9I96QD<W0$    \nM9G9A;%):    =!T   !V:65W36%N:7!$:7-P;&%Y3&ES=%1H<F5S:&]L9%)<\nM    *      H     \'-B    0SI<4\')O9W)A;2!&:6QE<UQ$87-S875L="!3\nM>7-T96UE<UQ3:6UU;&%T:6]N4V5R=FEC97-<5C92,C Q.\'A<=VEN7V(V-%Q3\nM34%<<VET95QG<F%P:&EC<T-O;F9I9RYE;G9T%    &]N0V%E1W)A<&AI8W-3\nM=&%R=\'5P 0   \',& P    ,, 0P!# $, 0P2!@$& 08!!@$2 08!!@$, 0P!\nM# $, 08!!@$& 08!!@$, 0P!# $, 0P!# $, 0P!!@0&! P!# $, @8#!@,&\nM 08#\' (0 2 !!@$/ @P"!@$) @P"# $& 0D!$ $& 0P$# (& A "!@$& 0D"\nM# (& @D#!@$) @P"!@$& @8"$@(0 @8!!@(& A4"$ (6 P8!"0(& @8#!@,0\nM 08!# $0 3 !# $6 1(! P$7 0,!"@(, 0T!# 42 0\\!"0(, 0T%# $& B0!\nM"0$, 08!# $, 08! P$, 1H!!@$6 0,!!@$B Q(!!@$& @D"$@$2 @8"#P$/\nM! 8\'!@(/ 08!"0(/ 0D"#P$) P\\!"00/ 0P"$@$) A(!!@$&!!(!"04/ 08"\nM#P$, 1 "#P$& 1L"#P(/ @8"# ,& 14$!@$) @\\!#P(& 0P!%0(5 @\\!!@$,\nM @\\#!@$2"PP"!@$, @P"# (, @P"# (& 08!# (, @P!!@$) @P"!@$, @P"\nM# $& 0D"# (& 08$"0(, @8!!@,) @P$!@$, @P"!@02 @P"!@$) @P!$@$)\nM A !#P$, @P"$ (, @P"!@(, 1 "# 0& A !$ $0 1 !$ $) A ""0(0 @D"\nM$ $& 0D"$ $& 08"# $, 0P"$ (,!@8"$ $0 1 #!@(0 0D"$ $) A !!@$)\nM A ""0(0 0\\"$ (, 18!# (0 @D"$ (& 08!"0(0 @8"#P$& 08"# $&!0P!\nM"0(, P8!!@,&! D"# ,& 0D"# $& 0D"$ ,& 08## (, A !!@$& 0D## $)\nM PP!# (,!@8*# (& 14"# $, 0P!# $) PP!$@$& 0D"!@(& 0P#$ $/ 0P!\nM$@$) 2 %#0(, 0P#( $4 0\\!!@$, @P"!@$& 1(%# $0 @P!"0,& 08#$ $0\nM 2 !!@$/ P8!!@(, 0D%# $) @P!!@$& 08!!@$& 08!!@$& 08!!@$& 08!\n8!@$& 08!!@$& 08!!@$& 08!!@$& 08!\n \nend\n'),
    'onCaeStartup':driverUtils.decodeFunction('begin 666 -\nM8P          !P   $,   !S4    \'0  &H! &H" &0! \'0# (,  0%T  !J\nM! !D @ 9:@4 :@( 9 , = 8 @P ! 70  &H\' &H" &0$ \'0( &0% &0& &0\'\nM &0( (,  P%D  !3* D   !.= \\   !R96-O=F5R1V5O;65T<GES"P   %9I\nM97=P;W)T.B Q= 4   !T:71L970/    8F%C:V=R;W5N9%-T>6QE= \\   !B\nM86-K9W)O=6YD0V]L;W)S!P   "-&1D9&1D9T$    \'1R86YS;\'5C96YC>4UO\nM9&5I @   "@)    = <   !S97-S:6]N= X   !J;W5R;F%L3W!T:6]N<W0)\nM    <V5T5F%L=65S= H   !#3T]21$E.051%= D   !V:65W<&]R=\'-T&0  \nM \'9I97=P;W)T06YN;W1A=&EO;D]P=&EO;G-T P   $]&1G0/    9W)A<&AI\nM8W-/<\'1I;VYS= 4   !33TQ)1"@     *      H     \',6    1#I<06)Q\nM5T1<86)A<75S7W8V+F5N=G0,    ;VY#8653=&%R=\'5P)P   \',(      (3\n% AH!$@$ \n \nend\n'),
    'outputKeywords':ON,
    'parameterized':OFF,
    'partsAndAssemblies':ON,
    'parval':OFF,
    'plugin_central_dir':'C:\\SIMULIA\\CAE\\plugins\\2018',
    'postOutput':OFF,
    'preDecomposition':ON,
    'restart':OFF,
    'restartEndStep':OFF,
    'restartIncrement':0,
    'restartStep':0,
    'restartWrite':OFF,
    'resultsFormat':ODB,
    'rezone':OFF,
    'rigidbody':OFF,
    'runCalculator':OFF,
    'soils':OFF,
    'soliter':OFF,
    'solverTypes':['DIRECT'],
    'standard_parallel':ALL,
    'staticNonlinear':ON,
    'steadyStateTransport':OFF,
    'step':ON,
    'subGen':OFF,
    'subGenLibs':[],
    'subGenTypes':[],
    'submodel':OFF,
    'substrLibDefs':OFF,
    'substructure':OFF,
    'symmetricModelGeneration':OFF,
    'thermal':OFF,
    'tmpdir':'C:\\Users\\prafu\\AppData\\Local\\Temp',
    'tracer':OFF,
    'unsymm':OFF,
    'visco':OFF,
    'xplSelect':OFF,
}
analysis = StandardAnalysis(options)
status = analysis.run()
sys.exit(status)