import os
import commands

MN_version_dir = 'MultiNest_v3.7/'


# compiling with MPI support?
out = commands.getstatusoutput('mpirun --version')
if out[0] == 0:
	MPI = True
else:
	MPI = False
	print 'Compiling MultiNest without MPI'


## get path to mpif90 executable
if MPI:
	out = commands.getstatusoutput('which mpif90')
	if out[0] == 0:
		mpif90_exec = out[1]
		mpif90_version = commands.getoutput(mpif90_exec + ' -v').split('\n')[0]
		print 'Found %s at %s' % (mpif90_version, mpif90_exec)
	else:
		raise RuntimeError('Compiling with MPI but mpif90 does not seem to be installed. Aborting!')

# # Configure the f90 compiler
FFLAGS = '-O3 -w -Wno-unused-parameter -ffree-line-length-none'
if MPI: FFLAGS += ' -DMPI -lmpi'
env = DefaultEnvironment()
env = env.Clone(tools=['gfortran'],
	            F90=mpif90_exec if MPI else 'gfortran',
	            LINK=mpif90_exec if MPI else 'gfortran',
	            # LINKFLAGS='-g',
	            F90FLAGS=FFLAGS,
	            FORTRANMODDIRPREFIX = '-J',  # option used by gfortran to specify where to put .mod files for compiled modules
	            FORTRANMODDIR = MN_version_dir
	            )

if GetOption('gfortran') is None:
	pass
else:
	if MPI:
		env['F90FLAGS'] += ' -f90=' + GetOption('gfortran')
		# print env['F90FLAGS']
	else:
		env['F90'] = GetOption('gfortran')
		# print env['F90']

# if ARGUMENTS.get('VERBOSE') != '1':
#     env['F90COMSTR'] = "Compiling $TARGET"
#     env['LINKCOMSTR'] = "Linking $TARGET"

source_files = Glob(MN_version_dir + '*.f90') + Glob(MN_version_dir + '*.F90')
source_files = [f for f in source_files if 'cwrapper.f90' not in str(f)]

# explicitly compile the source files to object+mod files
comps = env.SharedObject(source_files)

# only object files without .mod
objs = [obj for obj in comps if obj.get_suffix() in (".o", ".os")]

# create static library
static_lib = env.Library(target=MN_version_dir+'nest3', source=objs)
shared_lib = env.SharedLibrary(target=MN_version_dir+'nest3', source=objs, VariantDir=MN_version_dir)