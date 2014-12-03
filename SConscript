import os
import commands
import sys

Import('MPI')
if MPI: Import('mpif90_exec')

lib_prefix = "/usr/local"
MN_version_dir = 'MultiNest_v3.7/'




# # Configure the f90 compiler
FFLAGS = '-O3 -w -Wno-unused-parameter -ffree-line-length-none'
if MPI: FFLAGS += ' -DMPI -lmpi'

env = Environment()
env = env.Clone(F90=mpif90_exec if MPI else 'gfortran',
	            LINK=mpif90_exec if MPI else 'gfortran',
	            # LINKFLAGS='-g',
	            F90COMSTR='Compiling $SOURCE',
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
def print_cmd_line(s, target, src, env):
    """s is the original command line, target and src are lists of target
    and source nodes respectively, and env is the environment."""
    sys.stdout.write(" Making %s...\n"% (' and '.join([str(x) for x in target])))

env['PRINT_CMD_LINE_FUNC'] = print_cmd_line

source_files = Glob(MN_version_dir + '*.f90') + Glob(MN_version_dir + '*.F90')
source_files = [f for f in source_files if 'cwrapper.f90' not in str(f)]

# explicitly compile the source files to object+mod files
comps = env.SharedObject(source_files)

# only object files without .mod
objs = [obj for obj in comps if obj.get_suffix() in (".o", ".os")]

# create static library
static_lib = env.Library(target=MN_version_dir+'nest3', source=objs)
shared_lib = env.SharedLibrary(target=MN_version_dir+'nest3', source=objs, VariantDir=MN_version_dir)

env.Install(os.path.join(lib_prefix, "lib"), shared_lib)
env.Alias('install', [os.path.join(lib_prefix, "lib")])