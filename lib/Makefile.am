AUTOMAKE_OPTIONS=foreign

bin_PROGRAMS=main
src_dir=../src
main_SOURCES=$(src_dir)/main.cpp \
 	     $(src_dir)/mcmodel.cpp \
 	     $(src_dir)/utils.cpp \
	     $(src_dir)/TetHit.cpp \
             $(src_dir)/set_user_defined_data.cpp  $(src_dir)/PoissonSolver.cpp \
             $(src_dir)/ParticleUtils.cpp \
             $(src_dir)/ParSurfaceScatter.cpp \
             $(src_dir)/MultipleRefresh.cpp \
             $(src_dir)/HolePhononScatter.cpp \
             $(src_dir)/GetSurfRoughnessPhononScRate.cpp $(src_dir)/erf.cpp \
             $(src_dir)/ElectronPhononScatter.cpp \
             $(src_dir)/ElectronImpurityScatter.cpp \
	     $(src_dir)/Constant_GlobalVariable_MathFunction.cpp \
	     $(src_dir)/CellHit.cpp \
	     $(src_dir)/Band.cpp \
	     $(src_dir)/heat.cpp \
             $(src_dir)/myStatistic.cpp \
             $(src_dir)/fermi_int.f

INCLUDES=$(INCLUDE) -I../include

