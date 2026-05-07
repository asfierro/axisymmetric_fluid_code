export TRILINOS_LIB=/dir/to/Trilinos_16_hopper/lib64
export TRILINOS_INC=/dir/to/Trilinos_16_hopper/include
export KOKKOS_INC=/dir/to/Trilinos_16_hopper/include/kokkos
export FFTW3_INC=/dir/to/fftw3/include
export FFTW3_LIB=/dir/to/fftw3/lib

mpic++ -Wno-deprecated-declarations \
		main.cpp field_solver.cpp utilities.cpp species.cpp \
		reaction.cpp transport_data.cpp photon_solver.cpp input_deck.cpp \
		photoionization.cpp \
		-I$TRILINOS_INC -L$TRILINOS_LIB \
		-I$FFTW3_INC -L$FFTW3_LIB \
		-I$NETLIB_LAPACK_INC -L$NETLIB_LAPACK_LIB \
		-I$MATIO_INC -L$MATIO_LIB \
		-I$KOKKOS_INC \
		-lamesos2  -lbelos -lbelostpetra \
		-lbelosxpetra -lgaleri-xpetra -lifpack2  \
		-lkokkoscontainers \
		-lkokkoscore -lkokkoskernels \
		-lml -lmuelu -lmuelu-adapters \
		-lteuchoscomm \
		-lteuchoscore -lteuchoskokkoscomm \
		-lteuchoskokkoscompat -lteuchosnumerics \
		-lteuchosparameterlist -lteuchosparser \
		-lteuchosremainder -ltpetra \
		-ltpetraclassic \
		-ltpetraext \
		-ltpetrainout -ltrilinosss \
		-lxpetra \
		-lfftw3_mpi -lblas -llapack \
		-ldl -lfftw3 \
		-fopenmp -ffast-math -o main
