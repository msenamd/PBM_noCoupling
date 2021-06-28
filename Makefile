PBM_solver: PBM.cpp \
			./particle_1D/particle_1D.h \
			./particle_1D/TDMA.h \
			./particle_1D/calcConvection.h \
			./particle_1D/VegReaction.C \
			./particle_1D/particle_1D.C \
			./particle_1D/TDMA.C \
			./particle_1D/calcConvection.C \
			./particle_1D/VegReaction.C
			
		g++ -omp -o PBM_solver \
					PBM.cpp \
					./particle_1D/particle_1D.C \
					./particle_1D/TDMA.C \
					./particle_1D/calcConvection.C \
					./particle_1D/VegReaction.C
					
					

clean:

	rm *.o PBM_solver