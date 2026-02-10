NVCC = nvcc
CXXFLAGS = -O3 -lm
NVCCFLAGS = -gencode arch=compute_35,code=sm_35 \
            -gencode arch=compute_37,code=sm_37 \
            -gencode arch=compute_61,code=sm_61 \
            -Dnumber_of_myocyte=$(NUMBER_OF_MYOCYTE) \
            -Dnumber_of_1D_cables=$(NUMBER_OF_1D_CABLES)
LIBS = -std=c++14 -D_GLIBCXX_USE_CXX11_ABI=0 -lcusolver -lcuda -lcublas 
NUMBER_OF_MYOCYTE ?= 200
NUMBER_OF_1D_CABLES ?= 1
EXECUTABLE = X50_50%Na90%K_m$(NUMBER_OF_MYOCYTE)_c$(NUMBER_OF_1D_CABLES).run
.PHONY: all clean
all: $(EXECUTABLE)
$(EXECUTABLE): main.cu
	$(NVCC) $(NVCCFLAGS) $(CXXFLAGS) $(LIBS) $< -o $@ 
clean:
	rm -f $(EXECUTABLE)
