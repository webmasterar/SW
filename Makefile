ifdef SystemRoot
	RM = del /Q
	EXT = .exe
else
	ifeq ($(shell uname), Linux)
		RM = rm -f
		EXT =
	endif
endif

all:
	g++ -std=c++11 -D_USE_64 -msse3 -O3 -fomit-frame-pointer -funroll-loops -DNDEBUG -lm -I . sw.cpp -o sw$(EXT)

clean:
	$(RM) sw.o sw
