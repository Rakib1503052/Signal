#------------------------------------------------------------------------------#
# This makefile was generated by 'cbp2make' tool rev.147                       #
#------------------------------------------------------------------------------#


WORKDIR = %cd%

CC = gcc.exe
CXX = g++.exe
AR = ar.exe
LD = g++.exe
WINDRES = windres.exe

INC = 
CFLAGS = -O3
RESINC = 
LIBDIR = 
LIB = libs\\libMatrix_alg.a
LDFLAGS = -static-libstdc++ -static-libgcc -static

INC_DEBUG = $(INC)
CFLAGS_DEBUG = $(CFLAGS) -Wall -g
RESINC_DEBUG = $(RESINC)
RCFLAGS_DEBUG = $(RCFLAGS)
LIBDIR_DEBUG = $(LIBDIR)
LIB_DEBUG = $(LIB) 
LDFLAGS_DEBUG = $(LDFLAGS)
OBJDIR_DEBUG = obj\\Debug
DEP_DEBUG = 
OUT_DEBUG = bin\\Debug\\libSignal.a

INC_RELEASE = $(INC)
CFLAGS_RELEASE = $(CFLAGS) -O3 -Wall
RESINC_RELEASE = $(RESINC)
RCFLAGS_RELEASE = $(RCFLAGS)
LIBDIR_RELEASE = $(LIBDIR)
LIB_RELEASE = $(LIB) 
LDFLAGS_RELEASE = $(LDFLAGS) -O3
OBJDIR_RELEASE = obj\\Release
DEP_RELEASE = 
OUT_RELEASE = bin\\Release\\libSignal.a

OBJ_DEBUG = $(OBJDIR_DEBUG)\\fft.o $(OBJDIR_DEBUG)\\filter.o $(OBJDIR_DEBUG)\\getWindow.o $(OBJDIR_DEBUG)\\padding.o $(OBJDIR_DEBUG)\\welch.o $(OBJDIR_DEBUG)\\windowF.o

OBJ_RELEASE = $(OBJDIR_RELEASE)\\fft.o $(OBJDIR_RELEASE)\\filter.o $(OBJDIR_RELEASE)\\getWindow.o $(OBJDIR_RELEASE)\\padding.o $(OBJDIR_RELEASE)\\welch.o $(OBJDIR_RELEASE)\\windowF.o

all: debug release

clean: clean_debug clean_release

before_debug: 
	cmd /c if not exist bin\\Debug md bin\\Debug
	cmd /c if not exist $(OBJDIR_DEBUG) md $(OBJDIR_DEBUG)

after_debug: 

debug: before_debug out_debug after_debug

out_debug: before_debug $(OBJ_DEBUG) $(DEP_DEBUG)
	$(AR) rcs $(OUT_DEBUG) $(OBJ_DEBUG)

$(OBJDIR_DEBUG)\\fft.o: fft.cpp
	$(CXX) $(CFLAGS_DEBUG) $(INC_DEBUG) -c fft.cpp -o $(OBJDIR_DEBUG)\\fft.o

$(OBJDIR_DEBUG)\\filter.o: filter.cpp
	$(CXX) $(CFLAGS_DEBUG) $(INC_DEBUG) -c filter.cpp -o $(OBJDIR_DEBUG)\\filter.o

$(OBJDIR_DEBUG)\\getWindow.o: getWindow.cpp
	$(CXX) $(CFLAGS_DEBUG) $(INC_DEBUG) -c getWindow.cpp -o $(OBJDIR_DEBUG)\\getWindow.o

$(OBJDIR_DEBUG)\\padding.o: padding.cpp
	$(CXX) $(CFLAGS_DEBUG) $(INC_DEBUG) -c padding.cpp -o $(OBJDIR_DEBUG)\\padding.o

$(OBJDIR_DEBUG)\\welch.o: welch.cpp
	$(CXX) $(CFLAGS_DEBUG) $(INC_DEBUG) -c welch.cpp -o $(OBJDIR_DEBUG)\\welch.o

$(OBJDIR_DEBUG)\\windowF.o: windowF.cpp
	$(CXX) $(CFLAGS_DEBUG) $(INC_DEBUG) -c windowF.cpp -o $(OBJDIR_DEBUG)\\windowF.o

clean_debug: 
	cmd /c del /f $(OBJ_DEBUG) $(OUT_DEBUG)
	cmd /c rd bin\\Debug
	cmd /c rd $(OBJDIR_DEBUG)

before_release: 
	cmd /c if not exist bin\\Release md bin\\Release
	cmd /c if not exist $(OBJDIR_RELEASE) md $(OBJDIR_RELEASE)

after_release: 

release: before_release out_release after_release

out_release: before_release $(OBJ_RELEASE) $(DEP_RELEASE)
	$(AR) rcs $(OUT_RELEASE) $(OBJ_RELEASE)

$(OBJDIR_RELEASE)\\fft.o: fft.cpp
	$(CXX) $(CFLAGS_RELEASE) $(INC_RELEASE) -c fft.cpp -o $(OBJDIR_RELEASE)\\fft.o

$(OBJDIR_RELEASE)\\filter.o: filter.cpp
	$(CXX) $(CFLAGS_RELEASE) $(INC_RELEASE) -c filter.cpp -o $(OBJDIR_RELEASE)\\filter.o

$(OBJDIR_RELEASE)\\getWindow.o: getWindow.cpp
	$(CXX) $(CFLAGS_RELEASE) $(INC_RELEASE) -c getWindow.cpp -o $(OBJDIR_RELEASE)\\getWindow.o

$(OBJDIR_RELEASE)\\padding.o: padding.cpp
	$(CXX) $(CFLAGS_RELEASE) $(INC_RELEASE) -c padding.cpp -o $(OBJDIR_RELEASE)\\padding.o

$(OBJDIR_RELEASE)\\welch.o: welch.cpp
	$(CXX) $(CFLAGS_RELEASE) $(INC_RELEASE) -c welch.cpp -o $(OBJDIR_RELEASE)\\welch.o

$(OBJDIR_RELEASE)\\windowF.o: windowF.cpp
	$(CXX) $(CFLAGS_RELEASE) $(INC_RELEASE) -c windowF.cpp -o $(OBJDIR_RELEASE)\\windowF.o

clean_release: 
	cmd /c del /f $(OBJ_RELEASE) $(OUT_RELEASE)
	cmd /c rd bin\\Release
	cmd /c rd $(OBJDIR_RELEASE)

.PHONY: before_debug after_debug clean_debug before_release after_release clean_release

