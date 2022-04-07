OBJDIR = $(GARFIELD_HOME)/Object
SRCDIR = $(GARFIELD_HOME)/Source
INCDIR = $(GARFIELD_HOME)/Include
HEEDDIR = $(GARFIELD_HOME)/Heed
LIBDIR = $(GARFIELD_HOME)/Library
TARNAME = Cell_v00

# Compiler flags
CFLAGS = -Wall -Wextra -Wno-long-long \
	`root-config --cflags` \
	-O3 -fno-common -c \
	-I$(INCDIR) -I$(HEEDDIR)

# Debug flags
#CFLAGS += -g

LDFLAGS = -L$(LIBDIR) -lGarfield
LDFLAGS += `root-config --glibs` `root-config --ldflags` -lGeom \
	-lgfortran -lm

#LDFLAGS += -g

target: $(TARNAME).o AvalancheMC_My.o Mat_Vd_const.o
	$(CXX) $(TARNAME).o AvalancheMC_My.o Mat_Vd_const.o -o $(TARNAME) $(LDFLAGS)
	rm -f $(TARNAME).o AvalancheMC_My.o
	rm -f $(TARNAME).o Mat_Vd_const.o

$(TARNAME).o: $(TARNAME).cc AvalancheMC_My.hh Mat_Vd_const.hh
	$(CXX) $(CFLAGS) $(TARNAME).cc

Mat_Vd_const.o: Mat_Vd_const.cc Mat_Vd_const.hh $(INCDIR)/Medium.hh
	$(CXX) $(CFLAGS) Mat_Vd_const.cc

AvalancheMC_My.o: AvalancheMC_My.cc AvalancheMC_My.hh
	$(CXX) $(CFLAGS) AvalancheMC_My.cc

clean:
	rm $(TARNAME).o AvalancheMC_My.o $(TARNAME)
	rm $(TARNAME).o Mat_Vd_const.o $(TARNAME)
