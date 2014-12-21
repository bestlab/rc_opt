BINDIR = $(HOME)/bin
TARGETS = rc_opt qcalc qcalc_gen assign_tp assign_tp2 constN_proj \
	  pTP_crd qcalc_pbc xtr_pvec ave_pvec dcdinfo trunctraj crossval \
	  extract_tp drms ave_contacts qcontig qphi drms_pbc \
	  ferguson dcd2dcd dij qcalc_discrete quick_epert distance vector fret
SCRIPTS = qlist2pvec.py
OBJECTS = random_gen.o TrajFile.o swapbytes.o
CFLAGS = -O3 

all: $(TARGETS)

assign_tp: assign_tp.cc
	$(CXX) $(CFLAGS) -o assign_tp assign_tp.cc

extract_tp: extract_tp.cc TrajFile.o swapbytes.o
	$(CXX) $(CFLAGS) -o extract_tp extract_tp.cc TrajFile.o swapbytes.o

drms: drms.cc TrajFile.o swapbytes.o
	$(CXX) $(CFLAGS) -o drms drms.cc TrajFile.o swapbytes.o

dij: dij.cc TrajFile.o swapbytes.o
	$(CXX) $(CFLAGS) -o dij dij.cc TrajFile.o swapbytes.o

drms_pbc: drms_pbc.cc TrajFile.o swapbytes.o
	$(CXX) $(CFLAGS) -o drms_pbc drms_pbc.cc TrajFile.o swapbytes.o

assign_tp2: assign_tp2.cc
	$(CXX) $(CFLAGS) -o assign_tp2 assign_tp2.cc

rc_opt: random_gen.o common.o traj.o swapbytes.o rc_opt.cc config.h
	$(CXX) $(CFLAGS) -o rc_opt rc_opt.cc random_gen.o traj.o swapbytes.o common.o

ferguson: random_gen.o TrajFile.o swapbytes.o ferguson.cc config.h
	$(CXX) $(CFLAGS) -o ferguson ferguson.cc random_gen.o TrajFile.o swapbytes.o

ave_contacts: TrajFile.o swapbytes.o ave_contacts.cc config.h
	$(CXX) $(CFLAGS) -o ave_contacts ave_contacts.cc TrajFile.o swapbytes.o

ncdiff_contacts: TrajFile.o swapbytes.o ncdiff_contacts.cc config.h
	$(CXX) $(CFLAGS) -o ncdiff_contacts ncdiff_contacts.cc TrajFile.o swapbytes.o

ave_cmap: TrajFile.o swapbytes.o ave_cmap.cc config.h
	$(CXX) $(CFLAGS) -o ave_cmap ave_cmap.cc TrajFile.o swapbytes.o

crossval: crossval.cc config.h
	$(CXX) $(CFLAGS) -o crossval crossval.cc 

pTP_crd: pTP_crd.cc
	$(CXX) $(CFLAGS) -o pTP_crd pTP_crd.cc 

constN_proj: TrajFile.o swapbytes.o constN_proj.cc
	$(CXX) $(CFLAGS) -o constN_proj constN_proj.cc TrajFile.o swapbytes.o

qphi: qphi.cc TrajFile.o swapbytes.o
	$(CXX) $(CFLAGS) -o qphi TrajFile.o swapbytes.o qphi.cc

distance: distance.cc common.o traj.o swapbytes.o
	$(CXX) $(CFLAGS) -o distance common.o traj.o swapbytes.o distance.cc

fret: fret.cc common.o traj.o swapbytes.o
	$(CXX) $(CFLAGS) -o fret common.o traj.o swapbytes.o fret.cc

vector: vector.cc common.o traj.o swapbytes.o
	$(CXX) $(CFLAGS) -o vector common.o traj.o swapbytes.o vector.cc

qcalc_gen: qcalc_gen.cc TrajFile.o swapbytes.o
	$(CXX) $(CFLAGS) -o qcalc_gen TrajFile.o swapbytes.o qcalc_gen.cc

#go_correct: go_correct.cc TrajFile.o swapbytes.o
#	$(CXX) $(CFLAGS) -o go_correct TrajFile.o swapbytes.o go_correct.cc

qcalc_discrete: qcalc_discrete.cc TrajFile.o swapbytes.o
	$(CXX) $(CFLAGS) -o qcalc_discrete TrajFile.o swapbytes.o qcalc_discrete.cc

quick_epert: quick_epert.cc TrajFile.o swapbytes.o
	$(CXX) $(CFLAGS) -o quick_epert TrajFile.o swapbytes.o quick_epert.cc

qcontig: qcontig.cc TrajFile.o
	$(CXX) $(CFLAGS) -o qcontig TrajFile.o swapbytes.o qcontig.cc

lqcalc: lqcalc.cc TrajFile.o swapbytes.o
	$(CXX) $(CFLAGS) -o lqcalc TrajFile.o swapbytes.o lqcalc.cc

ldiff: ldiff.cc TrajFile.o swapbytes.o
	$(CXX) $(CFLAGS) -o ldiff TrajFile.o swapbytes.o ldiff.cc

qcalc: qcalc.cc TrajFile.o swapbytes.o
	$(CXX) $(CFLAGS) -o qcalc TrajFile.o swapbytes.o qcalc.cc

qcalc_debug: qcalc_debug.cc TrajFile.o swapbytes.o
	$(CXX) $(CFLAGS) -o qcalc_debug TrajFile.o swapbytes.o qcalc_debug.cc

qcalc_cmap: qcalc_cmap.cc TrajFile.o swapbytes.o
	$(CXX) $(CFLAGS) -o qcalc_cmap TrajFile.o swapbytes.o qcalc_cmap.cc

kmt: kmt.cc TrajFile.o swapbytes.o
	$(CXX) $(CFLAGS) -o kmt TrajFile.o swapbytes.o kmt.cc

qcalc_pbc: qcalc_pbc.cc TrajFile.o swapbytes.o
	$(CXX) $(CFLAGS) -o qcalc_pbc TrajFile.o swapbytes.o qcalc_pbc.cc

xtr_pvec: xtr_pvec.cc config.h
	$(CXX) $(CFLAGS) -o xtr_pvec xtr_pvec.cc

ave_pvec: ave_pvec.cc config.h
	$(CXX) $(CFLAGS) -o ave_pvec ave_pvec.cc

dcdinfo: dcdinfo.cc config.h TrajFile.o swapbytes.o
	$(CXX) $(CFLAGS) -o dcdinfo dcdinfo.cc TrajFile.o swapbytes.o

trunctraj: trunctraj.cc config.h TrajFile.o swapbytes.o
	$(CXX) $(CFLAGS) -o trunctraj trunctraj.cc TrajFile.o swapbytes.o

random_gen.o: random_gen.cc random_gen.h
	$(CXX) $(CFLAGS) -c random_gen.cc

TrajFile.o: TrajFile.cc TrajFile.h
	$(CXX) $(CFLAGS) -c TrajFile.cc

traj.o: traj.cc traj.h
	$(CXX) $(CFLAGS) -c traj.cc

common.o: common.cc common.h
	$(CXX) $(CFLAGS) -c common.cc

dcd2dcd: TrajFile.o swapbytes.o dcd2dcd.cc
	$(CXX) $(CFLAGS) -o dcd2dcd dcd2dcd.cc TrajFile.o swapbytes.o

fixaaqaa: TrajFile.o swapbytes.o fixaaqaa.cc
	$(CXX) $(CFLAGS) -o fixaaqaa fixaaqaa.cc TrajFile.o swapbytes.o

swapbytes.o: swapbytes.cc swapbytes.h
	$(CXX) $(CFLAGS) -c swapbytes.cc

install: $(TARGETS) $(SCRIPTS)
	-cp $(TARGETS) $(SCRIPTS) $(BINDIR)
clean:
	-rm $(OBJECTS)

distclean: clean
	-rm $(TARGETS)
	-rm *.o
