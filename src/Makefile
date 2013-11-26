OBJS = main.o Autotune.o SpExMCMC.o SpExBranchEvent.o Node.o Tree.o Settings.o SettingsParameter.o TraitModel.o TraitMCMC.o MbRandom.o TraitBranchEvent.o SpExModel.o Log.o
CC = g++
DEBUG = -g
GIT_COMMIT_ID = $(shell git rev-parse HEAD)
CFLAGS = -Wall -Wextra -c -m64 -O3 -funroll-loops $(DEBUG)
LFLAGS = -Wall -Wextra -m64 $(DEBUG)

bamm: $(OBJS)
	$(CC) $(LFLAGS) $(OBJS) -o bamm

main.o: main.cpp
	$(CC) $(CFLAGS) -DGIT_COMMIT_ID=\"$(GIT_COMMIT_ID)\" main.cpp

Autotune.o: Autotune.cpp MbRandom.h SpExModel.h TraitModel.h Settings.h
	$(CC) $(CFLAGS) Autotune.cpp

Log.o: Log.cpp Log.h
	$(CC) $(CFLAGS) Log.cpp
	
SpExMCMC.o: SpExMCMC.cpp SpExMCMC.h
	$(CC) $(CFLAGS) SpExMCMC.cpp

SpExBranchEvent.o: SpExBranchEvent.cpp SpExBranchEvent.h
	$(CC) $(CFLAGS) SpExBranchEvent.cpp

Node.o: Node.cpp Node.h
	$(CC) $(CFLAGS) Node.cpp

Tree.o: Tree.cpp Tree.h
	$(CC) $(CFLAGS) Tree.cpp

Settings.o: Settings.cpp Settings.h
	$(CC) $(CFLAGS) Settings.cpp

SettingsParameter.o: SettingsParameter.cpp SettingsParameter.h
	$(CC) $(CFLAGS) SettingsParameter.cpp

TraitModel.o: TraitModel.cpp TraitModel.h TraitBranchEvent.h
	$(CC) $(CFLAGS) TraitModel.cpp

TraitMCMC.o: TraitMCMC.cpp TraitMCMC.h
	$(CC) $(CFLAGS) TraitMCMC.cpp

MbRandom.o: MbRandom.cpp MbRandom.h
	$(CC) $(CFLAGS) MbRandom.cpp

TraitBranchEvent.o: TraitBranchEvent.cpp TraitBranchEvent.h
	$(CC) $(CFLAGS) TraitBranchEvent.cpp

SpExModel.o: SpExModel.cpp SpExModel.h SpExBranchEvent.h
	$(CC) $(CFLAGS) SpExModel.cpp

clean:
	rm -rf *.o bamm
