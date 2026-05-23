SUBDIRS1 := $(wildcard Mod/*)
SUBDIRS2 := $(wildcard Lib/)

.PHONY : lib $(SUBDIRS2)
lib : $(SUBDIRS2)

$(SUBDIRS2) :
	$(MAKE) build -C $@
	
.PHONY : mods $(SUBDIRS1) 
mods : $(SUBDIRS1) 
	
$(SUBDIRS1) :
	$(MAKE) build -C $@


