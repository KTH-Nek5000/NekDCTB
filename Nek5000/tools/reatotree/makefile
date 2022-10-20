# general makefile to compile 2D and 3D nek-p4est converters
.PHONY: all lib rea2Dtree rea3Dtree clean_lib clean

all: lib rea2Dtree rea3Dtree

lib:
	@cd ./p4est/; env CC="$(CC)" CXX="$(CXX)" FC="$(FC)" ./install

rea2Dtree:
	$(eval NDIM:=2)
	$(MAKE) --file=makefile_reatotree clean
	$(MAKE) --file=makefile_reatotree all

rea3Dtree:
	$(eval NDIM:=3)
	$(MAKE) --file=makefile_reatotree clean
	$(MAKE) --file=makefile_reatotree all

clean_lib:
	@cd ./p4est/; ./install clean

clean: clean_lib
	$(MAKE) --file=makefile_reatotree clean
	@rm -f $(bin_nek_tools)/rea2Dtree $(bin_nek_tools)/rea3Dtree
