#
# This Makefile can be used to automatically build the entire package.  
#
# If you make changes in the "Makefile" under any subdirectory, you can
# rebuild the system with "make clean" followed by "make all".
#
#
# Generic. On most systems, this should handle everything.
#
all:
	cd src; make libfnlsdp.a
	cd solver; make fnlsdp
	
#
# Clean out all of the directories.
# 
clean:
	cd src; make clean
	cd solver; make clean
	







