ARCH=$(shell uname)
ifeq ($(ARCH),Darwin)
	include Makefile.osx
else ifeq ($(ARCH),Linux)
	include Makefile.linux
else ifeq ($(ARCH),CYGWIN_NT-6.1-WOW64)
	include Makefile.cygwin
else
	$(warning Unsupported architecture. Trying Makefile.linux)
	include Makefile.linux
endif
