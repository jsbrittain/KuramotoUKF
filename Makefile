#SUBDIRS = src tests
all:
	make -C ./src
	make -C ./tests    
install:
	make -C ./src install
uninstall:
	make -C ./src uninstall
clean:
	make -C ./src clean
test:
	make -C ./tests
	make -C ./tests run
