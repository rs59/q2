all:
	mkdir bin
	cd bin && g++ -std=c++17 ../src/2023partition.cpp -o 2023partition
	cd bin && g++ -std=c++17 ../src/test.cpp -lcppunit -o test

test:
	cd bin && chmod +x 2023partition
	cd bin && ./test

clean:
	cd bin && $(RM) 2023partition test
	$(RM) ./bin