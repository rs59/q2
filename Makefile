all:
	mkdir bin
	cd bin && g++ -std=c++17 ../src/sequential.cpp -o sequential
	cd bin && g++ -std=c++17 ../src/test.cpp -lcppunit -o test

test:
	cd bin && chmod +x sequential
	cd bin && ./test

clean:
	cd bin && $(RM) sequential test
	$(RM) ./bin