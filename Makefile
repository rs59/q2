all:
	cd src && g++ -std=c++17 sequential.cpp -o sequential
	cd src && g++ -std=c++17 test.cpp -lcppunit -o test

test:
	cd src && chmod +x sequential
	cd src && ./test

clean:
	cd src && $(RM) sequential test