#include <iostream>
#include <string>

using namespace std;

class CustomException : public std::exception {
    private:
    string message;

    public:
    CustomException(string msg) : message(msg) {}
    string what () {
        return message;
    }
};


int main(int argc, char** argv)
{
	try {
	    cout << "Hello, World!" << endl;
	    if(argc!=6) {
	        throw CustomException("Wrong number of arguments provided");
	    }
	} catch (CustomException ce) {
        cout << ce.what();
        exit(1);
    }
}