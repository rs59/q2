#include <iostream>
#include <cstdio>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <array>

#include <cppunit/TestRunner.h>
#include <cppunit/TestResult.h>
#include <cppunit/TestResultCollector.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/BriefTestProgressListener.h>
#include <cppunit/extensions/TestFactoryRegistry.h>

// #include <2023partition.cpp>


// Custom execution function to get exit code, stderr, and stdout
// Returns value to test if equal to EXIT_SUCCESS, as well as output
std::pair<int, std::string> exec(const char* cmd)
{
    std::array<char, 128> buffer;
    std::string result;
    auto pipe = popen(cmd, "r");
    
    if (!pipe) throw std::runtime_error("popen() failed!");
    
    while (!feof(pipe))
    {
        if (fgets(buffer.data(), 128, pipe) != nullptr)
            result += buffer.data();
    }
    
    auto rc = pclose(pipe);
    
    return std::make_pair(rc, result);
}


class Execution_Tests_2023partition : public CPPUNIT_NS::TestCase
{
  CPPUNIT_TEST_SUITE(Execution_Tests_2023partition);
  CPPUNIT_TEST(allparamsgood_succeeds);
  CPPUNIT_TEST(noarguments_fails);
  CPPUNIT_TEST(test_1stparamwrong_fails);
  CPPUNIT_TEST(test_2ndparamwrong_fails);
  CPPUNIT_TEST(test_3rdparamwrong_fails);
  CPPUNIT_TEST(test_4thparamwrong_fails);
  CPPUNIT_TEST(test_5thparamwrong_fails);
  // TODO: Test that incorrectly structured input file is rejected
  CPPUNIT_TEST_SUITE_END();


public:
  void setUp(void) {}
  void tearDown(void) {}

protected:


  // Test that ./2023partition 1 2 1.10 ../resources/testModelA.gv ../resources/testModelA-2part.gv (proper arguments) succeeds
  void allparamsgood_succeeds(void) {
    bool DESIRED_SUCCESS = true; // Expect success
    std::pair<int, std::string> p = exec("./2023partition 1 2 1.10 ../resources/testModelA.gv ../resources/testModelA-2part.gv");
    std::cout << " expecting_" << (DESIRED_SUCCESS ? "SUCCESS" : "FAILURE") << " exit_code:" << p.first << " stdoutl1:" << p.second.substr(0, p.second.find('\n'));
    if(DESIRED_SUCCESS ? p.first!=EXIT_SUCCESS : p.first==EXIT_SUCCESS) { 
      exit(1);
    }
  }


  // Test that ./2023partition (no arguments) fails AND returns an error message
  void noarguments_fails(void) {
    bool DESIRED_SUCCESS = false; // Expect failure
    std::pair<int, std::string> p = exec("./2023partition");
    std::cout << " expecting_" << (DESIRED_SUCCESS ? "SUCCESS" : "FAILURE") << " exit_code:" << p.first << " stdoutl1:" << p.second.substr(0, p.second.find('\n'));
    if(DESIRED_SUCCESS ? p.first!=EXIT_SUCCESS : p.first==EXIT_SUCCESS) { 
      exit(1);
    }
    if(p.second=="") {
      exit(1);
    }
  }

  // Test that ./2023partition 3.1 2 1.10 ../resources/testModelA.gv ../resources/testModelA-2part.gv (wrong arguments) fails
  // Concretely: 1st parameter is not a positive int
  void test_1stparam_fails(void) {
    bool DESIRED_SUCCESS = false; // Expect failure
    std::pair<int, std::string> p = exec("./2023partition 3.1 2 1.10 ../resources/testModelA.gv ../resources/testModelA-2part.gv");
    std::cout << " expecting_" << (DESIRED_SUCCESS ? "SUCCESS" : "FAILURE") << " exit_code:" << p.first << " stdoutl1:" << p.second.substr(0, p.second.find('\n'));
    if(DESIRED_SUCCESS ? p.first!=EXIT_SUCCESS : p.first==EXIT_SUCCESS) { 
      exit(1);
    }
  }

  // Test that ./2023partition 1 A 1.10 ../resources/testModelA.gv ../resources/testModelA-2part.gv (wrong arguments) fails
  // Concretely: 2nd parameter is not a positive int
  void test_2ndparamwrong_fails(void) {
    bool DESIRED_SUCCESS = false; // Expect failure
    std::pair<int, std::string> p = exec("./2023partition 1 A 1.10 ../resources/testModelA.gv ../resources/testModelA-2part.gv");
    std::cout << " expecting_" << (DESIRED_SUCCESS ? "SUCCESS" : "FAILURE") << " exit_code:" << p.first << " stdoutl1:" << p.second.substr(0, p.second.find('\n'));
    if(DESIRED_SUCCESS ? p.first!=EXIT_SUCCESS : p.first==EXIT_SUCCESS) { 
      exit(1);
    }
  }

  // Test that ./2023partition 1 2 0.99 ../resources/testModelA.gv ../resources/testModelA-2part.gv (wrong arguments) fails
  // Concretely: 3rd parameter is not a float greater than 1
  void test_3rdparamwrong_fails(void) {
    bool DESIRED_SUCCESS = false; // Expect failure
    std::pair<int, std::string> p = exec("./2023partition 1 2 0.95 ../resources/testModelA.gv ../resources/testModelA-2part.gv");
    std::cout << " expecting_" << (DESIRED_SUCCESS ? "SUCCESS" : "FAILURE") << " exit_code:" << p.first << " stdoutl1:" << p.second.substr(0, p.second.find('\n'));
    if(DESIRED_SUCCESS ? p.first!=EXIT_SUCCESS : p.first==EXIT_SUCCESS) { 
      exit(1);
    }
  }

  // Test that ./2023partition 1 2 1.10 ../resources/testModelaaaaA.gv ../resources/testModelA-2part.gv (wrong arguments) fails
  // Concretely: 4th parameter file does not exist
  void test_4thparamwrong_fails(void) {
    bool DESIRED_SUCCESS = false; // Expect failure
    std::pair<int, std::string> p = exec("./2023partition 1 2 1.10 ../resources/testModelaaaA.gv ../resources/testModelA-2part.gv");
    std::cout << " expecting_" << (DESIRED_SUCCESS ? "SUCCESS" : "FAILURE") << " exit_code:" << p.first << " stdoutl1:" << p.second.substr(0, p.second.find('\n'));
    if(DESIRED_SUCCESS ? p.first!=EXIT_SUCCESS : p.first==EXIT_SUCCESS) { 
      exit(1);
    }
  }

  // Test that ./2023partition 1 2 1.10 ../resources/testModelA.gv ../asbdsad12/qq.gv (wrong arguments) fails
  // Concretely: 5th parameter file is not writable
  void test_5thparamwrong_fails(void) {
    bool DESIRED_SUCCESS = false; // Expect failure
    std::pair<int, std::string> p = exec("./2023partition 1 2 1.10 ../resources/testModelA.gv ../asbdsad12/qq.gv");
    std::cout << " expecting_" << (DESIRED_SUCCESS ? "SUCCESS" : "FAILURE") << " exit_code:" << p.first << " stdoutl1:" << p.second.substr(0, p.second.find('\n'));
    if(DESIRED_SUCCESS ? p.first!=EXIT_SUCCESS : p.first==EXIT_SUCCESS) { 
      exit(1);
    }
  }





};


CPPUNIT_TEST_SUITE_REGISTRATION(Execution_Tests_2023partition);

int main()

{
  CPPUNIT_NS::TestResult controller;

  CPPUNIT_NS::TestResultCollector result;
  controller.addListener(&result);

  CPPUNIT_NS::BriefTestProgressListener progress;
  controller.addListener(&progress);

  CPPUNIT_NS::TestRunner runner;
  runner.addTest(CPPUNIT_NS::TestFactoryRegistry::getRegistry().makeTest());
  runner.run(controller);

  return result.wasSuccessful() ? 0 : 1;
}