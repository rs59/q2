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


class Test : public CPPUNIT_NS::TestCase
{
  CPPUNIT_TEST_SUITE(Test);
  CPPUNIT_TEST(test2023Partition_execution);
  CPPUNIT_TEST(test2023Partition_execution2);
  CPPUNIT_TEST(test2023Partition_execution3);
  CPPUNIT_TEST(test2023Partition_execution4);
  CPPUNIT_TEST(test2023Partition_execution5);
  CPPUNIT_TEST(test2023Partition_execution6);
  CPPUNIT_TEST(test2023Partition_execution7);
  CPPUNIT_TEST_SUITE_END();

public:
  void setUp(void) {}
  void tearDown(void) {}

protected:

  // Test that ./2023partition (no arguments) fails AND returns an error message
  void test2023Partition_execution(void) {
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

  // Test that ./2023partition 1 2 1.10 ../resources/testModelA.gv ../resources/testModelA-2part.gv (proper arguments) succeeds
  void test2023Partition_execution2(void) {
    bool DESIRED_SUCCESS = true; // Expect success
    std::pair<int, std::string> p = exec("./2023partition 1 2 1.10 ../resources/testModelA.gv ../resources/testModelA-2part.gv");
    std::cout << " expecting_" << (DESIRED_SUCCESS ? "SUCCESS" : "FAILURE") << " exit_code:" << p.first << " stdoutl1:" << p.second.substr(0, p.second.find('\n'));
    if(DESIRED_SUCCESS ? p.first!=EXIT_SUCCESS : p.first==EXIT_SUCCESS) { 
      exit(1);
    }
  }

  // Test that ./2023partition 3.1 2 1.10 ../resources/testModelA.gv ../resources/testModelA-2part.gv (wrong arguments) fails
  // Concretely: 1st parameter is not a positive int
  void test2023Partition_execution3(void) {
    bool DESIRED_SUCCESS = false; // Expect failure
    std::pair<int, std::string> p = exec("./2023partition 3.1 2 1.10 ../resources/testModelA.gv ../resources/testModelA-2part.gv");
    std::cout << " expecting_" << (DESIRED_SUCCESS ? "SUCCESS" : "FAILURE") << " exit_code:" << p.first << " stdoutl1:" << p.second.substr(0, p.second.find('\n'));
    if(DESIRED_SUCCESS ? p.first!=EXIT_SUCCESS : p.first==EXIT_SUCCESS) { 
      exit(1);
    }
  }

  // Test that ./2023partition 1 A 1.10 ../resources/testModelA.gv ../resources/testModelA-2part.gv (wrong arguments) fails
  // Concretely: 2nd parameter is not a positive int
  void test2023Partition_execution4(void) {
    bool DESIRED_SUCCESS = false; // Expect failure
    std::pair<int, std::string> p = exec("./2023partition 1 A 1.10 ../resources/testModelA.gv ../resources/testModelA-2part.gv");
    std::cout << " expecting_" << (DESIRED_SUCCESS ? "SUCCESS" : "FAILURE") << " exit_code:" << p.first << " stdoutl1:" << p.second.substr(0, p.second.find('\n'));
    if(DESIRED_SUCCESS ? p.first!=EXIT_SUCCESS : p.first==EXIT_SUCCESS) { 
      exit(1);
    }
  }

  // Test that ./2023partition 1 2 0.99 ../resources/testModelA.gv ../resources/testModelA-2part.gv (wrong arguments) fails
  // Concretely: 3rd parameter is not a float greater than 1
  void test2023Partition_execution5(void) {
    bool DESIRED_SUCCESS = false; // Expect failure
    std::pair<int, std::string> p = exec("./2023partition 1 2 0.95 ../resources/testModelA.gv ../resources/testModelA-2part.gv");
    std::cout << " expecting_" << (DESIRED_SUCCESS ? "SUCCESS" : "FAILURE") << " exit_code:" << p.first << " stdoutl1:" << p.second.substr(0, p.second.find('\n'));
    if(DESIRED_SUCCESS ? p.first!=EXIT_SUCCESS : p.first==EXIT_SUCCESS) { 
      exit(1);
    }
  }

  // Test that ./2023partition 1 2 1.10 ../resources/testModelaaaaA.gv ../resources/testModelA-2part.gv (wrong arguments) fails
  // Concretely: 4th parameter file does not exist
  void test2023Partition_execution6(void) {
    bool DESIRED_SUCCESS = false; // Expect failure
    std::pair<int, std::string> p = exec("./2023partition 1 2 1.10 ../resources/testModelaaaA.gv ../resources/testModelA-2part.gv");
    std::cout << " expecting_" << (DESIRED_SUCCESS ? "SUCCESS" : "FAILURE") << " exit_code:" << p.first << " stdoutl1:" << p.second.substr(0, p.second.find('\n'));
    if(DESIRED_SUCCESS ? p.first!=EXIT_SUCCESS : p.first==EXIT_SUCCESS) { 
      exit(1);
    }
  }

  // Test that ./2023partition 1 2 1.10 ../resources/testModelA.gv ../asbdsad12/qq.gv (wrong arguments) fails
  // Concretely: 5th parameter file is not writable
  void test2023Partition_execution7(void) {
    bool DESIRED_SUCCESS = false; // Expect failure
    std::pair<int, std::string> p = exec("./2023partition 1 2 1.10 ../resources/testModelA.gv ../asbdsad12/qq.gv");
    std::cout << " expecting_" << (DESIRED_SUCCESS ? "SUCCESS" : "FAILURE") << " exit_code:" << p.first << " stdoutl1:" << p.second.substr(0, p.second.find('\n'));
    if(DESIRED_SUCCESS ? p.first!=EXIT_SUCCESS : p.first==EXIT_SUCCESS) { 
      exit(1);
    }
  }





};


CPPUNIT_TEST_SUITE_REGISTRATION(Test);

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