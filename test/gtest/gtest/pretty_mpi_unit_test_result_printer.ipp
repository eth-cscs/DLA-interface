// Copyright (c) 2016, Eidgenössische Technische Hochschule Zürich and
// Forschungszentrum Jülich GmbH.
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without modification,
// are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice, this
//    list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice,
//    this list of conditions and the following disclaimer in the documentation
//    and/or other materials provided with the distribution.
//
// 3. Neither the name of the copyright holder nor the names of its contributors
//    may be used to endorse or promote products derived from this software without
//    specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
// ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
// ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

// Source: Arbor library
// Note: the file has the following modifications:
//       - plain MPI functions
//       - rank 0 uses PrettyUnitTestResultPrinter output style
//         when printing to standard output.

#ifndef GTEST_PRETTY_MPI_UNIT_TEST_RESULT_PRINTER_H
#define GTEST_PRETTY_MPI_UNIT_TEST_RESULT_PRINTER_H

#include <cstdio>
#include <fstream>
#include <stdexcept>

#include <mpi.h>

#include "gtest/gtest.h"

/// A specialized listener designed for printing test results with MPI.
///
/// When tests are run with MPI, one instance of each test is run on
/// each rank. The default behavior of Google Test is for each test
/// instance to print to stdout. With more than one MPI rank, this creates
/// the usual MPI mess of output.
///
/// This specialization has the first rank (rank 0) print to stdout, and all MPI
/// ranks print their output to separate text files.
/// For each test a message is printed showing
///     - detailed messages about errors on rank 0
///     - a head count of errors that occured on other MPI ranks
#define skip_printing_failure "**** should not be printed ****"

namespace testing {
  namespace internal {
    AssertionResult AssertTestPassedGlobally(const char*, const std::vector<int>& failures) {
      bool flag = true;
      for (size_t i = 0; i < failures.size(); ++i) {
        if (failures[i] > 0) {
          std::cout << "    Rank " << i << " FAILED " << failures[i] << " tests\n" << std::endl;
          flag = false;
        }
      }
      return flag ? AssertionSuccess() : AssertionFailure();
    }
  }

  class PrettyMpiUnitTestResultPrinter : public internal::PrettyUnitTestResultPrinter {
    private:
    int rank_;
    int size_;
    std::ofstream fid_;
    char buffer_[1024];
    int test_case_failures_;
    int test_case_tests_;
    int test_failures_;

    bool does_print() const {
      return rank_ == 0;
    }

    void print(const char* s) {
      if (fid_) {
        fid_ << s;
      }
    }

    void print(const std::string& s) {
      print(s.c_str());
    }

    /// convenience function that handles the logic of using snprintf
    /// and forwarding the results to file and/or stdout.
    template <typename... Args>
    void printf_helper(const char* s, Args&&... args) {
      std::snprintf(buffer_, sizeof(buffer_), s, std::forward<Args>(args)...);
      print(buffer_);
    }

    public:
    PrettyMpiUnitTestResultPrinter(std::string f_base = "") {
      MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
      MPI_Comm_size(MPI_COMM_WORLD, &size_);

      if (f_base.empty()) {
        return;
      }
      std::string fname = f_base + "_" + std::to_string(rank_) + ".txt";
      fid_.open(fname);
      if (!fid_) {
        throw std::runtime_error("could not open file " + fname + " for test output");
      }
    }

    // Returns the file stream (can be used to output extra informations durint the test)
    std::ofstream& outStream() {
      return fid_;
    }

    /// Messages that are printed at the start and end of the test program.
    /// i.e. once only.
    virtual void OnTestProgramStart(const UnitTest&) override {
      printf_helper("*** test output for rank %d of %d\n\n", rank_, size_);
    }
    virtual void OnTestProgramEnd(const UnitTest&) override {
      printf_helper("*** end test output for rank %d of %d\n", rank_, size_);
    }

    virtual void OnEnvironmentsSetUpStart(const UnitTest& unit_test) override {
      if (does_print()) {
        PrettyUnitTestResultPrinter::OnEnvironmentsSetUpStart(unit_test);
      }
    }
    /// Messages that are printed at the start and end of each test case.
    /// On startup a counter that counts the number of tests that fail in
    /// this test case is initialized to zero, and will be incremented for each
    /// test that fails.
    virtual void OnTestCaseStart(const TestCase& test_case) override {
      test_case_failures_ = 0;
      test_case_tests_ = 0;
      if (does_print()) {
        PrettyUnitTestResultPrinter::OnTestCaseStart(test_case);
      }
    }
    virtual void OnTestCaseEnd(const TestCase& test_case) override {
      printf_helper("    PASSED %d of %d tests in %s\n", test_case_tests_ - test_case_failures_,
                    test_case_tests_, test_case.name());
      if (test_case_failures_ > 0) {
        printf_helper("    FAILED %d of %d tests in %s\n", test_case_failures_, test_case_tests_,
                      test_case.name());
      }
      print("\n");
      if (does_print()) {
        PrettyUnitTestResultPrinter::OnTestCaseEnd(test_case);
      }
    }

    // Called before a test starts.
    virtual void OnTestStart(const TestInfo& test_info) override {
      printf_helper("TEST:  %s::%s\n", test_info.test_case_name(), test_info.name());
      test_failures_ = 0;
      if (does_print()) {
        PrettyUnitTestResultPrinter::OnTestStart(test_info);
      }
    }

    virtual void OnTestIterationStart(const UnitTest& unit_test, int iteration) override {
      if (does_print()) {
        PrettyUnitTestResultPrinter::OnTestIterationStart(unit_test, iteration);
      }
    }

    // Called after a failed assertion or a SUCCEED() invocation.
    virtual void OnTestPartResult(const TestPartResult& test_part_result) override {
      if (strcmp(test_part_result.message(), skip_printing_failure) == 0)
        return;

      // indent all lines in the summary by 4 spaces
      std::string summary = "    " + std::string(test_part_result.summary());
      auto pos = summary.find("\n");
      while (pos != summary.size() && pos != std::string::npos) {
        summary.replace(pos, 1, "\n    ");
        pos = summary.find("\n", pos + 1);
      }

      printf_helper("  LOCAL_%s\n    %s:%d\n%s\n", test_part_result.failed() ? "FAIL" : "SUCCESS",
                    test_part_result.file_name(), test_part_result.line_number(), summary.c_str());

      // note that there was a failure in this test case
      if (test_part_result.failed()) {
        test_failures_++;
      }

      if (does_print()) {
        PrettyUnitTestResultPrinter::OnTestPartResult(test_part_result);
      }
    }

    virtual void OnTestIterationEnd(const UnitTest& unit_test, int iteration) override {
      if (does_print()) {
        PrettyUnitTestResultPrinter::OnTestIterationEnd(unit_test, iteration);
      }
    }

    // Called after a test ends.
    virtual void OnTestEnd(const TestInfo& test_info) override {
      test_case_tests_++;

      // gather the number of failures to report them
      std::vector<int> failures(size_);
      MPI_Allgather(&test_failures_, 1, MPI_INT, &failures[0], 1, MPI_INT, MPI_COMM_WORLD);

      // count the number of ranks that had errors
      int local_error = test_failures_ > 0;
      int global_errors = 0;
      MPI_Allreduce(&local_error, &global_errors, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

      if (global_errors > 0) {
        test_case_failures_++;
        printf_helper("  GLOBAL_FAIL on %d ranks\n", global_errors);
        for (size_t i = 0; i < failures.size(); ++i) {
          if (failures[i] > 0)
            printf_helper("    Rank %d FAILED %d tests\n", i, failures[i]);
        }
      }
      if (does_print()) {
        for (size_t i = 0; i < failures.size(); ++i) {
          if (failures[i] > 0) {
            std::cout << "    Rank " << i << " FAILED " << failures[i] << " tests" << std::endl;
            GTEST_NONFATAL_FAILURE_(skip_printing_failure);
          }
        }
        PrettyUnitTestResultPrinter::OnTestEnd(test_info);
      }
    }

    virtual void OnEnvironmentsTearDownStart(const UnitTest& unit_test) override {
      if (does_print()) {
        PrettyUnitTestResultPrinter::OnEnvironmentsTearDownStart(unit_test);
      }
    }
  };

  // Set the unit test listener to PrettyMpiUnitTestResultPrinter.
  // The name of the output files is f_base + "_" + mpi_rank.txt
  std::ofstream& setMPIListener(std::string f_base) {
    auto& listeners = testing::UnitTest::GetInstance()->listeners();
    // first delete the original printer
    delete listeners.Release(listeners.default_result_printer());
    // now add our custom printer
    auto mpi_listener = new PrettyMpiUnitTestResultPrinter(f_base);
    listeners.Append(mpi_listener);

    return mpi_listener->outStream();
  }
}

#endif
