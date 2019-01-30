#include <cstring>
#include <fstream>
#include <iostream>
#include <string>
#include "gtest/gtest.h"
#include "mpi_listener.h"
#include "communicator_grid.h"
#include "communicator_manager.h"
#include "test_ftn_dlai_main.h"
#include "test_ftn_setup.h"

using namespace dla_interface;
using namespace testing;

std::vector<comm::Communicator2DGrid*> comms;
std::ofstream* outstream;

int main(int argc, char** argv) {
  int set_comms;
  int icomms[8];
  test_ftn_dlai_setup(icomms, &set_comms);

  for (int i = 0; i < set_comms; ++i)
    comms.push_back(&comm::CommunicatorManager::getCommunicator2DGridFromBlacsContext(icomms[i]));

  comm::CommunicatorManager::getFallbackInfo().setFullReport(true);

  ::testing::InitGoogleTest(&argc, argv);

  const char* name = strrchr(argv[0], '/');
  if (name == nullptr)
    name = argv[0];
  else
    ++name;
  outstream = &::testing::setMPIListener(std::string("results_") + name);

  auto ret = RUN_ALL_TESTS();

  test_ftn_dlai_end();

  return ret;
}
