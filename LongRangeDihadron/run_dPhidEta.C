#ifdef __CINT__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#endif

#include "include/BasicForDihadron.h"
#include "Process_dPhidEta.cxx"

void run_dPhidEta() {
    Process_dPhidEta();
    std::cout << "\nProcess_dPhidEta completed successfully!" << std::endl;
}
