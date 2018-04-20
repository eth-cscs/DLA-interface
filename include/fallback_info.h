#ifndef DLA_INTERFACE_FALLBACK_INFO_H
#define DLA_INTERFACE_FALLBACK_INFO_H

#include <iostream>
#include <string>
#include "types.h"
#include "util_types.h"

namespace dla_interface {
  class FallbackInfo {
    public:
    void setFullReport(bool flag) {
      full_report_ = flag;
    }
    bool getFullReport() const {
      return full_report_;
    }
    void report(const char* func, SolverType old_solver, SolverType new_solver, std::string reason) {
      if (full_report_) {
        std::cerr << "WARNING (DLA_interface): \"" << func << "\" executed with " << new_solver
                  << " instead of " << old_solver << ". (" << reason << ")" << std::endl;
        ++status[func].first;
      }
      else {
        ++status[func].second;
      }
    }

    void finalReport(int rank, std::ostream& out) {
      for (auto val : status) {
        out << "DLA_interface report (Rank " << rank << "): \"" << val.first
            << "\" was executed with a different solver " << val.second.first
            << " (with full report) + " << val.second.second << " (without full report) times.\n";
      }
      status.clear();
    }

    private:
    std::map<std::string, std::pair<unsigned, unsigned>> status;
    bool full_report_;
  };
}

#endif  // DLA_INTERFACE_FALLBACK_INFO_H
