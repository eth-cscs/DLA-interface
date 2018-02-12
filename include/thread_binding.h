#ifndef DLA_INTERFACE_THREAD_BINDING_H
#define DLA_INTERFACE_THREAD_BINDING_H

#ifdef DLA_HAVE_HWLOC
#include <hwloc.h>
#endif
#include <cstdlib>
#include <iostream>
#include <string>

namespace dla_interface {
  namespace thread {

#ifdef DLA_HAVE_HWLOC
    class CpuSet {
      public:
      CpuSet() : cpuset_(hwloc_bitmap_alloc()) {}
      CpuSet(const CpuSet& rhs) : cpuset_(hwloc_bitmap_dup(rhs.cpuset_)) {}
      CpuSet(CpuSet&& rhs) : cpuset_(rhs.cpuset_) {
        rhs.cpuset_ = nullptr;
      }

      ~CpuSet() {
        if (cpuset_) {
          hwloc_bitmap_free(cpuset_);
        }
      }

      CpuSet& operator=(const CpuSet& rhs) {
        hwloc_bitmap_copy(cpuset_, rhs.cpuset_);
        return *this;
      }
      CpuSet& operator=(CpuSet&& rhs) {
        std::swap(cpuset_, rhs.cpuset_);
        return *this;
      }

      bool operator==(const CpuSet& rhs) const {
        return hwloc_bitmap_compare(cpuset_, rhs.cpuset_) == 0;
      }
      bool operator!=(const CpuSet& rhs) const {
        return !operator==(rhs);
      }

      void add(unsigned int id) {
        hwloc_bitmap_set(cpuset_, id);
      }

      void remove(unsigned int id) {
        hwloc_bitmap_clr(cpuset_, id);
      }

      std::string str() const {
        char* tmp;
        hwloc_bitmap_asprintf(&tmp, cpuset_);
        std::string str(tmp);
        std::free(tmp);
        return str;
      }

      const hwloc_cpuset_t hwlocPtr() const {
        return cpuset_;
      }

      hwloc_cpuset_t hwlocPtr() {
        return cpuset_;
      }

      private:
      hwloc_cpuset_t cpuset_;
    };

    class SystemTopology {
      public:
      SystemTopology() {
        hwloc_topology_init(&topo_);
        hwloc_topology_load(topo_);
      }
      SystemTopology(const SystemTopology& rhs) : topo_(nullptr) {
        hwloc_topology_dup(&topo_, rhs.topo_);
      }
      SystemTopology(SystemTopology&& rhs) : topo_(rhs.topo_) {
        rhs.topo_ = nullptr;
      }

      ~SystemTopology() {
        if (topo_) {
          hwloc_topology_destroy(topo_);
        }
      }

      SystemTopology& operator=(const SystemTopology& rhs) {
        hwloc_topology_destroy(topo_);
        hwloc_topology_dup(&topo_, rhs.topo_);
        return *this;
      }
      SystemTopology& operator=(SystemTopology&& rhs) {
        std::swap(topo_, rhs.topo_);
        return *this;
      }

      CpuSet getCpuBind() {
        CpuSet cpuset;
        hwloc_get_cpubind(topo_, cpuset.hwlocPtr(), HWLOC_CPUBIND_THREAD);
        return cpuset;
      }

      void setCpuBind(const CpuSet& cpuset) {
        printThreadBindingDebugInfo("setCpuBind called:\n  Previous value:");
        hwloc_set_cpubind(topo_, cpuset.hwlocPtr(), HWLOC_CPUBIND_THREAD);
        printThreadBindingDebugInfo("  New      value:");
      }

      hwloc_topology_t hwlocPtr() {
        return topo_;
      }

      private:
      void printThreadBindingDebugInfo(const char* str) {
#ifdef DLA_THREAD_DEBUG_INFO
        auto tmp = getCpuBind();
        std::cout << str << " " << tmp.str() << std::endl;
#endif
      }

      hwloc_topology_t topo_;
    };
#else
    class CpuSet {
      public:
      bool operator==(const CpuSet& /*rhs*/) const {
        return true;
      }
      bool operator!=(const CpuSet& /*rhs*/) const {
        return false;
      }

      void add(unsigned int /*id*/) {}
      void remove(unsigned int /*id*/) {}

      std::string str() const {
        return "<Warning: Cpu binding not Enabled>";
      }
    };
    class SystemTopology {
      public:
      CpuSet getCpuBind() {
        return CpuSet();
      }

      void setCpuBind(const CpuSet& /*cpuset*/) {}
    };
#endif
  }
}

#endif  // DLA_INTERFACE_THREAD_BINDING_H
