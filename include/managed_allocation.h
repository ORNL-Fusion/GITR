#ifndef _MANAGED_
#define _MANAGED_
class ManagedAllocation {
    public:
    void *operator new(size_t len) { void *ptr; cudaMallocManaged(&ptr, len); cudaDeviceSynchronize(); return ptr;
    }
    void operator delete(void *ptr) { cudaDeviceSynchronize(); cudaFree(ptr);
    } 
};
#endif
