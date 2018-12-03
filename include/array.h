#ifndef _ARRAY_
#define _ARRAY_
#ifdef __CUDACC__
#define CUDA_CALLABLE_MEMBER_DEVICE __device__
#define CUDA_CALLABLE_MEMBER_HOST __host__
#else
#define CUDA_CALLABLE_MEMBER_DEVICE
#define CUDA_CALLABLE_MEMBER_HOST
#endif

#include "managed_allocation.h"
namespace sim {
template <class T>
class Array : public ManagedAllocation {
    size_t n;
    T* data_;
    // @todo constexpr in c++14
    void free_data() {
    #if defined(CUDA)
    cudaFree((void*)data_);
    #else
    delete[] data_;
    #endif
    }
    public:
    //Array (const Array &a) {
    //    n = a.n; 
    //    cudaMallocManaged(&data, n); 
    //    memcpy(data, a.data, n);
    //}
    Array (const size_t nn) {
        n=nn; 
        cudaMallocManaged(&data_, n*sizeof(T)); 
    }
    Array (const size_t nn, T initial_value) {
        n=nn; 
        cudaMallocManaged(&data_, n*sizeof(T)); 
        for (std::size_t i = 0; i < n; i++) 
            data_[i] = initial_value;
    }
    // Also have to implement operator[], for example
    // // ...
    ~Array() {
       free_data();
     }
     __host__ __device__
     T *data() const {
              return this->data_;
                  }
     __host__ __device__
     T &operator[](const std::size_t index) {
       return data_[index];
     }
     
     __host__ __device__
     std::size_t size() const {
       return n;
     }
     __host__ __device__
                 T *begin() const {
                           return this->data();
                               }

     __host__ __device__
                     T *end() const {
                               return this->data() + this->size();
                                   }
    T &front() {
         return data_[0];
        }     
 };
}
#endif
