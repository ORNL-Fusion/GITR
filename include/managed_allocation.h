/*
The MIT License (MIT)

Copyright (c) 2016 Adam Simpson

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
*/

#pragma once

#if USE_CUDA > 0
#include "cuda_runtime.h"
#endif

#include <iostream>
#include <new>

/*! Inheritable class which overload new/delete operators to use cuda managed memory if CUDA is defined
 */
class ManagedAllocation {
public:
  ManagedAllocation() = default;

  ~ManagedAllocation() = default;

  ManagedAllocation(const ManagedAllocation &) = delete;

  ManagedAllocation &operator=(const ManagedAllocation &) = delete;

  ManagedAllocation(ManagedAllocation &&) noexcept = delete;

  ManagedAllocation &operator=(ManagedAllocation &&)      = delete;

  /*! new operator
   */
  static void *operator new(std::size_t size) {
#if USE_CUDA > 0
    void* data;
    auto err = cudaMallocManaged(&data, size);
    if(err != cudaSuccess) {
      //throw std::runtime_error("error allocating managed memory");
    }
    return data;
#else
    return ::operator new(size);
#endif
  }

  /*! new array operator
   */
  static void *operator new[](std::size_t size) {
    return ManagedAllocation::operator new(size);
  }

  /*! delete operator
   */
  static void operator delete(void *block) {
#if USE_CUDA > 0
    cudaFree((void*)block);
#else
    ::operator delete(block);
#endif
  }

  /*! delete array operator
   */
  static void operator delete[](void *block) {
    ManagedAllocation::operator delete(block);
  }

};
