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

#include "managed_allocation.h"
#include <memory>
#include "device.h"
#include <cstddef>
#include <stdexcept>
#include "Boundary.h"
#include <iomanip>
#include <iostream>
#include <vector>
namespace sim {


/*! Array like managed memory structure
 * Fixed capacity but variable size array structure allocated using cudaMallocManaged if CUDA is defined,
 * else system allocator If CUDA is defined both the data structure and underlying pointed to memory
 * are allocated using cudaMallocManaged
 */
  template<typename T>
  class Array : public ManagedAllocation {

    // @todo constexpr in c++14
    T *alloc_data() {
#if USE_CUDA > 0
      T* data;
      auto err = cudaMallocManaged(&data, sizeof(T)*capacity_);
      if(err != cudaSuccess){
        throw std::runtime_error("error allocating managed memory");
      }
      /* Captain! What are the contents of the "bin" variables in the CUDA vs OpenMP case?? */
      return data;
#else
      return new T[capacity_];
#endif
    }

    // @todo constexpr in c++14
    void free_data() {
#if USE_CUDA > 0
      cudaFree((void*)data_);
#else
      delete[] data_;
#endif
    }

  public:
    /*! Construct an array of fixed capacity
     */
    Array() : capacity_{1}, size_{1}, data_{alloc_data()} {}
    
    Array(const std::size_t capacity) : capacity_{capacity}, size_{capacity}, data_{alloc_data()} {}

    /*! Construct an array of fixed capacity and initialize values
     */
    Array(const std::size_t capacity, T initial_value) : capacity_{capacity}, size_{capacity}, data_{alloc_data()} {
        for (std::size_t i = 0; i < size_; i++) {
        data_[i] = initial_value;
      }
    }
    
    Array(const std::vector<T> initial_vector) : capacity_{initial_vector.size()}, size_{initial_vector.size()}, data_{alloc_data()} {
        for (std::size_t i = 0; i < size_; i++) {
        data_[i] = initial_vector[i];
      }
    }

    //Array(const std::size_t capacity, T initial_value, int on) : capacity_{capacity}, size_{capacity}, data_{alloc_data()} {
 
    //    std::cout << "capacity " << capacity << std::endl;
    //    std::cout << "size_ " << size_ << std::endl;
    //    std::cout << "inside array creation int " << on << std::endl;
    //    std::size_t major_size = size_;
    //    for (std::size_t i = 0; i < major_size; i++) {
    //    std::cout << "inside array loop int "<< i << std::endl;
    //    data_[i] = initial_value;
    //  }
    //}

    /*! Destruct array memory
     */
    ~Array() {
      free_data();
    }

    /* Captain! Is this needed? Delete the copy assignment and copy constructor */
    /*! Copy constructor
     */
    Array(const Array &source) = delete;
    /*
    {
      capacity_ = source.capacity_;
      size_ = source.size_;
      data_ = alloc_data();

      if (data_) {
        for (std::size_t i = 0; i < size_; i++) {
          data_[i] = source[i];
        }
      }
    }
    */

    Array &operator=(const Array &source) = delete;
    /*
    {
        for(int i=0;i<source.size();i++)
        {
            data_[i] = source[i];
        }
        return *this;
    }
    */
    
    Array &operator=(const std::vector<T> &source)// = delete;
    {
        for(int i=0;i<source.size();i++)
        {
            data_[i] = source[i];
        }
        return *this;
    }
    Array(Array &&) noexcept = delete;

    Array &operator=(Array &&)      = delete;

    /*! Array size getter
     * @return the number of in use elements in the array
     */
    DEVICE_CALLABLE
    std::size_t size() const {
      return size_;
    }

    /*! Array capacity getter
     * @return the maximum number of elements in the array
     */
    std::size_t capacity() const {
      return capacity_;
    }

    /*! Array capacity getter
     * @return the maximum number of elements in the array
     */
    std::size_t available() const {
      return capacity() - size();
    }

    /*! Return reference to first element of Array
     * @return reference to first element
     */
    T &front() {
      return data_[0];
    }
    /*       
    *  This function will %resize the %vector to the specified
               *         *  number of elements.  If the number is smaller than the
               *                *  %vector's current size the %vector is truncated, otherwise
               *                       *  default constructed elements are appended.
               *                              */
    void resize(const T __new_size)
    {
      if (__new_size > size())
      {
          free_data();
          capacity_ = __new_size;
          size_ = __new_size;
          data_ = alloc_data();
          
      }      
      else if (__new_size < size())
      {

      }
    }
    /*! Add element to end of the array
     * Copies argument to the back of array and increased the size by one
     */
/*
    void push_back(const T &value) {
      if (size_ + 1 > capacity_)
        throw std::runtime_error("Not enough capacity to push_back");
      else
        data_[size_++] = value;
    }
*/
    /*! Add multiple elements of a single value to end of the array
     * Copies argument value push_count times and increased the size by push_count
     */
/*
    void push_back(const T &value, const size_t push_count) {
      if (size_ + 1 > capacity_)
        throw std::runtime_error("Not enough capacity to push_back");
      else {
        for (std::size_t i = 0; i < push_count; i++) {
          this->push_back(value);
        }
      }
    }
*/
    /*! Add multiple elements to end of the array
     * Copies elements starting at argument values_ptr to the back of array and increased the size by push_count
     */
/*
    void push_back(const T *values_ptr, const size_t push_count) {
      if (size_ + 1 > capacity_)
        throw std::runtime_error("Not enough capacity to push_back");
      else {
        for (std::size_t i = 0; i < push_count; i++) {
          this->push_back(values_ptr[i]);
        }
      }
    }
*/
    /*! Remove element from end of the array
     * Remove element from end of array by reducing size by 1, element is not destructed
     */
    /*
    void pop_back() {
      if (size_ == 0)
        throw std::runtime_error("Array popped_back with 0 size");
      else
        size_--;
    }
    */

    /*! Remove multiple elements from end of the array
     * Remove element from end of array by reducing size by 1, element is not destructed
     */
     /*
    void pop_back(std::size_t pop_count) {
      if (size_ == 0 && pop_count != 0)
        throw std::runtime_error("Array popped_back with 0 size");
      else
        size_ -= pop_count;
    }
    */

    /*! Getter for pointer to underlying data
     */
    DEVICE_CALLABLE
    T *data() {
      return this->data_;
    }

    /*! const getter for pointer to underlying data
     */
    DEVICE_CALLABLE
    T *data() const {
      return this->data_;
    }

    /*! Subscript operator, []
     * Retrieve reference to element using subscript notation
     */
    DEVICE_CALLABLE
    T &operator[](const std::size_t index) {
      return data_[index];
    }

    /*! const subscript operator, []
     *  Retrieve const reference to element using subscript notation
     */
    DEVICE_CALLABLE
    const T &operator[](const std::size_t index) const {
      return data_[index];
    }
/*
    DEVICE_CALLABLE
    const T* begin() {
      return this->data();
    }

    DEVICE_CALLABLE
    const T* end() {
      return this->data() + this->size();
    }
*/
    DEVICE_CALLABLE
    T *begin() const {
      return this->data();
    }

    DEVICE_CALLABLE
    T *end() const {
      return this->data() + this->size();
    }


//private:
// @todo DEVICE_CALLABLE cant use private member variables
  public:
    std::size_t capacity_;
    std::size_t size_;
    T *data_;

void print(std::string const label)
{
  std::cout << label << '\n';
  //if constexpr (std::is_floating_point<P>::value)
  //{
    for (auto i = 0; i < size_; ++i)
      std::cout << std::setw(12) << std::setprecision(4) << std::scientific
                << std::right << data_[i];
  //}
  //else
  //{
  //  for (auto i = 0; i < size(); ++i)
  //    std::cout << std::right << (*this)(i) << " ";
  //}
  std::cout << '\n';
}
  };

  /*! begin iterator for range based for loops
   */
  /*
  template<typename T>
  DEVICE_CALLABLE
  const T *begin(const Array<T> &array) {
    return array.data();
  }
  */

  /*! end iterator for range based for loops
   */
   /*
  template<typename T>
  DEVICE_CALLABLE
  const T *end(const Array<T> &array) {
    return array.data() + array.size();
  }
  */




}
