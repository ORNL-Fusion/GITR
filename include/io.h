#pragma once

#include <algorithm>
#include <array>
#include <experimental/iterator>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <span>

#include "csv.h"

template <std::size_t N>
class csv_row_stacker {
public:
  explicit csv_row_stacker(std::filesystem::path file_name);

  ~csv_row_stacker();

  /* stack a row in the csv: time followed by n_columns of data */
  template <typename T>
  bool stack(double t, std::span<T, N> p);
  template <typename T, typename = std::enable_if_t<std::is_arithmetic_v<T>>>
  bool stack(double t, T p);

private:
  std::filesystem::path file_path;
  std::fstream stream;
};

template <std::size_t N>
csv_row_stacker<N>::csv_row_stacker(std::filesystem::path file_name) : file_path(std::move(file_name)) {
  stream.open(file_path, std::ios_base::out);

  if (stream.is_open() == false) {
    std::cout << "file error creating: " << file_path << '\n';
    exit(1);
  }
}

template <std::size_t N>
csv_row_stacker<N>::~csv_row_stacker() = default;

template <std::size_t N>
template <typename T>
bool csv_row_stacker<N>::stack(double t, std::span<T, N> p) {
  stream << std::setprecision(20) << t << ',';
  std::copy(std::begin(p), std::end(p), std::experimental::make_ostream_joiner(stream, ','));
  stream << '\n';
  return stream.good();
}

template <std::size_t N>
template <typename T, typename>
bool csv_row_stacker<N>::stack(double t, T p) {
  stream << std::setprecision(20) << t << ',' << p << '\n';
  return stream.good();
}
