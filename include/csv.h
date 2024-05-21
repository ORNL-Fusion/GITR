#pragma once

#include <algorithm>
#include <charconv>
#include <filesystem>
#include <fstream>
#include <string>
#include <string_view>
#include <vector>

template <typename InputIt, typename ForwardIt, typename BinOp>
void for_each_token(InputIt first, InputIt last, ForwardIt s_first, ForwardIt s_last, BinOp binary_op) {
  while (first != last) {
    const auto pos = std::find_first_of(first, last, s_first, s_last);

    binary_op(first, pos);

    if (pos == last)
      break;

    first = std::next(pos);
  }
}

template <typename T>
std::vector<T> split(std::string_view str) {
  constexpr char delims[] = ", \n\t\r\f";

  std::vector<T> row;

  for_each_token(std::cbegin(str), std::cend(str), std::cbegin(delims), std::cend(delims),
                 [&row](auto first, auto second) {
                   if (first != second) {
                     T value;
                     std::from_chars(first, second, value);
                     row.push_back(value);
                   }
                 });

  return row;
}

template <typename T>
std::vector<std::vector<T>> transpose_matrix(const std::vector<std::vector<T>> &rows) {
  int r = rows.size();

  int const c = rows[0].size();

  std::vector<std::vector<T>> columns(c);

  for (int i = 0; i < c; i++) {
    columns[i].resize(r);
  }

  /* change r x c to c x r */
  for (int i = 0; i < r; i++) {
    for (int j = 0; j < c; j++) {
      columns[j][i] = rows[i][j];
    }
  }
  return columns;
}

template <typename T>
std::vector<std::vector<T>> csv_columns(std::filesystem::path const &file_name, bool skip_first, bool transpose) {
  std::vector<std::vector<T>> rows;

  /* open the file */
  std::fstream file(file_name, std::ios::in);

  /* you should really assert things like this... */
  if (!file) {
    std::stringstream message;
    message << "csv file " << file_name << " failed to open!";
    throw std::runtime_error(message.str());
  }

  std::string line;

  if (skip_first) {
    getline(file, line);
  }

  while (getline(file, line)) {
    std::vector<T> row = split<T>(line);
    if (!row.empty()) {
      rows.push_back(std::move(row));
    }
  }

  auto number_columns = rows[0].size();
  if (!std::all_of(std::next(std::cbegin(rows)), std::cend(rows),
                   [number_columns](auto row) { return row.size() == number_columns; })) {
    std::stringstream message;
    message << "Error loading csv file " << file_name << ". Number of columns not identical!";
    throw std::runtime_error(message.str());
  }

  /* transpose and return columns */
  if (transpose) {
    return transpose_matrix(rows);
  } else {
    return rows;
  }
}
