#pragma once

#include <array>

template<typename T, std::size_t N>
using array_1d = std::array<T, N>;

template<typename T, std::size_t N1, std::size_t N2>
using array_2d = std::array<std::array<T, N2>, N1>;

template<typename T, std::size_t N1, std::size_t N2, std::size_t N3>
using array_3d = std::array<std::array<std::array<T, N3>, N2>, N1>;

template<typename T, std::size_t N1, std::size_t N2, std::size_t N3, std::size_t N4>
using array_4d = std::array<std::array<std::array<std::array<T, N4>, N3>, N2>, N1>;

template<typename T, std::size_t N1, std::size_t N2, std::size_t N3, std::size_t N4, std::size_t N5>
using array_5d = std::array<std::array<std::array<std::array<std::array<T, N5>, N4>, N3>, N2>, N1>;

template<typename T, std::size_t N1, std::size_t N2, std::size_t N3, std::size_t N4, std::size_t N5, std::size_t N6>
using array_6d = std::array<std::array<std::array<std::array<std::array<std::array<T, N6>, N5>, N4>, N3>, N2>, N1>;